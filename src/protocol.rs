//! Boundary protocol: how the cell is driven over time.
//!
//! Every variant ultimately produces a pair of electrode potentials each
//! step (returned as `BoundaryPotentials`). Galvanostatic variants use a PI
//! controller that adjusts the applied voltage based on the previous step's
//! measured current — the underlying solver remains potentiostatic, but
//! the outer loop tracks a current target the way a real galvanostat does.

use crate::poisson::BoundaryPotentials;

/// PI current controller with low-pass filtered input. Drives applied
/// voltage so the *filtered* measured current converges on `target`. The
/// filter is essential because per-step current is dominated by thermal
/// fluctuations (~10x the typical signal); without smoothing, the PI
/// chases noise. Anti-windup is handled by clamping the integrator to
/// `voltage_max`.
#[derive(Clone, Debug)]
pub struct CurrentController {
    pub target: f32,
    pub kp: f32,
    pub ki: f32,
    /// Exponential smoothing on input current. `1.0` = no smoothing,
    /// `0.05` = ~20-step time constant.
    pub filter_alpha: f32,
    pub voltage_max: f32,
    pub voltage: f32,
    pub integral: f32,
    pub filtered_current: f32,
}

impl CurrentController {
    pub fn new(target: f32, kp: f32, ki: f32, filter_alpha: f32, voltage_max: f32) -> Self {
        Self {
            target,
            kp,
            ki,
            filter_alpha,
            voltage_max,
            voltage: 0.0,
            integral: 0.0,
            filtered_current: 0.0,
        }
    }

    pub fn update(&mut self, measured_current: f32, dt: f32) -> f32 {
        let alpha = self.filter_alpha.clamp(0.0, 1.0);
        self.filtered_current =
            (1.0 - alpha) * self.filtered_current + alpha * measured_current;

        let err = self.target - self.filtered_current;
        self.integral =
            (self.integral + self.ki * err * dt).clamp(-self.voltage_max, self.voltage_max);
        self.voltage =
            (self.kp * err + self.integral).clamp(-self.voltage_max, self.voltage_max);
        self.voltage
    }
}

#[derive(Clone, Debug)]
pub enum Protocol {
    /// Hold a fixed terminal voltage forever.
    Potentiostatic { voltage: f32 },
    /// Hold a fixed terminal current forever (PI controller drives V).
    Galvanostatic { controller: CurrentController },
    /// Periodic voltage step: relax at 0 then pulse at `pulse_voltage`.
    StepVoltage {
        relax_steps: usize,
        pulse_steps: usize,
        pulse_voltage: f32,
    },
    /// Periodic current step: relax at 0 then pulse at `pulse_current` via PI.
    StepCurrent {
        relax_steps: usize,
        pulse_steps: usize,
        pulse_current: f32,
        controller: CurrentController,
    },
}

#[derive(Clone, Debug)]
pub struct ProtocolState {
    pub protocol: Protocol,
    pub step_index: usize,
}

impl ProtocolState {
    pub fn new(protocol: Protocol) -> Self {
        Self {
            protocol,
            step_index: 0,
        }
    }

    /// Compute the BCs to apply this step. `last_current` is the current
    /// measured at the end of the previous step (used by galvanostatic
    /// variants). Increments `step_index`.
    pub fn tick(&mut self, dt: f32, last_current: f32) -> BoundaryPotentials {
        let v = match &mut self.protocol {
            Protocol::Potentiostatic { voltage } => *voltage,
            Protocol::Galvanostatic { controller } => controller.update(last_current, dt),
            Protocol::StepVoltage {
                relax_steps,
                pulse_steps,
                pulse_voltage,
            } => {
                let cycle = *relax_steps + *pulse_steps;
                let phase = self.step_index % cycle.max(1);
                if phase < *relax_steps {
                    0.0
                } else {
                    *pulse_voltage
                }
            }
            Protocol::StepCurrent {
                relax_steps,
                pulse_steps,
                pulse_current,
                controller,
            } => {
                let cycle = *relax_steps + *pulse_steps;
                let phase = self.step_index % cycle.max(1);
                controller.target = if phase < *relax_steps {
                    0.0
                } else {
                    *pulse_current
                };
                controller.update(last_current, dt)
            }
        };
        self.step_index += 1;
        BoundaryPotentials {
            left: 0.5 * v,
            right: -0.5 * v,
        }
    }

    /// Currently-targeted current, if any. Useful for measurement reporting
    /// in galvanostatic modes.
    pub fn target_current(&self) -> Option<f32> {
        match &self.protocol {
            Protocol::Galvanostatic { controller } => Some(controller.target),
            Protocol::StepCurrent { controller, .. } => Some(controller.target),
            _ => None,
        }
    }
}
