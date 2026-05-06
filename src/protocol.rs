//! Boundary protocol: how the cell is driven over time.
//!
//! Galvanostatic: prescribe total current; let the BV machinery distribute it
//! across the surface. Terminal voltage is measured.
//! Potentiostatic: prescribe terminal voltage; current emerges.
//! StepPulse: piecewise-constant galvanostatic schedule for DCR-style tests.

#[derive(Clone, Debug)]
pub enum Protocol {
    /// Hold a fixed current forever.
    Galvanostatic { current: f32 },
    /// Hold a fixed terminal voltage forever.
    Potentiostatic { voltage: f32 },
    /// (relax_steps, pulse_steps, pulse_current) loop.
    StepPulse {
        relax_steps: usize,
        pulse_steps: usize,
        pulse_current: f32,
    },
}

#[derive(Clone, Copy, Debug)]
pub enum Drive {
    Current(f32),
    Voltage(f32),
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

    /// Returns the drive to apply this step and advances internal state.
    pub fn tick(&mut self) -> Drive {
        let drive = match &self.protocol {
            Protocol::Galvanostatic { current } => Drive::Current(*current),
            Protocol::Potentiostatic { voltage } => Drive::Voltage(*voltage),
            Protocol::StepPulse {
                relax_steps,
                pulse_steps,
                pulse_current,
            } => {
                let cycle = relax_steps + pulse_steps;
                let phase = self.step_index % cycle;
                if phase < *relax_steps {
                    Drive::Current(0.0)
                } else {
                    Drive::Current(*pulse_current)
                }
            }
        };
        self.step_index += 1;
        drive
    }
}
