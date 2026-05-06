//! Top-level cell state. Owns the particles, the grid, the protocol state,
//! and orchestrates one step of the pipeline.

use crate::domain::Domain;
use crate::grid::Grid;
use crate::langevin;
use crate::measure::Measurement;
use crate::particle::Particle;
use crate::poisson::{self, BoundaryPotentials};
use crate::protocol::{Drive, ProtocolState};
use crate::reactions::{self, BvParams};

#[derive(Clone, Debug)]
pub struct CellParams {
    pub dt: f32,
    pub kt: f32,
    pub epsilon0_eff: f32,
    pub poisson_tol: f32,
    pub poisson_max_iters: usize,
    pub bv_deposition: BvParams,
    pub bv_sei: BvParams,
}

pub struct Cell {
    pub domain: Domain,
    pub grid: Grid,
    pub particles: Vec<Particle>,
    pub protocol: ProtocolState,
    pub params: CellParams,
    pub step_index: usize,
    /// Total left-half charge from previous step; used to compute current as
    /// the rate of net charge flow across the midplane (x = 0).
    pub charge_left_prev: f32,
}

impl Cell {
    pub fn new(
        domain: Domain,
        grid: Grid,
        particles: Vec<Particle>,
        protocol: ProtocolState,
        params: CellParams,
    ) -> Self {
        let charge_left_prev = particles
            .iter()
            .filter(|p| p.pos.x < 0.0)
            .map(|p| p.charge())
            .sum();
        Self {
            domain,
            grid,
            particles,
            protocol,
            params,
            step_index: 0,
            charge_left_prev,
        }
    }

    /// One simulation step:
    ///   1. clear / deposit rho
    ///   2. set BCs from current protocol drive
    ///   3. solve Poisson
    ///   4. Langevin step on mobile particles
    ///   5. surface reactions (deposition + SEI)
    ///   6. record measurement
    pub fn step(&mut self) -> Measurement {
        let drive = self.protocol.tick();
        let bv_dep = self.params.bv_deposition;
        let bv_sei = self.params.bv_sei;
        let dt = self.params.dt;
        let kt = self.params.kt;

        self.grid.clear_rho();
        // Collect (pos, charge) pairs into a temporary; deposit_charges
        // borrows grid mutably.
        {
            let charges: Vec<_> = self
                .particles
                .iter()
                .filter_map(|p| {
                    let q = p.charge();
                    if q == 0.0 {
                        None
                    } else {
                        Some((p.pos, q))
                    }
                })
                .collect();
            self.grid.deposit_charges(charges);
        }

        let bcs = boundary_potentials_from_drive(drive);
        poisson::solve(
            &mut self.grid,
            bcs,
            self.params.epsilon0_eff,
            self.params.poisson_tol,
            self.params.poisson_max_iters,
        );

        langevin::step(&mut self.particles, &self.grid, &self.domain, kt, dt);

        reactions::deposition(self, bv_dep, dt);
        reactions::sei_formation(self, bv_sei, dt);

        let m = Measurement::sample(self, drive);
        self.step_index += 1;
        m
    }

    /// Net charge in the left half of the domain (x < 0). Used with
    /// `charge_left_prev` to derive a midplane current.
    pub fn charge_in_left_half(&self) -> f32 {
        self.particles
            .iter()
            .filter(|p| p.pos.x < 0.0)
            .map(|p| p.charge())
            .sum()
    }
}

/// Map a Drive (potentiostatic V, or galvanostatic I) into electrode BCs.
/// Day 2 supports potentiostatic only; Drive::Current produces open-circuit
/// (zero) BCs and the caller is expected not to use it yet.
fn boundary_potentials_from_drive(drive: Drive) -> BoundaryPotentials {
    match drive {
        Drive::Voltage(v) => BoundaryPotentials {
            left: 0.5 * v,
            right: -0.5 * v,
        },
        Drive::Current(_) => BoundaryPotentials {
            left: 0.0,
            right: 0.0,
        },
    }
}
