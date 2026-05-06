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
}

impl Cell {
    pub fn new(
        domain: Domain,
        grid: Grid,
        particles: Vec<Particle>,
        protocol: ProtocolState,
        params: CellParams,
    ) -> Self {
        Self {
            domain,
            grid,
            particles,
            protocol,
            params,
            step_index: 0,
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

        self.grid.clear_rho();
        self.grid.deposit_charges(
            self.particles
                .iter()
                .filter_map(|p| {
                    let q = p.charge();
                    if q == 0.0 {
                        None
                    } else {
                        Some((p.pos, q))
                    }
                }),
        );

        let bcs = boundary_potentials_from_drive(drive);
        poisson::solve(
            &mut self.grid,
            bcs,
            self.params.epsilon0_eff,
            self.params.poisson_tol,
            self.params.poisson_max_iters,
        );

        langevin::step(&mut self.particles, &self.grid, self.params.dt);

        reactions::deposition(self, self.params.bv_deposition, self.params.dt);
        reactions::sei_formation(self, self.params.bv_sei, self.params.dt);

        let m = Measurement::sample(self, drive);
        self.step_index += 1;
        m
    }
}

fn boundary_potentials_from_drive(_drive: Drive) -> BoundaryPotentials {
    todo!("translate Drive::Current into BCs that produce the requested net flux, or pass Drive::Voltage straight through")
}
