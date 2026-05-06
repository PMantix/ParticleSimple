//! Top-level cell state. Owns the particles, the grid, the protocol state,
//! and orchestrates one step of the pipeline.

use crate::domain::Domain;
use crate::grid::Grid;
use crate::langevin;
use crate::measure::Measurement;
use crate::particle::Particle;
use crate::poisson::{self, BoundaryPotentials};
use crate::protocol::ProtocolState;
use crate::reactions::{self, BvParams, ReactionCounts};

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
    /// Total left-half charge from the previous step; used to compute the
    /// midplane current as the rate of net charge flow across x = 0.
    pub charge_left_prev: f32,
    /// Most recent measured boundary current (= the cell's external current
    /// at the left electrode). Fed back to galvanostatic controllers.
    pub last_current: f32,
    /// Most recent boundary potentials applied; used by Measurement to
    /// report the applied terminal voltage.
    pub last_bcs: BoundaryPotentials,
    /// Reaction events counted in the most recent step.
    pub last_reaction_counts: ReactionCounts,
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
            last_current: 0.0,
            last_bcs: BoundaryPotentials {
                left: 0.0,
                right: 0.0,
            },
            last_reaction_counts: ReactionCounts::default(),
        }
    }

    pub fn step(&mut self) -> Measurement {
        let dt = self.params.dt;
        let kt = self.params.kt;
        let bv_dep = self.params.bv_deposition;
        let bv_sei = self.params.bv_sei;

        let bcs = self.protocol.tick(dt, self.last_current);
        self.last_bcs = bcs;

        self.grid.clear_rho();
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

        poisson::solve(
            &mut self.grid,
            bcs,
            self.params.epsilon0_eff,
            self.params.poisson_tol,
            self.params.poisson_max_iters,
        );

        langevin::step(&mut self.particles, &self.grid, &self.domain, kt, dt);

        let mut counts = reactions::deposition(self, bv_dep, dt);
        counts.sei_formed = reactions::sei_formation(self, bv_sei, dt);
        self.last_reaction_counts = counts;

        let m = Measurement::sample(self);
        // Roll forward the previous-charge baseline so next step's midplane
        // current is per-step, not cumulative.
        self.charge_left_prev = self.charge_in_left_half();
        // Boundary current feeds back to the PI controller.
        self.last_current = m.current;
        self.step_index += 1;
        m
    }

    /// Net charge in the left half of the domain (x < 0).
    pub fn charge_in_left_half(&self) -> f32 {
        self.particles
            .iter()
            .filter(|p| p.pos.x < 0.0)
            .map(|p| p.charge())
            .sum()
    }
}
