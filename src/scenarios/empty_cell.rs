//! Bare electrolyte between two electrodes, small applied voltage.
//!
//! Tuned so the EDL formation transient is visible inside ~5k steps:
//!   - Domain L = 20, dx = 0.625 (32 cells across)
//!   - D = 1, dt = 0.05  =>  sigma_diff = 0.32 < dx/2 (stable)
//!   - n_pairs = 500 (n = 5/area, lambda_D ~ 0.5 with epsilon = 100)
//!   - Conductivity sigma ~ 2*n*q^2*D/kT = 400  =>  R0 ~ L/(sigma*A) ~ 0.005
//!   - C_dl ~ epsilon/lambda_D ~ 200/area  =>  C_total ~ 1000
//!   - tau_RC ~ R0 * C_total ~ 5 sim time = 100 steps
//!
//! With no surface reactions, V is held constant and the current decays
//! after the EDL has charged.

use glam::Vec2;
use rand::Rng;

use crate::cell::{Cell, CellParams};
use crate::domain::Domain;
use crate::grid::Grid;
use crate::measure::CsvSink;
use crate::particle::Particle;
use crate::protocol::{Protocol, ProtocolState};
use crate::reactions::BvParams;
use crate::species::Species;

pub fn run() {
    let domain = Domain::new(10.0, 5.0);
    let grid = Grid::new(&domain, 32, 16);

    let n_pairs = 500;
    let mut particles = Vec::with_capacity(2 * n_pairs);
    let mut rng = rand::rng();
    for _ in 0..n_pairs {
        let cation_pos = Vec2::new(
            rng.random_range(-domain.half_width..domain.half_width),
            rng.random_range(-domain.half_height..domain.half_height),
        );
        particles.push(Particle::new(cation_pos, Species::Cation));
        let anion_pos = Vec2::new(
            rng.random_range(-domain.half_width..domain.half_width),
            rng.random_range(-domain.half_height..domain.half_height),
        );
        particles.push(Particle::new(anion_pos, Species::Anion));
    }

    let protocol = ProtocolState::new(Protocol::Potentiostatic { voltage: 0.1 });

    let bv_zero = BvParams::default();
    let params = CellParams {
        dt: 0.05,
        kt: 0.025,
        epsilon0_eff: 100.0,
        poisson_tol: 1e-4,
        poisson_max_iters: 200,
        bv_deposition: bv_zero,
        bv_sei: bv_zero,
    };

    let mut cell = Cell::new(domain, grid, particles, protocol, params);

    let path = "empty_cell.csv";
    let mut sink = CsvSink::new(path).expect("create CSV");
    let n_steps = 5000;
    for _ in 0..n_steps {
        let m = cell.step();
        sink.write(&m).expect("write CSV");
    }
    println!(
        "empty_cell: ran {} steps, {} particles, wrote {}",
        n_steps,
        cell.particles.len(),
        path
    );
}
