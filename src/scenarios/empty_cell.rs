//! Bare electrolyte between two electrodes, small applied voltage.
//! Sanity check: the cell should produce a charge-redistribution transient
//! that decays into a steady state once the double layers are formed. With
//! no surface reactions, the long-time current is small (zero in the limit).

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
    let domain = Domain::new(50.0, 25.0);
    let grid = Grid::new(&domain, 64, 32);

    let n_pairs = 2000;
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

    let bv_zero = BvParams {
        i0: 0.0,
        alpha: 0.5,
        kt: 0.025,
        eq_potential: 0.0,
    };
    let params = CellParams {
        dt: 0.05,
        kt: 0.025,
        epsilon0_eff: 20.0,
        poisson_tol: 1e-4,
        poisson_max_iters: 200,
        bv_deposition: bv_zero,
        bv_sei: bv_zero,
    };

    let mut cell = Cell::new(domain, grid, particles, protocol, params);

    let path = "empty_cell.csv";
    let mut sink = CsvSink::new(path).expect("create CSV");
    let n_steps = 1000;
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
