//! Headline scenario: voltage step pulse train with relaxation, run long
//! enough that morphology-driven and SEI-driven trends are visible over
//! cycles.
//!
//! Apply 0 -> +V -> 0 pulses. With BV deposition, V > 0 strips at left
//! and plates at right; the V_bulk(t) signature shows R0 / R1 / R2
//! separation in each pulse. Solvent reduces irreversibly to SEI on the
//! reducing-side metal surface, blocking ion transport in those grid
//! cells.
//!
//! Over many cycles we expect:
//!   - per-cycle metal_count drifts upward at the right (deposit growth)
//!   - sei_fraction grows monotonically as solvent gets consumed
//!   - voltage_bulk during pulse drops further over cycles (R0/R2 rising
//!     from SEI; R1 falling from roughened surface area)

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

    let n_pairs = 1000;
    let n_solvent = 200;
    let n_metal_rows = 4;
    let n_metal_per_row = 16;
    let mut particles = Vec::with_capacity(
        2 * n_pairs + n_solvent + 2 * n_metal_rows * n_metal_per_row,
    );
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
    for _ in 0..n_solvent {
        let pos = Vec2::new(
            rng.random_range(-domain.half_width..domain.half_width),
            rng.random_range(-domain.half_height..domain.half_height),
        );
        particles.push(Particle::new(pos, Species::Solvent));
    }

    let dy = (2.0 * domain.half_height) / n_metal_per_row as f32;
    for row in 0..n_metal_rows {
        let depth = (row as f32 + 0.5) * grid.dx;
        let x_left = -domain.half_width + depth;
        let x_right = domain.half_width - depth;
        for i in 0..n_metal_per_row {
            let y = -domain.half_height + (i as f32 + 0.5) * dy;
            particles.push(Particle::new(Vec2::new(x_left, y), Species::Metal));
            particles.push(Particle::new(Vec2::new(x_right, y), Species::Metal));
        }
    }

    // 10 cycles of relax(500) + pulse(1000) = 15000 steps total.
    let protocol = ProtocolState::new(Protocol::StepVoltage {
        relax_steps: 500,
        pulse_steps: 1000,
        pulse_voltage: 0.1,
    });

    let bv_deposition = BvParams {
        i0: 0.3,
        alpha: 0.5,
        kt: 0.025,
        eq_potential: 0.0,
        reaction_dx_factor: 1.5,
    };
    // SEI: lower i0 (irreversible passivation is slow), eq_potential set
    // such that only the reducing-side electrolyte forms it.
    let bv_sei = BvParams {
        i0: 0.05,
        alpha: 0.5,
        kt: 0.025,
        eq_potential: 0.0,
        reaction_dx_factor: 2.0,
    };

    let params = CellParams {
        dt: 0.05,
        kt: 0.025,
        epsilon0_eff: 100.0,
        poisson_tol: 1e-4,
        poisson_max_iters: 200,
        bv_deposition,
        bv_sei,
    };

    let mut cell = Cell::new(domain, grid, particles, protocol, params);

    let path = "pulse_dcr.csv";
    let mut sink = CsvSink::new(path).expect("create CSV");
    let n_steps = 15000;
    for _ in 0..n_steps {
        let m = cell.step();
        sink.write(&m).expect("write CSV");
    }

    let metal_total = cell
        .particles
        .iter()
        .filter(|p| p.species == Species::Metal)
        .count();
    let sei_total = cell
        .particles
        .iter()
        .filter(|p| p.species == Species::Sei)
        .count();
    let cation_total = cell
        .particles
        .iter()
        .filter(|p| p.species == Species::Cation)
        .count();
    println!(
        "pulse_dcr: ran {} steps. End: {} cations, {} metal, {} sei. wrote {}",
        n_steps, cation_total, metal_total, sei_total, path
    );
}
