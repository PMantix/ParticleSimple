//! Headline scenario: voltage step pulse train with relaxation, run long
//! enough that morphology-driven and SEI-driven trends are visible over
//! cycles.
//!
//! Apply 0 -> +V -> 0 pulses. With BV deposition, V > 0 strips at left
//! and plates at right; the V_bulk(t) signature shows R0 / R1 / R2
//! separation in each pulse. Solvent reduces irreversibly to SEI on the
//! reducing-side metal surface, blocking ion transport in those grid
//! cells.

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

pub const RELAX_STEPS: usize = 500;
pub const PULSE_STEPS: usize = 1000;
pub const PULSE_VOLTAGE: f32 = 0.1;

/// Construct a freshly-initialized Cell wired up for the pulse_dcr
/// experiment. Used by both the headless `run()` driver and the GUI bin.
pub fn setup() -> Cell {
    let domain = Domain::new(10.0, 5.0);
    let grid = Grid::new(&domain, 32, 16);

    let n_pairs = 1000;
    let n_solvent = 200;
    let n_metal_rows = 4;
    let n_metal_per_row = 16;
    let mut particles =
        Vec::with_capacity(2 * n_pairs + n_solvent + 2 * n_metal_rows * n_metal_per_row);
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

    let protocol = ProtocolState::new(Protocol::StepVoltage {
        relax_steps: RELAX_STEPS,
        pulse_steps: PULSE_STEPS,
        pulse_voltage: PULSE_VOLTAGE,
    });

    let bv_deposition = BvParams {
        i0: 0.3,
        alpha: 0.5,
        kt: 0.025,
        eq_potential: 0.0,
        reaction_dx_factor: 1.5,
    };
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

    Cell::new(domain, grid, particles, protocol, params)
}

pub fn run() {
    let mut cell = setup();
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
