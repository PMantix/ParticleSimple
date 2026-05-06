//! Headline scenario: voltage step pulse with relaxation.
//!
//! Apply a 0 -> +V -> 0 pulse train and record the current response.
//! With BV deposition active, V > 0 drives stripping at the left and
//! plating at the right; the transient I(t) signature exposes:
//!
//!   R0  - immediate ohmic conductance (large initial current spike)
//!   R1  - charge-transfer relaxation as overpotential at the surface
//!         catches up with the imposed V (current decays over tau1)
//!   R2  - longer-time concentration polarization (current keeps drifting
//!         down as the reaction zones deplete)
//!
//! The voltage step is the dual of the galvanostatic current step the
//! real DCR data uses. Same R/tau information, more numerically stable
//! in this toy because we don't fight a controller into saturation.

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
    // Seed several rows of Metal at each electrode so stripping has a
    // sustained reservoir. Single-row seeding depletes within ~hundreds
    // of pulse steps and stalls the response.
    let n_metal_rows = 4;
    let n_metal_per_row = 16;
    let mut particles =
        Vec::with_capacity(2 * n_pairs + 2 * n_metal_rows * n_metal_per_row);
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

    // 4 cycles of relax(500) + pulse(1000) = 6000 steps.
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
    let bv_sei = BvParams::default();

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
    let n_steps = 6000;
    for _ in 0..n_steps {
        let m = cell.step();
        sink.write(&m).expect("write CSV");
    }
    println!(
        "pulse_dcr: ran {} steps, {} particles total, {} metal at end, wrote {}",
        n_steps,
        cell.particles.len(),
        cell.particles
            .iter()
            .filter(|p| p.species == Species::Metal)
            .count(),
        path
    );
}
