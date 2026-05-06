//! PI controller sanity check.
//!
//! Hold a constant target current. The controller adjusts applied voltage
//! each step to drive measured current toward the target. After settling,
//! `current` ~ target and `voltage_applied` ~ target * R_ohmic.

use glam::Vec2;
use rand::Rng;

use crate::cell::{Cell, CellParams};
use crate::domain::Domain;
use crate::grid::Grid;
use crate::measure::CsvSink;
use crate::particle::Particle;
use crate::protocol::{CurrentController, Protocol, ProtocolState};
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

    // NULL SANITY TEST. With no surface reactions, the cell has no DC
    // path -- any applied V just charges the EDL until current stops, so
    // a non-zero target is physically impossible and the controller will
    // saturate trying. Set target=0 instead; we verify only that the PI
    // is wired correctly (voltage stays near 0, current fluctuates around
    // 0 from thermal noise). The real galvanostatic test happens once
    // deposition reactions exist next round.
    let controller = CurrentController::new(
        /* target */ 0.0,
        /* kp */ 0.0005,
        /* ki */ 0.05,
        /* filter_alpha */ 0.05,
        /* voltage_max */ 0.5,
    );
    let protocol = ProtocolState::new(Protocol::Galvanostatic { controller });

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

    let path = "galvanostatic_check.csv";
    let mut sink = CsvSink::new(path).expect("create CSV");
    let n_steps = 5000;
    for _ in 0..n_steps {
        let m = cell.step();
        sink.write(&m).expect("write CSV");
    }
    println!(
        "galvanostatic_check: ran {} steps, target current 0.0, wrote {}",
        n_steps, path
    );
}
