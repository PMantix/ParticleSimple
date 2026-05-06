//! Overdamped Langevin integrator.
//!
//! Update rule for each mobile particle:
//!   x_{n+1} = x_n + (q * D / kT) * E(x_n) * dt + sqrt(2 * D * dt) * xi
//! where xi is a 2D unit-variance Gaussian. This is the right dynamical
//! regime for ions in a viscous electrolyte and lets dt scale as dx^2 / D
//! rather than at the inertial-MD timescale.

use glam::Vec2;
use rand_distr::{Distribution, Normal};

use crate::domain::Domain;
use crate::grid::Grid;
use crate::particle::Particle;

/// Advance all mobile particles by one Langevin step. Particles are reflected
/// at all four domain walls; surface reactions (deposition / SEI) are handled
/// in a separate pass after this one.
pub fn step(particles: &mut [Particle], grid: &Grid, domain: &Domain, kt: f32, dt: f32) {
    let normal = Normal::<f64>::new(0.0, 1.0).expect("valid normal");
    let mut rng = rand::rng();
    let two_dt = 2.0 * dt;
    let hw = domain.half_width;
    let hh = domain.half_height;

    for p in particles.iter_mut() {
        let props = p.species.props();
        if !props.mobile {
            continue;
        }

        let e = grid.field_at(p.pos);
        let drift_coeff = props.charge * props.d * dt / kt;
        let drift = e * drift_coeff;

        let sigma = (two_dt * props.d).sqrt();
        let xi = Vec2::new(
            normal.sample(&mut rng) as f32,
            normal.sample(&mut rng) as f32,
        );
        let diffuse = xi * sigma;

        let mut new_pos = p.pos + drift + diffuse;

        // Reflect across walls (one bounce; clamp afterwards in case of overshoot).
        if new_pos.x < -hw {
            new_pos.x = -2.0 * hw - new_pos.x;
        } else if new_pos.x > hw {
            new_pos.x = 2.0 * hw - new_pos.x;
        }
        if new_pos.y < -hh {
            new_pos.y = -2.0 * hh - new_pos.y;
        } else if new_pos.y > hh {
            new_pos.y = 2.0 * hh - new_pos.y;
        }
        new_pos.x = new_pos.x.clamp(-hw, hw);
        new_pos.y = new_pos.y.clamp(-hh, hh);

        p.pos = new_pos;
    }
}
