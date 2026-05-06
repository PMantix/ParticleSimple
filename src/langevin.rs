//! Overdamped Langevin integrator.
//!
//! Update rule for each mobile particle:
//!   x_{n+1} = x_n + (q * D_eff / kT) * E(x_n) * dt + sqrt(2 * D_eff * dt) * xi
//! where xi is a 2D unit-variance Gaussian and D_eff = D * mobility(cell).
//! Reading mobility per cell lets SEI-blocked regions slow ions passing
//! through them.

use glam::Vec2;
use rand_distr::{Distribution, Normal};

use crate::domain::Domain;
use crate::grid::Grid;
use crate::particle::Particle;

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

        let (ix, iy) = grid.cell_of(p.pos);
        let mobility = grid.mobility[grid.idx(ix, iy)].max(0.0);
        let d_eff = props.d * mobility;

        let e = grid.field_at(p.pos);
        let drift_coeff = props.charge * d_eff * dt / kt;
        let drift = e * drift_coeff;

        let sigma = (two_dt * d_eff).sqrt();
        let xi = Vec2::new(
            normal.sample(&mut rng) as f32,
            normal.sample(&mut rng) as f32,
        );
        let diffuse = xi * sigma;

        let mut new_pos = p.pos + drift + diffuse;

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
