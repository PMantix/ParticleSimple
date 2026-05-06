//! Overdamped Langevin integrator.
//!
//! Update rule for each mobile particle:
//!   x_{n+1} = x_n + (q / gamma) * E(x_n) * dt + sqrt(2 * D * dt) * xi
//! where xi is a 2D unit-variance Gaussian. This is the right dynamical
//! regime for ions in a viscous electrolyte and lets dt scale as dx^2 / D
//! rather than at the inertial-MD timescale.

use crate::grid::Grid;
use crate::particle::Particle;

/// Advance all mobile particles by one Langevin step.
///
/// `mobility_factor(species)` lets per-species drag differ from the bulk
/// (e.g. ions inside SEI cells get a multiplicative slowdown).
pub fn step(_particles: &mut [Particle], _grid: &Grid, _dt: f32) {
    todo!("overdamped Langevin step with CIC field sampling and Gaussian noise")
}
