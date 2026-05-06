//! Stochastic surface reactions: deposition / stripping (Butler-Volmer)
//! and (next round) solvent reduction to SEI.
//!
//! Per-particle Bernoulli events with Bernoulli probability
//!   p_event = 1 - exp(-r * dt)
//! where r is the BV rate at the particle's grid cell.
//!
//! Plating branch (cathodic):    r_p = i0 * exp(-(1 - alpha) * eta / kT)
//! Stripping branch (anodic):    r_s = i0 * exp( alpha       * eta / kT)
//! eta = phi_local - eq_potential.
//!
//! With phi_left > 0 (left electrode held positive): eta > 0 there, so
//! stripping dominates at the left and plating at the right. The boundary
//! current is positive when net stripping happens at the left electrode
//! (Li -> Li+ + e-, electrons flow from electrode out through the wire,
//! conventional current flows from wire into the electrode).

use rand::Rng;

use crate::cell::Cell;
use crate::species::Species;

#[derive(Clone, Copy, Debug)]
pub struct BvParams {
    /// Exchange current (per-particle base rate).
    pub i0: f32,
    /// Symmetry factor (0..1). Use 0.5 unless asymmetry is needed.
    pub alpha: f32,
    /// Thermal scale kT in sim units (overpotential normalization).
    pub kt: f32,
    /// Equilibrium reduction potential of this species at the electrode.
    pub eq_potential: f32,
    /// How far from the electrode boundary a Cation can react (in dx units).
    pub reaction_dx_factor: f32,
}

impl Default for BvParams {
    fn default() -> Self {
        Self {
            i0: 0.0,
            alpha: 0.5,
            kt: 0.025,
            eq_potential: 0.0,
            reaction_dx_factor: 1.5,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ReactionCounts {
    pub plate_left: u32,
    pub strip_left: u32,
    pub plate_right: u32,
    pub strip_right: u32,
    pub sei_formed: u32,
}

/// Run BV deposition / stripping at both electrodes. Returns event counts.
pub fn deposition(cell: &mut Cell, bv: BvParams, dt: f32) -> ReactionCounts {
    if bv.i0 <= 0.0 {
        return ReactionCounts::default();
    }
    let mut counts = ReactionCounts::default();
    let reaction_distance = cell.grid.dx * bv.reaction_dx_factor;
    let half_w = cell.domain.half_width;
    let mut rng = rand::rng();

    let mut to_metal: Vec<usize> = Vec::new();
    let mut to_cation_left: Vec<usize> = Vec::new();
    let mut to_cation_right: Vec<usize> = Vec::new();

    for (i, p) in cell.particles.iter().enumerate() {
        let phi = cell.grid.sample_phi(p.pos);
        let eta = phi - bv.eq_potential;

        match p.species {
            Species::Cation => {
                let near_left = (p.pos.x + half_w) < reaction_distance;
                let near_right = (half_w - p.pos.x) < reaction_distance;
                if !near_left && !near_right {
                    continue;
                }
                let r = bv.i0 * f32::exp(-(1.0 - bv.alpha) * eta / bv.kt);
                let prob = 1.0 - f32::exp(-r * dt);
                if rng.random::<f32>() < prob {
                    to_metal.push(i);
                    if near_left {
                        counts.plate_left += 1;
                    } else {
                        counts.plate_right += 1;
                    }
                }
            }
            Species::Metal => {
                let r = bv.i0 * f32::exp(bv.alpha * eta / bv.kt);
                let prob = 1.0 - f32::exp(-r * dt);
                if rng.random::<f32>() < prob {
                    if p.pos.x < 0.0 {
                        to_cation_left.push(i);
                        counts.strip_left += 1;
                    } else {
                        to_cation_right.push(i);
                        counts.strip_right += 1;
                    }
                }
            }
            _ => {}
        }
    }

    // Apply species changes. Plated metals stay where the cation was (lets
    // morphology emerge from where ions actually arrive). Stripped cations
    // are nudged a fraction of a cell into the bulk to avoid immediate
    // re-plating in the same step.
    for i in to_metal {
        cell.particles[i].species = Species::Metal;
    }
    let nudge = cell.grid.dx * 0.5;
    for i in to_cation_left {
        cell.particles[i].species = Species::Cation;
        cell.particles[i].pos.x = (-half_w + nudge).max(cell.particles[i].pos.x);
    }
    for i in to_cation_right {
        cell.particles[i].species = Species::Cation;
        cell.particles[i].pos.x = (half_w - nudge).min(cell.particles[i].pos.x);
    }

    counts
}

/// Solvent -> SEI conversion. Stub for Day 4a; lands next round.
pub fn sei_formation(_cell: &mut Cell, _bv: BvParams, _dt: f32) -> u32 {
    0
}
