//! Stochastic surface reactions: deposition / stripping (Butler-Volmer)
//! and solvent reduction to SEI.
//!
//! Per-particle Bernoulli events: p_event = 1 - exp(-r * dt), where r is
//! the BV rate at the particle's grid cell.
//!
//! Plating branch (cathodic):  r_p = i0 * exp(-(1 - alpha) * eta / kT)
//! Stripping branch (anodic):  r_s = i0 * exp( alpha       * eta / kT)
//! eta = phi_local - eq_potential.
//!
//! Stripping is restricted to *surface* metals: a metal is strippable only
//! if its grid cell has no metal occupancy in the bulk-side neighbor cell.
//! Without this rule, deposits "evaporate" from the interior, defeating
//! the morphology-evolution goal.

use rand::Rng;

use crate::cell::Cell;
use crate::species::Species;

#[derive(Clone, Copy, Debug)]
pub struct BvParams {
    pub i0: f32,
    pub alpha: f32,
    pub kt: f32,
    pub eq_potential: f32,
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

/// Per-cell metal occupancy count, used by the surface-strip rule and by
/// SEI's "near a metal" gate. Built once per step.
fn build_metal_count(cell: &Cell) -> Vec<u16> {
    let n = cell.grid.nx * cell.grid.ny;
    let mut counts = vec![0u16; n];
    for p in &cell.particles {
        if p.species == Species::Metal {
            let (ix, iy) = cell.grid.cell_of(p.pos);
            counts[cell.grid.idx(ix, iy)] += 1;
        }
    }
    counts
}

/// True iff the cell adjacent to (ix, iy) in the bulk direction has no
/// metal -- i.e. this metal is at the deposit surface and exposed to
/// electrolyte.
fn is_surface_metal(
    cell: &Cell,
    metal_count: &[u16],
    ix: usize,
    iy: usize,
    x: f32,
) -> bool {
    let bulk_ix = if x < 0.0 {
        ix.saturating_add(1)
    } else if ix > 0 {
        ix - 1
    } else {
        return true; // Edge case at right boundary's leftward direction.
    };
    if bulk_ix >= cell.grid.nx {
        return true;
    }
    metal_count[cell.grid.idx(bulk_ix, iy)] == 0
}

pub fn deposition(cell: &mut Cell, bv: BvParams, dt: f32) -> ReactionCounts {
    if bv.i0 <= 0.0 {
        return ReactionCounts::default();
    }
    let mut counts = ReactionCounts::default();
    let reaction_distance = cell.grid.dx * bv.reaction_dx_factor;
    let half_w = cell.domain.half_width;
    let mut rng = rand::rng();

    let metal_count = build_metal_count(cell);

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
                let (ix, iy) = cell.grid.cell_of(p.pos);
                if !is_surface_metal(cell, &metal_count, ix, iy, p.pos.x) {
                    continue;
                }
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

/// Solvent -> Sei (irreversible, cathodic branch only). Forms only where
/// eta is sufficiently negative -- in our toy this is the reducing
/// electrode side at positive applied V. Reduces grid mobility in the
/// cell where SEI is formed, which slows ion transport through that
/// region (the R0/R2-impedance mechanism).
pub fn sei_formation(cell: &mut Cell, bv: BvParams, dt: f32) -> u32 {
    if bv.i0 <= 0.0 {
        return 0;
    }
    let metal_count = build_metal_count(cell);
    let reaction_distance = cell.grid.dx * bv.reaction_dx_factor;
    let r_cells = (reaction_distance / cell.grid.dx).ceil() as i32;
    let mut rng = rand::rng();

    let mut to_sei: Vec<usize> = Vec::new();

    for (i, p) in cell.particles.iter().enumerate() {
        if p.species != Species::Solvent {
            continue;
        }
        // Must be near a metal surface (within r_cells in any direction).
        let (ix, iy) = cell.grid.cell_of(p.pos);
        let mut near_metal = false;
        'search: for di in -r_cells..=r_cells {
            for dj in -r_cells..=r_cells {
                let nx = ix as i32 + di;
                let ny = iy as i32 + dj;
                if nx < 0 || ny < 0 || nx >= cell.grid.nx as i32 || ny >= cell.grid.ny as i32 {
                    continue;
                }
                if metal_count[cell.grid.idx(nx as usize, ny as usize)] > 0 {
                    near_metal = true;
                    break 'search;
                }
            }
        }
        if !near_metal {
            continue;
        }

        let phi = cell.grid.sample_phi(p.pos);
        let eta = phi - bv.eq_potential;
        // Only cathodic branch; only when eta is favorable (negative).
        if eta > 0.0 {
            continue;
        }
        let r = bv.i0 * f32::exp(-(1.0 - bv.alpha) * eta / bv.kt);
        let prob = 1.0 - f32::exp(-r * dt);
        if rng.random::<f32>() < prob {
            to_sei.push(i);
        }
    }

    let count = to_sei.len() as u32;
    for i in to_sei {
        cell.particles[i].species = Species::Sei;
        let (ix, iy) = cell.grid.cell_of(cell.particles[i].pos);
        let idx = cell.grid.idx(ix, iy);
        // Each SEI particle reduces local mobility multiplicatively, with
        // a floor so transport never fully stops.
        cell.grid.mobility[idx] = (cell.grid.mobility[idx] * 0.5).max(0.05);
    }

    count
}
