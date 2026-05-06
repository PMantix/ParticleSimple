//! Stochastic surface reactions: deposition / stripping (Butler-Volmer) and
//! solvent reduction to SEI.
//!
//! Day 2 stub: both functions are no-ops. The Cell::step pipeline calls them
//! every step so the plumbing exists, but the BV machinery is filled in next
//! round once empty_cell shows clean ohmic behavior.

use crate::cell::Cell;

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
}

/// Deposition / stripping at the metal anode: Cation <-> Metal.
pub fn deposition(_cell: &mut Cell, _bv: BvParams, _dt: f32) {
    // TODO(reactions): BV-driven plate/strip at the anode surface.
}

/// Solvent reduction at the metal anode: Solvent -> Sei (irreversible).
pub fn sei_formation(_cell: &mut Cell, _bv: BvParams, _dt: f32) {
    // TODO(reactions): BV-driven irreversible solvent->SEI conversion.
}
