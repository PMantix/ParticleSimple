//! Stochastic surface reactions: deposition / stripping (Butler-Volmer) and
//! solvent reduction to SEI.
//!
//! All reactions fire as per-particle Bernoulli events with rate r:
//!   p_event = 1 - exp(-r * dt)
//! and r is set by Butler-Volmer:
//!   r = i0 * (exp(alpha * eta / kT) - exp(-(1 - alpha) * eta / kT))
//! where eta is the local overpotential (potential at the particle's grid
//! cell minus the equilibrium reduction potential of the species).

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
/// Cations within reaction range of an exposed metal site can plate; metal
/// atoms at the surface can strip back to cations under reverse polarization.
pub fn deposition(_cell: &mut Cell, _bv: BvParams, _dt: f32) {
    todo!("BV-driven plate/strip at the anode surface")
}

/// Solvent reduction at the metal anode: Solvent -> Sei.
/// Same BV form, much smaller i0, irreversible (no backward branch).
pub fn sei_formation(_cell: &mut Cell, _bv: BvParams, _dt: f32) {
    todo!("BV-driven irreversible solvent->SEI conversion at the anode")
}
