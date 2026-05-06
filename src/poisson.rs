//! Mean-field Poisson solver on the grid.
//!
//! Day 1: Jacobi / SOR with electrode Dirichlet BCs (phi prescribed at the
//! left and right electrodes) and reflecting BCs in y. This is correct but
//! slow for large grids; FFT or multigrid is the upgrade path.

use crate::grid::Grid;

#[derive(Clone, Copy, Debug)]
pub struct BoundaryPotentials {
    pub left: f32,
    pub right: f32,
}

/// Solve -Laplacian(phi) = rho / epsilon0_eff with Dirichlet BCs on
/// the left/right electrodes. Returns once the relative residual drops below
/// `tol` or `max_iters` is reached.
pub fn solve(
    _grid: &mut Grid,
    _bcs: BoundaryPotentials,
    _epsilon0_eff: f32,
    _tol: f32,
    _max_iters: usize,
) -> usize {
    todo!("Jacobi/SOR Poisson solve with Dirichlet x-boundaries, reflecting y")
}
