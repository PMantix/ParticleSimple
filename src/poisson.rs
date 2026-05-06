//! Mean-field Poisson solver on the grid.
//!
//! Day 1: Jacobi with Dirichlet BCs at the left/right electrodes and
//! reflecting (zero-flux) BCs in y. Warm-started from the previous step's
//! phi, which makes typical iteration counts well below `max_iters`. FFT or
//! multigrid is the upgrade path when the grid grows.

use crate::grid::Grid;

#[derive(Clone, Copy, Debug)]
pub struct BoundaryPotentials {
    pub left: f32,
    pub right: f32,
}

/// Solve `-Laplacian(phi) = rho / epsilon0_eff` with Dirichlet BCs at
/// ix = 0 and ix = nx-1, reflecting in y. Returns the iteration count used.
pub fn solve(
    grid: &mut Grid,
    bcs: BoundaryPotentials,
    epsilon0_eff: f32,
    tol: f32,
    max_iters: usize,
) -> usize {
    let nx = grid.nx;
    let ny = grid.ny;

    // Apply Dirichlet BCs in-place (don't iterate them).
    for iy in 0..ny {
        let i_left = grid.idx(0, iy);
        let i_right = grid.idx(nx - 1, iy);
        grid.phi[i_left] = bcs.left;
        grid.phi[i_right] = bcs.right;
    }

    let inv_dx2 = 1.0 / (grid.dx * grid.dx);
    let inv_dy2 = 1.0 / (grid.dy * grid.dy);
    let denom = 2.0 * (inv_dx2 + inv_dy2);

    let mut phi_new = grid.phi.clone();

    for iter in 0..max_iters {
        let mut max_diff: f32 = 0.0;
        for iy in 0..ny {
            // Reflecting (zero-flux) BCs in y: ghost cell mirrors the interior.
            let iy_prev = if iy == 0 { 1 } else { iy - 1 };
            let iy_next = if iy == ny - 1 { ny - 2 } else { iy + 1 };
            for ix in 1..(nx - 1) {
                let p_left = grid.phi[grid.idx(ix - 1, iy)];
                let p_right = grid.phi[grid.idx(ix + 1, iy)];
                let p_down = grid.phi[grid.idx(ix, iy_prev)];
                let p_up = grid.phi[grid.idx(ix, iy_next)];
                let r = grid.rho[grid.idx(ix, iy)];
                let new_val = ((p_left + p_right) * inv_dx2
                    + (p_down + p_up) * inv_dy2
                    + r / epsilon0_eff)
                    / denom;
                let i = grid.idx(ix, iy);
                let d = (new_val - grid.phi[i]).abs();
                if d > max_diff {
                    max_diff = d;
                }
                phi_new[i] = new_val;
            }
        }
        // Preserve Dirichlet BCs in the swap buffer too.
        for iy in 0..ny {
            phi_new[grid.idx(0, iy)] = bcs.left;
            phi_new[grid.idx(nx - 1, iy)] = bcs.right;
        }
        std::mem::swap(&mut grid.phi, &mut phi_new);
        if max_diff < tol {
            return iter + 1;
        }
    }
    max_iters
}
