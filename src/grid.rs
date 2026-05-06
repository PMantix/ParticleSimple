//! Coarse 2D grid for charge density, electric potential, and SEI blocking.
//!
//! Particles deposit charge onto the grid via cloud-in-cell (CIC); the field
//! is sampled back via bilinear interpolation of phi for force evaluation.
//! Grid resolution is chosen so that one cell ~ a few Debye lengths under
//! the chosen unit system.

use glam::Vec2;

use crate::domain::Domain;

#[derive(Clone, Debug)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    pub dx: f32,
    pub dy: f32,
    /// Origin (lower-left corner) in domain coordinates.
    pub origin: Vec2,
    /// Charge density rho[ix + iy*nx].
    pub rho: Vec<f32>,
    /// Electric potential phi[ix + iy*nx]. Updated by Poisson solver.
    pub phi: Vec<f32>,
    /// Per-cell mobility mask (1.0 = free, < 1.0 = blocked by SEI).
    pub mobility: Vec<f32>,
}

impl Grid {
    pub fn new(domain: &Domain, nx: usize, ny: usize) -> Self {
        assert!(nx >= 2 && ny >= 2, "grid must be at least 2x2");
        let dx = (2.0 * domain.half_width) / nx as f32;
        let dy = (2.0 * domain.half_height) / ny as f32;
        let n = nx * ny;
        Self {
            nx,
            ny,
            dx,
            dy,
            origin: Vec2::new(-domain.half_width, -domain.half_height),
            rho: vec![0.0; n],
            phi: vec![0.0; n],
            mobility: vec![1.0; n],
        }
    }

    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        ix + iy * self.nx
    }

    pub fn clear_rho(&mut self) {
        self.rho.iter_mut().for_each(|v| *v = 0.0);
    }

    /// Cloud-in-cell charge deposition. Each point charge contributes to its
    /// four nearest grid nodes weighted by bilinear overlap, normalized by
    /// cell area so rho has units of charge / area.
    pub fn deposit_charges<I: IntoIterator<Item = (Vec2, f32)>>(&mut self, items: I) {
        let area = self.dx * self.dy;
        let nx_f = (self.nx - 1) as f32;
        let ny_f = (self.ny - 1) as f32;
        for (p, q) in items {
            let gx = ((p.x - self.origin.x) / self.dx).clamp(0.0, nx_f);
            let gy = ((p.y - self.origin.y) / self.dy).clamp(0.0, ny_f);
            let ix = (gx as usize).min(self.nx - 2);
            let iy = (gy as usize).min(self.ny - 2);
            let wx = gx - ix as f32;
            let wy = gy - iy as f32;
            let qa = q / area;
            let i00 = self.idx(ix, iy);
            let i10 = self.idx(ix + 1, iy);
            let i01 = self.idx(ix, iy + 1);
            let i11 = self.idx(ix + 1, iy + 1);
            self.rho[i00] += qa * (1.0 - wx) * (1.0 - wy);
            self.rho[i10] += qa * wx * (1.0 - wy);
            self.rho[i01] += qa * (1.0 - wx) * wy;
            self.rho[i11] += qa * wx * wy;
        }
    }

    /// Sample E = -grad phi at a point. Computes the gradient of the bilinear
    /// interpolant over the cell containing the point.
    pub fn field_at(&self, p: Vec2) -> Vec2 {
        let nx_f = (self.nx - 1) as f32;
        let ny_f = (self.ny - 1) as f32;
        let gx = ((p.x - self.origin.x) / self.dx).clamp(0.0, nx_f);
        let gy = ((p.y - self.origin.y) / self.dy).clamp(0.0, ny_f);
        let ix = (gx as usize).min(self.nx - 2);
        let iy = (gy as usize).min(self.ny - 2);
        let wx = gx - ix as f32;
        let wy = gy - iy as f32;
        let p00 = self.phi[self.idx(ix, iy)];
        let p10 = self.phi[self.idx(ix + 1, iy)];
        let p01 = self.phi[self.idx(ix, iy + 1)];
        let p11 = self.phi[self.idx(ix + 1, iy + 1)];
        let ex = -((1.0 - wy) * (p10 - p00) + wy * (p11 - p01)) / self.dx;
        let ey = -((1.0 - wx) * (p01 - p00) + wx * (p11 - p10)) / self.dy;
        Vec2::new(ex, ey)
    }

    /// Locate the grid cell containing a point.
    pub fn cell_of(&self, p: Vec2) -> (usize, usize) {
        let nx_f = (self.nx - 1) as f32;
        let ny_f = (self.ny - 1) as f32;
        let gx = ((p.x - self.origin.x) / self.dx).clamp(0.0, nx_f);
        let gy = ((p.y - self.origin.y) / self.dy).clamp(0.0, ny_f);
        let ix = (gx as usize).min(self.nx - 1);
        let iy = (gy as usize).min(self.ny - 1);
        (ix, iy)
    }

    /// Sample phi at a point via bilinear interpolation. Used by reactions
    /// to read local Galvani potential at a particle's position.
    pub fn sample_phi(&self, p: Vec2) -> f32 {
        let nx_f = (self.nx - 1) as f32;
        let ny_f = (self.ny - 1) as f32;
        let gx = ((p.x - self.origin.x) / self.dx).clamp(0.0, nx_f);
        let gy = ((p.y - self.origin.y) / self.dy).clamp(0.0, ny_f);
        let ix = (gx as usize).min(self.nx - 2);
        let iy = (gy as usize).min(self.ny - 2);
        let wx = gx - ix as f32;
        let wy = gy - iy as f32;
        let p00 = self.phi[self.idx(ix, iy)];
        let p10 = self.phi[self.idx(ix + 1, iy)];
        let p01 = self.phi[self.idx(ix, iy + 1)];
        let p11 = self.phi[self.idx(ix + 1, iy + 1)];
        (1.0 - wx) * (1.0 - wy) * p00
            + wx * (1.0 - wy) * p10
            + (1.0 - wx) * wy * p01
            + wx * wy * p11
    }

    /// Average phi over cells whose centers fall in [x_min, x_max].
    /// Used for terminal-voltage measurement next to electrode boundaries.
    pub fn slab_potential(&self, x_min: f32, x_max: f32) -> f32 {
        let mut sum = 0.0;
        let mut count = 0;
        for ix in 0..self.nx {
            let cx = self.origin.x + (ix as f32 + 0.5) * self.dx;
            if cx >= x_min && cx <= x_max {
                for iy in 0..self.ny {
                    sum += self.phi[self.idx(ix, iy)];
                    count += 1;
                }
            }
        }
        if count == 0 {
            0.0
        } else {
            sum / count as f32
        }
    }
}
