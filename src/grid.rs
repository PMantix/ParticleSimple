//! Coarse 2D grid for charge density, electric potential, and SEI blocking.
//!
//! Particles deposit charge onto the grid via cloud-in-cell (CIC); the field
//! is sampled back via the same kernel for force evaluation. Grid resolution
//! is chosen so that one cell ~ a few Debye lengths under the chosen unit
//! system.

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

    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        ix + iy * self.nx
    }

    pub fn cell_of(&self, p: Vec2) -> (usize, usize) {
        let fx = ((p.x - self.origin.x) / self.dx).max(0.0);
        let fy = ((p.y - self.origin.y) / self.dy).max(0.0);
        let ix = (fx as usize).min(self.nx - 1);
        let iy = (fy as usize).min(self.ny - 1);
        (ix, iy)
    }

    pub fn clear_rho(&mut self) {
        self.rho.iter_mut().for_each(|v| *v = 0.0);
    }

    /// Deposit point charges onto rho via CIC.
    pub fn deposit_charges<I: IntoIterator<Item = (Vec2, f32)>>(&mut self, _items: I) {
        todo!("CIC deposition")
    }

    /// Sample the field E = -grad phi at a point via CIC.
    pub fn field_at(&self, _p: Vec2) -> Vec2 {
        todo!("CIC field sampling")
    }

    /// Volume-weighted potential difference between two x-slabs.
    /// Used for terminal-voltage measurement.
    pub fn slab_potential(&self, _x_min: f32, _x_max: f32) -> f32 {
        todo!("slab potential averaging")
    }
}
