//! Rectangular simulation domain with two electrode boundaries (left, right)
//! and reflecting / periodic boundaries on top and bottom.

use glam::Vec2;

#[derive(Clone, Copy, Debug)]
pub struct Domain {
    /// Half-width (electrodes sit at x = ±half_width).
    pub half_width: f32,
    /// Half-height.
    pub half_height: f32,
}

impl Domain {
    pub fn new(half_width: f32, half_height: f32) -> Self {
        Self {
            half_width,
            half_height,
        }
    }

    pub fn contains(&self, p: Vec2) -> bool {
        p.x.abs() <= self.half_width && p.y.abs() <= self.half_height
    }

    pub fn area(&self) -> f32 {
        4.0 * self.half_width * self.half_height
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Electrode {
    Left,
    Right,
}

impl Electrode {
    pub fn x_position(self, domain: &Domain) -> f32 {
        match self {
            Electrode::Left => -domain.half_width,
            Electrode::Right => domain.half_width,
        }
    }

    /// Outward normal direction for this electrode.
    pub fn normal(self) -> Vec2 {
        match self {
            Electrode::Left => Vec2::new(-1.0, 0.0),
            Electrode::Right => Vec2::new(1.0, 0.0),
        }
    }
}
