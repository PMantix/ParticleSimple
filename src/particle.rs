use glam::Vec2;

use crate::species::Species;

/// A single particle. Overdamped, so no velocity field is stored.
#[derive(Clone, Debug)]
pub struct Particle {
    pub pos: Vec2,
    pub species: Species,
}

impl Particle {
    pub fn new(pos: Vec2, species: Species) -> Self {
        Self { pos, species }
    }

    pub fn charge(&self) -> f32 {
        self.species.props().charge
    }
}
