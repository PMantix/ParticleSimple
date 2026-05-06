//! Particle species and their dimensionless properties.
//!
//! All numeric values are in simulation units (no fs / Å / amu). The model
//! relies on similitude: only the ratios between species' diffusivities and
//! charges matter. Pick units so that the diffusion timescale across the
//! domain is `O(100)` steps and the reaction timescale is `O(1-10)` steps.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Species {
    /// Mobile cation (Li+ analog).
    Cation,
    /// Mobile counter-anion (PF6- analog).
    Anion,
    /// Mobile neutral solvent. Reducible at the metal anode -> SEI.
    Solvent,
    /// Deposited metal. Immobile, neutral, conductive (extends anode geometry).
    Metal,
    /// Solid electrolyte interphase. Immobile, neutral, blocks transport.
    Sei,
}

#[derive(Clone, Copy, Debug)]
pub struct Props {
    /// Charge in elementary units (sign included).
    pub charge: f32,
    /// Diffusion coefficient in sim units (length^2 / step).
    pub d: f32,
    /// Excluded-volume radius (sim length).
    pub sigma: f32,
    /// Whether the particle moves under Langevin dynamics each step.
    pub mobile: bool,
}

impl Species {
    pub fn props(self) -> Props {
        match self {
            Species::Cation => Props {
                charge: 1.0,
                d: 1.0,
                sigma: 0.5,
                mobile: true,
            },
            Species::Anion => Props {
                charge: -1.0,
                d: 1.0,
                sigma: 0.5,
                mobile: true,
            },
            Species::Solvent => Props {
                charge: 0.0,
                d: 0.5,
                sigma: 0.6,
                mobile: true,
            },
            Species::Metal => Props {
                charge: 0.0,
                d: 0.0,
                sigma: 0.5,
                mobile: false,
            },
            Species::Sei => Props {
                charge: 0.0,
                d: 0.0,
                sigma: 0.7,
                mobile: false,
            },
        }
    }
}
