//! Per-step diagnostics + CSV output.
//!
//! Records terminal voltage, current (rate of midplane charge flow), exposed
//! metal-electrolyte interface area, and SEI coverage. From a recorded stream
//! we can fit an R0 / (R1, tau1) / (R2, tau2) ECM after a step pulse.

use std::fs::File;
use std::path::Path;

use crate::cell::Cell;
use crate::protocol::Drive;
use crate::species::Species;

#[derive(Clone, Copy, Debug)]
pub struct Measurement {
    pub step: usize,
    /// Terminal V_left - V_right (sim potential units).
    pub voltage: f32,
    /// Net current crossing the midplane (sim charge / step).
    pub current: f32,
    /// Number of immobile Metal particles (proxy for A_eff while we don't yet
    /// distinguish surface from interior).
    pub exposed_metal_sites: usize,
    /// Fraction of grid cells flagged as SEI-blocked.
    pub sei_fraction: f32,
}

impl Measurement {
    pub fn sample(cell: &Cell, _drive: Drive) -> Self {
        let half_w = cell.domain.half_width;
        let probe = 0.05 * half_w;
        let v_left = cell.grid.slab_potential(-half_w, -half_w + probe);
        let v_right = cell.grid.slab_potential(half_w - probe, half_w);
        let voltage = v_left - v_right;

        let charge_left = cell.charge_in_left_half();
        let dt = cell.params.dt;
        let current = -(charge_left - cell.charge_left_prev) / dt;

        let exposed_metal_sites = cell
            .particles
            .iter()
            .filter(|p| p.species == Species::Metal)
            .count();

        let blocked = cell.grid.mobility.iter().filter(|&&m| m < 1.0).count() as f32;
        let total = cell.grid.mobility.len() as f32;
        let sei_fraction = if total > 0.0 { blocked / total } else { 0.0 };

        Self {
            step: cell.step_index,
            voltage,
            current,
            exposed_metal_sites,
            sei_fraction,
        }
    }
}

/// Stream measurements to a CSV file. Closes automatically on drop.
pub struct CsvSink {
    writer: csv::Writer<File>,
}

impl CsvSink {
    pub fn new<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let mut writer = csv::Writer::from_path(path)?;
        writer
            .write_record(["step", "voltage", "current", "exposed_metal", "sei_fraction"])
            .map_err(std::io::Error::other)?;
        Ok(Self { writer })
    }

    pub fn write(&mut self, m: &Measurement) -> std::io::Result<()> {
        self.writer
            .write_record(&[
                m.step.to_string(),
                m.voltage.to_string(),
                m.current.to_string(),
                m.exposed_metal_sites.to_string(),
                m.sei_fraction.to_string(),
            ])
            .map_err(std::io::Error::other)?;
        Ok(())
    }
}
