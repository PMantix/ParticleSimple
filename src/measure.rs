//! Per-step diagnostics + CSV output.
//!
//! `current` is the boundary current at the left electrode, in conventional
//! sign: positive when net stripping happens there (electrons leave the
//! left electrode through the external wire). This is the current a real
//! galvanostat would track.
//!
//! `midplane_current` is the older measurement (rate of net charge flow
//! across x=0), kept as a diagnostic. In steady state the two should agree;
//! during EDL transients they differ.

use std::fs::File;
use std::path::Path;

use crate::cell::Cell;
use crate::species::Species;

#[derive(Clone, Copy, Debug)]
pub struct Measurement {
    pub step: usize,
    pub voltage_applied: f32,
    pub voltage_bulk: f32,
    /// Boundary current at the left electrode (= external circuit current).
    pub current: f32,
    /// Rate of net charge flow across x=0 (legacy diagnostic).
    pub midplane_current: f32,
    /// Reaction event counts this step.
    pub plate_left: u32,
    pub strip_left: u32,
    pub plate_right: u32,
    pub strip_right: u32,
    pub sei_formed: u32,
    /// Number of immobile Metal particles (proxy for total deposit mass).
    pub metal_count: usize,
    /// Fraction of grid cells flagged as SEI-blocked.
    pub sei_fraction: f32,
}

impl Measurement {
    pub fn sample(cell: &Cell) -> Self {
        let voltage_applied = cell.last_bcs.left - cell.last_bcs.right;

        let v_left = column_mean(&cell.grid, 1);
        let v_right = column_mean(&cell.grid, cell.grid.nx - 2);
        let voltage_bulk = v_left - v_right;

        let charge_left = cell.charge_in_left_half();
        let dt = cell.params.dt;
        let midplane_current = -(charge_left - cell.charge_left_prev) / dt;

        let counts = cell.last_reaction_counts;
        // Boundary current at left, conventional sign:
        //   positive when net stripping at left (Li -> Li+ + e-).
        let current = (counts.strip_left as f32 - counts.plate_left as f32) / dt;

        let metal_count = cell
            .particles
            .iter()
            .filter(|p| p.species == Species::Metal)
            .count();

        let blocked = cell.grid.mobility.iter().filter(|&&m| m < 1.0).count() as f32;
        let total = cell.grid.mobility.len() as f32;
        let sei_fraction = if total > 0.0 { blocked / total } else { 0.0 };

        Self {
            step: cell.step_index,
            voltage_applied,
            voltage_bulk,
            current,
            midplane_current,
            plate_left: counts.plate_left,
            strip_left: counts.strip_left,
            plate_right: counts.plate_right,
            strip_right: counts.strip_right,
            sei_formed: counts.sei_formed,
            metal_count,
            sei_fraction,
        }
    }
}

fn column_mean(grid: &crate::grid::Grid, ix: usize) -> f32 {
    let mut sum = 0.0;
    for iy in 0..grid.ny {
        sum += grid.phi[grid.idx(ix, iy)];
    }
    sum / grid.ny as f32
}

pub struct CsvSink {
    writer: csv::Writer<File>,
}

impl CsvSink {
    pub fn new<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let mut writer = csv::Writer::from_path(path)?;
        writer
            .write_record([
                "step",
                "voltage_applied",
                "voltage_bulk",
                "current",
                "midplane_current",
                "plate_left",
                "strip_left",
                "plate_right",
                "strip_right",
                "metal_count",
                "sei_fraction",
            ])
            .map_err(std::io::Error::other)?;
        Ok(Self { writer })
    }

    pub fn write(&mut self, m: &Measurement) -> std::io::Result<()> {
        self.writer
            .write_record(&[
                m.step.to_string(),
                m.voltage_applied.to_string(),
                m.voltage_bulk.to_string(),
                m.current.to_string(),
                m.midplane_current.to_string(),
                m.plate_left.to_string(),
                m.strip_left.to_string(),
                m.plate_right.to_string(),
                m.strip_right.to_string(),
                m.metal_count.to_string(),
                m.sei_fraction.to_string(),
            ])
            .map_err(std::io::Error::other)?;
        Ok(())
    }
}
