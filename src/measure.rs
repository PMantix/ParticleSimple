//! Per-step diagnostics + CSV output.

use std::fs::File;
use std::path::Path;

use crate::cell::Cell;
use crate::species::Species;

#[derive(Clone, Copy, Debug)]
pub struct Measurement {
    pub step: usize,
    /// Applied terminal voltage = last_bcs.left - last_bcs.right.
    pub voltage_applied: f32,
    /// Terminal voltage measured one cell inside the electrodes (captures
    /// EDL screening relative to the applied BC).
    pub voltage_bulk: f32,
    /// Net current crossing the midplane this step (sim charge / time).
    pub current: f32,
    /// Number of immobile Metal particles (proxy for A_eff while we don't
    /// yet distinguish surface from interior).
    pub exposed_metal_sites: usize,
    /// Fraction of grid cells flagged as SEI-blocked.
    pub sei_fraction: f32,
}

impl Measurement {
    pub fn sample(cell: &Cell) -> Self {
        let voltage_applied = cell.last_bcs.left - cell.last_bcs.right;

        // Average phi over the column just inside each Dirichlet boundary
        // (ix=1 on the left, ix=nx-2 on the right). Captures EDL screening.
        let v_left = column_mean(&cell.grid, 1);
        let v_right = column_mean(&cell.grid, cell.grid.nx - 2);
        let voltage_bulk = v_left - v_right;

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
            voltage_applied,
            voltage_bulk,
            current,
            exposed_metal_sites,
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

/// Stream measurements to a CSV file. Closes automatically on drop.
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
                "exposed_metal",
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
                m.exposed_metal_sites.to_string(),
                m.sei_fraction.to_string(),
            ])
            .map_err(std::io::Error::other)?;
        Ok(())
    }
}
