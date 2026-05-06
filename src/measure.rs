//! Per-step diagnostics + CSV output.
//!
//! Records terminal voltage, applied current, exposed metal-electrolyte
//! interface area, and SEI coverage. From a recorded `Measurement` stream we
//! can fit an R0 / (R1, tau1) / (R2, tau2) ECM after a step.

use std::fs::File;
use std::path::Path;

use crate::cell::Cell;
use crate::protocol::Drive;

#[derive(Clone, Copy, Debug)]
pub struct Measurement {
    pub step: usize,
    /// Terminal V_left - V_right (sim potential units).
    pub voltage: f32,
    /// Net current crossing the cell (sim charge / step).
    pub current: f32,
    /// Number of exposed Metal/Cation interface contacts (proxy for A_eff).
    pub exposed_metal_sites: usize,
    /// Fraction of anode surface cells flagged as SEI-blocked.
    pub sei_fraction: f32,
}

impl Measurement {
    pub fn sample(_cell: &Cell, _drive: Drive) -> Self {
        todo!("read terminal voltage from grid, count surface contacts, compute SEI fraction")
    }
}

/// Stream measurements to a CSV file. Closed automatically on drop.
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
