//! Live macroquad GUI for the pulse_dcr scenario.
//!
//! Build:
//!   cargo build --release --bin psimple_gui --features gui
//! Run:
//!   cargo run   --release --bin psimple_gui --features gui
//!
//! Controls
//!   space  pause / resume
//!   r      reset to step 0
//!   + / -  steps per render frame (1..64)
//!   q/esc  quit

use macroquad::prelude::*;

use particle_simple::domain::Domain;
use particle_simple::measure::Measurement;
use particle_simple::scenarios;
use particle_simple::species::Species;

const WINDOW_W: f32 = 1200.0;
const WINDOW_H: f32 = 720.0;

fn window_conf() -> Conf {
    Conf {
        window_title: "ParticleSimple — pulse_dcr".to_string(),
        window_width: WINDOW_W as i32,
        window_height: WINDOW_H as i32,
        high_dpi: true,
        ..Default::default()
    }
}

struct Viewport {
    x: f32,
    y: f32,
    w: f32,
    h: f32,
    scale: f32,
}

fn make_viewport(domain: &Domain) -> Viewport {
    let dw = 2.0 * domain.half_width;
    let dh = 2.0 * domain.half_height;
    let max_w = WINDOW_W - 240.0;
    let max_h = WINDOW_H - 220.0;
    let scale = (max_w / dw).min(max_h / dh);
    let w = dw * scale;
    let h = dh * scale;
    let x = (WINDOW_W - w) / 2.0;
    let y = 110.0 + (max_h - h) / 2.0;
    Viewport { x, y, w, h, scale }
}

fn sim_to_screen(px: f32, py: f32, domain: &Domain, vp: &Viewport) -> (f32, f32) {
    let sx = vp.x + (px + domain.half_width) * vp.scale;
    let sy = vp.y + (domain.half_height - py) * vp.scale;
    (sx, sy)
}

fn species_style(s: Species) -> (Color, f32) {
    match s {
        Species::Cation => (Color::new(0.95, 0.30, 0.30, 1.0), 2.5),
        Species::Anion => (Color::new(0.30, 0.55, 0.95, 1.0), 2.5),
        Species::Solvent => (Color::new(0.55, 0.85, 0.55, 0.55), 2.5),
        Species::Metal => (Color::new(0.97, 0.70, 0.15, 1.0), 5.0),
        Species::Sei => (Color::new(0.20, 0.60, 0.25, 1.0), 4.0),
    }
}

fn legend_item(x: f32, y: f32, label: &str, style: (Color, f32)) {
    draw_circle(x, y, style.1, style.0);
    draw_text(label, x + 12.0, y + 5.0, 16.0, Color::new(0.9, 0.9, 0.9, 1.0));
}

#[macroquad::main(window_conf)]
async fn main() {
    let mut cell = scenarios::pulse_dcr::setup();
    let mut paused = false;
    let mut steps_per_frame: u32 = 4;
    let mut last: Option<Measurement> = None;

    let cycle = scenarios::pulse_dcr::RELAX_STEPS + scenarios::pulse_dcr::PULSE_STEPS;
    let relax = scenarios::pulse_dcr::RELAX_STEPS;

    loop {
        // ---- input ----
        if is_key_pressed(KeyCode::Space) {
            paused = !paused;
        }
        if is_key_pressed(KeyCode::R) {
            cell = scenarios::pulse_dcr::setup();
            last = None;
        }
        if is_key_pressed(KeyCode::Equal) || is_key_pressed(KeyCode::KpAdd) {
            steps_per_frame = (steps_per_frame * 2).min(64);
        }
        if is_key_pressed(KeyCode::Minus) || is_key_pressed(KeyCode::KpSubtract) {
            steps_per_frame = (steps_per_frame / 2).max(1);
        }
        if is_key_pressed(KeyCode::Q) || is_key_pressed(KeyCode::Escape) {
            break;
        }

        // ---- step ----
        if !paused {
            for _ in 0..steps_per_frame {
                last = Some(cell.step());
            }
        }

        // ---- draw ----
        clear_background(Color::new(0.10, 0.10, 0.12, 1.0));
        let vp = make_viewport(&cell.domain);

        // cell outline + electrodes
        draw_rectangle_lines(vp.x, vp.y, vp.w, vp.h, 2.0, Color::new(0.45, 0.45, 0.50, 1.0));
        let elec = Color::new(0.75, 0.75, 0.78, 1.0);
        draw_rectangle(vp.x - 8.0, vp.y, 8.0, vp.h, elec);
        draw_rectangle(vp.x + vp.w, vp.y, 8.0, vp.h, elec);
        draw_text("(+) anode side", vp.x - 24.0, vp.y - 10.0, 16.0, elec);
        draw_text(
            "(-) cathode side",
            vp.x + vp.w - 60.0,
            vp.y - 10.0,
            16.0,
            elec,
        );

        // particles (draw immobile species first so mobile ions render on top)
        for p in &cell.particles {
            if matches!(p.species, Species::Metal | Species::Sei) {
                let (c, r) = species_style(p.species);
                let (sx, sy) = sim_to_screen(p.pos.x, p.pos.y, &cell.domain, &vp);
                draw_circle(sx, sy, r, c);
            }
        }
        for p in &cell.particles {
            if !matches!(p.species, Species::Metal | Species::Sei) {
                let (c, r) = species_style(p.species);
                let (sx, sy) = sim_to_screen(p.pos.x, p.pos.y, &cell.domain, &vp);
                draw_circle(sx, sy, r, c);
            }
        }

        // pulse phase progress bar under the cell
        let step = last.as_ref().map(|m| m.step).unwrap_or(0);
        let phase = step % cycle;
        let in_pulse = phase >= relax;
        let phase_color = if in_pulse {
            Color::new(0.95, 0.40, 0.30, 1.0)
        } else {
            Color::new(0.40, 0.55, 0.95, 1.0)
        };
        let phase_y = vp.y + vp.h + 14.0;
        draw_rectangle(vp.x, phase_y, vp.w, 6.0, Color::new(0.18, 0.18, 0.20, 1.0));
        let frac = phase as f32 / cycle as f32;
        draw_rectangle(vp.x, phase_y, vp.w * frac, 6.0, phase_color);
        let phase_label = if in_pulse { "pulse" } else { "relax" };
        draw_text(
            &format!("cycle phase: {}  ({}/{})", phase_label, phase, cycle),
            vp.x,
            phase_y + 24.0,
            14.0,
            Color::new(0.7, 0.7, 0.7, 1.0),
        );

        // ---- header ----
        draw_text("ParticleSimple — pulse_dcr (live)", 20.0, 32.0, 24.0, WHITE);

        let info = if let Some(m) = &last {
            format!(
                "step {:>5}   V_app {:+.3}   V_bulk {:+.3}   I {:+.0}   metal {:>3}   sei {:.2}",
                m.step,
                m.voltage_applied,
                m.voltage_bulk,
                m.current,
                m.metal_count,
                m.sei_fraction
            )
        } else {
            "ready".to_string()
        };
        draw_text(&info, 20.0, 60.0, 18.0, Color::new(0.92, 0.92, 0.92, 1.0));

        let pause_str = if paused { "PAUSED" } else { "running" };
        draw_text(
            &format!("{}  |  {} step(s) / frame", pause_str, steps_per_frame),
            20.0,
            82.0,
            16.0,
            Color::new(0.7, 0.7, 0.7, 1.0),
        );

        // ---- footer: legend + controls ----
        let legend_y = WINDOW_H - 50.0;
        legend_item(30.0, legend_y, "Cation", species_style(Species::Cation));
        legend_item(130.0, legend_y, "Anion", species_style(Species::Anion));
        legend_item(230.0, legend_y, "Solvent", species_style(Species::Solvent));
        legend_item(340.0, legend_y, "Metal", species_style(Species::Metal));
        legend_item(440.0, legend_y, "SEI", species_style(Species::Sei));
        draw_text(
            "[space] pause   [r] reset   [+/-] speed   [q/esc] quit",
            30.0,
            WINDOW_H - 22.0,
            16.0,
            Color::new(0.6, 0.6, 0.6, 1.0),
        );

        next_frame().await;
    }
}
