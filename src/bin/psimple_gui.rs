//! Live macroquad GUI for the pulse_dcr scenario.
//!
//! Build:
//!   cargo build --release --bin psimple_gui --features gui
//! Run:
//!   cargo run   --release --bin psimple_gui --features gui
//!
//! Visual layers, back to front:
//!   1. dark backdrop
//!   2. phi (electric potential) field rendered from the grid as a low-res
//!      red/blue divergent texture, linearly upsampled to the viewport
//!   3. electrode plates and the cell rectangle frame
//!   4. SEI and metal particles, drawn with rim/fill/highlight to read as
//!      shaded discs rather than flat dots
//!   5. solvent (faint), then cations and anions on top
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

const WINDOW_W: f32 = 1280.0;
const WINDOW_H: f32 = 760.0;

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
    scale: f32, // pixels per sim length unit
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

/// Map a phi value in [-phi_max, +phi_max] to a divergent red/blue background
/// color. Negative goes blue (cathode side), positive goes red (anode side),
/// near zero goes near-black so subtle structure stays readable.
fn phi_to_bg(phi: f32, phi_max: f32) -> Color {
    let t = (phi / phi_max).clamp(-1.0, 1.0);
    let mag = t.abs();
    if t > 0.0 {
        Color::new(0.10 + mag * 0.55, 0.06 + mag * 0.05, 0.10, 1.0)
    } else {
        Color::new(0.08, 0.08 + mag * 0.05, 0.10 + mag * 0.55, 1.0)
    }
}

#[derive(Clone, Copy)]
struct ParticleStyle {
    rim: Color,
    fill: Color,
    highlight: Color,
    radius_px: f32,
}

fn species_style(s: Species) -> ParticleStyle {
    match s {
        Species::Cation => ParticleStyle {
            rim: Color::new(0.55, 0.20, 0.05, 1.0),
            fill: Color::new(0.98, 0.55, 0.25, 0.95),
            highlight: Color::new(1.0, 0.90, 0.65, 0.85),
            radius_px: 4.5,
        },
        Species::Anion => ParticleStyle {
            rim: Color::new(0.05, 0.20, 0.55, 1.0),
            fill: Color::new(0.30, 0.65, 0.98, 0.95),
            highlight: Color::new(0.75, 0.90, 1.0, 0.85),
            radius_px: 4.5,
        },
        Species::Solvent => ParticleStyle {
            rim: Color::new(0.35, 0.50, 0.35, 0.40),
            fill: Color::new(0.70, 0.88, 0.70, 0.35),
            highlight: Color::new(0.95, 1.00, 0.95, 0.20),
            radius_px: 5.0,
        },
        Species::Metal => ParticleStyle {
            rim: Color::new(0.40, 0.28, 0.05, 1.0),
            fill: Color::new(0.95, 0.78, 0.30, 1.0),
            highlight: Color::new(1.0, 0.95, 0.70, 0.95),
            radius_px: 11.0,
        },
        Species::Sei => ParticleStyle {
            rim: Color::new(0.12, 0.22, 0.12, 1.0),
            fill: Color::new(0.40, 0.60, 0.32, 0.92),
            highlight: Color::new(0.60, 0.85, 0.55, 0.70),
            radius_px: 9.0,
        },
    }
}

fn draw_particle(x: f32, y: f32, style: ParticleStyle) {
    let r = style.radius_px;
    draw_circle(x, y, r, style.rim);
    draw_circle(x, y, r * 0.75, style.fill);
    draw_circle(x - r * 0.30, y - r * 0.30, r * 0.30, style.highlight);
}

fn legend_chip(x: f32, y: f32, label: &str, style: ParticleStyle) {
    draw_particle(x, y, style);
    draw_text(label, x + 14.0, y + 5.0, 16.0, Color::new(0.92, 0.92, 0.92, 1.0));
}

#[macroquad::main(window_conf)]
async fn main() {
    let mut cell = scenarios::pulse_dcr::setup();
    let mut paused = false;
    let mut steps_per_frame: u32 = 4;
    let mut last: Option<Measurement> = None;

    let cycle = scenarios::pulse_dcr::RELAX_STEPS + scenarios::pulse_dcr::PULSE_STEPS;
    let relax = scenarios::pulse_dcr::RELAX_STEPS;
    let pulse_v = scenarios::pulse_dcr::PULSE_VOLTAGE;

    // Preallocate the phi backdrop texture (grid resolution; linearly
    // upsampled when blitted into the viewport).
    let nx = cell.grid.nx;
    let ny = cell.grid.ny;
    let mut phi_image = Image::gen_image_color(nx as u16, ny as u16, BLACK);
    let phi_texture = Texture2D::from_image(&phi_image);
    phi_texture.set_filter(FilterMode::Linear);

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

        // ---- update phi backdrop ----
        let phi_max = pulse_v * 0.7; // saturates the colormap a bit before BC clamp
        for iy in 0..ny {
            for ix in 0..nx {
                let phi = cell.grid.phi[cell.grid.idx(ix, iy)];
                let c = phi_to_bg(phi, phi_max);
                phi_image.set_pixel(ix as u32, (ny - 1 - iy) as u32, c);
            }
        }
        phi_texture.update(&phi_image);

        // ---- draw ----
        clear_background(Color::new(0.05, 0.05, 0.07, 1.0));
        let vp = make_viewport(&cell.domain);

        // backdrop: phi field, linearly filtered
        draw_texture_ex(
            &phi_texture,
            vp.x,
            vp.y,
            WHITE,
            DrawTextureParams {
                dest_size: Some(vec2(vp.w, vp.h)),
                ..Default::default()
            },
        );

        // cell frame and electrode plates
        draw_rectangle_lines(vp.x, vp.y, vp.w, vp.h, 2.0, Color::new(0.55, 0.55, 0.60, 1.0));
        let elec_left = Color::new(0.85, 0.55, 0.45, 1.0);
        let elec_right = Color::new(0.45, 0.65, 0.85, 1.0);
        draw_rectangle(vp.x - 8.0, vp.y, 8.0, vp.h, elec_left);
        draw_rectangle(vp.x + vp.w, vp.y, 8.0, vp.h, elec_right);
        draw_text(
            "(+) anode",
            vp.x - 26.0,
            vp.y - 10.0,
            16.0,
            elec_left,
        );
        draw_text(
            "(−) cathode",
            vp.x + vp.w - 30.0,
            vp.y - 10.0,
            16.0,
            elec_right,
        );

        // particles, drawn in z-like order (deposits and SEI under the ions)
        for p in &cell.particles {
            if !matches!(p.species, Species::Metal | Species::Sei) {
                continue;
            }
            let style = species_style(p.species);
            let (sx, sy) = sim_to_screen(p.pos.x, p.pos.y, &cell.domain, &vp);
            draw_particle(sx, sy, style);
        }
        for p in &cell.particles {
            if !matches!(p.species, Species::Solvent) {
                continue;
            }
            let style = species_style(p.species);
            let (sx, sy) = sim_to_screen(p.pos.x, p.pos.y, &cell.domain, &vp);
            draw_particle(sx, sy, style);
        }
        for p in &cell.particles {
            if !matches!(p.species, Species::Cation | Species::Anion) {
                continue;
            }
            let style = species_style(p.species);
            let (sx, sy) = sim_to_screen(p.pos.x, p.pos.y, &cell.domain, &vp);
            draw_particle(sx, sy, style);
        }

        // pulse phase progress bar
        let step = last.as_ref().map(|m| m.step).unwrap_or(0);
        let phase = step % cycle;
        let in_pulse = phase >= relax;
        let phase_color = if in_pulse {
            Color::new(0.95, 0.40, 0.30, 1.0)
        } else {
            Color::new(0.40, 0.55, 0.95, 1.0)
        };
        let phase_y = vp.y + vp.h + 18.0;
        draw_rectangle(vp.x, phase_y, vp.w, 6.0, Color::new(0.18, 0.18, 0.22, 1.0));
        let frac = phase as f32 / cycle as f32;
        draw_rectangle(vp.x, phase_y, vp.w * frac, 6.0, phase_color);
        let phase_label = if in_pulse { "pulse" } else { "relax" };
        draw_text(
            &format!("{}  ({}/{})", phase_label, phase, cycle),
            vp.x,
            phase_y + 24.0,
            14.0,
            Color::new(0.7, 0.7, 0.72, 1.0),
        );

        // header
        draw_text(
            "ParticleSimple — pulse_dcr (live)",
            20.0,
            34.0,
            24.0,
            WHITE,
        );

        let info = if let Some(m) = &last {
            format!(
                "step {:>5}    V_app {:+.3}    V_bulk {:+.3}    I {:+.0}    metal {:>3}    sei {:.2}",
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
        draw_text(&info, 20.0, 64.0, 18.0, Color::new(0.92, 0.92, 0.93, 1.0));

        let pause_str = if paused { "PAUSED" } else { "running" };
        draw_text(
            &format!("{}    {} step(s) / frame", pause_str, steps_per_frame),
            20.0,
            86.0,
            16.0,
            Color::new(0.7, 0.7, 0.72, 1.0),
        );

        // legend + controls
        let legend_y = WINDOW_H - 50.0;
        legend_chip(30.0, legend_y, "Cation", species_style(Species::Cation));
        legend_chip(140.0, legend_y, "Anion", species_style(Species::Anion));
        legend_chip(250.0, legend_y, "Solvent", species_style(Species::Solvent));
        legend_chip(370.0, legend_y, "Metal", species_style(Species::Metal));
        legend_chip(470.0, legend_y, "SEI", species_style(Species::Sei));
        draw_text(
            "[space] pause   [r] reset   [+/-] speed   [q/esc] quit",
            30.0,
            WINDOW_H - 22.0,
            16.0,
            Color::new(0.6, 0.6, 0.62, 1.0),
        );

        next_frame().await;
    }
}
