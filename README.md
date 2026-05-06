# ParticleSimple

A qualitative 2D particle simulator for studying electrochemical impulse-response
signatures driven by morphology. Forked in spirit from ParticleSim but rebuilt
to operate in a dimensionless, similitude-scaled regime so that the impulse
window (ohmic / charge-transfer / diffusion timescales) is reachable in seconds
of compute.

## Why this exists

The full ParticleSim runs inertial MD at a 5 fs timestep, which makes 10 s of
physical pulse data unreachable by ~12 orders of magnitude. ParticleSimple
discards inertia, dipole solvent, pairwise Coulomb, and the rescaling
thermostat, replacing them with:

- Overdamped Langevin dynamics
- Mean-field Poisson on a coarse grid
- Stochastic Butler-Volmer reactions at electrode boundaries
- Deposition / SEI rules that feed back on field and transport

The model is not a quantitative cell simulator. It is a qualitative generator
of morphology and a probe of the resulting impulse response.

## Goals

1. Reproduce the three-timescale shape of the DCR pulse
   (R0 ohmic + R1 charge-transfer + R2 diffusion).
2. Show R1 dropping as the deposition front roughens.
3. Show R0 and R2 rising as SEI accumulates.

## Run

### Headless (CSV out, fast)

```
cargo run --release --bin psimple -- --scenario pulse_dcr
cargo run --release --bin psimple -- --scenario empty_cell
cargo run --release --bin psimple -- --scenario galvanostatic_check
```

Each scenario writes a CSV in the current working directory. Plot with:

```
python3 scripts/plot_pulse_dcr.py
```

### Live GUI (macroquad, watch the run unfold)

```
cargo run --release --bin psimple_gui --features gui
```

Window opens with the cell viewport, particle species color-coded
(cation / anion / solvent / metal / SEI), and live readouts for the
applied voltage, bulk voltage, boundary current, metal count, and SEI
fraction. Controls:

| key | action |
| --- | --- |
| `space` | pause / resume |
| `r` | reset to step 0 |
| `+` / `-` | steps per render frame (1..64) |
| `q` / `esc` | quit |

## Layout

```
src/
  particle.rs   species.rs   domain.rs       data types
  grid.rs       poisson.rs   langevin.rs     physics core
  reactions.rs  protocol.rs  cell.rs         orchestration
  measure.rs                                 CSV out
  scenarios/    *.rs                         named experiments
  bin/psimple_gui.rs                         macroquad live window
scripts/
  plot_pulse_dcr.py                          matplotlib plots from CSV
```
