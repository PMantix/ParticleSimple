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

## Goals (in order)

1. Reproduce the three-timescale shape of the DCR pulse
   (R0 ohmic + R1 charge-transfer + R2 diffusion).
2. Show R1 dropping as the deposition front roughens.
3. Show R0 and R2 rising as SEI accumulates.

## Run

```
cargo run --release --bin psimple -- --scenario pulse_dcr
```

Scenarios live under `src/scenarios/` and dump CSV for external plotting.

## Status

Day 1 sketch: types and module boundaries are in place; physics bodies are
`todo!()` placeholders. See the audit notes in the parent ParticleSim repo
for design rationale.
