//! Headline experiment: galvanostatic step pulse with relaxation.
//! Apply 0 -> I -> 0 step train, record V(t), fit R0 / (R1,tau1) / (R2,tau2)
//! against the recorded trace. Repeat across cycles to watch R1 fall as the
//! anode roughens and R0/R2 rise as SEI accumulates.

pub fn run() {
    todo!(
        "spawn cell with cation/anion/solvent populations, run StepPulse protocol \
         for ~10 pulses, write CSV including exposed_metal and sei_fraction columns"
    )
}
