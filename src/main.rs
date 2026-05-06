use std::env;

use particle_simple::scenarios;

fn main() {
    let args: Vec<String> = env::args().collect();
    let scenario = args
        .iter()
        .position(|a| a == "--scenario")
        .and_then(|i| args.get(i + 1))
        .map(String::as_str)
        .unwrap_or("empty_cell");

    match scenario {
        "empty_cell" => scenarios::empty_cell::run(),
        "pulse_dcr" => scenarios::pulse_dcr::run(),
        other => {
            eprintln!("unknown scenario: {other}");
            eprintln!("available: empty_cell, pulse_dcr");
            std::process::exit(1);
        }
    }
}
