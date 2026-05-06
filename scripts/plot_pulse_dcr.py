#!/usr/bin/env python3
"""Plot ParticleSimple pulse_dcr CSV output.

Produces two figures:
  pulse_dcr_timeseries.png  - full-run time-series overview
  pulse_dcr_overlay.png     - pulse 1 vs pulse 5 vs pulse 10 V_bulk overlay
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def smooth(y: np.ndarray, window: int = 50) -> np.ndarray:
    return (
        pd.Series(y).rolling(window=window, min_periods=1, center=True).mean().values
    )


def shade_pulses(ax, n_steps: int, relax: int, pulse: int) -> None:
    cycle = relax + pulse
    for cycle_start in range(0, n_steps, cycle):
        ps = cycle_start + relax
        pe = ps + pulse
        if ps < n_steps:
            ax.axvspan(ps, min(pe, n_steps), alpha=0.10, color="red", linewidth=0)


def plot_timeseries(df: pd.DataFrame, relax: int, pulse: int, out: Path) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    n = len(df)
    for ax in axes:
        shade_pulses(ax, n, relax, pulse)

    # Voltage
    ax = axes[0]
    ax.plot(df["step"], df["voltage_applied"], color="black", lw=0.6, label="V applied")
    ax.plot(df["step"], df["voltage_bulk"], color="C0", lw=0.6, label="V bulk (slab)")
    ax.set_ylabel("Voltage")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_title("pulse_dcr: 10 voltage-step cycles, BV deposition + SEI active")

    # Current (raw + smoothed)
    ax = axes[1]
    ax.plot(df["step"], df["current"], color="gray", lw=0.3, alpha=0.5, label="raw")
    ax.plot(
        df["step"],
        smooth(df["current"].values, window=50),
        color="C3",
        lw=1.2,
        label="smoothed (50)",
    )
    ax.set_ylabel("Boundary current\n(strip_L − plate_L)/dt")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)

    # Metal + SEI
    ax = axes[2]
    ax.plot(df["step"], df["metal_count"], color="C1", lw=1.0, label="metal_count")
    ax.set_ylabel("Metal count", color="C1")
    ax.tick_params(axis="y", labelcolor="C1")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(
        df["step"], df["sei_fraction"], color="C2", lw=1.0, label="sei_fraction"
    )
    ax2.set_ylabel("SEI fraction", color="C2")
    ax2.tick_params(axis="y", labelcolor="C2")

    # Cumulative reaction events
    ax = axes[3]
    ax.plot(df["step"], df["plate_right"].cumsum(), color="C3", lw=1.0, label="plate right")
    ax.plot(df["step"], df["strip_left"].cumsum(), color="C0", lw=1.0, label="strip left")
    ax.plot(
        df["step"],
        df["plate_left"].cumsum(),
        color="lightcoral",
        lw=0.8,
        ls="--",
        label="plate left",
    )
    ax.plot(
        df["step"],
        df["strip_right"].cumsum(),
        color="lightblue",
        lw=0.8,
        ls="--",
        label="strip right",
    )
    ax.set_ylabel("Cumulative events")
    ax.set_xlabel("Step (sim time / dt)")
    ax.legend(loc="upper left", fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_pulse_overlay(
    df: pd.DataFrame, relax: int, pulse: int, out: Path
) -> None:
    cycle = relax + pulse
    n = len(df)
    n_cycles = n // cycle
    pulses_to_show = [1, n_cycles // 2, n_cycles]
    pulses_to_show = [p for p in pulses_to_show if 1 <= p <= n_cycles]
    pulses_to_show = sorted(set(pulses_to_show))

    fig, (ax_v, ax_i) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(pulses_to_show)))

    for i, p in enumerate(pulses_to_show):
        ps = (p - 1) * cycle + relax
        pe = ps + pulse
        seg = df[(df["step"] >= ps) & (df["step"] < pe)].copy()
        seg["t_rel"] = seg["step"] - ps
        ax_v.plot(
            seg["t_rel"],
            seg["voltage_bulk"],
            color=colors[i],
            lw=0.7,
            label=f"pulse {p}",
        )
        ax_i.plot(
            seg["t_rel"],
            smooth(seg["current"].values, window=30),
            color=colors[i],
            lw=1.0,
            label=f"pulse {p}",
        )

    ax_v.axhline(0.1, color="black", ls=":", lw=0.6, label="V applied")
    ax_v.set_ylabel("V_bulk")
    ax_v.legend(loc="upper right", fontsize=9)
    ax_v.grid(True, alpha=0.3)
    ax_v.set_title("Pulse comparison: V_bulk and current evolve as morphology + SEI grow")

    ax_i.set_xlabel("Step from pulse start")
    ax_i.set_ylabel("Boundary current (smoothed)")
    ax_i.legend(loc="upper right", fontsize=9)
    ax_i.grid(True, alpha=0.3)
    ax_i.set_xlim(0, pulse)

    plt.tight_layout()
    plt.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="pulse_dcr.csv")
    p.add_argument("--relax", type=int, default=500)
    p.add_argument("--pulse", type=int, default=1000)
    p.add_argument("--out-dir", default=".")
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    out_dir = Path(args.out_dir)
    plot_timeseries(df, args.relax, args.pulse, out_dir / "pulse_dcr_timeseries.png")
    plot_pulse_overlay(df, args.relax, args.pulse, out_dir / "pulse_dcr_overlay.png")


if __name__ == "__main__":
    main()
