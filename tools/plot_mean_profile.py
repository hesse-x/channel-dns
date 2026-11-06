#!/usr/bin/env python3
"""Plot the single-rank channel mean profile in friction units."""
import argparse
import json
import os
from pathlib import Path
import tomllib

os.environ.setdefault("MPLCONFIGDIR", "/tmp/channel-dns-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(run_dir: Path, expected_samples=None, sample_start_step=None, sample_end_step=None,
         reference=None, separate_walls=False, caption=None):
    config = tomllib.loads((run_dir / "config.toml").read_text())
    if config["parm"]["INDEX_GRAD_PRE"] != 1:
        raise ValueError("this plot expects fixed pressure-gradient forcing")
    data = np.loadtxt(run_dir / "Result/1/V_Pre_mean_Field_ 0.dat")
    gradient = (run_dir / "Result/0/DPDX_DPDZ.dat").read_text().split()
    samples = int(gradient[0])
    if expected_samples is not None and samples != expected_samples:
        raise ValueError(f"expected {expected_samples} samples, found {samples}")
    if not np.isfinite(data).all():
        raise ValueError("non-finite mean field")
    dpdx, dpdz = map(float, gradient[1:3])
    utau = (dpdx**2 + dpdz**2)**0.25  # h=1, rho=1
    reynolds = config["parm"]["Re"]
    retau = reynolds * utau
    y, u = data[:, 0], data[:, 1]
    n = len(y) // 2 + 1
    distance = 1.0 - np.abs(y[:n])
    # In a developing flow, prescribed pressure gradient does not equal wall drag.
    # Use the measured mean wall shear for + units, and retain forcing units separately.
    coefficients = np.polynomial.chebyshev.chebfit(y, u, len(y) - 1)
    wall_derivative = np.polynomial.chebyshev.chebval(
        [-1., 1.], np.polynomial.chebyshev.chebder(coefficients))
    wall_utau = np.sqrt(np.abs(wall_derivative) / reynolds)
    mean_wall_utau = float(np.sqrt(np.mean(np.abs(wall_derivative)) / reynolds))
    if not np.isfinite(mean_wall_utau) or mean_wall_utau <= 0:
        raise ValueError("invalid mean wall shear for wall-unit normalization")
    wall_retau = reynolds * mean_wall_utau
    reference_data = np.loadtxt(reference, comments="#") if reference is not None else None
    reference_retau = float(reference_data[-1, 1]) if reference_data is not None else None
    for suffix, scale, actual_wall_units in [("", mean_wall_utau, True),
                                            ("_pressure_units", utau, False)]:
        yp = distance * reynolds * scale
        upper_scale = wall_utau[1] if separate_walls and actual_wall_units else scale
        lower_scale = wall_utau[0] if separate_walls and actual_wall_units else scale
        upper_y = distance * reynolds * upper_scale
        lower_y = distance * reynolds * lower_scale
        upper, lower = u[:n] / upper_scale, u[::-1][:n] / lower_scale
        symmetric = (upper + lower) / 2
        if separate_walls:
            columns = np.column_stack([upper_y, upper, lower_y, lower])
            header = ("y_plus_upper,u_plus_upper,y_plus_lower,u_plus_lower" if actual_wall_units
                      else "y_forcing_upper,u_forcing_upper,y_forcing_lower,u_forcing_lower")
        else:
            columns = np.column_stack([yp, upper, lower, symmetric])
            header = ("y_plus,u_plus_upper,u_plus_lower,u_plus_symmetric" if actual_wall_units
                      else "y_forcing,u_forcing_upper,u_forcing_lower,u_forcing_symmetric")
        np.savetxt(run_dir / f"mean_profile{suffix}.csv", columns, delimiter=",",
                   header=header, comments="")
        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        fig.subplots_adjust(left=.11, right=.98, bottom=.18, top=.87)
        upper_label, lower_label = "Upper wall", "Lower wall"
        if separate_walls and actual_wall_units:
            upper_label += rf", $Re_\tau={wall_utau[1]*reynolds:.2f}$"
            lower_label += rf", $Re_\tau={wall_utau[0]*reynolds:.2f}$"
        ax.semilogx(upper_y[1:], upper[1:], "o-" if separate_walls else "--",
                    label=upper_label, markersize=2.5, alpha=.85)
        ax.semilogx(lower_y[1:], lower[1:], "--" if separate_walls else ":",
                    label=lower_label, alpha=.9)
        if not separate_walls:
            ax.semilogx(yp[1:], symmetric[1:], "o-", markersize=3,
                        color="black", label="Both-wall average")
        if actual_wall_units:
            if reference_data is not None:
                ax.semilogx(reference_data[1:, 1], reference_data[1:, 2], color="#cf3333",
                            linewidth=1.8, label=rf"MKM 1999, $Re_\tau={reference_retau:.2f}$")
            viscous = np.geomspace(max(yp[1], .1), 10, 100)
            ax.semilogx(viscous, viscous, color="gray", label=r"$u^+=y^+$")
            ax.set(xlabel=r"$y^+$ (mean wall shear)", ylabel=r"$\langle u\rangle^+$",
                   title=(rf"Separate wall normalization; forcing $Re_\tau={retau:.0f}$" if separate_walls
                          else rf"Sampled mean: $Re_{{\tau,wall}}={wall_retau:.1f}$, forcing $Re_\tau={retau:.0f}$"))
        else:
            ax.set(xlabel=r"$d u_{\tau,forcing}/\nu$", ylabel=r"$\langle u\rangle/u_{\tau,forcing}$",
                   title=rf"Pressure-forcing units: $Re_{{\tau,forcing}}={retau:.0f}$, wall $Re_\tau={wall_retau:.1f}$")
        ax.grid(True, which="both", alpha=.25)
        ax.legend()
        window = (f"Steps {sample_start_step:,}–{sample_end_step:,}; "
                  if sample_start_step is not None and sample_end_step is not None else "")
        footer = caption if caption is not None else (
            f"{window}{samples:,} samples\nStatistical stationarity has not been established")
        fig.text(.5, .02, footer,
                 ha="center", va="bottom", fontsize=8)
        fig.savefig(run_dir / f"mean_profile{suffix}.png", dpi=180)
        fig.savefig(run_dir / f"mean_profile{suffix}.svg")
        plt.close(fig)
    orders = np.arange(0, len(coefficients), 2)
    bulk_u = float(np.sum(coefficients[::2] / (1-orders**2)))
    report = {"samples": samples, "Re_tau_pressure": retau,
              "sample_start_step": sample_start_step, "sample_end_step": sample_end_step,
              "u_tau_pressure": utau,
              "Re_tau_wall_mean_profile": (wall_utau * reynolds).tolist(),
              "Re_tau_wall_symmetric": wall_retau,
              "u_tau_wall_symmetric": mean_wall_utau,
              "centerline_u_plus": float(u[n-1] / mean_wall_utau),
              "centerline_u_forcing_units": float(u[n-1] / utau),
              "bulk_u": bulk_u,
              "mean_wall_stress_over_pressure_forcing": (mean_wall_utau / utau)**2,
              "normalization": "mean_profile wall shear, averaged over both walls",
              "note": "Finite-duration statistics; stationarity and grid convergence are not established."}
    if separate_walls:
        for key in ["Re_tau_wall_symmetric", "u_tau_wall_symmetric", "centerline_u_plus",
                    "mean_wall_stress_over_pressure_forcing"]:
            report.pop(key)
        report.update(normalization="each wall uses its own mean wall shear; no wall averaging",
                      u_tau_wall_lower_upper=wall_utau.tolist(),
                      centerline_u_plus_lower_upper=(u[n-1]/wall_utau).tolist())
    if caption is not None:
        report['figure_caption'] = caption
    if reference_data is not None:
        yp = distance * reynolds * mean_wall_utau
        symmetric = (u[:n] + u[::-1][:n]) / (2 * mean_wall_utau)
        overlap = (yp > 0) & (yp <= reference_data[-1, 1])
        ref_u = np.interp(yp[overlap], reference_data[:, 1], reference_data[:, 2])
        if separate_walls:
            rows = []
            for side, values, wall_scale in [("upper", u[:n], wall_utau[1]),
                                             ("lower", u[::-1][:n], wall_utau[0])]:
                wall_y = distance * reynolds * wall_scale
                wall_u = values / wall_scale
                valid = (wall_y > 0) & (wall_y <= reference_data[-1, 1])
                target_u = np.interp(wall_y[valid], reference_data[:, 1], reference_data[:, 2])
                rows.extend((side, x, value, target, value-target)
                            for x, value, target in zip(wall_y[valid], wall_u[valid], target_u))
            np.savetxt(run_dir / "mkm_comparison.csv", np.asarray(rows, dtype=object),
                       delimiter=",", fmt="%s",
                       header="wall,y_plus,u_plus_simulation,u_plus_MKM,difference", comments="")
        else:
            np.savetxt(run_dir / "mkm_comparison.csv",
                       np.column_stack([yp[overlap], symmetric[overlap], ref_u,
                                        symmetric[overlap]-ref_u]), delimiter=",",
                       header="y_plus,u_plus_simulation,u_plus_MKM,difference", comments="")
        report["reference"] = {"file": str(Path(reference).resolve()),
                               "Re_tau": reference_retau,
                               "centerline_u_plus": float(reference_data[-1, 2]),
                               "note": "Different wall Reynolds numbers; differences are not a discretization error estimate."}
    (run_dir / "profile_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--expected-samples", type=int)
    parser.add_argument("--sample-start-step", type=int)
    parser.add_argument("--sample-end-step", type=int)
    parser.add_argument("--reference", type=Path, help="MKM mean table: columns y/h, y+, U+")
    parser.add_argument("--separate-walls", action="store_true", help="normalize each wall independently; no averaged curve")
    parser.add_argument("--caption", help="custom figure description; use actual run provenance")
    args = parser.parse_args()
    print(json.dumps(plot(args.run_dir, args.expected_samples,
                          args.sample_start_step, args.sample_end_step, args.reference,
                          args.separate_walls, args.caption), indent=2))
