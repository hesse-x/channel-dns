#!/usr/bin/env python3
"""Paired full-timestep comparisons of executable variants on one MPI/thread layout."""
import argparse
import csv
import pathlib
import re
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", action="append", required=True, help="label=/absolute/executable")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--grid", type=int, default=128)
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--work-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    variants = [(v.split("=", 1)[0], pathlib.Path(v.split("=", 1)[1]).resolve()) for v in args.variant]
    if (min(args.ranks, args.threads, args.steps, args.repeats) < 1 or args.grid < 8 or args.grid % 2
            or len({label for label, _ in variants}) != len(variants)
            or any(not re.fullmatch(r"[A-Za-z0-9_-]+", label) for label, _ in variants)):
        parser.error("invalid layout, grid, repeats, or variant labels")
    root = pathlib.Path(__file__).resolve().parents[2]
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    template = (root / "examples/config.toml").read_text()

    def configure(directory, steps, output=False):
        values = dict(NX=args.grid, NY=args.grid, NZ=args.grid, NX2=args.grid, NY2=args.grid,
                      NZ2=args.grid, NXS=args.grid//2, NYS=args.grid-2, NZS=args.grid//2,
                      INDEX_Serial=int(args.ranks == 1), NSTEPS=steps,
                      NDUMFRQ=steps if output else args.steps+1,
                      NPRNFRQ=steps if output else args.steps+1,
                      EPS_init="0.0", EPS_disturb="0.0", NCOUNT_disturb=0)
        text = template.replace('format = "ascii"', 'format = "binary"')
        for key, value in values.items():
            text, count = re.subn(rf"^{key} = .*?$", f"{key} = {value}", text, flags=re.M)
            assert count == 1
        (directory / "config.toml").write_text(text)

    def run(label, executable):
        result = subprocess.run(["mpiexec", "--bind-to", "core", "--map-by", f"slot:PE={args.threads}",
                                 "-n", str(args.ranks), str(executable), "--threads", str(args.threads), "--timing", str(work/label)],
                                capture_output=True, text=True, check=True)
        matches = re.findall(r"timestep wall time = ([\d.eE+-]+) s", result.stdout)
        return float(matches[-1]), result.stdout + result.stderr

    for label, executable in variants:
        directory = work/label
        directory.mkdir()
        configure(directory, 0)
        print(f"Warm {label}", flush=True)
        _, log = run(label, executable)
        (directory / "warm.log").write_text(log)
        configure(directory, args.steps)
    rows = []
    for repeat in range(args.repeats):
        for label, executable in (variants if repeat % 2 == 0 else variants[::-1]):
            elapsed, log = run(label, executable)
            (work / label / f"run-{repeat}.log").write_text(log)
            rows.append(dict(variant=label, grid=args.grid, steps=args.steps, ranks=args.ranks,
                             threads=args.threads, repeat=repeat, timestep_seconds=elapsed))
            print(f"{label} repeat={repeat}: {elapsed:.6f}s", flush=True)
    with (work / "results.csv").open("w") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for label, _ in variants:
        print(f"median {label}: {statistics.median(r['timestep_seconds'] for r in rows if r['variant']==label):.6f}s", flush=True)
    for label, executable in variants:
        configure(work/label, min(args.steps, 10), output=True)
        _, log = run(label, executable)
        (work / label / "precision.log").write_text(log)
    for label, _ in variants[1:]:
        subprocess.run(["python3", str(root/"tests/tools/compare_fields.py"),
                        str(work/variants[0][0]), str(work/label), "--tolerance", "1e-10"], check=True)
    print(f"Results: {work}")


if __name__ == "__main__":
    main()
