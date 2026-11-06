#!/usr/bin/env python3
"""Sequential, fixed-core-budget comparisons of MPI and singleton-thread-pool layouts."""
import argparse
import csv
import pathlib
import re
import statistics
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=pathlib.Path)
    parser.add_argument("--pure-mpi-executable", type=pathlib.Path,
                        help="use a pure MPI build for layouts with one thread per rank")
    parser.add_argument("--grid", type=int, default=64)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--layouts", default="8x1,4x2,2x4,1x8")
    parser.add_argument("--serial-baseline", action="store_true",
                        help="also measure 1 rank x 1 thread as a separate speedup baseline")
    parser.add_argument("--work-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    layouts = [tuple(map(int, layout.split("x"))) for layout in args.layouts.split(",")]
    if (args.grid < 8 or args.grid % 2 or args.steps < 1 or args.repeats < 1
            or any(r < 1 or t < 1 for r, t in layouts)
            or len({r * t for r, t in layouts}) != 1):
        parser.error("use an even grid >=8, positive steps/repeats and equal core budgets")
    if args.serial_baseline and (1, 1) not in layouts:
        layouts = [(1, 1)] + layouts
    executable = args.executable.resolve()
    pure_mpi = args.pure_mpi_executable.resolve() if args.pure_mpi_executable else None
    root = pathlib.Path(__file__).resolve().parents[2]
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    template = (root / "examples/config.toml").read_text()

    def config(ranks, steps, output=False):
        n = args.grid
        values = dict(NX=n, NY=n, NZ=n, NXS=n//2, NYS=n-2, NZS=n//2,
                      NX2=n, NY2=n, NZ2=n, INDEX_Serial=int(ranks == 1),
                      NSTEPS=steps, NDUMFRQ=steps if output else args.steps+1,
                      NPRNFRQ=steps if output else args.steps+1,
                      EPS_init="0.0", EPS_disturb="0.0", NCOUNT_disturb=0)
        text = template.replace('format = "ascii"', 'format = "binary"')
        for key, value in values.items():
            text, count = re.subn(rf"^{key} = .*?$", f"{key} = {value}", text, flags=re.M)
            assert count == 1, key
        return text

    def run(ranks, threads, directory, timing=True):
        binary = pure_mpi if pure_mpi and threads == 1 else executable
        command = ["mpiexec", "--bind-to", "core", "--map-by", f"slot:PE={threads}",
                   "-n", str(ranks), str(binary), "--threads", str(threads)]
        if timing:
            command += ["--timing"]
        command += [str(directory)]
        start = time.perf_counter()
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        wall = time.perf_counter() - start
        matches = re.findall(r"timestep wall time = ([\d.eE+-]+) s", result.stdout)
        return float(matches[-1]) if timing else 0, wall, result.stdout + result.stderr

    samples = {layout: [] for layout in layouts}
    for ranks, threads in layouts:
        directory = work / f"np{ranks}-t{threads}"
        directory.mkdir()
        (directory / "config.toml").write_text(config(ranks, 0))
        print(f"Warm {ranks}x{threads}", flush=True)
        _, _, log = run(ranks, threads, directory)
        (directory / "warm.log").write_text(log)
        (directory / "config.toml").write_text(config(ranks, args.steps))
    for repeat in range(args.repeats):
        for ranks, threads in (layouts if repeat % 2 == 0 else layouts[::-1]):
            directory = work / f"np{ranks}-t{threads}"
            step, wall, log = run(ranks, threads, directory)
            samples[ranks, threads].append((step, wall))
            (directory / f"run-{repeat}.log").write_text(log)
            print(f"{ranks}x{threads} repeat={repeat} steps={step:.6f}s process={wall:.6f}s", flush=True)
    with (work / "results.csv").open("w") as output:
        writer = csv.writer(output)
        writer.writerow(["grid", "steps", "ranks", "threads", "repeat", "timestep_seconds", "process_seconds"])
        for (ranks, threads), times in samples.items():
            for repeat, (step, wall) in enumerate(times):
                writer.writerow([args.grid, args.steps, ranks, threads, repeat, step, wall])
    for (ranks, threads), times in samples.items():
        print(f"median {ranks}x{threads}: {statistics.median(t[0] for t in times):.6f}s", flush=True)
    # Compare identical completed steps across decompositions and thread counts.
    precision_steps = min(args.steps, 10)
    directories = []
    for ranks, threads in layouts:
        directory = work / f"np{ranks}-t{threads}"
        (directory / "config.toml").write_text(config(ranks, precision_steps, output=True))
        _, _, log = run(ranks, threads, directory)
        (directory / "precision.log").write_text(log)
        directories.append(directory)
    for actual in directories[1:]:
        subprocess.run(["python3", str(root / "tests/tools/compare_fields.py"),
                        str(directories[0]), str(actual), "--tolerance", "1e-10"], check=True)
    print(f"Results and logs: {work}")


if __name__ == "__main__":
    main()
