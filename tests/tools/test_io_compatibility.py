#!/usr/bin/env python3
"""End-to-end TOML, text, binary, restart and MPI compatibility checks."""
import argparse
import importlib.util
import pathlib
import shutil
import struct
import subprocess
import tempfile
import sys

sys.dont_write_bytecode = True

spec = importlib.util.spec_from_file_location("compare", pathlib.Path(__file__).with_name("compare_fields.py"))
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=pathlib.Path)
    parser.add_argument("source", type=pathlib.Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    config = (args.source / "examples/config.toml").read_text().replace('format = "ascii"', 'format = "both"')
    config = config.replace("NSTEPS = 5", "NSTEPS = 2").replace("Nstat = 50000", "Nstat = 0")
    prefix = [args.mpiexec, "-n", "2"] if args.mpiexec else []
    if prefix:
        config = config.replace("INDEX_Serial = 1", "INDEX_Serial = 0")

    def run(directory, text):
        (directory / "config.toml").write_text(text)
        subprocess.run(prefix + [str(args.executable), "--threads", str(args.threads), str(directory)], check=True, capture_output=True, text=True)

    with tempfile.TemporaryDirectory(prefix="channel-io-integration-") as temporary:
        root = pathlib.Path(temporary)
        both = root / "both"
        both.mkdir()
        run(both, config)
        text_header, text_fields = compare.read_text(both)
        binary_header, binary_fields = compare.read_cpp(both)
        assert text_header == binary_header
        assert text_fields == binary_fields, "physical ordering or round-trip precision mismatch"
        for path in both.glob("field_S_*.bin"):
            data = path.read_bytes()
            header = compare.HEADER.unpack_from(data)
            values = compare.nums(path.with_suffix(".dat").read_text())
            assert values[:2] == [header[3], header[4]]
            assert values[2:] == list(struct.unpack_from(f"<{(len(data)-128)//8}d", data, 128))

        # Each format resumes the identical checkpoint; physical output must agree.
        continuation = config.replace("INDEX_CONT = 0", "INDEX_CONT = 1").replace("NSTEPS = 2", "NSTEPS = 1")
        for format_name in ("ascii", "binary"):
            directory = root / format_name
            shutil.copytree(both, directory)
            run(directory, continuation.replace('format = "both"', f'format = "{format_name}"'))
        th, tf = compare.read_text(root / "ascii")
        bh, bf = compare.read_cpp(root / "binary")
        assert th == bh and tf == bf, "text and binary restart diverged"

        if not prefix:
            # The time-dependent boundary update branch must remain usable via TOML.
            run(root / "ascii", config.replace("INDEX_UNSTEADY_BC = 0", "INDEX_UNSTEADY_BC = 1"))
        print("TOML, text/binary layouts and continuation agree" + (" (MPI)" if prefix else ""))


if __name__ == "__main__":
    main()
