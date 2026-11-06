#!/usr/bin/env python3
"""Compare text or binary channel-dns physical fields and EKT."""

import argparse
import math
import pathlib
import re
import struct
import sys
import tomllib

NUMBER = re.compile(
    r"[-+]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][-+]?\d+)?|nan|inf(?:inity)?)",
    re.IGNORECASE,
)
HEADER = struct.Struct("<8sIIqd12q")
MAGIC = b"CHDNS\0\0\0"


def nums(text):
    return [float(v.replace("D", "E").replace("d", "e")) for v in NUMBER.findall(text)]


def mesh_size(directory):
    if (directory / "config.toml").exists():
        mesh = tomllib.loads((directory / "config.toml").read_text())["mesh"]
        return mesh["NX"], mesh["NY"] + 1, mesh["NZ"]
    lines = (directory / "meshparm_channel_with_Microblow.txt").read_text().splitlines()
    values = [int(v) for line in lines[1::2] for v in NUMBER.findall(line)]
    return values[0], values[1] + 1, values[2]


def find_one(directory, pattern):
    files = sorted(directory.glob(pattern))
    if len(files) != 1:
        raise ValueError(f"expected one {pattern} in {directory}, found {len(files)}")
    return files[0]


def read_text(directory):
    nx, global_nyp, nz = mesh_size(directory)
    paths = sorted(directory.glob("field_P_*.dat"))
    if not paths:
        raise ValueError(f"no text physical field in {directory}")
    fields = [[0.0] * (nx * global_nyp * nz) for _ in range(4)]
    blocks = []
    checkpoint = None
    ranks = set()
    for path in paths:
        match = re.fullmatch(r"field_P_\s*(\d+)\.dat", path.name)
        if not match:
            raise ValueError(f"cannot determine text rank from {path.name}")
        rank = int(match.group(1))
        if rank in ranks:
            raise ValueError(f"duplicate text rank {rank}")
        values = nums(path.read_text())
        payload_count = len(values) - 2
        denominator = 4 * nx * nz
        if payload_count <= 0 or payload_count % denominator:
            raise ValueError(f"text physical payload length mismatch: {path}")
        local_nyp = payload_count // denominator
        current = (int(values[0]), values[1])
        if checkpoint is None:
            checkpoint = current
        elif checkpoint != current:
            raise ValueError("mixed text physical checkpoint metadata")
        payload = values[2:]
        if any(not math.isfinite(value) for value in payload):
            raise ValueError(f"text physical payload contains non-finite values: {path}")
        blocks.append((rank, local_nyp, payload))
        ranks.add(rank)

    blocks.sort()
    if ranks != set(range(len(blocks))):
        raise ValueError("missing text physical rank file")
    counts = [local_nyp for _, local_nyp, _ in blocks]
    if sum(counts) != global_nyp:
        raise ValueError(f"text y slices {counts} do not cover {global_nyp} points")
    base, remainder = divmod(global_nyp, len(blocks))
    expected_counts = [base + (rank >= len(blocks) - remainder) for rank in range(len(blocks))]
    if counts != expected_counts:
        raise ValueError(f"unexpected text y decomposition: {counts} != {expected_counts}")

    ystart = 0
    for _, local_nyp, payload in blocks:
        local_count = nx * local_nyp * nz
        cursor = 0
        for iz in range(nz):
            for iy in range(local_nyp):
                for ix in range(nx):
                    logical = (iz * global_nyp + ystart + iy) * nx + ix
                    for component in range(3):
                        fields[component][logical] = payload[cursor]
                        cursor += 1
        pressure = payload[3 * local_count:]
        for iz in range(nz):
            for iy in range(local_nyp):
                src = (iz * local_nyp + iy) * nx
                dst = (iz * global_nyp + ystart + iy) * nx
                fields[3][dst:dst + nx] = pressure[src:src + nx]
        ystart += local_nyp
    return (*checkpoint, nx, global_nyp, nz), fields


def read_cpp(directory):
    paths = sorted(directory.glob("field_P_*.bin")) or sorted(directory.glob("field_P_*.dat"))
    if not paths:
        raise ValueError(f"no field_P_*.dat in {directory}")
    global_nx, global_nyp, global_nz = mesh_size(directory)
    fields = [[0.0] * (global_nx * global_nyp * global_nz) for _ in range(4)]
    covered = [False] * global_nyp
    ranks = set()
    checkpoint = None
    expected_nproc = None
    for path in paths:
        data = path.read_bytes()
        if len(data) < HEADER.size:
            raise ValueError(f"truncated C++ physical header: {path}")
        header = HEADER.unpack_from(data)
        (magic, version, kind, ntpit, time, nx, nyp, nz, components, zmodes, xmodes,
         rank, nproc, xstart, xcount, ystart, ycount) = header
        if magic != MAGIC or version != 1 or kind != 2 or components != 4:
            raise ValueError(f"invalid C++ physical header: {path}")
        if nx != global_nx or nz != global_nz or nyp != ycount or zmodes != 0 or xmodes != 0:
            raise ValueError(f"physical dimensions do not match mesh: {path}")
        if rank in ranks or rank < 0 or rank >= nproc:
            raise ValueError(f"duplicate or invalid physical rank: {rank}")
        if xstart != 0 or xcount != global_nx or ystart < 0 or ycount <= 0 or ystart + ycount > global_nyp:
            raise ValueError(f"invalid physical block metadata: {path}")
        if checkpoint is None:
            checkpoint = (ntpit, time)
            expected_nproc = nproc
        elif checkpoint != (ntpit, time) or expected_nproc != nproc:
            raise ValueError("mixed physical checkpoint metadata")
        local_count = nx * nyp * nz
        expected = HEADER.size + 4 * local_count * 8
        if len(data) != expected:
            raise ValueError(f"C++ physical payload length mismatch: {len(data)} != {expected}")
        values = struct.unpack_from(f"<{4 * local_count}d", data, HEADER.size)
        if any(not math.isfinite(value) for value in values):
            raise ValueError("C++ physical payload contains non-finite values")
        for iy in range(ystart, ystart + ycount):
            if covered[iy]:
                raise ValueError(f"overlapping physical y block at {iy}")
            covered[iy] = True
        for component in range(4):
            local = values[component * local_count:(component + 1) * local_count]
            for iz in range(nz):
                for iy in range(nyp):
                    src = (iz * nyp + iy) * nx
                    dst = (iz * global_nyp + ystart + iy) * global_nx
                    fields[component][dst:dst + nx] = local[src:src + nx]
        ranks.add(rank)
    if expected_nproc != len(paths) or ranks != set(range(expected_nproc)):
        raise ValueError("missing physical rank file")
    if not all(covered):
        raise ValueError("physical y blocks do not cover the mesh")
    return (*checkpoint, global_nx, global_nyp, global_nz), fields


def read_fields(directory):
    paths = sorted(directory.glob("field_P_*.bin")) or sorted(directory.glob("field_P_*.dat"))
    if not paths:
        raise ValueError(f"no physical field in {directory}")
    return read_cpp(directory) if paths[0].read_bytes()[:8] == MAGIC else read_text(directory)


def compare_field(name, reference, actual, shape, tolerance):
    if len(reference) != len(actual):
        raise ValueError(f"{name}: element count differs")
    errors = [abs(a - r) for r, a in zip(reference, actual)]
    index = max(range(len(errors)), key=errors.__getitem__)
    max_abs = errors[index]
    ref_norm = max(abs(value) for value in reference)
    relative = max_abs / max(ref_norm, 1e-300)
    nx, nyp, _ = shape
    iz, rem = divmod(index, nx * nyp)
    iy, ix = divmod(rem, nx)
    print(f"{name}: n={len(reference)} ref_max={ref_norm:.17e} max_abs={max_abs:.17e} "
          f"relative={relative:.17e} at (ix={ix},iy={iy},iz={iz})")
    return relative <= tolerance


def read_ekt(directory):
    values = nums(find_one(directory, "EKT*.dat").read_text())
    if len(values) % 2:
        raise ValueError("EKT has a partial trailing row")
    return list(zip(values[0::2], values[1::2]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=pathlib.Path)
    parser.add_argument("actual", type=pathlib.Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    reference_header, reference_fields = read_fields(args.reference)
    actual_header, actual_fields = read_fields(args.actual)
    if reference_header != actual_header:
        raise ValueError(f"physical headers differ: {reference_header} != {actual_header}")
    nx, nyp, nz = reference_header[2:]
    results = [compare_field(name, ref, got, (nx, nyp, nz), args.tolerance)
               for name, ref, got in zip(("u", "v", "w", "p"), reference_fields, actual_fields)]
    passed = all(results)
    reference_ekt, actual_ekt = read_ekt(args.reference), read_ekt(args.actual)
    if len(reference_ekt) != len(actual_ekt):
        raise ValueError("EKT row count differs")
    for i, ((rt, rv), (ct, cv)) in enumerate(zip(reference_ekt, actual_ekt), 1):
        if rt != ct:
            raise ValueError(f"EKT time differs at row {i}: {rt} != {ct}")
    passed &= compare_field("EKT", [v for _, v in reference_ekt],
                            [v for _, v in actual_ekt], (len(reference_ekt), 1, 1), args.tolerance)
    if not passed:
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, struct.error) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(2)
