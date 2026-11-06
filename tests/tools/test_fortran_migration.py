#!/usr/bin/env python3
"""Deterministic nonzero-mode checkpoints: rank invariance and restart validation."""
import argparse
import importlib.util
import math
import pathlib
import struct
import subprocess
import tempfile
import sys

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('compare', pathlib.Path(__file__).with_name('compare_fields.py'))
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


def partition(n, ranks):
    q, r = divmod(n, ranks)
    counts = [q + (rank >= ranks-r) for rank in range(ranks)]
    starts = [sum(counts[:rank]) for rank in range(ranks)]
    return starts, counts


def checkpoint(directory, ranks, config, nxs=6, nzs=6, nys=12):
    directory.mkdir()
    (directory / 'config.toml').write_text(config.replace('INDEX_Serial = 1', f'INDEX_Serial = {int(ranks == 1)}'))
    xb, xc = partition(nxs+1, ranks)
    _, yc = partition(17, ranks)
    for rank in range(ranks):
        header = compare.HEADER.pack(compare.MAGIC, 1, 1, 0, 0., 16, yc[rank], 16, 3,
                                     2*nzs+1, xc[rank], rank, ranks, xb[rank], xc[rank], 0, nys+1)
        payload = []
        for z in range(-nzs, nzs+1):
            for y in range(nys+1):
                for x in range(xb[rank], xb[rank]+xc[rank]):
                    for c in range(3):
                        value = 0j
                        if x == z == c == 0:
                            if y == 0:
                                value = .75
                            elif y == 2:
                                value = -.75
                        if x == z == y == 0 and c == 1:
                            value = .001
                        # Both signs, retained high modes, and complex phases exercise all transforms.
                        if x in (1, nxs) and abs(z) in (1, nzs) and y <= 4:
                            value += complex(math.sin(13*x + 7*z + 3*y + c), math.cos(5*x - z + y + 11*c)) * 1.e-5
                        payload.extend((value.real, value.imag))
        (directory / f'field_S_{rank:2d}.bin').write_bytes(header + struct.pack(f'<{len(payload)}d', *payload))


def compare_runs(reference, actual):
    rh, rf = compare.read_fields(reference)
    ah, af = compare.read_fields(actual)
    assert rh == ah, (rh, ah)
    worst = 0.
    for name, left, right in zip(('u', 'v', 'w', 'p'), rf, af):
        error = max(abs(x-y) for x, y in zip(left, right))
        scale = max(abs(x) for x in left)
        assert error <= 2.e-11 + 2.e-10*scale, (name, error, scale)
        worst = max(worst, error)
    # Plane statistics must also agree after the same evolution.
    for path in (reference / 'Result').glob('*/*.dat'):
        if path.name == 'DPDX_DPDZ.dat' or len(list((reference / 'Result/1').glob('*.dat'))) == 1 == len(list((actual / 'Result/1').glob('*.dat'))):
            other = actual / path.relative_to(reference)
            left, right = compare.nums(path.read_text()), compare.nums(other.read_text())
            assert len(left) == len(right)
            assert max(abs(x-y) for x,y in zip(left,right)) <= 2.e-8
    return worst


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=pathlib.Path)
    parser.add_argument('source', type=pathlib.Path)
    parser.add_argument('--mpiexec', required=True)
    args = parser.parse_args()
    base = (args.source / 'examples/config.toml').read_text()
    replacements = {
        'NX = 32': 'NX = 16', 'NY = 32': 'NY = 16', 'NZ = 32': 'NZ = 16',
        'NXS = 16': 'NXS = 6', 'NYS = 30': 'NYS = 12', 'NZS = 16': 'NZS = 6',
        'NX2 = 32': 'NX2 = 24', 'NY2 = 32': 'NY2 = 24', 'NZ2 = 32': 'NZ2 = 24',
        'NSTEPS = 5': 'NSTEPS = 3', 'INDEX_CONT = 0': 'INDEX_CONT = 1',
        'Nstat = 50000': 'Nstat = 0', 'format = "ascii"': 'format = "both"',
    }
    for old,new in replacements.items():
        base = base.replace(old,new)
    def run(executable, directory, ranks, success=True):
        proc = subprocess.run([args.mpiexec, '-n', str(ranks), str(executable), str(directory)],
                              text=True, capture_output=True, timeout=60)
        if success and proc.returncode:
            raise AssertionError(proc.stdout + proc.stderr)
        if not success:
            assert proc.returncode != 0, 'invalid input unexpectedly succeeded'
    with tempfile.TemporaryDirectory(prefix='dns-fortran-regression-') as tmp:
        root = pathlib.Path(tmp)
        for case in ('dealiased', 'full', 'two-dimensional', 'pressure-gradient'):
            config = base
            nxs, nzs = 6,6
            if case == 'full':
                config = config.replace('NXS = 6','NXS = 8').replace('NZS = 6','NZS = 8')
                config = config.replace('NX2 = 24','NX2 = 16').replace('NZ2 = 24','NZ2 = 16')
                nxs,nzs = 8,8
            elif case == 'two-dimensional':
                config = config.replace('NZS = 6','NZS = 0')
                nzs = 0
            elif case == 'pressure-gradient':
                config = config.replace('INDEX_GRAD_PRE = 0','INDEX_GRAD_PRE = 1').replace('THETA = 0.0e0','THETA = 30.0e0')
            serial = root / f'{case}-np1'
            checkpoint(serial, 1, config, nxs, nzs)
            run(args.executable, serial, 1)
            for ranks in (2,3,4):
                parallel = root / f'{case}-np{ranks}'
                checkpoint(parallel, ranks, config, nxs, nzs)
                run(args.executable, parallel, ranks)
                print(case, f'np1/np{ranks}', 'max_abs', compare_runs(serial, parallel))
        # Invalid checkpoints must fail before progressing into a collective transform.
        for failure in ('truncated','mixed-time','nan'):
            directory = root / failure
            checkpoint(directory, 2, base)
            path = directory / 'field_S_ 1.bin'
            data = bytearray(path.read_bytes())
            if failure == 'truncated':
                data = data[:-8]
            elif failure == 'mixed-time':
                struct.pack_into('<d', data, 24, .1)
            else:
                struct.pack_into('<d', data, 128, math.nan)
            path.write_bytes(data)
            run(args.executable, directory, 2, success=False)
        for name,old,new in (
            ('fractional-grid','NX = 16','NX = 16.5'),
            ('float-grid','NX = 16','NX = 16.0'),
            ('string-grid','NX = 16','NX = \"16\"'),
            ('boolean-real','XL = 1.0e1','XL = true'),
            ('integer-overflow','NSTEPS = 3','NSTEPS = 2147483648'),
            ('missing-key','NX = 16','# omitted NX'),
            ('unknown-key','NX = 16','NX = 16\nunexpected = 1'),
            ('unknown-table','[mesh]','[unused]'),
            ('duplicate-table-case','[parm]','[parm]\n[PARM]'),
            ('extra-number','NX = 16','NX = 16, 32'),
            ('space-separated-number','NX = 16','NX = 16 32'),
            ('list-directed-terminator','NX = 16','NX = 16/'),
            ('missing-exponent','DT = 3.0e-3','DT = 3.0e'),
            ('hex-integer','NX = 16','NX = 0x10'),
            ('numeric-separator','NX = 16','NX = 1_6'),
            ('nested-table','[mesh]','[mesh.nested]'),
            ('quoted-key','NX = 16','\"NX\" = 16'),
            ('multiline-array','a_V = [1.0e0, 1.0e0, 1.0e0]','a_V = [\n1, 1, 1]'),
            ('array-trailing-comma','a_V = [1.0e0, 1.0e0, 1.0e0]','a_V = [1, 1, 1,]'),
            ('array-null-item','a_V = [1.0e0, 1.0e0, 1.0e0]','a_V = [1,,1]'),
            ('unterminated-string','format = \"both\"','format = \"both'),
            ('escaped-string','format = "both"',r'format = "bi\u006eary"'),
            ('long-line','NX = 16','NX = 16' + ' '*4100 + 'oops'),
            ('bad-format','format = "both"','format = "binary-extra"'),
            ('array-short','a_V = [1.0e0, 1.0e0, 1.0e0]','a_V = [1, 1]'),
            ('array-long','a_V = [1.0e0, 1.0e0, 1.0e0]','a_V = [1, 1, 1, 1]'),
            ('nonfinite','DT = 3.0e-3','DT = nan'),
            ('duplicate-key','NX = 16','NX = 16\nnx = 16'),
        ):
            directory = root / name
            directory.mkdir()
            text = base.replace(old,new)
            assert text != base, name
            (directory / 'config.toml').write_text(text)
            run(args.executable, directory, 1, success=False)
        print('PASS: nonzero modes, rank invariance, zero mode, restart validation')


if __name__ == '__main__':
    main()
