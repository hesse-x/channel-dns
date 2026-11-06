#!/usr/bin/env python3
"""Check worker selection, invalid CLI arguments and numerical equivalence."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('source', type=Path)
    parser.add_argument('--pure-mpi', action='store_true')
    args = parser.parse_args()
    executable = str(args.executable.resolve())
    with tempfile.TemporaryDirectory(prefix='channel-thread-args-') as temporary:
        root = Path(temporary)
        template = (args.source / 'examples/config.toml').read_text()
        for label, flags, expected in [
            ('default', [], os.cpu_count() or 1),
            ('one', ['--threads', '1'], 1),
            ('eight', ['-t', '8'], 8),
        ]:
            directory = root / label
            directory.mkdir()
            (directory / 'config.toml').write_text(template)
            result = subprocess.run([executable, *flags, str(directory)],
                                    capture_output=True, text=True, check=True)
            if args.pure_mpi:
                assert 'inline computation' in result.stdout, result.stdout
            else:
                match = re.search(r'(\d+) compute worker\(s\)', result.stdout)
                assert match and int(match[1]) == expected, result.stdout
        for flags in [['--threads'], ['--threads', '0'], ['--threads', '-1'],
                      ['--threads', 'abc'], ['--threads', '2.5'],
                      ['--threads', '999999999999999999999999']]:
            result = subprocess.run([executable, *flags], capture_output=True, text=True)
            assert result.returncode == 2, (flags, result)
            assert 'positive integer' in result.stderr, (flags, result.stderr)
        for label in ['default', 'eight']:
            subprocess.run(['python3', str(args.source / 'tests/tools/compare_fields.py'),
                            str(root / 'one'), str(root / label), '--tolerance', '1e-10'],
                           check=True)
    print('Thread arguments: PASS')


if __name__ == '__main__':
    main()
