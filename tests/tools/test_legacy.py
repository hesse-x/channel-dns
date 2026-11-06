#!/usr/bin/env python3
"""Check preserved source hashes and the FFTPACK/FFTW numerical parity."""
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile

modern,legacy,aligned,launcher,source = sys.argv[1:]
root=pathlib.Path(source)
manifest=json.loads((root/'legacy/SOURCE_MANIFEST.json').read_text())
for name,expected in manifest['sha256'].items():
    assert hashlib.sha256((root/'legacy'/name).read_bytes()).hexdigest()==expected, name
with tempfile.TemporaryDirectory(prefix='dns-legacy-test-') as tmp:
    work=pathlib.Path(tmp)/'comparison'
    subprocess.run([sys.executable,str(root/'tools/compare_legacy.py'),
        '--modern',modern,'--legacy',legacy,'--aligned-legacy',aligned,'--mpiexec',launcher,
        '--work-dir',str(work),'--grids','16','--ranks','1','3','--check-only'],check=True)
    results=json.loads((work/'results.json').read_text())
    for row in results['precision']:
        for component in ('u','v','w','p'):
            assert row['aligned_errors'][component]['max_abs']<1.e-9, (row['case'],component)
    print('PASS: preserved legacy sources and aligned zero-mode numerical parity')
