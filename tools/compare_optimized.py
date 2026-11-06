#!/usr/bin/env python3
"""Matched whole-loop comparison of the pre-optimization, optimized and optional C++ executables."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import time

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('legacy_compare',ROOT/'tools/compare_legacy.py')
helper=importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before',type=Path,required=True)
    parser.add_argument('--after',type=Path,required=True)
    parser.add_argument('--cpp',type=Path)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--grids',type=int,nargs='+',default=[32,64])
    parser.add_argument('--ranks',type=int,nargs='+',default=[1,2,4])
    parser.add_argument('--steps',type=int,default=60)
    parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args()
    if min(args.grids)<8 or any(g%2 for g in args.grids) or min(args.ranks)<1 or min(args.steps,args.repeats)<1:
        parser.error('invalid benchmark layout')
    variants=[('before',args.before.resolve()),('after',args.after.resolve())]
    if args.cpp: variants.append(('cpp',args.cpp.resolve()))
    work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=False)
    result=dict(method='MPI maximum whole-loop wall time; includes first two setup steps and last diagnostic; no field output during timed runs',
                steps=args.steps,repeats=args.repeats,timing=[],precision=[],summary=[],
                executable_sha256={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in variants},
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('src/*') if p.suffix in ('.f90','.cpp')})
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    def save(): (work/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    def run(executable,directory,ranks):
        cmd=['mpiexec','--bind-to','core','--map-by','core','-n',str(ranks),str(executable),'--timing',str(directory)]
        start=time.perf_counter()
        proc=subprocess.run(cmd,capture_output=True,text=True,env=env,timeout=600)
        text=proc.stdout+proc.stderr
        (directory/'run.log').write_text(text)
        if proc.returncode: raise RuntimeError(text[-8000:])
        match=re.search(r'timestep wall time =\s*([0-9.eE+-]+) s',text)
        if not match: raise AssertionError(text[-2000:])
        return dict(loop_seconds=float(match[1]),process_seconds=time.perf_counter()-start,command=cmd)
    # Accuracy uses identical nonzero-mode checkpoints, including full Nyquist spectra.
    for case,full in [('dealiased',False),('full',True)]:
        for ranks in args.ranks:
            config=helper.parameters(min(args.grids),ranks,5,full=full)
            base=work/f'precision-{case}-np{ranks}'
            for label,exe in variants:
                directory=base/label
                helper.prepare(directory,config,ranks)
                run(exe,directory,ranks)
            for label,_ in variants[1:]:
                errors=helper.compare_outputs(base/'before',base/label,require_pressure=True)
                result['precision'].append(dict(case=case,ranks=ranks,variant=label,errors=errors))
            save()
            print(f'precision {case} np={ranks}: PASS',flush=True)
    for grid in args.grids:
        for ranks in args.ranks:
            config=helper.parameters(grid,ranks,args.steps)
            config['control'].update(NDUMFRQ=args.steps+1,NPRNFRQ=args.steps+1)
            config['output'].update(physical_interval=0)
            # Untimed per-executable warm run creates wisdom. Timed runs retain the
            # same input state and wisdom, and rotate execution order each repeat.
            directories={}
            for label,exe in variants:
                directory=work/f'g{grid}-np{ranks}'/label
                helper.prepare(directory,config,ranks)
                run(exe,directory,ranks)
                directories[label]=directory
            for repeat in range(args.repeats):
                ordered=variants[repeat%len(variants):]+variants[:repeat%len(variants)]
                for label,exe in ordered:
                    directory=directories[label]
                    row=run(exe,directory,ranks)
                    row.update(variant=label,grid=grid,ranks=ranks,repeat=repeat,
                               ms_per_step=1000*row['loop_seconds']/args.steps)
                    (directory/f'run-{repeat}.log').write_text((directory/'run.log').read_text())
                    result['timing'].append(row);save()
                    print(f'{label} grid={grid} np={ranks} repeat={repeat}: {row["ms_per_step"]:.3f} ms/step',flush=True)
            summary=dict(grid=grid,ranks=ranks)
            for label,_ in variants:
                samples=[r['ms_per_step'] for r in result['timing'] if r['grid']==grid and r['ranks']==ranks and r['variant']==label]
                summary[label]=dict(median_ms=statistics.median(samples),min_ms=min(samples),max_ms=max(samples))
            summary['speedup']=summary['before']['median_ms']/summary['after']['median_ms']
            if args.cpp: summary['after_over_cpp']=summary['after']['median_ms']/summary['cpp']['median_ms']
            result['summary'].append(summary);save()
    print(work/'results.json',flush=True)


if __name__=='__main__': main()
