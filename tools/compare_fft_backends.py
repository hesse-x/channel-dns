#!/usr/bin/env python3
"""Isolate FFT backend cost in the same MPI solver; 100 warmup + 200 measured steps."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import platform
import subprocess
import compare_legacy as bench


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fftw',type=Path,required=True)
    p.add_argument('--fftpack',type=Path,required=True)
    p.add_argument('--work-dir',type=Path,required=True)
    p.add_argument('--grids',type=int,nargs='+',default=[32,64,128])
    p.add_argument('--ranks',type=int,default=8)
    p.add_argument('--warmup',type=int,default=100)
    p.add_argument('--steps',type=int,default=200)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--skip-physical-check',action='store_true',help='skip full physical field comparison; always stream final spectral checks')
    p.add_argument('--measure-memory',action='store_true')
    p.add_argument('--timeout',type=int,default=600,help='per MPI run timeout in seconds')
    a=p.parse_args()
    if min(a.grids)<8 or any(g%2 for g in a.grids) or a.ranks<1 or a.steps<1 or a.warmup<0 or a.repeats<1:
        p.error("invalid grid, ranks, steps, warmup or repeats")
    work=a.work_dir.resolve();work.mkdir(parents=True,exist_ok=False)
    variants={'fftw':a.fftw.resolve(),'fftpack':a.fftpack.resolve()}
    result={'host':platform.platform(),
            'cpu':subprocess.check_output(['lscpu'],text=True),
            'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
            'options':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
            'binary_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in variants.items()},
            'source_sha256':{str(v.relative_to(bench.ROOT)):hashlib.sha256(v.read_bytes()).hexdigest() for v in sorted((bench.ROOT/'src').glob('*')) if v.suffix in ('.f90','.cpp')},
            'precision':[],'timing':[],'summary':[]}
    def save(): (work/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    # Keep regression tests and other MPI jobs out of this timing session.
    for grid in a.grids:
        if not a.skip_physical_check:
            config=bench.parameters(grid,a.ranks,5)
            base=work/f'precision-g{grid}'
            hashes=[]
            for label,exe in variants.items():
                directory=base/label
                hashes.append(bench.prepare(directory,config,a.ranks))
                bench.run(exe,directory,a.ranks,0,'mpiexec',True,True)
            assert hashes[0]==hashes[1]
            errors=bench.compare_outputs(base/'fftw',base/'fftpack',require_pressure=True)
            result['precision'].append(dict(grid=grid,errors=errors))
            print(f'precision grid={grid}: {errors}',flush=True);save()
        config=bench.parameters(grid,a.ranks,a.warmup+a.steps)
        config['control']['NPRNFRQ']=a.warmup+a.steps+1
        config['output']['physical_interval']=0
        for repeat in range(a.repeats):
            labels=list(variants)
            if repeat%2: labels.reverse()
            base=work/f'timing-g{grid}-r{repeat}'
            hashes=[]
            for label in labels:
                directory=base/label
                hashes.append(bench.prepare(directory,config,a.ranks))
                row=bench.run(variants[label],directory,a.ranks,a.warmup,'mpiexec',True,True,timeout=a.timeout,measure_memory=a.measure_memory)
                assert row['measured_steps']==a.steps
                row.update(variant=label,grid=grid,ranks=a.ranks,repeat=repeat)
                result['timing'].append(row);save()
                print(f'{label} grid={grid} repeat={repeat}: {row["compute_ms_per_step"]:.3f} ms/step',flush=True)
            assert hashes[0]==hashes[1]
            errors=bench.compare_restart(base/'fftw',base/'fftpack',a.ranks)
            for row in result['timing'][-2:]: row['final_spectral_max_abs']=errors
            save()
        times={label:[r['compute_ms_per_step'] for r in result['timing'] if r['grid']==grid and r['variant']==label] for label in variants}
        row={'grid':grid,'ranks':a.ranks}
        for label,values in times.items():
            row[label+'_median_ms']=statistics.median(values)
            row[label+'_range_ms']=[min(values),max(values)]
        row['fftpack_over_fftw']=row['fftpack_median_ms']/row['fftw_median_ms']
        result['summary'].append(row);save()
        print(row,flush=True)


if __name__=='__main__': main()
