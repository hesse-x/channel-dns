#!/usr/bin/env python3
"""Compare the FFTW solver and byte-preserved legacy/FFTPACK on identical TOML inputs."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import os
import pathlib
import platform
import re
import statistics
import subprocess
import sys
import time
import tomllib

sys.dont_write_bytecode = True
ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('compare_fields', ROOT/'tests/tools/compare_fields.py')
fields = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fields)


def partition(n, ranks):
    q, r = divmod(n, ranks)
    counts = [q + (rank >= ranks-r) for rank in range(ranks)]
    return [sum(counts[:rank]) for rank in range(ranks)], counts


def parameters(grid, ranks, steps, full=False, laminar=False):
    config = tomllib.loads((ROOT/'examples/config.toml').read_text())
    m, c, p = config['mesh'], config['control'], config['parm']
    m.update(NX=grid, NY=grid, NZ=grid, NXS=grid//2 if full else grid//3,
             NYS=grid-2, NZS=grid//2 if full else grid//3,
             NX2=grid if full else 3*grid//2, NY2=3*grid//2, NZ2=grid if full else 3*grid//2,
             INDEX_Serial=int(ranks==1))
    c.update(INDEX_CONT=1, NSTEPS=steps, NDUMFRQ=steps, NPRNFRQ=steps,
             EPS_init=0., NCOUNT_disturb=0, EPS_disturb=0., Nstat=2_000_000_000)
    config['output'] = dict(format='ascii', physical_interval=steps)
    return config


def write_toml(directory,config):
    lines=[]
    for section, values in config.items():
        lines.append(f'[{section}]')
        for key, value in values.items():
            lines.append(f'{key} = {json.dumps(value)}')
        lines.append('')
    (directory/'config.toml').write_text('\n'.join(lines))


def write_legacy_parameters(directory, config):
    m,p,c=(config[key] for key in ('mesh','parm','control'))
    records = {
        'meshparm_channel_with_Microblow.txt': [
            ('NX,NY,NZ,INDEX_disturb_Z', [m[k] for k in ('NX','NY','NZ','INDEX_disturb_Z')]),
            ('NXS,NYS,NZS',[m[k] for k in ('NXS','NYS','NZS')]),
            ('NX2,NY2,NZ2',[m[k] for k in ('NX2','NY2','NZ2')]),
            ('INDEX_Serial',[m['INDEX_Serial']])],
        'parm_channel_with_Microblow.txt': [
            ('XL,ZL,THETA',[p[k] for k in ('XL','ZL','THETA')]),
            ('INDEX_GRAD_PRE,Re',[p[k] for k in ('INDEX_GRAD_PRE','Re')]),
            ('INDEX_UNSTEADY_BC',[p['INDEX_UNSTEADY_BC']]),
            ('G_pre,a_pre,b_pre',[p[k] for k in ('G_pre','a_pre0','b_pre0')])],
        'control_channel_with_Microblow.txt': [
            ('INDEX_CONT,EPS_init',[c[k] for k in ('INDEX_CONT','EPS_init')]),
            ('NSTEPS,NDUMFRQ,NPRNFRQ,DT',[c[k] for k in ('NSTEPS','NDUMFRQ','NPRNFRQ','DT')]),
            ('NCOUNT_disturb,NDISTURBFRQ,EPS_disturb',[c[k] for k in ('NCOUNT_disturb','NDISTURBFRQ','EPS_disturb')]),
            ('V_Microblow0',[c['V_Microblow0']]),('Nstat',[c['Nstat']])],
    }
    for name, pairs in records.items():
        lines=[]
        for title,values in pairs:
            lines.extend([title,', '.join(str(v) for v in values)])
        if name.startswith('parm_'):
            lines.append('a_V,b_V')
            lines.extend(f'{a}, {b}' for a,b in zip(p['a_V'],p['b_V']))
        (directory/name).write_text('\n'.join(lines)+'\n')


def initial_checkpoint(directory,config,ranks,laminar=False):
    """Smooth deterministic divergence-free perturbation; both solvers read the same text."""
    m,p,c = (config[key] for key in ('mesh','parm','control'))
    xb,xc=partition(m['NXS']+1,ranks)
    hashes={}
    for rank in range(ranks):
        path=directory/f'field_S_{rank:2d}.dat'
        with path.open('w') as out:
            out.write('0 0.0\n')
            for z in range(-m['NZS'],m['NZS']+1):
                for y in range(m['NYS']+1):
                    for x in range(xb[rank],xb[rank]+xc[rank]):
                        values=[0j,0j,0j]
                        if x==z==0:
                            if y==0: values[0]=.75
                            if y==2: values[0]=-.75
                            if y==0: values[1]=c['V_Microblow0']
                        if not laminar and x in (1,2) and abs(z)==1:
                            # psi = epsilon*(1-y^2)^2*cos(kx*x)*cos(kz*z)
                            # u=dy(psi), v=-dx(psi) satisfies divergence and no-slip walls.
                            amplitude=1.e-3/(4*x)
                            if y in (1,3): values[0]+=amplitude*(-1 if y==1 else 1)
                            phi={0:3/8,2:-1/2,4:1/8}.get(y,0.)
                            values[1]+=-1j*(2*math.pi*x/p['XL'])*amplitude*phi
                        for value in values:
                            out.write(f'({value.real:.17e},{value.imag:.17e})\n')
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def prepare(directory,config,ranks,laminar=False):
    directory.mkdir(parents=True,exist_ok=False)
    write_toml(directory,config)
    write_legacy_parameters(directory,config)
    return initial_checkpoint(directory,config,ranks,laminar)


def compare_outputs(legacy,modern,require_pressure=False):
    lh,lf=fields.read_fields(legacy)
    mh,mf=fields.read_fields(modern)
    if lh != mh:
        raise AssertionError(f'checkpoint metadata differs: {lh} != {mh}')
    errors={}
    for name, left, right in zip(('u','v','w','p'),lf,mf):
        if not all(math.isfinite(v) for v in left+right):
            raise AssertionError('nonfinite output')
        delta=[b-a for a,b in zip(left,right)]
        max_abs=max(abs(v) for v in delta)
        norm=max(abs(v) for v in left)
        errors[name]=dict(max_abs=max_abs,relative=max_abs/max(norm,1.e-300),
                          rms=math.sqrt(sum(v*v for v in delta)/len(delta)))
        if name=='p':
            offset=statistics.fmean(delta)
            errors[name]['gauge_offset']=offset
            errors[name]['gauge_aligned_max_abs']=max(abs(v-offset) for v in delta)
            nx,ny,nz=lh[2:]
            plane_delta=[sum(delta[(z*ny+y)*nx+x] for z in range(nz) for x in range(nx))/(nx*nz) for y in range(ny)]
            errors[name]['nonzero_horizontal_mode_max_abs']=max(abs(delta[(z*ny+y)*nx+x]-plane_delta[y]) for z in range(nz) for y in range(ny) for x in range(nx))
        if (name!='p' or require_pressure) and max_abs>2.e-9+2.e-9*norm:
            raise AssertionError(f'{name} differs: {errors[name]}')
    le,me=fields.read_ekt(legacy),fields.read_ekt(modern)
    if len(le)!=len(me) or any(a[0]!=b[0] for a,b in zip(le,me)):
        raise AssertionError('EKT times differ')
    errors['ekt']={'max_abs':max(abs(a[1]-b[1]) for a,b in zip(le,me))}
    return errors



def compare_restart(legacy,modern,ranks):
    """Stream final spectral fields, so 128-grid validation does not need full copies."""
    maximum=[0.,0.,0.]
    norms=[0.,0.,0.]
    for rank in range(ranks):
        def values(path):
            with path.open() as stream:
                header=fields.nums(next(stream))
                yield header
                pending=None
                for line in stream:
                    for number in fields.nums(line):
                        if pending is None:
                            pending=number
                        else:
                            yield complex(pending,number)
                            pending=None
                if pending is not None: raise AssertionError('partial complex checkpoint')
        left=values(legacy/f'field_S_{rank:2d}.dat')
        right=values(modern/f'field_S_{rank:2d}.dat')
        if next(left)!=next(right): raise AssertionError('final checkpoint metadata differs')
        count=0
        for count,(a,b) in enumerate(itertools.zip_longest(left,right),1):
            if a is None or b is None: raise AssertionError('final checkpoint size differs')
            if not all(math.isfinite(v) for v in (a.real,a.imag,b.real,b.imag)):
                raise AssertionError('nonfinite final spectral field')
            c=(count-1)%3
            maximum[c]=max(maximum[c],abs(a-b));norms[c]=max(norms[c],abs(a))
        if count==0: raise AssertionError('empty final checkpoint')
    if any(error>2.e-9+2.e-9*norm for error,norm in zip(maximum,norms)):
        raise AssertionError(f'long-run spectral field mismatch: {maximum}')
    return dict(zip(('u','v','w'),maximum))


def run(executable,directory,ranks,warmup,launcher,bind,modern,timeout=600,measure_memory=False):
    command=[launcher]
    if bind: command+=['--bind-to','core','--map-by','core']
    command+=['-n',str(ranks)]
    if measure_memory:
        # Separate rank logs: MPI may interleave individual stderr writes.
        command+=['sh','-c',
                  'out="$1"; shift; exec /usr/bin/time -v -o "$out/memory-rank${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}.log" "$@"',
                  'dns-memory',str(directory)]
    command+=[str(executable),'--benchmark-warmup',str(warmup)]
    if modern: command+=[str(directory)]
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',LC_ALL='C',
             OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    start=time.perf_counter()
    proc=subprocess.run(command,cwd=directory,env=env,capture_output=True,text=True,timeout=timeout)
    wall=time.perf_counter()-start
    log=proc.stdout+proc.stderr
    (directory/'run.log').write_text(log)
    if proc.returncode:
        raise RuntimeError(f'{command}:\n{log[-6000:]}')
    result={'process_seconds':wall}
    if measure_memory:
        memory_logs=sorted(directory.glob('memory-rank*.log'))
        peaks=[int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',path.read_text())[1]) for path in memory_logs]
        if len(peaks)!=ranks: raise AssertionError('missing per-rank peak memory')
        result['rank_peak_rss_kib']=peaks
        result['sum_rank_peak_rss_gib']=sum(peaks)/1024**2
    for key,label in [('compute_seconds','compute wall time'),('loop_seconds','loop wall time')]:
        match=re.search(r'benchmark '+label+r' =\s*([0-9.eE+-]+) s',log)
        if not match: raise AssertionError(f'missing {label}: {log[-2000:]}')
        result[key]=float(match[1])
    measured=int(re.search(r'benchmark measured steps =\s*(\d+)',log)[1])
    if measured<=0: raise AssertionError('no measured steps')
    result['measured_steps']=measured
    result['compute_ms_per_step']=1000*result['compute_seconds']/measured
    result['command']=command
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modern',type=pathlib.Path,required=True)
    parser.add_argument('--legacy',type=pathlib.Path,required=True)
    parser.add_argument('--aligned-legacy',type=pathlib.Path)
    parser.add_argument('--work-dir',type=pathlib.Path,required=True)
    parser.add_argument('--grids',type=int,nargs='+',default=[32,64,128])
    parser.add_argument('--ranks',type=int,nargs='+',default=[1,4,8])
    parser.add_argument('--steps',type=int,default=200)
    parser.add_argument('--warmup',type=int,default=100)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--mpiexec',default='mpiexec')
    parser.add_argument('--no-bind',action='store_true')
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args()
    if min(args.grids)<8 or any(g%2 for g in args.grids) or min(args.ranks)<1 or args.steps<1 or args.warmup<0 or args.repeats<1:
        parser.error('invalid grid, rank count, steps, warmup or repeats')
    work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=False)
    modern=args.modern.resolve();legacy=args.legacy.resolve()
    result={'host':platform.platform(),'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
            'cpu':subprocess.check_output(['lscpu'],text=True),'options':vars(args).copy(),
            'precision':[],'timing':[],'summary':[]}
    result['source_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(ROOT.glob('src/*')) if path.suffix in ('.f90','.cpp')}
    result['options']={key:str(v) if isinstance(v,pathlib.Path) else v for key,v in result['options'].items()}
    def save(): (work/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    # Accuracy is evaluated separately with short, equal-length trajectories and final text output.
    for case,full,laminar in [('dealiased',False,False),('full',True,False),('laminar',False,True)]:
        for ranks in args.ranks:
            grid=min(args.grids)
            config=parameters(grid,ranks,5,full,laminar)
            base=work/f'precision-{case}-np{ranks}'
            hashes=[]
            for label,exe in [('legacy',legacy),('modern',modern)]:
                directory=base/label
                hashes.append(prepare(directory,config,ranks,laminar))
                run(exe,directory,ranks,0,args.mpiexec,not args.no_bind,label=='modern')
            assert hashes[0]==hashes[1], 'different initial checkpoints'
            errors=compare_outputs(base/'legacy',base/'modern')
            row=dict(case=case,grid=grid,ranks=ranks,errors=errors,initial_sha256=hashes[0])
            if args.aligned_legacy:
                directory=base/'aligned-legacy'
                prepare(directory,config,ranks,laminar)
                run(args.aligned_legacy.resolve(),directory,ranks,0,args.mpiexec,not args.no_bind,False)
                row['aligned_errors']=compare_outputs(directory,base/'modern',require_pressure=True)
            result['precision'].append(row)
            save()
            print(f'precision {case} grid={grid} np={ranks}: max_abs={max(v["max_abs"] for v in errors.values()):.3e}',flush=True)
    if not args.check_only:
        for grid in args.grids:
            for ranks in args.ranks:
                config=parameters(grid,ranks,args.steps+args.warmup)
                # Only a final spectral checkpoint is written, after compute timing.
                # Full physical fields are reserved for the short precision cases.
                config['control']['NPRNFRQ']=args.steps+args.warmup+1
                config['output']['physical_interval']=0
                # Alternate order to reduce thermal/frequency bias; every run starts at the same state.
                for repeat in range(args.repeats):
                    variants=[('legacy',legacy),('modern',modern)]
                    if repeat%2: variants.reverse()
                    for label,exe in variants:
                        directory=work/f'timing-g{grid}-np{ranks}-r{repeat}'/label
                        initial_hashes=prepare(directory,config,ranks)
                        row=run(exe,directory,ranks,args.warmup,args.mpiexec,not args.no_bind,label=='modern')
                        if row['measured_steps']!=args.steps: raise AssertionError('incorrect measured step count')
                        row.update(variant=label,grid=grid,ranks=ranks,repeat=repeat,initial_sha256=initial_hashes)
                        result['timing'].append(row);save()
                        print(f'{label} grid={grid} np={ranks} repeat={repeat}: {row["compute_ms_per_step"]:.3f} ms/step',flush=True)
                    final_errors=compare_restart(work/f'timing-g{grid}-np{ranks}-r{repeat}'/'legacy',
                                                 work/f'timing-g{grid}-np{ranks}-r{repeat}'/'modern',ranks)
                    for row in result['timing'][-2:]: row['final_spectral_max_abs']=final_errors
                    save()
                times={label:[row['compute_ms_per_step'] for row in result['timing'] if row['grid']==grid and row['ranks']==ranks and row['variant']==label] for label in ('legacy','modern')}
                summary=dict(grid=grid,ranks=ranks,legacy_median_ms=statistics.median(times['legacy']),
                             modern_median_ms=statistics.median(times['modern']),
                             legacy_range_ms=[min(times['legacy']),max(times['legacy'])],modern_range_ms=[min(times['modern']),max(times['modern'])])
                summary['speedup']=summary['legacy_median_ms']/summary['modern_median_ms']
                result['summary'].append(summary);save()
    print(f'Results: {work}/results.json',flush=True)


if __name__=='__main__':
    main()
