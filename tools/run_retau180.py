#!/usr/bin/env python3
"""Run or resume a single-rank case, optionally pin workers, then plot its means."""
import argparse
from datetime import datetime, timezone
import json
import os
import resource
from pathlib import Path
import struct
import subprocess
import sys
import time
import tomllib

from plot_mean_profile import plot


def cpu_list(value):
    cpus = [int(cpu) for cpu in value.split(',')]
    if len(set(cpus)) != len(cpus) or min(cpus) < 0:
        raise argparse.ArgumentTypeError('CPU ids must be unique and non-negative')
    return cpus


def pin(child, workers, controls):
    deadline = time.monotonic() + 15
    expected = len(workers) + 3  # main thread, two OpenMPI helpers, compute workers
    while time.monotonic() < deadline:
        tids = sorted(int(p.name) for p in Path(f'/proc/{child.pid}/task').iterdir())
        if len(tids) >= expected:
            break
        if child.poll() is not None:
            raise RuntimeError('solver exited before thread pinning')
        time.sleep(.002)
    if len(tids) != expected:
        raise RuntimeError(f'expected {expected} threads for pinning, found {len(tids)}')
    mapping = {}
    for tid, cpu in zip(tids[-len(workers):], workers):
        os.sched_setaffinity(tid, {cpu})
        mapping[tid] = dict(role='worker', cpus=sorted(os.sched_getaffinity(tid)))
    if controls:
        for index, tid in enumerate(tids[:-len(workers)]):
            os.sched_setaffinity(tid, {controls[index % len(controls)]})
            mapping[tid] = dict(role='control_or_MPI', cpus=sorted(os.sched_getaffinity(tid)))
    return mapping


def rewind_diagnostics(directory, step, time_value, stamp):
    # Drop diagnostic rows from discarded progress beyond the restart, preserving originals.
    paths = [(directory/'EKT 0.dat', lambda line: float(line.split()[0]) <= time_value + 1e-12),
             (directory/'diagnostics.csv', lambda line: line.startswith('step,') or
                                                       int(line.split(',')[0]) <= step)]
    for path, keep in paths:
        if path.exists():
            text = path.read_text()
            filtered = ''.join(line for line in text.splitlines(keepends=True) if keep(line))
            if filtered != text:
                path.rename(path.with_name(path.name + '.before_resume.' + stamp))
                path.write_text(filtered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_dir', type=Path)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--worker-cpus', type=cpu_list)
    parser.add_argument('--control-cpus', type=cpu_list)
    parser.add_argument('--no-plot', action='store_true', help='save statistics without drawing figures')
    args = parser.parse_args()
    if args.threads <= 0 or (args.worker_cpus and len(args.worker_cpus) != args.threads):
        parser.error('positive thread count and one CPU per pinned worker required')
    root = Path(__file__).resolve().parents[1]
    directory = args.run_dir.resolve()
    config = tomllib.loads((directory/'config.toml').read_text())
    control = config['control']
    state_path = directory/'run_state.json'
    start_step, start_time = 0, 0.
    resumed = control['INDEX_CONT'] == 1
    if resumed:
        with (directory/'field_S_ 0.bin').open('rb') as stream:
            magic, version, kind, start_step, start_time = struct.unpack('<8sIIqd', stream.read(32))
        if magic != b'CHDNS\0\0\0' or version != 1 or kind != 1:
            raise ValueError('runner expects a native binary restart')
        if start_step > control['Nstat']:
            raise ValueError('statistics restart is not supported after averaging has started')
    elif state_path.exists():
        raise FileExistsError(state_path)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    if state_path.exists():
        previous = json.loads(state_path.read_text())
        if previous['status'] not in ('stopped', 'failed', 'complete'):
            raise RuntimeError(f"cannot resume state {previous['status']}")
        state_path.rename(directory/f'run_state.before_resume.{stamp}.json')
    if resumed:
        rewind_diagnostics(directory, start_step, start_time, stamp)
    state = dict(status='starting', controller_pid=os.getpid(),
                 started_utc=datetime.now(timezone.utc).isoformat(), python=sys.executable,
                 threads=args.threads, ranks=1, start_step=start_step,
                 target_step=start_step+control['NSTEPS'])

    def save():
        temporary = state_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(state, indent=2)+'\n')
        temporary.replace(state_path)

    started = time.monotonic()
    # Children inherit this limit. On this host Apport writes un-packaged binary
    # cores under /var/lib/apport/coredump when the requested limit is nonzero.
    _, hard_core = resource.getrlimit(resource.RLIMIT_CORE)
    resource.setrlimit(resource.RLIMIT_CORE, (hard_core, hard_core))
    state['core_limit_bytes'] = hard_core
    state['core_pattern'] = Path('/proc/sys/kernel/core_pattern').read_text().strip()
    child = None
    save()
    try:
        command = [str(root/'build-channel-dns/bin/channel_dns'),
                   '--threads', str(args.threads), '--timing', str(directory)]
        with (directory/'solver.log').open('a' if resumed else 'w') as log:
            if resumed:
                log.write(f'\nResuming at step {start_step} with {args.threads} workers ({stamp})\n')
                log.flush()
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, cwd=root)
            state.update(status='running', solver_pid=child.pid, command=command)
            save()
            if args.worker_cpus:
                state['affinity'] = pin(child, args.worker_cpus, args.control_cpus)
                save()
            result = child.wait()
        if result == -15:
            state.update(status='stopped', stop_reason='solver_received_SIGTERM')
            return
        if result:
            raise RuntimeError(f'solver exited with code {result}; see solver.log')
        expected = max(0, start_step+control['NSTEPS']-max(start_step, control['Nstat']))
        if args.no_plot:
            actual = int((directory/'Result/0/DPDX_DPDZ.dat').read_text().split()[0])
            if actual != expected:
                raise ValueError(f'expected {expected} samples, found {actual}')
            state.update(samples=actual, plot_generated=False)
        else:
            state.update(status='plotting')
            save()
            state['profile'] = plot(directory, expected)
        state.update(status='complete')
    except BaseException as error:
        if child is not None and child.poll() is None:
            child.terminate()
            child.wait()
        state.update(status='failed', error=str(error))
        raise
    finally:
        state.update(elapsed_seconds=time.monotonic()-started,
                     finished_utc=datetime.now(timezone.utc).isoformat())
        save()


if __name__ == '__main__':
    main()
