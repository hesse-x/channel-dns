#!/usr/bin/env python3
"""Summarize channel_dns on/off-CPU intervals from perf --switch-events cycles data.
Non-preempted off-CPU is descriptive: it does not measure recoverable queue cost.
"""
import subprocess,re,json,collections,argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("data")
parser.add_argument("--output", required=True)
args=parser.parse_args()
cmd=['perf','script','-i',args.data,'--show-switch-events','-F','comm,pid,tid,time,event,period,ip,sym,dso']
p=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True,stderr=subprocess.DEVNULL)
header=re.compile(r'^\s*(\S+)\s+(\d+)/(\d+)\s+(\d+\.\d+):\s+(.*)$')
switches=collections.defaultdict(list);periods=collections.Counter();workers=set();owners=set();pid=None;current=None;tmin=None;tmax=None
for line in p.stdout:
    m=header.match(line)
    if m:
        comm,proc,tid,t,rest=m.groups();proc=int(proc);tid=int(tid);t=float(t);current=None
        if comm!='channel_dns':continue
        if 'PERF_RECORD_SWITCH' in rest:
            switches[tid].append((t,' IN' in rest,'preempt' in rest));continue
        if 'cpu_core/cycles/u:' in rest:
            pid=proc;owners.add(proc);current=tid;periods[tid]+=int(rest.split()[0])
            tmin=t if tmin is None else min(tmin,t);tmax=t if tmax is None else max(tmax,t)
    elif current is not None and any(s in line for s in ('ThreadPool::worker()', 'CpuDevice::Impl::worker()')):workers.add(current)
if p.wait() != 0:
    raise SystemExit("perf script failed")
if tmin is None or tmax is None:
    raise SystemExit("No channel_dns cpu_core/cycles/u samples found")
result=[]
for tid in sorted(workers|owners):
    ev=sorted(switches[tid]);prior=[x for x in ev if x[0]<=tmin]
    if prior:state=prior[-1][1];preempt=prior[-1][2]
    else:
        after=next((x for x in ev if x[0]>=tmin),None)
        state=not after[1] if after else True;preempt=False
    last=tmin;on=blocked=preempted=0.;nsw=0
    for t,newstate,flag in ev:
        if t<=tmin:continue
        stop=min(t,tmax);dt=max(0,stop-last)
        if state:on+=dt
        elif preempt:preempted+=dt
        else:blocked+=dt
        last=stop
        if t>=tmax:break
        state=newstate;preempt=flag;nsw+=1
    if last<tmax:
        dt=tmax-last
        if state:on+=dt
        elif preempt:preempted+=dt
        else:blocked+=dt
    result.append(dict(tid=tid,role='owner' if tid in owners else 'worker',oncpu=on,blocked=blocked,preempted=preempted,switch_records=nsw,cycles_share=periods[tid]/sum(periods.values())))
output=dict(window_seconds=tmax-tmin,owner_pid=pid,threads=result)
open(args.output,'w').write(json.dumps(output,indent=2))
print(json.dumps(output,indent=2))
