"""Summarize completed alternating controls; validate actual worlds and outputs."""
from pathlib import Path
import json
import math
import re
import statistics

root=Path(__file__).resolve().parent
rows=[]
for ranks,threads in ((1,16),(4,4)):
    for size in (32,64):
        path=root/'paired'/f'L{size}-r{ranks}-t{threads}'
        if not (path/'run.log').exists():
            continue
        text=(path/'run.log').read_text()
        records=re.findall(r'^PAIRED (baseline|candidate) (\d+) ([0-9.]+) (\S+)$',text,re.M)
        if len(records)!=6:
            continue
        for label,value in (('WORLD',ranks),('THREADS',threads),('BLAS_THREADS',1),('NATIVE_BLAS_THREADS',1)):
            assert sorted((int(a),int(b)) for a,b in re.findall(rf'^{label} (\d+) (\d+)$',text,re.M))==[(r,value) for r in range(ranks)]
        times={version:{int(rep):float(seconds) for v,rep,seconds,energy in records if v==version}
               for version in ('baseline','candidate')}
        assert all(sorted(values)==[1,2,3] for values in times.values())
        assert all(math.isfinite(float(seconds)) and float(seconds)>0 and math.isfinite(float(energy))
                   for v,rep,seconds,energy in records)
        delta=0.0
        for rep in (1,2,3):
            arrays=[]
            for version in ('baseline','candidate'):
                array=[[float(x) for x in line.split()] for line in (path/f'{version}-{rep}/zvo_out.dat').read_text().splitlines()]
                assert len(array)==300 and all(len(row)==6 and all(math.isfinite(x) for x in row) for row in array)
                arrays.append(array)
            delta=max(delta,max(abs(x-y) for a,b in zip(*arrays) for x,y in zip(a,b)))
        medians={version:statistics.median(values.values()) for version,values in times.items()}
        rows.append(dict(sites=size,ranks=ranks,threads=threads,times=times,medians=medians,
                         candidate_percent_change=100*(medians['candidate']/medians['baseline']-1),
                         candidate_pair_wins=sum(times['candidate'][rep]<times['baseline'][rep] for rep in (1,2,3)),
                         output_max_abs_observed=delta))
report=dict(cases=rows,complete=len(rows)==4,
            protocol='Each process warmed both paths, then ran three alternating pairs. A typed Ref gate selector was installed before all warmups; no timed eval, JIT or Profile. Production API max-rank timing, Opt300, total320 samples. Selector overhead is common; separate unmodified-production source pilots and exact audits are preserved.')
(root/'paired-summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
