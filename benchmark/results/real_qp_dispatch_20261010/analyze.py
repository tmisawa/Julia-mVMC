from pathlib import Path
import json, math, re, statistics, sys

root=Path(__file__).parent
cases=[]
for folder in sorted(root.glob('paired-L*-r*-t*-s*')):
    match=re.fullmatch(r'paired-L(\d+)-r(\d+)-t(\d+)-s(\d+)',folder.name)
    size,ranks,threads,steps=map(int,match.groups())
    log=(folder/'run.log').read_text()
    records=re.findall(r'^PAIRED (baseline|candidate) ([123]) ([0-9.e+-]+) ([0-9.e+-]+)$',log,re.M)
    assert len(records)==6, folder
    times={name:{} for name in ['baseline','candidate']}
    output=[]
    for name,rep,seconds,energy in records:
        assert math.isfinite(float(energy)) and float(seconds)>0
        assert rep not in times[name]
        times[name][rep]=float(seconds)
        rows=[[float(x) for x in line.split()] for line in (folder/f'{name}-{rep}'/'zvo_out.dat').read_text().splitlines() if line.strip()]
        assert len(rows)==steps and all(len(row)==6 for row in rows)
        assert all(math.isfinite(x) for row in rows for x in row)
        output.append(rows)
    reference=output[0]
    # Existing long-run repeatability budget from docs/NUMERICAL_COMPARISONS.md.
    # This is supplementary to exact RNG/control audits and inverse residual tests.
    for rows in output:
        for row, ref in zip(rows, reference):
            for actual, expected in zip(row, ref):
                assert abs(actual-expected) <= 1e-11 + 1e-11*max(abs(actual),abs(expected)), (folder,actual,expected)
    maximum=max(abs(a-b) for rows in output for row,ref in zip(rows,reference) for a,b in zip(row,ref))
    medians={name:statistics.median(values.values()) for name,values in times.items()}
    cases.append(dict(sites=size,ranks=ranks,threads=threads,steps=steps,times=times,medians=medians,
        candidate_percent_change=100*(medians['candidate']/medians['baseline']-1),
        pair_wins=sum(times['candidate'][str(i)]<times['baseline'][str(i)] for i in range(1,4)),
        output_max_abs_observed=maximum))
assert cases
(root/'paired-summary.json').write_text(json.dumps({'cases':cases,'complete':True},indent=2)+'\n')
print(json.dumps(cases,indent=2))
