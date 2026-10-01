"""python -m revision2 --input revision2/example.json --output results/v2"""
import argparse
import csv
from dataclasses import fields
import json
from pathlib import Path
import sys
from .model import Crust,Planet,Event
from .evolution import evolve,sample_events


def construct(cls,data):
    if not isinstance(data,dict): raise ValueError(f'{cls.__name__} must be an object')
    unknown=data.keys()-{f.name for f in fields(cls)}
    if unknown: raise ValueError(f'Unknown {cls.__name__} fields: {sorted(unknown)}')
    return cls(**data)


def load(path):
    data=json.loads(Path(path).read_text())
    if set(data)-{'schema_version','planet','events','random_events'} or data.get('schema_version')!=2:
        raise ValueError('Expected revision2 schema_version 2')
    raw=dict(data['planet']);raw['crust']=construct(Crust,raw.get('crust',{}))
    p=construct(Planet,raw);p.validate()
    events=[construct(Event,e) for e in data.get('events',[])]
    events+=sample_events(p,**data.get('random_events',{}))
    for e in events:e.validate(p.age_myr)
    return p,events


def write_outputs(output,result):
    output.mkdir(parents=True,exist_ok=True)
    (output/'evolution.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    rows=[]
    for q in result['epochs']:
        rows.append({'epoch':q['epoch'],'age_myr':q['age_myr'],
            'temperature_k':q['equilibrium']['temperature_k'],
            'pressure_bar':q['equilibrium']['pressure_pa']/1e5,
            'ocean_crust_contact':q['ocean_crust_contact'],'tectonics_active':q['tectonics_active'],
            'oxidation_capacity_mol_o_m2':q['oxidation_capacity_mol_o_m2'],
            'energy_residual_w_m2':q['energy_residual_w_m2'],**q['surface_fractions']})
    with (output/'history.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    final=result['epochs'][-1]
    with (output/'layers.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(final['layers'][0]));w.writeheader();w.writerows(final['layers'])


def main(argv=None):
    parser=argparse.ArgumentParser(description='Conserved crust-up equilibrium epochs. Epoch count × duration sets age.')
    parser.add_argument('--input',type=Path,default=Path(__file__).with_name('example.json'))
    parser.add_argument('--output',type=Path,default=Path('results/revision2'))
    args=parser.parse_args(argv)
    try:
        p,events=load(args.input)
        def progress(i,n,q):
            print(f'Epoch {i}/{n}: {q.age_myr:g} Myr, {q.equilibrium.temperature_k:.2f} K, '
                  f'{q.equilibrium.pressure_pa/1e5:.4g} bar, tectonics={q.tectonics_active}',flush=True)
        result=evolve(p,events,progress)
        write_outputs(args.output,result)
    except (ValueError,TypeError,KeyError,RuntimeError,OSError) as exc:
        print(f'Error: {exc}',file=sys.stderr);return 1
    print(f'Wrote {args.output.resolve()}');return 0


if __name__=='__main__': raise SystemExit(main())
