"""A physical-time convergence experiment, independent of the demo's outputs.

Run: python -m revision2.validate_resolution
The printed JSON can be redirected to a report file. This takes ~1 minute.
"""
import json
from .model import Planet,Crust,Event
from .evolution import evolve


def main():
    rows=[]
    for fraction in (.05,.025,.0125):
        p=Planet(epochs=2,epoch_duration_myr=2.,max_transport_step_myr=1.,
                 max_escape_fraction_per_step=fraction,stellar_flux_w_m2=1800.,
                 crust=Crust(cao_mol_m2=1000.))
        result=evolve(p,[Event(1.5,'comet',elements={'H':20000.,'O':10000.}),
                        Event(3.,'impact',strip_fraction=.2,
                              mantle_release_fraction=.001,resurface_fraction=.01)])
        final=result['epochs'][-1]
        rows.append({'maximum_step_myr':1.,'maximum_fractional_change':fraction,
                     'temperature_k':final['equilibrium']['temperature_k'],
                     'pressure_pa':final['equilibrium']['pressure_pa'],
                     'escaped_h_mol_m2':result['final_state']['escaped']['H'],
                     'maximum_budget_error_mol_m2':max(abs(v) for v in final['element_budget_error_mol_m2'].values())})
    print(json.dumps(rows,indent=2,allow_nan=False))


if __name__=='__main__':main()
