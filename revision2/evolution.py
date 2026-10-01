"""Only reservoir transport and scheduled events advance geological time."""
from dataclasses import asdict
import copy
import math
import random
from .model import State, Event, inventory, add, transfer
from .thermo import ELEMENTS, ATOMIC_MASS, MOLAR_MASS, elements_of, water_psat
from .stack import solve_stack

SECONDS_MYR = 365.25*86400*1e6
KB = 1.380649e-23
NA = 6.02214076e23
G = 6.67430e-11


def resurface(p,state,snapshot,fraction):
    """Bury old reactive solids, expose a finite fresh veneer. Never reset O."""
    actual = min(fraction,state.fresh_veneer_units)
    if actual <= 0: return 0.
    old = elements_of({s:n*actual for s,n in snapshot.equilibrium.condensed.items()
                       if s not in ('water','ice')})
    add(state.accessible,old,-1); add(state.buried,old)
    fresh = inventory(); add(fresh,p.crust.fresh_elements,actual)
    add(state.fresh_crust,fresh,-1); add(state.accessible,fresh)
    state.fresh_veneer_units -= actual
    state.buried_veneer_units += actual
    state.assert_conserved()
    return actual


def escape_rates(p,snapshot):
    """Elemental loss rates (mol/m²/Myr) at this equilibrium boundary.

    Phase exchange is instantaneous, so the cold trap limits water delivery,
    not the amount of vapor found in one arbitrarily long integration step.
    """
    lost = inventory()
    eq = snapshot.equilibrium
    if not p.escape_enabled or not eq.gas_mass: return lost
    total_moles = sum(eq.gas.values())
    mu = eq.gas_mass/total_moles
    tx = p.exosphere_temperature_k
    height = 8.314462618*tx/(mu*p.gravity)
    nex = 1/(p.collision_cross_section_m2*height)
    nsurf = eq.pressure_pa/(KB*tx)
    altitude = height*max(0.,math.log(max(nsurf/nex,1.)))
    if altitude > p.radius_m:
        from .stack import UnsupportedState
        raise UnsupportedState('Exobase altitude exceeds planet radius; hydrodynamic escape required')
    radius = p.radius_m+altitude
    area_ratio = (radius/p.radius_m)**2
    requests = {}
    for species,n in eq.gas.items():
        if n <= 0 or species == 'H2O': continue
        molecular_mass = MOLAR_MASS[species]/NA
        lam = G*p.mass_kg*molecular_mass/(KB*tx*radius)
        flux = nex*(n/total_moles)*math.sqrt(KB*tx/(2*math.pi*molecular_mass))*(1+lam)*math.exp(-min(lam,745))*area_ratio
        requests[species] = flux*SECONDS_MYR/NA
    energy = p.escape_efficiency*p.xuv_flux_w_m2*p.radius_m/(4*G*p.mass_kg)*SECONDS_MYR
    requested_mass = sum(MOLAR_MASS[s]*n for s,n in requests.items())
    factor = min(1.,energy/requested_mass) if requested_mass else 1.
    for species,n in requests.items(): add(lost,elements_of({species:n*factor}))
    energy -= requested_mass*factor
    water = eq.gas.get('H2O',0.)
    if water and energy > 0:
        skin = max(120.,(snapshot.absorbed_w_m2/(2*5.670374419e-8))**.25)
        cold_pressure = max(eq.pressure_pa*.1,1e-12)
        xwater = min(water/total_moles,water_psat(skin)/cold_pressure)
        diffusion_mol = 2.5e17*min(1.,2*xwater)*SECONDS_MYR/NA
        lost['H'] += min(diffusion_mol,energy/ATOMIC_MASS['H'])
    return lost


def escape(p,state,snapshot,dt_myr,start_accessible=None,start_snapshot=None):
    if dt_myr < 0: raise ValueError('Negative escape interval')
    end_rates = escape_rates(p,snapshot)
    start = start_accessible if start_accessible is not None else state.accessible
    start_rates = escape_rates(p,start_snapshot) if start_snapshot is not None else end_rates
    lost = inventory()
    for e in ELEMENTS:
        available = state.accessible[e]
        old = start[e]
        if e == 'O':
            available=max(0.,available-state.accessible['Ca'])
            old=max(0.,old-start['Ca'])
        if available <= 0: continue
        # Exponential integration handles stiff trace-gas loss. Positive supply
        # is distributed across the interval, not dumped in before escape.
        rate,reference = (start_rates[e],old) if old > 1e-12 else (end_rates[e],available)
        if rate <= 0 or dt_myr == 0: continue
        z=rate/reference*dt_myr
        retained_fraction=math.exp(-min(z,745.))
        supply_fraction=-math.expm1(-z)/z if z else 1.
        base=min(old,available)
        supplied=available-base
        retained=base*retained_fraction+supplied*supply_fraction
        lost[e]=max(0.,min(available,available-retained))
    energy=p.escape_efficiency*p.xuv_flux_w_m2*p.radius_m/(4*G*p.mass_kg)*SECONDS_MYR*dt_myr
    requested_mass=sum(ATOMIC_MASS[e]*n for e,n in lost.items())
    if requested_mass > energy:
        lost={e:n*energy/requested_mass for e,n in lost.items()}
    add(state.accessible,lost,-1); add(state.escaped,lost)
    state.assert_conserved()
    return lost


def apply_event(p,state,snapshot,event):
    """Impact order: strip existing gas, deliver, release mantle, replace crust."""
    before = snapshot.equilibrium
    escaped_before = dict(state.escaped)
    moved = inventory()
    if event.kind == 'volcanic_rate':
        state.volcanic_multiplier = event.volcanic_multiplier
    else:
        if event.strip_fraction:
            stripped = elements_of({s:n*event.strip_fraction for s,n in before.gas.items()})
            add(state.accessible,stripped,-1); add(state.escaped,stripped)
        add(state.accessible,event.elements); add(state.delivered,event.elements)
        moved = transfer(state.mantle,state.accessible,event.mantle_release_fraction)
    renewed = resurface(p,state,snapshot,event.resurface_fraction)
    state.assert_conserved()
    after = solve_stack(p,state,before.temperature_k)
    record = {'event':asdict(event),'actual_resurfaced_fraction':renewed,
              'mantle_released_elements':moved,
              'escaped_elements':{e:state.escaped[e]-escaped_before[e] for e in ELEMENTS},
              'pressure_before_pa':before.pressure_pa,'pressure_after_pa':after.equilibrium.pressure_pa,
              'temperature_before_k':before.temperature_k,'temperature_after_k':after.equilibrium.temperature_k}
    return after,record


def sample_events(p,comets_per_myr=0.,impacts_per_myr=0.):
    """Poisson event times in physical time, independent of epoch/substep count.

    Event magnitudes are explicit gameplay priors, not an accretion simulation.
    Users may supply a hand-authored schedule instead. Same-time events retain
    their listed order, and events exactly on an epoch boundary occur first.
    """
    rng = random.Random(p.seed)
    events=[]
    for kind,rate in [('comet',comets_per_myr),('impact',impacts_per_myr)]:
        if not math.isfinite(rate) or rate < 0: raise ValueError('Invalid event rate')
        if not rate: continue
        t=rng.expovariate(rate)
        while t <= p.age_myr:
            water_mol=10**rng.uniform(3,5)
            elements={'H':2*water_mol,'O':water_mol,'C':water_mol*.01,'N':water_mol*.02,'S':water_mol*.001}
            events.append(Event(t,kind,f'{kind}-{len(events)+1}',elements,
                strip_fraction=rng.uniform(.1,.7) if kind=='impact' else 0.,
                mantle_release_fraction=rng.uniform(.001,.02) if kind=='impact' else 0.,
                resurface_fraction=rng.uniform(.01,.2) if kind=='impact' else 0.))
            t+=rng.expovariate(rate)
    return sorted(events,key=lambda e:e.time_myr)


def evolve(planet,events=(),progress=None):
    p=copy.deepcopy(planet)
    state=State.initial_state(p)
    schedule=copy.deepcopy(list(events))
    for event in schedule: event.validate(p.age_myr)
    schedule.sort(key=lambda e:e.time_myr)
    snapshot=solve_stack(p,state)
    initial=snapshot.to_dict()
    history=[]; event_log=[]; transport_log=[]; cursor=0
    for epoch in range(1,p.epochs+1):
        end=epoch*p.epoch_duration_myr
        # Process t=0 and endpoint events even if no time remains in the epoch.
        while True:
            if cursor < len(schedule) and schedule[cursor].time_myr <= state.age_myr+1e-10:
                snapshot,record=apply_event(p,state,snapshot,schedule[cursor])
                event_log.append(record);cursor+=1
                continue
            if state.age_myr >= end-1e-10: break
            tectonic=snapshot.tectonics_active
            vr=p.volcanic_fraction_per_myr*state.volcanic_multiplier
            rr=p.resurfacing_fraction_per_myr*state.volcanic_multiplier
            if tectonic:
                vr+=p.tectonic_volcanic_fraction_per_myr
                rr+=p.tectonic_resurfacing_fraction_per_myr
            event_time=schedule[cursor].time_myr if cursor<len(schedule) else math.inf
            target=min(end,event_time,state.age_myr+p.max_transport_step_myr)
            rates=escape_rates(p,snapshot)
            solids=elements_of({s:n for s,n in snapshot.equilibrium.condensed.items() if s not in ('water','ice')})
            for e,rate in rates.items():
                abundance=state.accessible[e]
                source=state.mantle[e]*vr
                if state.fresh_veneer_units > 0:
                    source+=rr*(p.crust.fresh_elements[e]-solids[e])
                net=abs(rate-source)
                significant=max(1e-8,(state.initial[e]+state.delivered[e])*1e-10)
                if net > 0 and abundance > significant:
                    target=min(target,state.age_myr+p.max_escape_fraction_per_step*abundance/net)
            dt=target-state.age_myr
            if dt <= 1e-12:
                raise RuntimeError('Transport timestep collapsed; refine the escape closure')
            start_accessible=dict(state.accessible)
            moved=transfer(state.mantle,state.accessible,-math.expm1(-vr*dt))
            renewed=resurface(p,state,snapshot,-math.expm1(-rr*dt))
            state.age_myr=target
            supplied=solve_stack(p,state,snapshot.equilibrium.temperature_k)
            loss=escape(p,state,supplied,dt,start_accessible,snapshot)
            snapshot=solve_stack(p,state,supplied.equilibrium.temperature_k) if sum(loss.values()) else supplied
            state.assert_conserved()
            transport_log.append({'start_myr':target-dt,'end_myr':target,'tectonics_active':tectonic,
                'volcanic_rate_per_myr':vr,'resurfacing_rate_per_myr':rr,
                'mantle_released_elements':moved,'escaped_elements':loss,'actual_resurfaced_fraction':renewed})
        row=snapshot.to_dict();row['epoch']=epoch
        row['reservoirs']={k:dict(getattr(state,k)) for k in ('accessible','mantle','fresh_crust','buried','escaped','delivered')}
        row['element_budget_error_mol_m2']=state.budget_error()
        history.append(row)
        if progress: progress(epoch,p.epochs,snapshot)
    return {'schema_version':2,'planet':asdict(p),'initial':initial,'epochs':history,
            'events':event_log,'transport':transport_log,'final_state':asdict(state)}
