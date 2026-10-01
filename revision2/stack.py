"""One energy/chemistry/pressure solution, followed by a read-only layer stack."""
from dataclasses import asdict, dataclass
import math
from scipy.optimize import brentq
from .thermo import equilibrium_column, MOLAR_MASS, DENSITY, R, EquilibriumError

SIGMA = 5.670374419e-8
MINERAL_NAMES = ('Fe','FeO','Fe3O4','Fe2O3','CaO')


class UnsupportedState(RuntimeError): pass


@dataclass
class Layer:
    material: str
    mass_kg_m2: float
    thickness_m: float
    bottom_pressure_pa: float
    top_pressure_pa: float
    bottom_temperature_k: float
    top_temperature_k: float


@dataclass
class Snapshot:
    age_myr: float
    equilibrium: object
    layers: list
    olr_w_m2: float
    absorbed_w_m2: float
    energy_residual_w_m2: float
    ocean_crust_contact: bool
    tectonics_active: bool
    oxidation_capacity_mol_o_m2: float
    surface_fractions: dict
    limitations: list

    def to_dict(self):
        out = asdict(self)
        out['equilibrium']['partial_pressures_pa'] = self.equilibrium.partial_pressures
        return out


def optical_depth(p, eq):
    return sum(p.opacity_m2_kg.get(s,0.)*MOLAR_MASS[s]*n for s,n in eq.gas.items())


def solve_stack(p, state, previous_temperature=None):
    absorbed = p.stellar_flux_w_m2*(1-p.albedo)/4 + p.geothermal_flux_w_m2
    cache = {}
    def evaluate(t):
        if t not in cache:
            eq = equilibrium_column(state.accessible,t,p.gravity)
            olr = SIGMA*t**4/(1+.75*optical_depth(p,eq))
            cache[t] = eq,olr,olr-absorbed
        return cache[t]
    # Continue the previous stable branch locally. An equilibrium-only repeat
    # is idempotent; no geological clock or arbitrary relaxation iterations.
    chosen = None
    if previous_temperature is not None:
        eq,olr,residual = evaluate(previous_temperature)
        if abs(residual) < 1e-6:
            chosen = eq,olr,residual
        else:
            direction = 1 if residual < 0 else -1
            last = previous_temperature
            for distance in (.05,.1,.25,.5,1.,2.,5.,10.):
                candidate = max(200.,min(1200.,previous_temperature+direction*distance))
                a,b = sorted((last,candidate))
                if a != b and evaluate(a)[2] <= 0 <= evaluate(b)[2]:
                    t = brentq(lambda t:evaluate(t)[2],a,b,xtol=1e-7)
                    eq,olr,error = evaluate(t)
                    if abs(error) <= .01:
                        chosen = eq,olr,error
                        break
                last = candidate
    if chosen is None:
        reference = previous_temperature if previous_temperature is not None else 300.
        grid = sorted(set([float(t) for t in range(200,1201,10)] +
                          [273.15,273.17,647.095,647.097,reference]))
        brackets = [(a,b) for a,b in zip(grid,grid[1:])
                    if evaluate(a)[2] <= 0 <= evaluate(b)[2]]
        brackets.sort(key=lambda ab: abs((ab[0]+ab[1])/2-reference))
        for a,b in brackets:
            t = brentq(lambda t:evaluate(t)[2],a,b,xtol=1e-7)
            eq,olr,residual = evaluate(t)
            if abs(residual) <= .01:
                chosen = eq,olr,residual
                break
    if chosen is None:
        raise UnsupportedState('No resolved stable energy-balance root in the supported 200–1200 K range')
    eq,olr,residual = chosen
    if eq.pressure_pa > 1e7:
        raise UnsupportedState('Pressure exceeds 100 bar: a nonideal atmosphere EOS is required')
    layers,limits = build_layers(p,eq)
    water_mass = eq.condensed.get('water',0)*MOLAR_MASS['water']
    contact = water_mass >= p.minimum_contact_water_kg_m2 and water_mass > 0
    # This revision's allowed condensed shell has no HP ice. Deep shells are
    # rejected in build_layers, not misidentified as ocean-rock contact.
    iron = sum(eq.condensed.get(s,0)*{'Fe':1,'FeO':1,'Fe3O4':3,'Fe2O3':2,'FeS':1,'FeS2':1}[s]
               for s in ('Fe','FeO','Fe3O4','Fe2O3','FeS','FeS2'))
    bound_o = sum(eq.condensed.get(s,0)*{'Fe':0,'FeO':1,'Fe3O4':4,'Fe2O3':3,'FeS':0,'FeS2':0}[s]
                  for s in ('Fe','FeO','Fe3O4','Fe2O3','FeS','FeS2'))
    sulfide_s = eq.condensed.get('FeS',0)+2*eq.condensed.get('FeS2',0)
    # Oxygen uptake to ferric iron and S(VI); bound sulfate is already oxidized.
    capacity = max(0.,1.5*iron-bound_o+3*sulfide_s)
    carbonate = eq.condensed.get('carbonate',0)*MOLAR_MASS['carbonate']
    silicate = p.crust.inert_veneer_mass+sum(eq.condensed.get(s,0)*MOLAR_MASS[s] for s in MINERAL_NAMES)
    sulfide = sum(eq.condensed.get(s,0)*MOLAR_MASS[s] for s in ('FeS','FeS2'))
    sulfate = eq.condensed.get('sulfate',0)*MOLAR_MASS['sulfate']
    total = silicate+carbonate+sulfide+sulfate
    fractions = {s:getattr(p.crust,s)*silicate/total for s in ('felsic','mafic','ultramafic')}
    fractions.update(carbonate=carbonate/total,sulfide=sulfide/total,sulfate=sulfate/total)
    return Snapshot(state.age_myr,eq,layers,olr,absorbed,residual,contact,
                    p.tectonics_enabled and contact,capacity,fractions,limits)


def build_layers(p, eq):
    """Bottom-to-top shell. All interface P/T agree; no reservoir mutation.

    Condensed shells are isothermal in this first closure, as is the diagnostic
    atmosphere. Thermal transport/EOS through a deep shell is not invented here.
    """
    t = eq.temperature_k
    upper_layers = []
    for name in ('carbonate','FeS','FeS2','sulfate','graphite','sulfur','water','ice'):
        m = eq.condensed.get(name,0)*MOLAR_MASS[name]
        if m > 1e-9: upper_layers.append((name,m,m/DENSITY[name]))
    basal_pressure = eq.pressure_pa+p.gravity*sum(m for _,m,_ in upper_layers)
    if basal_pressure > 2e8:
        raise UnsupportedState('Condensed-shell base exceeds 0.2 GPa: high-pressure water/mineral EOS required')
    structural_mass = p.crust.thickness_m*p.crust.density_kg_m3
    # Structural crust is a fixed inert scaffold; the active veneer is a small
    # separately accounted chemical reservoir at its top.
    active_mass = p.crust.inert_veneer_mass+sum(eq.condensed.get(s,0)*MOLAR_MASS[s] for s in MINERAL_NAMES)
    layers = [Layer('silicate crust (inert scaffold)',structural_mass,p.crust.thickness_m,
                    basal_pressure+p.gravity*(structural_mass+active_mass),
                    basal_pressure+p.gravity*active_mass,t,t),
              Layer('reactive silicate veneer',active_mass,active_mass/p.crust.density_kg_m3,
                    basal_pressure+p.gravity*active_mass,basal_pressure,t,t)]
    pressure = basal_pressure
    for name,m,depth in upper_layers:
        top = pressure-p.gravity*m
        layers.append(Layer(name,m,depth,pressure,top,t,t)); pressure=top
    if eq.gas_mass > 0:
        mean_mu = eq.gas_mass/sum(eq.gas.values())
        height = R*t/(mean_mu*p.gravity)
        # Last bin includes the exponentially thin tail in its column mass;
        # reported thickness is the altitude enclosing 99.9999% of the column.
        for i in range(8):
            bottom = eq.pressure_pa*math.exp(-14*i/8)
            top = eq.pressure_pa*math.exp(-14*(i+1)/8) if i < 7 else 0.
            layers.append(Layer('atmosphere', (bottom-top)/p.gravity,14*height/8,bottom,top,t,t))
    limits = ['Grey radiation with prescribed albedo and absorption coefficients.',
              'Isothermal shell and hydrostatic atmospheric geometry; no clouds or circulation.',
              'Pure-phase equilibrium: no aqueous ions, salinity or weathering kinetics.',
              'Reactive CaO is a finite alkalinity proxy, not a full silicate mineral assemblage.',
              'Surface fractions are mass fractions of the active substrate, not mapped area coverage.']
    if t < 300 and any(eq.elements[e] > 0 for e in ('Fe','Ca','S')):
        limits.append('Some mineral/sulfur thermochemical fits are extrapolated below their ~300 K lower bounds.')
    if eq.pressure_pa > 1e6: limits.append('Ideal-gas approximation above 10 bar is qualitative.')
    return layers, limits
