"""Restricted Gibbs equilibrium, using Cantera's GRI gas data and pure phases.

All public inventories are mol/m², never molecular mass fractions. Pressure is
Pa. Pure water chemical potentials are matched to saturation vapor pressure.
The iron and calcium phases are finite buffers, not a complete mineral assemblage.
"""
from dataclasses import dataclass
import math
import cantera as ct
import numpy as np
from scipy.optimize import brentq

R = ct.gas_constant / 1000
BAR = 100000.0
ELEMENTS = ('H', 'C', 'N', 'O', 'Ar', 'He', 'Fe', 'Ca', 'S')
ATOMIC_MASS = {'H': .001008, 'C': .012011, 'N': .014007,
               'O': .015999, 'Ar': .03995, 'He': .004002602, 'Fe': .055845, 'Ca': .040078, 'S': .03206}
GAS_NAMES = ('H2', 'H', 'O2', 'O', 'OH', 'H2O', 'CO', 'CO2', 'C',
             'CH4', 'N2', 'N', 'NH3', 'AR', 'He', 'S', 'S2', 'H2S', 'SO2', 'SO3')
SOURCE = {s.name: s for s in ct.Species.list_from_file('gri30.yaml')}
SOURCE.update({s.name:s for s in ct.Species.list_from_file('nasa_gas.yaml')
               if s.name in ('S','S2','H2S','SO2','SO3')})
SOURCE['He'] = ct.Species.from_dict({'name': 'He', 'composition': {'He': 1},
    'thermo': {'model': 'constant-cp', 'T0': 298.15, 'h0': 0,
               's0': 126150., 'cp0': 20786.}})
COMPOSITION = {k: dict(SOURCE[k].composition) for k in GAS_NAMES}
COMPOSITION.update({'water': {'H': 2, 'O': 1}, 'ice': {'H': 2, 'O': 1},
                    'graphite': {'C': 1}, 'FeO': {'Fe': 1, 'O': 1},
                    'Fe2O3': {'Fe': 2, 'O': 3}, 'Fe3O4': {'Fe': 3, 'O': 4},
                    'Fe': {'Fe': 1}, 'CaO': {'Ca': 1, 'O': 1},
                    'carbonate': {'Ca': 1, 'C': 1, 'O': 3}, 'FeS': {'Fe':1,'S':1},
                    'FeS2': {'Fe':1,'S':2}, 'sulfate': {'Ca':1,'S':1,'O':4},
                    'sulfur': {'S':1}})
MOLAR_MASS = {s: sum(ATOMIC_MASS[e]*a for e,a in c.items())
              for s,c in COMPOSITION.items()}
DENSITY = {'water': 1000., 'ice': 917., 'graphite': 2200., 'FeO': 5745., 'Fe2O3': 5250., 'Fe3O4': 5170., 'Fe': 7874., 'CaO': 3340., 'carbonate': 2710., 'FeS':4840., 'FeS2':5010., 'sulfate':2960., 'sulfur':2000.}


def elements_of(species):
    result = {e: 0. for e in ELEMENTS}
    for s,n in species.items():
        for e,a in COMPOSITION[s].items(): result[e] += n*a
    return result


def mass(elements):
    return sum(ATOMIC_MASS[e]*n for e,n in elements.items())


def water_psat(t):
    """Ice: Murphy & Koop 2005; liquid: Wagner/IAPWS curve. Pa."""
    if t < 273.16:
        return math.exp(9.550426 - 5723.265/t + 3.53068*math.log(t) - .00728332*t)
    if t >= 647.096: return math.inf
    z = 1-t/647.096
    f = (-7.85951783*z + 1.84408259*z**1.5 - 11.7866497*z**3
         + 22.6807411*z**3.5 - 15.9618719*z**4 + 1.80122502*z**7.5)
    return 22.064e6*math.exp(647.096/t*f)


def pure_phase(name, t, g, density):
    # Match g exactly at this evaluation temperature; no extrapolation of this
    # constructed phase is used. Cantera handles the condensed pressure term.
    sp = ct.Species.from_dict({'name': name, 'composition': COMPOSITION[name],
        'thermo': {'model': 'constant-cp', 'T0': t, 'h0': g*1000,
                   's0': 0., 'cp0': 50000.},
        'equation-of-state': {'model': 'constant-volume', 'density': density}})
    return ct.Solution(thermo='fixed-stoichiometry', species=[sp])


# NASA condensed species are distributed with Cantera (nasa_condensed.yaml).
# Minerals have ~298–300 K lower bounds; colder evaluations explicitly report
# extrapolation in the snapshot. Never silently present them as calibrated.
SOLIDS = {s.name:s for s in ct.Species.list_from_file('nasa_condensed.yaml')}
MINERALS = {'FeO':'FeO(s)', 'Fe2O3':'Fe2O3(s)', 'Fe3O4':'Fe3O4(s)',
            'CaO':'CaO(s)', 'carbonate':'CaCO3(caL)', 'FeS2':'FeS2(s)', 'sulfate':'CaSO4(s)'}


def mineral_g(name, t):
    key = ('Fe(a)' if t < 1184 else 'Fe(c)') if name == 'Fe' else MINERALS.get(name,'')
    if name == 'FeS': key = 'FeS(a)' if t < 411 else ('FeS(b)' if t < 598 else 'FeS(c)')
    if name == 'sulfur': key = 'S(cr1)' if t < 368.3 else ('S(cr2)' if t < 388.36 else 'S(L)')
    thermo = SOLIDS[key].thermo
    return (thermo.h(t)-t*thermo.s(t))/1000


class EquilibriumError(RuntimeError): pass


@dataclass
class Equilibrium:
    temperature_k: float
    pressure_pa: float
    gas: dict
    condensed: dict
    pressure_relative_error: float = 0.

    @property
    def gas_mass(self): return sum(MOLAR_MASS[s]*n for s,n in self.gas.items())

    @property
    def partial_pressures(self):
        total = sum(self.gas.values())
        return {s: self.pressure_pa*n/total if total else 0. for s,n in self.gas.items()}

    @property
    def elements(self): return elements_of({**self.gas, **self.condensed})


def equilibrium_tp(inventory, t, pressure):
    """Global minimum of the restricted ideal-gas/pure-phase Gibbs model."""
    if not 200 <= t <= 1200: raise ValueError('Thermochemical range is 200–1200 K')
    if pressure <= 0: raise ValueError('Equilibrium pressure must be positive')
    if any(e not in ELEMENTS or not math.isfinite(n) or n < 0 for e,n in inventory.items()):
        raise ValueError('Invalid elemental inventory')
    if inventory.get('O',0) < inventory.get('Ca',0):
        raise ValueError('The calcium buffer requires at least one O per Ca')
    scale = sum(inventory.values())
    if not scale: return Equilibrium(t,0.,{}, {})
    gas = ct.Solution(thermo='ideal-gas', species=[SOURCE[k] for k in GAS_NAMES])
    gas.TP = t, pressure
    phases = [(gas, 1.)]
    if inventory.get('H',0) and inventory.get('O',0) and t < 647.096:
        name = 'ice' if t < 273.16 else 'water'
        # NASA standard chemical potentials refer to one atmosphere.
        gv = gas.standard_gibbs_RT[gas.species_index('H2O')]*R*t
        # standard_gibbs_RT already contains the ideal gas pressure correction.
        gw = gv + R*t*math.log(water_psat(t)/pressure)
        gw -= MOLAR_MASS[name]/DENSITY[name]*(water_psat(t)-ct.one_atm)
        phases.append((pure_phase(name,t,gw,DENSITY[name]),0.))
    if inventory.get('C',0): phases.append((ct.Solution('graphite.yaml'),0.))
    if inventory.get('Fe',0):
        for name in ('Fe','FeO','Fe3O4','Fe2O3') + (('FeS','FeS2') if inventory.get('S',0) else ()):
            phases.append((pure_phase(name,t,mineral_g(name,t),DENSITY[name]),0.))
    if inventory.get('Ca',0):
        for name in ('CaO','carbonate') + (('sulfate',) if inventory.get('S',0) else ()):
            phases.append((pure_phase(name,t,mineral_g(name,t),DENSITY[name]),0.))
    if inventory.get('S',0):
        phases.append((pure_phase('sulfur',t,mineral_g('sulfur',t),DENSITY['sulfur']),0.))
    mix = ct.Mixture(phases)
    mix.T, mix.P = t, pressure
    initial = np.zeros(mix.n_species)
    for e,n in inventory.items():
        if e in ('Fe','Ca'): continue
        if e == 'O': n -= inventory.get('Ca',0)
        atom = 'AR' if e == 'Ar' else e
        initial[mix.species_index(0,atom)] = n/scale
    if inventory.get('Fe',0):
        idx = next(i for i in range(mix.n_phases) if mix.phase(i).species_names == ['Fe'])
        initial[mix.species_index(idx,'Fe')] = inventory['Fe']/scale
    if inventory.get('Ca',0):
        idx = next(i for i in range(mix.n_phases) if mix.phase(i).species_names == ['CaO'])
        initial[mix.species_index(idx,'CaO')] = inventory['Ca']/scale
    # Equilibrate the gas before allowing condensed phases to appear. Starting
    # from atomic H/O can incorrectly eliminate the whole gas phase in VCS.
    if initial[:gas.n_species].sum() > 0:
        gas.TPX = t, pressure, initial[:gas.n_species]
        initial_mass = sum(initial[i]*MOLAR_MASS[k] for i,k in enumerate(gas.species_names))
        gas.equilibrate('TP', solver='auto', rtol=1e-11, max_steps=3000)
        initial[:gas.n_species] = gas.X * initial_mass / (gas.mean_molecular_weight/1000)
    mix.species_moles = initial
    try:
        mix.equilibrate('TP',solver='vcs',rtol=1e-10,max_steps=3000,log_level=0)
    except ct.CanteraError as exc:
        raise EquilibriumError(f'Gibbs solve failed at {t:.2f} K, {pressure:.4g} Pa') from exc
    values = mix.species_moles*scale
    gases = dict(zip(gas.species_names, values[:gas.n_species].tolist()))
    condensed = {}
    offset = gas.n_species
    for phase,_ in phases[1:]:
        name = 'graphite' if phase.species_names == ['C(gr)'] else phase.species_names[0]
        condensed[name] = float(values[offset]); offset += 1
    result = Equilibrium(t,pressure,gases,condensed)
    for e in ELEMENTS:
        if abs(result.elements[e]-inventory.get(e,0)) > 2e-8*max(inventory.get(e,0),1.):
            raise EquilibriumError(f'Element residual for {e}')
    return result


def equilibrium_column(inventory, t, gravity):
    """Solve surface P = g * atmospheric mass, including condensation.

    A pure condensable has a discontinuity in the fixed-TP phase amount at
    saturation. At that boundary the lever rule selects the conserved mixture
    whose atmospheric weight supplies precisely the coexistence pressure.
    """
    if gravity <= 0: raise ValueError('Gravity must be positive')
    upper = max(gravity*mass(inventory),1e-12)
    if not sum(inventory.values()): return Equilibrium(t,0.,{}, {})
    cache = {}
    def at(logp):
        if logp not in cache:
            p = math.exp(logp)
            q = equilibrium_tp(inventory,t,p)
            cache[logp] = (q,gravity*q.gas_mass-p)
        return cache[logp]
    lo,hi = math.log(1e-12),math.log(max(upper*1.000001,1e-11))
    qlo,flo = at(lo)
    if flo <= 0:
        qlo.pressure_pa = gravity*qlo.gas_mass
        return qlo
    root = brentq(lambda x: at(x)[1],lo,hi,xtol=1e-10,rtol=1e-12)
    q,residual = at(root)
    if abs(residual) > 1e-5*max(q.pressure_pa,1e-8):
        # Bracket arbitrarily closely in P; blend coexisting phases, not distinct
        # climates. Both endpoints conserve the identical elemental inventory.
        a,fa = at(root-1e-7); b,fb = at(root+1e-7)
        if fa*fb > 0: raise EquilibriumError('Unresolved pressure/phase boundary')
        w = -fb/(fa-fb)
        def blend(x,y): return {k:w*x.get(k,0)+(1-w)*y.get(k,0) for k in x.keys()|y.keys()}
        q = Equilibrium(t,math.exp(root),blend(a.gas,b.gas),blend(a.condensed,b.condensed))
    q.pressure_relative_error = abs(gravity*q.gas_mass-q.pressure_pa)/max(q.pressure_pa,1e-8)
    if q.pressure_relative_error > 2e-5: raise EquilibriumError('Hydrostatic pressure did not converge')
    return q
