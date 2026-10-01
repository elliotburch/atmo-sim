"""Configuration and conserved reservoirs. Inputs use SI and mol of atoms/m²."""
from dataclasses import dataclass, field
import math
from .thermo import ELEMENTS, elements_of, mass


def inventory(values=None):
    out = dict.fromkeys(ELEMENTS, 0.)
    for e,n in (values or {}).items():
        if e not in out or isinstance(n,bool) or not isinstance(n,(int,float)) or not math.isfinite(n) or n < 0:
            raise ValueError(f'Invalid elemental abundance: {e}={n}')
        out[e] = float(n)
    return out


def add(target, source, factor=1.):
    for e in ELEMENTS: target[e] += source.get(e,0.)*factor


def transfer(source, destination, fraction):
    if not 0 <= fraction <= 1: raise ValueError('Transfer fraction outside [0,1]')
    moved = {e:source[e]*fraction for e in ELEMENTS}
    add(source,moved,-1); add(destination,moved)
    return moved


@dataclass
class Crust:
    felsic: float = .2
    mafic: float = .65
    ultramafic: float = .15
    # Exchangeable surface veneer, not the whole structural crust.
    active_silicate_kg_m2: float = 1000.
    feo_mol_m2: float = 1000.
    cao_mol_m2: float = 100.
    fes_mol_m2: float = 100.
    fresh_veneer_replacements: float = 20.
    thickness_m: float = 35000.
    density_kg_m3: float = 3000.

    def validate(self):
        for name,value in vars(self).items():
            if not math.isfinite(value) or value < 0: raise ValueError(f'Invalid crust {name}')
        if abs(self.felsic+self.mafic+self.ultramafic-1) > 1e-9:
            raise ValueError('Felsic, mafic and ultramafic fractions must sum to one')
        if self.active_silicate_kg_m2 <= mass(self.fresh_elements):
            raise ValueError('Active silicate mass must exceed its reactive mineral mass')
        if self.density_kg_m3 <= 0 or self.thickness_m <= 0:
            raise ValueError('Crust thickness and density must be positive')

    @property
    def fresh_elements(self):
        return elements_of({'FeO':self.feo_mol_m2,'CaO':self.cao_mol_m2,'FeS':self.fes_mol_m2})

    @property
    def inert_veneer_mass(self):
        return self.active_silicate_kg_m2-mass(self.fresh_elements)


@dataclass
class Planet:
    name: str = 'Crust-up reference'
    mass_kg: float = 5.972e24
    radius_m: float = 6.371e6
    core_mass_fraction: float = .32
    stellar_flux_w_m2: float = 1600.
    albedo: float = .25
    geothermal_flux_w_m2: float = .08
    # Explicit grey absorption coefficients, m²/kg; NOT fitted planet classes.
    opacity_m2_kg: dict = field(default_factory=lambda:{'H2O':.002,'CO2':.0005,'CH4':.002,'H2':.00002,'NH3':.002})
    primordial_elements: dict = field(default_factory=lambda:{'H':4e5,'O':2.05e5,'C':1e3,'N':4e5,'Ar':1000.})
    mantle_volatiles: dict = field(default_factory=lambda:{'H':4e6,'O':2.02e6,'C':1e4,'N':2e5})
    crust: Crust = field(default_factory=Crust)
    epochs: int = 20
    epoch_duration_myr: float = 10.
    max_transport_step_myr: float = 10.
    max_escape_fraction_per_step: float = .05
    volcanic_fraction_per_myr: float = .0001
    resurfacing_fraction_per_myr: float = .0001
    tectonics_enabled: bool = True
    tectonic_volcanic_fraction_per_myr: float = .0002
    tectonic_resurfacing_fraction_per_myr: float = .0005
    minimum_contact_water_kg_m2: float = 1.
    escape_enabled: bool = True
    exosphere_temperature_k: float = 800.
    collision_cross_section_m2: float = 3e-19
    xuv_flux_w_m2: float = .0046
    escape_efficiency: float = .15
    seed: int = 1

    @property
    def gravity(self): return 6.67430e-11*self.mass_kg/self.radius_m**2

    @property
    def age_myr(self): return self.epochs*self.epoch_duration_myr

    def validate(self):
        self.crust.validate()
        if not self.name: raise ValueError('Planet requires a name')
        if type(self.epochs) is not int or self.epochs < 1: raise ValueError('Epochs must be a positive integer')
        if type(self.seed) is not int: raise ValueError('Seed must be an integer')
        if type(self.tectonics_enabled) is not bool or type(self.escape_enabled) is not bool:
            raise ValueError('Activity switches must be booleans')
        for name,value in vars(self).items():
            if isinstance(value,(float,int)) and not isinstance(value,bool):
                if not math.isfinite(value) or (name != 'seed' and value < 0): raise ValueError(f'Invalid {name}')
        for name in ('mass_kg','radius_m','epoch_duration_myr','max_transport_step_myr',
                     'collision_cross_section_m2','exosphere_temperature_k'):
            if getattr(self,name) <= 0: raise ValueError(f'{name} must be positive')
        if not 0 < self.max_escape_fraction_per_step <= .25: raise ValueError('Escape step fraction must be in (0,0.25]')
        if not 0 <= self.core_mass_fraction < 1: raise ValueError('Invalid differentiated core fraction')
        if not 0 <= self.albedo < 1 or not 0 <= self.escape_efficiency <= 1:
            raise ValueError('Albedo/escape efficiency outside range')
        from .thermo import GAS_NAMES
        for s,k in self.opacity_m2_kg.items():
            if s not in GAS_NAMES or not math.isfinite(k) or k < 0: raise ValueError('Invalid opacity')
        self.primordial_elements = inventory(self.primordial_elements)
        self.mantle_volatiles = inventory(self.mantle_volatiles)
        if any(self.primordial_elements[e] or self.mantle_volatiles[e] for e in ('Fe','Ca')):
            raise ValueError('Supply crustal Fe/Ca through the crust, not volatile reservoirs')


@dataclass
class Event:
    time_myr: float
    kind: str
    label: str = ''
    elements: dict = field(default_factory=dict)
    strip_fraction: float = 0.
    mantle_release_fraction: float = 0.
    resurface_fraction: float = 0.
    volcanic_multiplier: float = 1.

    def validate(self, age):
        if self.kind not in ('comet','impact','volcanic_rate','resurfacing'): raise ValueError('Unknown event kind')
        if not math.isfinite(self.time_myr) or not 0 <= self.time_myr <= age: raise ValueError('Event outside system age')
        self.elements = inventory(self.elements)
        if self.elements['Fe'] or self.elements['Ca']: raise ValueError('Impact delivery currently supports volatile elements only')
        for k in ('strip_fraction','mantle_release_fraction','resurface_fraction'):
            if not math.isfinite(getattr(self,k)) or not 0 <= getattr(self,k) <= 1: raise ValueError(f'Invalid {k}')
        if not math.isfinite(self.volcanic_multiplier) or self.volcanic_multiplier < 0: raise ValueError('Invalid volcanic multiplier')
        # Do not silently ignore supplied effects on the wrong event type.
        if self.kind == 'volcanic_rate' and (sum(self.elements.values()) or self.strip_fraction or self.mantle_release_fraction or self.resurface_fraction):
            raise ValueError('A volcanic_rate event only sets the rate multiplier')
        if self.kind == 'comet' and (self.strip_fraction or self.mantle_release_fraction or self.resurface_fraction):
            raise ValueError('Use impact for delivery plus stripping/outgassing/resurfacing')
        if self.kind == 'resurfacing' and (sum(self.elements.values()) or self.strip_fraction or self.mantle_release_fraction):
            raise ValueError('A resurfacing event only replaces crust')


@dataclass
class State:
    accessible: dict
    mantle: dict
    fresh_crust: dict
    buried: dict
    escaped: dict
    delivered: dict
    initial: dict
    fresh_veneer_units: float
    buried_veneer_units: float = 0.
    volcanic_multiplier: float = 1.
    age_myr: float = 0.

    @classmethod
    def initial_state(cls,p):
        p.validate()
        accessible = inventory(p.primordial_elements)
        add(accessible,p.crust.fresh_elements)
        fresh = inventory()
        add(fresh,p.crust.fresh_elements,p.crust.fresh_veneer_replacements)
        initial = inventory(accessible)
        add(initial,p.mantle_volatiles); add(initial,fresh)
        column_mass = p.mass_kg/(4*math.pi*p.radius_m**2)
        modeled_shell = (mass(initial)+p.crust.thickness_m*p.crust.density_kg_m3
                         +p.crust.inert_veneer_mass*(1+p.crust.fresh_veneer_replacements))
        if modeled_shell >= (1-p.core_mass_fraction)*column_mass:
            raise ValueError('Reservoirs and crust exceed the differentiated silicate-shell mass')
        return cls(accessible,inventory(p.mantle_volatiles),fresh,inventory(),inventory(),
                   inventory(),initial,p.crust.fresh_veneer_replacements)

    def budget_error(self):
        error = inventory()
        for e in ELEMENTS:
            total = sum(getattr(self,k)[e] for k in ('accessible','mantle','fresh_crust','buried','escaped'))
            expected = self.initial[e]+self.delivered[e]
            error[e] = total-expected
        return error

    def assert_conserved(self):
        for k in ('accessible','mantle','fresh_crust','buried','escaped','delivered'):
            for e,n in getattr(self,k).items():
                if not math.isfinite(n) or n < -1e-8: raise AssertionError(f'Invalid {k}/{e}: {n}')
        for e,error in self.budget_error().items():
            if abs(error) > 1e-8*max(1.,self.initial[e]+self.delivered[e]):
                raise AssertionError(f'Elemental budget failed: {e}, {error}')
