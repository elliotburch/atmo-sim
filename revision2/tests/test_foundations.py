import copy
import math
import unittest
from revision2.model import Planet,Crust,State,Event,inventory,add
from revision2.thermo import (equilibrium_column,equilibrium_tp,elements_of,ELEMENTS,
                              MOLAR_MASS,water_psat)
from revision2.stack import solve_stack,build_layers,UnsupportedState
from revision2.evolution import evolve,apply_event,resurface,escape,sample_events


def quiet(**kwargs):
    return Planet(escape_enabled=False,volcanic_fraction_per_myr=0.,
                  resurfacing_fraction_per_myr=0.,tectonic_volcanic_fraction_per_myr=0.,
                  tectonic_resurfacing_fraction_per_myr=0.,**kwargs)


class ChemistryTests(unittest.TestCase):
    def check_elements(self,q,inv):
        for e in ELEMENTS:
            self.assertAlmostEqual(q.elements[e]/max(inv.get(e,0),1),inv.get(e,0)/max(inv.get(e,0),1),places=7,msg=e)

    def test_conservation_across_temperatures_and_redox(self):
        for t in (220.,300.,600.,1100.):
            for oxygen in (1500.,20000.):
                inv={'H':10000.,'C':1000.,'O':oxygen,'N':3000.,'Fe':200.,'Ca':20.,'S':100.,'Ar':3.}
                with self.subTest(t=t,oxygen=oxygen):
                    q=equilibrium_column(inv,t,9.81)
                    self.check_elements(q,inv)
                    self.assertLess(q.pressure_relative_error,2e-5)
                    self.assertTrue(all(n>=0 for n in [*q.gas.values(),*q.condensed.values()]))

    def test_pure_water_coexistence(self):
        q=equilibrium_column({'H':2e6,'O':1e6},300.,9.81)
        self.assertGreater(q.condensed['water'],0)
        self.assertGreater(q.gas['H2O'],0)
        self.assertAlmostEqual(q.pressure_pa/water_psat(300.),1.,delta=.005)
        self.check_elements(q,{'H':2e6,'O':1e6})

    def test_partial_pressure_saturation_in_mixture(self):
        q=equilibrium_column({'H':2e6,'O':1e6,'N':1e6},300.,9.81)
        self.assertAlmostEqual(q.partial_pressures['H2O']/water_psat(300.),1.,delta=.005)
        self.assertGreater(q.pressure_pa,q.partial_pressures['H2O'])

    def test_zero_water_never_generates_hydrogen(self):
        q=equilibrium_column({'N':10000.,'C':1000.,'O':2000.},300.,9.81)
        self.assertEqual(q.elements['H'],0.)
        self.assertEqual(q.gas['H2O'],0.)
        self.assertNotIn('water',q.condensed)

    def test_carbonate_limited_by_calcium(self):
        q=equilibrium_column({'Ca':100.,'O':10000.,'C':1000.,'N':1000.},320.,9.81)
        self.assertAlmostEqual(q.condensed['carbonate'],100.,places=5)
        self.assertAlmostEqual(q.elements['C'],1000.,places=6)
        self.assertGreater(q.gas['CO2'],899.)

    def test_sulfides_form_without_inventing_iron_or_sulfur(self):
        inv={'Fe':100.,'S':100.,'N':1000.}
        q=equilibrium_column(inv,350.,9.81)
        self.check_elements(q,inv)
        self.assertGreater(q.condensed['FeS']+q.condensed['FeS2'],1.)

    def test_finite_iron_capacity_leaves_excess_oxygen(self):
        inv={'Fe':100.,'O':200.,'N':1000.}
        q=equilibrium_column(inv,320.,9.81)
        self.assertAlmostEqual(q.condensed['Fe2O3'],50.,places=5)
        self.assertAlmostEqual(q.gas['O2'],25.,places=5)
        self.check_elements(q,inv)

    def test_unsupported_chemistry_temperature_rejected(self):
        with self.assertRaises(ValueError):equilibrium_tp({'N':100.},150.,1000.)


class StackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=quiet(stellar_flux_w_m2=1600.)
        cls.state=State.initial_state(cls.p)
        cls.snapshot=solve_stack(cls.p,cls.state)

    def test_energy_hydrostatic_and_element_closure(self):
        q=self.snapshot
        self.assertLess(abs(q.energy_residual_w_m2),.01)
        self.assertAlmostEqual(q.equilibrium.pressure_pa/self.p.gravity,q.equilibrium.gas_mass,delta=1e-3)
        for e in ELEMENTS:
            self.assertAlmostEqual(q.equilibrium.elements[e],self.state.accessible[e],delta=1e-5)

    def test_layer_interfaces_and_unique_water(self):
        layers=self.snapshot.layers
        for a,b in zip(layers,layers[1:]):
            self.assertAlmostEqual(a.top_pressure_pa,b.bottom_pressure_pa,delta=1e-5)
            self.assertEqual(a.top_temperature_k,b.bottom_temperature_k)
        self.assertLessEqual(sum(l.material=='water' for l in layers),1)
        atmosphere=sum(l.mass_kg_m2 for l in layers if l.material=='atmosphere')
        self.assertAlmostEqual(atmosphere,self.snapshot.equilibrium.gas_mass,delta=.001)

    def test_stack_does_not_modify_inventories(self):
        before=copy.deepcopy(self.state)
        solve_stack(self.p,self.state,self.snapshot.equilibrium.temperature_k)
        self.assertEqual(self.state,before)

    def test_map_fractions_normalized_and_preserve_silicate_ratio(self):
        f=self.snapshot.surface_fractions
        self.assertAlmostEqual(sum(f.values()),1.)
        self.assertAlmostEqual(f['felsic']/f['mafic'],self.p.crust.felsic/self.p.crust.mafic)
        self.assertGreaterEqual(f['sulfide'],0.)
        self.assertGreaterEqual(f['carbonate'],0.)

    def test_ocean_gate_wet_dry_frozen_and_disabled(self):
        self.assertTrue(self.snapshot.tectonics_active)
        cases=[(dict(primordial_elements={'N':4e5,'O':1000.}),False),
               (dict(stellar_flux_w_m2=1000.),False),
               (dict(tectonics_enabled=False),False)]
        for options,expected in cases:
            p=quiet(**options);q=solve_stack(p,State.initial_state(p))
            self.assertEqual(q.tectonics_active,expected)

    def test_deep_shell_rejected_without_fake_contact(self):
        q=copy.deepcopy(self.snapshot.equilibrium)
        q.condensed['water']=1e12
        with self.assertRaises(UnsupportedState):build_layers(self.p,q)


class EvolutionTests(unittest.TestCase):
    def test_static_epochs_do_not_age_chemistry(self):
        p=quiet(epochs=2,epoch_duration_myr=3.,max_transport_step_myr=3.)
        result=evolve(p)
        self.assertEqual([q['age_myr'] for q in result['epochs']],[3.,6.])
        for q in result['epochs']:
            self.assertAlmostEqual(q['equilibrium']['temperature_k'],result['initial']['equilibrium']['temperature_k'],places=5)
            self.assertEqual(q['reservoirs']['escaped'],inventory())
            self.assertTrue(all(abs(v)<1e-6 for v in q['element_budget_error_mol_m2'].values()))

    def test_comet_then_impact_preserves_full_budget(self):
        p=quiet(epochs=1,epoch_duration_myr=2.,max_transport_step_myr=2.)
        events=[Event(0.,'comet',elements={'H':20000.,'O':10000.}),
                Event(2.,'impact',strip_fraction=.2,mantle_release_fraction=.001,resurface_fraction=.1)]
        result=evolve(p,events)
        self.assertEqual(len(result['events']),2)
        self.assertEqual(result['final_state']['delivered']['H'],20000.)
        self.assertGreater(result['final_state']['escaped']['N'],0)
        self.assertTrue(all(abs(v)<1e-6 for v in result['epochs'][0]['element_budget_error_mol_m2'].values()))
        self.assertEqual(result['events'][1]['event']['time_myr'],2.)

    def test_stripping_removes_only_existing_gas(self):
        p=quiet();state=State.initial_state(p);q=solve_stack(p,state)
        before=dict(state.accessible)
        gas_elements=elements_of(q.equilibrium.gas)
        q,record=apply_event(p,state,q,Event(0.,'impact',strip_fraction=.5,elements=inventory()))
        for e in ELEMENTS:
            self.assertAlmostEqual(state.accessible[e],before[e]-.5*gas_elements[e],delta=1e-6)

    def test_resurfacing_buries_old_oxygen_and_uses_finite_stock(self):
        p=quiet(crust=Crust(fresh_veneer_replacements=.2));state=State.initial_state(p);q=solve_stack(p,state)
        renewed=resurface(p,state,q,.5)
        self.assertEqual(renewed,.2)
        self.assertGreater(state.buried['O'],0.)
        self.assertEqual(state.fresh_veneer_units,0.)
        self.assertEqual(resurface(p,state,q,.5),0.)
        state.assert_conserved()

    def test_hydrogen_escape_retains_oxygen(self):
        p=quiet(stellar_flux_w_m2=1800.);p.escape_enabled=True
        state=State.initial_state(p);q=solve_stack(p,state)
        oxygen_before=state.accessible['O']
        lost=escape(p,state,q,.01)
        self.assertGreater(lost['H'],0.)
        self.assertLess(lost['O'],1e-10)
        self.assertAlmostEqual(state.accessible['O'],oxygen_before,places=7)
        state.assert_conserved()

    def test_continuous_escape_supply_obeys_energy_ceiling(self):
        from revision2.evolution import SECONDS_MYR,G
        from revision2.thermo import mass
        p=quiet();p.escape_enabled=True
        state=State.initial_state(p);q=solve_stack(p,state)
        before=dict(state.accessible)
        delivered=inventory({'H':1e7,'O':5e6})
        add(state.accessible,delivered);add(state.delivered,delivered)
        # Use a prescribed existing snapshot to isolate the transport ledger.
        lost=escape(p,state,q,.01,before,q)
        ceiling=p.escape_efficiency*p.xuv_flux_w_m2*p.radius_m/(4*G*p.mass_kg)*SECONDS_MYR*.01
        self.assertLessEqual(mass(lost),ceiling*(1+1e-12))
        state.assert_conserved()

    def test_seeded_event_timeline_independent_of_epoch_count(self):
        a=quiet(epochs=20,epoch_duration_myr=10.)
        b=quiet(epochs=40,epoch_duration_myr=5.)
        self.assertEqual(sample_events(a,.03,.01),sample_events(b,.03,.01))

    def test_tectonic_rates_are_gated_not_background_rates(self):
        p=Planet(epochs=1,epoch_duration_myr=.01,max_transport_step_myr=.01,escape_enabled=False,
                 primordial_elements={'N':4e5},volcanic_fraction_per_myr=.01,tectonic_volcanic_fraction_per_myr=.1)
        r=evolve(p)
        self.assertFalse(r['transport'][0]['tectonics_active'])
        self.assertEqual(r['transport'][0]['volcanic_rate_per_myr'],.01)

    def test_invalid_fractions_and_event_age(self):
        with self.assertRaises(ValueError):Crust(felsic=.9).validate()
        with self.assertRaises(ValueError):Event(11.,'comet').validate(10.)
        with self.assertRaises(ValueError):Planet(epochs=0).validate()


if __name__=='__main__':unittest.main()
