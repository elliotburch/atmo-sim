# Crust-up revision 2

A new implementation, independent of `atmo-sim.py`. It starts **after planetary
differentiation**, with a fixed inert core/scaffold, finite mantle volatiles,
fresh reactive rock, a surface veneer, and a specified primordial elemental
inventory. It does not simulate differentiation or choose an atmosphere by planet
class. The primordial inventory is allowed to condense and react immediately.

## Run

From the `atmo-sim` directory:

```sh
python3 -m venv .venv-v2
.venv-v2/bin/python -m pip install -r revision2/requirements.txt
.venv-v2/bin/python -m revision2 --input revision2/example.json --output results/revision2
.venv-v2/bin/python -m unittest discover -s revision2/tests -v
```

The old launcher's Matplotlib dependency is unnecessary. Outputs are
`evolution.json` (complete snapshots, reservoirs, events, and transport audit),
`history.csv` (one row per epoch), and `layers.csv` (final crust-up stack).
All computation succeeds before output writing starts. Existing named output
files are overwritten. No files are written on a simulation failure.

## Units and time

- Elemental inventories: **moles of atoms per square metre**. H=2, O=1 supplies
  one mole of water's elements, without requiring that they remain water.
- Species inventories in results: moles of molecules/formula units per m².
- Pressure: Pa internally; the history table also provides bar.
- Surface material proportions: fractions summing to one, not percentages.
- Rates: fractions per Myr. An interval removes `1-exp(-rate * duration)` of the
  available source reservoir, rather than drawing from an infinite supply.
- Age = `epochs * epoch_duration_myr`. Both the initial equilibrium and every
  completed epoch are exported. There is no independent hidden 4.5 Gyr age.
- `max_transport_step_myr` controls numerical resolution within epochs. It does
  not change the system age. Event times split transport steps exactly.

For example, 450 epochs of 10 Myr represent 4.5 Gyr. To compare integration
resolutions, keep that age fixed and reduce `max_transport_step_myr`. Chemistry
has no reaction-rate clock: it is re-equilibrated after reservoir changes.

## Conserved material

H, C, N, O, Ar, He, Fe, Ca and S are tracked independently in:

1. accessible atmosphere + surface condensates + reactive minerals;
2. mantle volatiles;
3. finite fresh-crust stock;
4. buried solids;
5. escaped elements.

For **each element**, after every transport operation and event:

```
accessible + mantle + fresh crust + buried + escaped = initial + delivered
```

The fixed structural core/crust and unreactive silicate mantle do not exchange
chemistry. Initial modeled material must fit inside the non-core planetary mass.
The inert portion of a surface veneer follows its own finite replacement count;
resurfacing transfers equal inert veneer masses from fresh stock to burial.
Planet radius/gravity remain fixed: late deliveries and losses are assumed small
relative to the bulk planet.

There are no arbitrary gas redox fractions. Cantera minimizes Gibbs energy at
specified T/P, using a restricted gas/mineral phase set. A separate pressure
root requires atmospheric weight to match surface pressure. Partial pressures
use mole fractions, not gas weight fractions. Pure-water vapor/condensate
coexistence uses a lever rule when the fixed-pressure phase amount is discontinuous.
Every equilibrium solve checks its elemental residuals.

## Minerals and map-facing material fractions

`crust.felsic`, `mafic`, and `ultramafic` sum to one **within the silicate
component**. They are retained independently of deposits; they do not silently
imply new chemical capacities. Set capacities explicitly:

- `feo_mol_m2`: fresh ferrous mineral capacity, represented as FeO;
- `cao_mol_m2`: fresh alkalinity, represented as CaO-equivalent;
- `fes_mol_m2`: fresh sulfide inventory;
- `active_silicate_kg_m2`: mass of the exchangeable veneer, not the whole crust;
- `fresh_veneer_replacements`: finite stock of additional identical veneers.

The equilibrium mineral set includes metallic Fe, FeO, magnetite, hematite,
FeS, pyrite, CaO, calcite, calcium sulfate, graphite and condensed sulfur.
Carbonate deposition **consumes finite Ca and C**; sulfides **consume finite Fe
and S**. Both can react away under a new equilibrium. Sulfate is included so
oxidized sulfur has a mineral destination rather than being forced into sulfides.

Each snapshot exports `surface_fractions` for felsic, mafic, ultramafic,
carbonate, sulfide and sulfate. These normalize the active substrate's mineral
masses; the silicate proportions retain their original ratio. They are useful
map-generator inputs, **not claimed area coverages**. Separate layer masses
provide graphite, sulfur, water and ice cover. Burial removes deposits from the
active substrate and records their elements permanently in `buried`.

`oxidation_capacity_mol_o_m2` is the additional oxygen needed to take reactive
iron to Fe(III) and sulfide sulfur to S(VI). It is a diagnostic of actual minerals,
not an inexhaustible scalar that is reset. Equilibrium may oxidize or reduce this
finite assemblage. Resurfacing buries its current oxygen and sulfur and exposes
fresh material by an explicit, conserved transfer. Actual fresh replacement can
be less than requested when stock is exhausted.

## Equilibrium stack and climate

At every temperature candidate:

1. Speciate the accessible elements among gas and allowed pure condensed phases.
2. Solve `surface pressure = gravity * gas mass column`.
3. Compute optical depth from explicit mass absorption coefficients.
4. Evaluate `OLR = sigma * T^4 / (1 + 3*tau/4)`.

Select a resolved stable root of `OLR = absorbed stellar + geothermal flux`,
continuing the nearby stable temperature branch. The root search follows the prior stable branch locally and falls back to a
10 K scan; very narrow branches
could still require a finer search. Failed solves and unresolved phase-boundary
roots raise errors. Returned energy residuals must be below 0.01 W/m².

Only then build the layers, **without moving mass**: inert crust, reactive
veneer, mineral deposits, condensed water/ice, atmosphere. Layer interfaces have
matching pressure and temperature. Atmospheric bins integrate the actual gas
column and include the tenuous tail in the final bin; the reported height ends
at 14 scale heights. The gas is well mixed and the geometry is isothermal.
These are deliberately reduced closures, not a resolved radiative-convective
column. There is no separate supercritical "liquid water" reservoir.

## Long-term processes and events

Background volcanism transfers a fraction of remaining mantle volatile atoms to
the accessible pool. Background resurfacing exposes fresh veneers. A
`volcanic_rate` event sets the multiplier of both background rates until another
such event changes it. Tectonic contributions are added separately.

Tectonics requires `tectonics_enabled` **and** surface liquid water exceeding
`minimum_contact_water_kg_m2`. In this supported shallow-shell model, deposits
are treated as permeable crust and liquid water contacts the crust. Dry and
ice-covered states disable tectonics; there is no subsurface ocean model yet.
Deep shells requiring high-pressure ice are rejected rather than falsely enabling
tectonics. Ocean contact is a requested model rule, not a claim that it is a
sufficient physical criterion for plate tectonics.

Explicit event kinds:

| Kind | Effect |
| --- | --- |
| `comet` | Deliver specified volatile elemental abundances |
| `impact` | Strip a fraction of existing gas; deliver material; release mantle volatiles; resurface |
| `volcanic_rate` | Set a background volcanic/resurfacing multiplier |
| `resurfacing` | Bury and replace a fraction of the reactive veneer |

Impact stripping applies **before** impact-driven outgassing. Condensed oceans
are not directly stripped. This first revision retains mass-transfer effects,
not transient shock heating, ejecta, impact chemistry kinetics or cooling.
Each event is logged with actual transferred elements and before/after pressure
and temperature. Events at an epoch endpoint occur before its snapshot;
equal-time events follow input order. Time-zero events are supported.

Optional `random_events` accepts `comets_per_myr` and `impacts_per_myr`. Seeded
Poisson arrival times are sampled in physical time independently of numerical
resolution. Magnitudes are documented gameplay priors in `sample_events`, not
an inferred bombardment history. Explicit event schedules are preferable for
experiments.

Escape combines a Jeans-flux estimate at a mean-free-path exobase with a
cold-trap-limited hydrogen-loss approximation. Exosphere temperature and XUV flux
are explicit boundary inputs. Losses share an energy ceiling and cannot exceed
available accessible elements. Equilibrium phase exchange can replenish vapor
continuously; loss is not capped at one instantaneous vapor column. Hydrogen removed from water leaves
oxygen behind. Adaptive substeps bound fractional net inventory changes; exponential loss
integration includes continuous volcanic supply and avoids a stiff trace-gas
timestep limit. `max_escape_fraction_per_step` defaults to 0.05. Reducing it, as
well as `max_transport_step_myr`, tests numerical convergence; it is not a calibrated stellar-wind or hydrodynamic escape model.

## Deliberate limits / next physics work

This revision establishes the contracts and conservation, not a full replacement
for every speculative feature in the old prototype:

- Supported solve range: 200–1200 K; no Titan/Triton cryochemistry or magma EOS.
- Minerals and sulfur gases below ~300 K extrapolate their thermochemical fits;
  cold snapshots flag this limitation.
- Gas pressure above 100 bar and condensed basal pressure above 0.2 GPa are
  rejected. Between 10 and 100 bar, ideal-gas results are flagged as qualitative.
- No deep water phase diagram, high-pressure ice, supercritical fluid structure,
  salinity, dissolved ions, carbonate pH model, clouds or regional climate.
- CaO is an alkalinity proxy; felsic/mafic/ultramafic fractions are preserved for
  maps, not yet connected to a complete silicate reaction network.
- The exposed reactive veneer is chemically accessible even beneath surface ice;
  ice isolation and subsurface exchange require a subsequent interface model.
- The isothermal stack does not resolve conductive geothermal profiles, a moist
  adiabat, vertical chemistry or latent-heat transport. Geothermal flux enters
  the global energy boundary only.
- Bolometric and XUV fluxes are fixed; stellar evolution is a future boundary
  function, not implicitly tied to a Sun-like age law.

The next useful extension is a pressure/temperature-resolved water shell with
explicit interface accessibility, then a vertical radiative-convective closure.
Keep the elemental ledger and event/transport contracts unchanged when adding
those solvers.

See [VALIDATION.md](VALIDATION.md) for the tests and a reproducible resolution experiment.

## Data and references

- [Cantera equilibrium API](https://www.cantera.org/stable/python/thermo.html):
  gas pre-equilibration and multiphase VCS Gibbs minimization. Thermochemical data
  are supplied by Cantera's `gri30.yaml`, `nasa_gas.yaml`,
  `nasa_condensed.yaml`, and `graphite.yaml`; their file headers retain provenance.
- [IAPWS water saturation formulation](https://www.iapws.org/relguide/Supp-sat.html).
- [Murphy & Koop ice saturation](https://doi.org/10.1256/qj.04.94).

No thermochemical table is fetched at simulation runtime.
