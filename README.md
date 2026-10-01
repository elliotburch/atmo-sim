# Terminal Burn atmosphere simulation

## New crust-up revision

The independent [revision 2 model](revision2/README.md) adds conserved elemental equilibrium, finite reactive crust, carbonate and sulfide reservoirs, ocean-gated tectonics, and scheduled geological events. See its [validation report](revision2/VALIDATION.md) and [example input](revision2/example.json). The launcher described below remains the original prototype.

Run the existing `atmo-sim.py` library through `run_planets.py`. Python 3.9+ and
Matplotlib are required for the per-planet plots. Install the dependency with
`python3 -m pip install -r requirements.txt`. The simulator and its physics are
unchanged.

```sh
python3 run_planets.py --help
python3 run_planets.py --list
python3 run_planets.py --planet "Earth" --output results/earth
python3 run_planets.py --all --output results/all
python3 run_planets.py --all --epochs 8 --output results/smoke
python3 run_planets.py --planets my_planets.json --planet "My Planet"
python3 -m unittest -v
```

The default is 128 evolution epochs across the library's fixed 4.5 Gyr interval.
Fewer epochs are useful for smoke checks, but change numerical results. Named
selection is exact except for letter case. Input and library defaults resolve
beside the launcher, so it can be invoked from another directory. The default
output directory is `results` in the current directory. Existing matching CSVs
are overwritten. All simulations finish before CSV writing begins; simulation
failures return a nonzero exit status. Output files from an earlier run are not
deleted on failure.

## Example inputs

`example_planets.json` contains 14 scenarios reconstructed from the user's
description, **not the recovered inputs of the earlier chat's 14-world grid**:
Earth, Mars, Venus, Moon, Titan, Triton, Super Titan (one Earth mass), Cold Water
World (30% water), Migrated Hot Water World (30% water, cold formation and hot
current orbit), Partial Melt World, Molten World, two tidally locked Earth-like
worlds, and Desert World. Seeds and otherwise unspecified quantities are explicit
choices for this set. The names describe intended scenarios, not calibrated
reproductions of the Solar System or guarantees of the final climate. In
particular, inspect the water reservoirs and regional temperatures to assess
supercritical water and desert equator/pole behavior.

The JSON object has `schema_version: 1` and a `planets` array. Each array entry
maps directly to `Planet` keyword arguments, with a nested `crust` object mapping
to `Crust`. Omitted optional fields use the library's dataclass defaults; omitted
or null bulk inventories use its seeded formation draws. Unknown dataclass fields,
duplicate names, invalid numeric types, and basic out-of-range values are rejected.
Top-level description/units metadata does not enter the simulator.

Minimal custom input:

```json
{
  "schema_version": 1,
  "planets": [{
    "name": "My Planet",
    "mass_e": 1.0,
    "radius_e": 1.0,
    "flux": 1361.0,
    "formation_flux": 1361.0,
    "seed": 1,
    "crust": {"oxidation": 0.2}
  }]
}
```

Mass/radius are in Earth units, fluxes in W/m², and bulk inventory values are
fractions of planet mass (0.30 means 30%). Water inventory is H₂O mass; carbon,
nitrogen, and sulfur are elemental mass. The library clamps inventories to its
own supported limits: H₂O 0.50, C 0.10, N 0.03, S 0.10. Crust fields retain the
units in their names. Redox/oxidation range from −1 to +1.

## CSV outputs

| File | Contents |
| --- | --- |
| `summary.csv` | One final row per planet, using `result_row`, plus formation erosion diagnostics |
| `history.csv` | The library's sparse epoch samples: age, temperature, pressure, crust phase and melt fraction |
| `columns.csv` | Final three-region climate diagnostics from `column_rows` |
| `layers.csv` | Final ordered material layers, pressures, temperatures, thickness and coupling |
| `reservoirs.csv` | Final atmospheric, condensed, crustal and escaped inventories, kg/m² |
| `vertical_profiles.csv` | Final regional altitude/pressure/temperature and condensable mole-fraction profiles |

Each run also writes `layer_cake.png` and `atmosphere_composition.png` under a
subdirectory named for each planet (for example, `results/Earth/`). The layer
cake shows the crust and condensed shells beneath each column's resolved
surface-to-tropopause profile; the hatched upper-atmosphere height is schematic.
Surface temperatures are labeled for each column. Atmospheric pie slices are
reported in ppm, with total surface pressure in the title.

History ages are the library's epoch midpoint labels; history is neither a
complete timestep log nor the final current-star solve. `summary.csv` describes
that final solve. Its species columns are lower-atmosphere mole fractions;
pressure includes lower and upper atmosphere. `initial_target_pressure_bar` is
the formation erosion target, not a post-equilibration measured pressure.
Reservoir `mixed` denotes the library's aggregate organic deposits.

The launcher calls `evolve(planet, epochs=...)` directly and consumes its
`(state, bulk, redox, history)` return tuple. It does not duplicate the evolution
loop, reorder initialization, or add chemistry. The library remains a
gameplay-oriented reduced-order model, not a precision climate model.
