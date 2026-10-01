# Initial verification

Test environment: Python 3.11, Cantera 3.2.0, NumPy 2.4.6, SciPy 1.17.1.

23 foundation tests pass. These test element conservation across temperatures
and reducing/oxidizing inventories, water coexistence and mixed-gas partial
pressure, finite carbonate and iron capacities, sulfide formation, no invented
water, energy/pressure closure, contiguous layer interfaces, map fractions,
ocean-gated tectonics, event ordering, finite fresh-rock replacement, retained
oxygen during H escape, the escape energy ceiling, fixed epoch age, and
resolution-independent random event schedules. Closed-system epochs do not
continue changing chemistry merely because the simulator was called again.

The 8-epoch demonstration completes, including all five explicit events. Its
outputs are examples of the model, not targets used to fit it.

## Resolution experiment

Run from the repository root:

```sh
.venv-v2/bin/python -m revision2.validate_resolution
```

This holds age at 4 Myr, keeps the same delivery and impact times, and halves
the maximum fractional inventory change allowed per internal step:

| Fractional limit | Final T (K) | Final pressure (Pa) | Escaped H (mol/m²) |
| --- | ---: | ---: | ---: |
| 0.05 | 289.5634 | 60,859.15 | 226,663.32 |
| 0.025 | 289.5573 | 60,952.22 | 228,047.79 |
| 0.0125 | 289.5541 | 61,001.23 | 228,780.62 |

The last refinement changes temperature by 0.0032 K, pressure by 0.080%, and
escaped hydrogen by 0.32%. Absolute per-element ledger errors stay below
2e-9 mol/m² in this experiment. This is evidence of improving numerical
resolution for this fixture, not a universal error bound or physical calibration.
Climate bifurcations and other inputs still require their own resolution checks.

## What these tests do not establish

They do not validate the grey opacity coefficients, escape boundary parameters,
impact magnitude priors, accessibility of minerals beneath ice, or the actual
occurrence of tectonics. The README lists the supported phase range and the
approximations. Thick water shells and magma worlds are explicitly outside this
revision; they are not assigned fabricated phase stacks.
