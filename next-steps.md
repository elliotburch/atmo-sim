# Terminal Burn Planetary Evolution Simulator
## Current State, Open Problems, and Roadmap

This report summarizes where the planetary evolution simulator currently stands, what the most important next steps are, and several related design conversations that should inform future development.

---

## 1. Current design goal

The simulator is intended to generate physically plausible starting planets for Terminal Burn rather than provide a perfect first-principles planetary science model.

The preferred approach is to generate planets from continuous physical and chemical variables rather than explicit planet classes. A planet should become Earth-like, Titan-like, carbon-rich, snowball, desert, magma-ocean, etc. because of its mass, inventories, redox state, climate, escape history, and surface phases.

The major emergent inputs are becoming:

- mass and radius
- stellar flux and formation flux
- bulk H2O, C, N, and S inventories
- mineral redox state
- primordial atmospheric retention
- surface pressure and temperature
- ocean/ice/crust geometry
- ocean-silicate coupling
- photochemistry and escape
- horizontal heat transport

---

## 2. Current initialization pipeline

The current intended initialization order is:

```text
bulk elemental / volatile inventories
        ↓
initial redox-dependent speciation
        ↓
deep vs accessible reservoirs
        ↓
primordial climate / phase equilibrium
        ↓
layer cake
        ↓
actual atmosphere remaining above the surface
        ↓
primordial impact erosion
        ↓
re-equilibrate at reduced pressure
        ↓
begin long-term evolution
```

---

## 3. Volatile initialization

Planets start with independent bulk inventories:

```python
water_fraction
carbon_fraction
nitrogen_fraction
sulfur_fraction
```

These can either be explicitly specified for designed regression worlds or randomly generated.

This fixes a major problem with water worlds. Previously, something like a 4% volatile fraction meant that 4% of the entire planet was divided among H2O, N2, CO2, CO, CH4, and other gases. That produced absurd nitrogen and carbon inventories.

A water world can now instead look like:

```text
4% H2O
0.1% C
0.01% N
0.2% S
```

which is much more useful.

---

## 4. Bulk vs accessible material

The model now conceptually distinguishes:

```text
bulk planetary inventory
        ↓
deep / inaccessible reservoir
        +
surface-accessible reservoir
```

Only the accessible fraction participates directly in the current atmosphere-ocean-crust evolution.

This distinction is crucial because a planet can contain enormous total volatile or carbon inventories without having those materials immediately available to the atmosphere.

Long term, the deep reservoir can remain a reduced-order bookkeeping system rather than a full mantle simulation.

---

## 5. Primordial atmosphere erosion

A recent addition is formation-era atmospheric erosion.

The idea is that giant impacts and violent late accretion remove a large fraction of the atmosphere after the first surface has already formed.

The retained atmosphere is assigned a stochastic target pressure, broadly around:

```text
0.1–100 bar
```

with the possibility of much thinner atmospheres for small or vulnerable bodies.

This is only an initialization prior. It does **not** cap later evolution.

A planet may start at 2 bar and later evolve to 100 bar, or start at 20 bar and later lose nearly everything.

This provides a useful separation between:

```text
candidate primordial atmosphere
retained post-impact atmosphere
final game-era atmosphere
```

---

## 6. Layer cake / surface structure

The simulator already supports a useful reduced-order vertical structure including:

- atmosphere
- organic deposits
- volatile frost
- liquid methane
- water ice
- liquid water
- high-pressure ice
- solid silicate crust
- partially molten silicate
- magma
- carbonaceous crust

This allows useful states such as:

```text
ocean directly on rock
ocean over high-pressure ice
ice shell over subsurface ocean
global ice shell
surface magma
```

A particularly important quantity is whether an ocean can chemically contact silicate rock.

That controls weathering, carbonate formation, redox buffering, and eventually ocean chemistry.

---

## 7. Two-atmosphere model

The atmosphere currently uses two conserved boxes:

```text
LOWER ATMOSPHERE
    ↓ mixing / cold trap
UPPER ATMOSPHERE
    ↓ escape
SPACE
```

The lower atmosphere handles greenhouse physics, condensation, clouds, and surface interaction.

The upper atmosphere handles photochemistry, methane destruction, hydrogen escape, XUV interaction, and solar-wind loss.

A diagnostic tropopause / cold trap limits how much water reaches the upper atmosphere.

This was a major improvement because warm lower atmospheres no longer automatically expose the entire ocean to photolysis.

---

## 8. Three-column climate

The climate currently uses three geographic columns.

For normal rotating planets:

```text
equator
midlatitude
pole
```

For tidally locked planets:

```text
substellar
terminator
antistellar
```

They share a global atmosphere but have different local:

- temperatures
- albedos
- cloud fractions
- ice fractions
- ocean fractions
- condensed-water inventories

This allows states such as:

- hot equator with polar seas
- equatorial ocean with polar ice
- snowball
- substellar ocean with nightside ice
- terminator habitable belt
- nightside atmospheric collapse

---

## 9. Heat transport

Horizontal heat transport was increased substantially.

The older transport law isolated latitude bands too strongly, even at Earth-like pressure.

The new pressure-dependent law rises much faster with atmospheric pressure and approaches strong redistribution for dense atmospheres.

This especially improved dense tidally locked worlds, reducing implausibly enormous day-night temperature contrasts.

The model may still need tuning, but the behavior is now in the right qualitative regime.

---

## 10. Vertical atmosphere and cloud model

The simulator now has a pressure-coordinate vertical profile for each climate column.

It currently treats:

```text
H2O
CH4
CO2
```

as condensable species.

The vertical solver estimates:

- dry and moist adiabatic profiles
- saturation
- cloud base
- cloud top
- condensate
- precipitable vapor
- cold-trap abundance

This was added because the earlier climate model treated too much water as though it were distributed through the whole atmosphere.

A remaining problem is that the climate solver does not yet iterate fully against the vertically diagnosed water column.

---

## 11. Standard regression worlds

The current standard grid contains:

- Ashen Vale — dry abiotic Earth analog
- Pelagos — ocean Earth analog
- Cythera — Venus analog
- Ares — Mars analog
- Selene — Moon analog
- Amber Veil — Titan analog
- Pale Lantern — Triton analog
- Nereid Deep — cold deep-water world
- Wandersea — migrated water world
- Obsidian Garden — carbon-rich Earth
- Sable Breath — reducing terrestrial
- Atlas — 5 Earth-mass super-Earth
- Emberglass — partial-melt test
- Caldera — magma-ocean test

These are best thought of as regression and behavior tests rather than exact Solar System reproductions.

---

## 12. Water-world progress

The water worlds were one of the largest earlier failures.

Before the inventory and erosion changes, they produced thousands of bars of atmosphere, often dominated by H2, N2, CO, and other non-water volatiles.

After separating H2O/C/N/S inventories and moving primordial erosion after the layer cake, their behavior became much more reasonable.

Nereid Deep approximately became:

```text
post-layer-cake atmosphere: ~412 bar
after primordial erosion:   ~2 bar
final atmosphere:           <1 bar
```

Wandersea became roughly:

```text
post-layer-cake atmosphere: ~930 bar
after primordial erosion:   ~2 bar
final atmosphere:           a few bar
```

The large pre-erosion numbers still indicate that the accessible fraction and initial equilibrium can be improved, but the final atmospheres are no longer pathological.

---

## 13. Biggest current geochemical problem: carbon

Carbon-rich worlds tend to rebuild very thick atmospheres after primordial erosion.

This now appears to be a genuine reservoir-model problem rather than an initialization bug.

The missing competition is approximately:

```text
bulk carbon
    ↓
mantle carbon
graphite / diamond
carbides
carbonate minerals
dissolved ocean carbon
organics
atmospheric CO2 / CO / CH4
```

At present, too much accessible carbon eventually finds its way into the atmosphere.

This is probably the highest-value next geochemical system to improve.

---

# Recommended Next Work

## Priority 1 — Carbon reservoir model

Add a reduced-order carbon cycle with at least:

```text
atmospheric carbon
dissolved ocean carbon
carbonate sediment / crust
refractory carbon
reduced carbon / organics
```

The main atmospheric pathway should become:

```text
CO2 atmosphere
    ↕
dissolved CO2
    ↕
carbonic acid / bicarbonate / carbonate
    ↕
carbonate rock
```

Exchange rates should depend on:

- water availability
- temperature
- ocean-silicate coupling
- exposed fresh rock
- redox state
- crust composition

This would help keep carbon-rich planets from automatically turning their entire carbon inventory into atmospheric CO2.

---

## Priority 2 — Condensed CO2 and CH4 must become real reservoirs

The current vertical model can diagnose CO2 or CH4 condensation, but condensate is not always transferred into persistent surface reservoirs.

That needs to become real mass movement:

```text
atmosphere
    ↓ condensation
cloud
    ↓ deposition
surface frost / liquid / clathrate
```

This is necessary for:

- atmospheric collapse
- permanent nightside cold traps
- seasonal volatile caps
- methane lakes / frost
- CO2 ice worlds

---

## Priority 3 — Elemental conservation

The chemistry should eventually operate on elemental budgets.

Track at least:

```text
H
C
N
O
S
Ar
```

and potentially Fe or an abstract crustal oxidation capacity.

Then enforce:

```text
current elemental inventory
+
escaped inventory
=
initial inventory
+
explicit sources
```

This would eliminate the current places where oxygen or hydrogen is implicitly created or consumed.

---

## Priority 4 — Fix oxidation chemistry

Specific chemistry that should eventually become stoichiometric:

### Water photolysis

```text
H2O
 ↓
H escapes
O remains
```

The retained oxygen should either:

- form O2
- oxidize CO
- oxidize CH4 products
- oxidize the crust

### CO oxidation

CO should consume a real oxidant rather than simply disappear.

### Methane photochemistry

Methane destruction should conserve C and H explicitly and route products into:

- H / H2 escape
- CO
- CO2
- hydrocarbons
- organic haze / deposits

---

## Priority 5 — Replace the current climate temperature iteration with an OLR solver

The existing climate formulation is still a reduced-order approximation.

A stronger architecture would solve:

```text
OLR(T, pressure, composition)
=
absorbed stellar flux
```

over candidate surface temperatures.

This would allow physically meaningful branches:

```text
snowball
cold stable
temperate
warm ocean
moist greenhouse
runaway greenhouse
```

and make runaway detection much more robust.

---

## Priority 6 — Fully couple vertical humidity to climate

Right now, the vertical atmosphere solver diagnoses a more realistic water-vapor column after the climate calculation.

Eventually the model should iterate:

```text
temperature
    ↓
vertical water profile
    ↓
greenhouse opacity
    ↓
new temperature
```

until convergence.

---

## Priority 7 — Surface water coverage should depend on terrain

A major conceptual lesson from the hot-desert / polar-sea test was:

Add a basin depth parameter to planets -- this would reflect how deep a body of water can become and remain localized. Once a bin exceeds its basin depth, the additional fluid is added to the nearest bin. 

```text
water mass ≠ surface coverage
```

A few thousand kg/m2 of water might form:

- scattered lakes
- deep isolated basins
- groundwater
- inland seas
- shallow global ocean

depending on topography.

Once the planetary terrain generator exists, local basin capacity should determine ocean coverage.

---

## User priority -- Better model transient events

Epochs should contain rare events including - deposition of volatiles by comet impacts. Massive impacts that can remove atmosphere. Volcanic activity changing intensity (could help with the carbon problem). 

# Side Conversations and Future Chemistry

## Redox vs pH

A useful distinction emerged:

### Redox

Redox is mainly an **initial mineralogical / geochemical state**.

It describes the preferred oxidation state of the crust, mantle, and outgassed species.

Examples:

```text
reducing:
H2
CH4
CO
H2S
NH3

oxidizing:
CO2
SO2
NOx
sulfates
nitrates
```

### pH

pH should be an **emergent property of a liquid reservoir**.

It depends on:

```text
dissolved atmospheric gases
mineral buffering
water-rock interaction
salts / ions
weathering
temperature
```

So a planet should generally not start with an arbitrary pH number. It should start with mineralogy and chemistry from which ocean pH eventually emerges.

---

## Nitrogen chemistry

Nitrogen has a useful redox sequence:

```text
NH3
 ↓
N2
 ↓
NO
 ↓
NO2
 ↓
HNO3 / nitrate
```

In that sense, NH3 and HNO3 are close to opposite ends of the nitrogen redox spectrum.

Potential environments include:

- ammonia-rich reducing atmospheres
- ordinary N2 atmospheres
- oxidizing NOx chemistry
- nitric-acid precipitation
- nitrate-rich surface deposits

---

## Sulfur chemistry

Sulfur has an especially useful planetary oxidation ladder:

```text
H2S
 ↓
S
 ↓
SO2
 ↓
SO3
 ↓
H2SO4 / sulfate
```

This could naturally generate:

- H2S volcanic atmospheres
- elemental sulfur deposits
- SO2-rich planets
- sulfuric acid clouds
- sulfate-rich crusts

Sulfur is probably the best trace chemistry system to add after carbon is stable.

---

## Sulfuric acid

Sulfuric acid is important because it can dominate clouds and albedo while remaining only a trace fraction of the total atmosphere.

A Venus-like pathway is:

```text
SO2
+
H2O
+
UV / oxidants
    ↓
H2SO4
    ↓
acid aerosol / cloud deck
```

This should probably be represented as photochemical aerosol formation rather than ordinary equilibrium gas condensation.

---

## Nitric acid

Possible pathway:

```text
N2
 ↓ lightning / energetic particles / hot chemistry
NO
 ↓
NO2
 ↓
HNO3
 ↓
rainout / nitrate deposits
```

This could create strongly oxidized acidic worlds.

---

## Carbonic acid

Carbon chemistry is likely the most important aqueous acid system.

```text
CO2(g)
 ↕
CO2(aq)
 ↕
H2CO3
 ↕
HCO3-
 ↕
CO3--
```

This directly connects atmosphere, climate, ocean pH, weathering, and carbonate deposition.

This should probably be implemented before detailed nitric- or sulfuric-acid ocean chemistry.

---

## Ammonia

NH3 is not currently simulated.

A future ammonia pathway could include:

```text
NH3 / NH4-bearing interior
      ↓
cryovolcanism
      ↓
atmospheric NH3
      ↓
clouds / condensation
      +
photolysis
      ↓
NH4+ in water
N2 + H products
```

Ammonia is especially relevant to cold, reducing outer-system planets.

It also pushes water chemistry alkaline:

```text
NH3 + H2O ↔ NH4+ + OH-
```

---

## Silicon and mineral buffering

Silicon is not expected to behave like carbon in normal terrestrial atmospheres.

Its main importance is that silicate minerals determine:

- Ca and Mg availability
- Fe2+/Fe3+ buffering
- Na and K availability
- weathering capacity
- carbonate precipitation potential
- redox capacity

Rather than track detailed mineral species immediately, a reduced-order model could track a few crustal buffering capacities.

At very high temperatures, however, silicon becomes atmospheric and interesting.

Magma planets could have rock-vapor atmospheres containing things such as:

- SiO
- Na
- K
- O
- O2
- metal vapor

This opens the door to a genuine rock-vapor weather cycle.

---

# Different Liquid / Condensable Cycles

One of the more interesting side directions was to generalize the idea of a hydrological cycle.

Instead of assuming water is always the active planetary liquid, the simulator could eventually support generic condensable cycles:

```text
surface reservoir
    ↓
evaporation / sublimation
    ↓
atmospheric vapor
    ↓
cloud
    ↓
precipitation
    ↓
surface reservoir
```

Possible active cycles include:

- H2O
- CH4 / C2H6
- NH3-H2O
- CO2
- N2
- SO2
- H2SO4/H2O droplets
- silicate / metal vapor

This would allow the dominant “weather” cycle to emerge from temperature, pressure, and composition.

Examples:

### Temperate terrestrial
Water performs the active weather cycle.

### Titan-like
Water behaves almost like rock while methane / ethane form rain, clouds, rivers, and lakes.

### Extremely cold world
N2 or CH4 may migrate between atmosphere and surface.

### Venus-like
H2SO4/H2O droplets form photochemical clouds and possibly precipitation that evaporates before reaching the ground. (or forms bodies of surface water if conditions are right!)

### Magma world
Rock evaporates on the hot side, travels through a silicate-vapor atmosphere, and condenses on cooler regions.

That would be a very interesting long-term extension because it produces exotic environments without requiring bespoke planet classes.

---

# Recommended Immediate Work Session

If development resumes from the current state, the best next focused session would probably be:

1. Add the simplified carbon reservoir / carbonate system.
2. Make CO2 and CH4 condensation move real mass into surface reservoirs.
3. Add elemental conservation diagnostics for C, H, O, and N.
4. Rerun the full standard grid and identify which worlds move substantially.
5. Only after that, revisit the OLR climate solver.

Once those systems are stable, sulfur / ammonia / acid chemistry becomes much safer to add because the underlying mass bookkeeping and ocean chemistry will have a stronger foundation.

---

# Long-Term Design Principle

The strongest design direction so far is:

> Generate reservoirs, chemistry, phase states, and energy balance; let planet types emerge.

The simulator should ideally never need to say:

```text
PlanetType = Titan
PlanetType = AcidWorld
PlanetType = CarbonWorld
```

Instead those worlds should emerge from combinations of:

```text
mass
gravity
stellar flux
H2O / C / N / S inventory
mineral redox
primordial atmospheric retention
surface pressure
condensable phase behavior
ocean-rock coupling
escape
photochemistry
aqueous chemistry
mineral buffering
```

That approach fits Terminal Burn particularly well because it allows unusual planetary environments to arise naturally from the same system rather than from hand-authored categories.