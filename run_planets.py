#!/usr/bin/env python3
"""Load Planet/Crust inputs and run the unchanged planetary evolution library."""

import argparse
import csv
from dataclasses import fields
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent


def load_library(path):
    spec = importlib.util.spec_from_file_location("planet_simulator", path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Cannot import simulator: {path}")
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolves annotations through the module registry.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def construct(cls, values, label):
    if not isinstance(values, dict):
        raise ValueError(f"{label} must be a JSON object")
    unknown = values.keys() - {f.name for f in fields(cls)}
    if unknown:
        raise ValueError(f"{label}: unknown fields: {', '.join(sorted(unknown))}")
    try:
        return cls(**values)
    except TypeError as exc:
        raise ValueError(f"{label}: {exc}") from exc


def validate_planet(p):
    if not isinstance(p.name, str) or not p.name.strip():
        raise ValueError("Every planet needs a nonempty name")
    for owner in (p, p.crust):
        for f in fields(owner):
            value = getattr(owner, f.name)
            if f.name in ("name", "crust"):
                continue
            if f.name == "tidally_locked":
                if type(value) is not bool:
                    raise ValueError(f"{p.name}: tidally_locked must be a boolean")
                continue
            if f.name == "seed":
                if type(value) is not int:
                    raise ValueError(f"{p.name}: seed must be an integer")
                continue
            if value is None and f.name in (
                "water_fraction", "carbon_fraction", "nitrogen_fraction",
                "sulfur_fraction", "formation_redox",
            ):
                continue
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError(f"{p.name}: {f.name} must be a finite number")
            if f.name in ("formation_redox", "oxidation"):
                valid = -1 <= value <= 1
            elif f.name in ("mass_e", "radius_e", "flux", "formation_flux"):
                valid = value > 0
            elif f.name.endswith("fraction") or f.name in (
                "carbon_richness", "albedo_rock", "mafic", "felsic", "ultramafic", "porosity",
            ):
                valid = 0 <= value <= 1
            else:
                valid = value >= 0
            if not valid:
                raise ValueError(f"{p.name}: {f.name} is outside its allowed range")


def load_planets(path, sim):
    with path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ValueError("Planet file must be an object with schema_version: 1")
    definitions = document.get("planets")
    if not isinstance(definitions, list) or not definitions:
        raise ValueError("Planet file must contain a nonempty planets array")
    planets, seen = [], set()
    for definition in definitions:
        if not isinstance(definition, dict):
            raise ValueError("Each planets entry must be an object")
        values = dict(definition)
        values["crust"] = construct(sim.Crust, values.get("crust", {}), "Crust")
        planet = construct(sim.Planet, values, "Planet")
        validate_planet(planet)
        key = planet.name.casefold()
        if key in seen:
            raise ValueError(f"Duplicate planet name: {planet.name}")
        seen.add(key)
        planets.append(planet)
    return planets


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def planet_directory_name(name):
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    return cleaned or "planet"


def write_planet_plots(planet, state, sim, pressure_bar, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    output.mkdir(parents=True, exist_ok=True)
    layers = list(reversed(state.layers[1:]))
    material_colors = {
        layer.material.value: plt.get_cmap("tab10")(index % 10)
        for index, layer in enumerate(layers)
    }
    columns = state.climate_columns
    figure, axis = plt.subplots(figsize=(9, 7))
    x_positions = list(range(len(columns)))
    maximum_height = 0.0

    for x_position, column in zip(x_positions, columns):
        height = 0.0
        for layer in layers:
            axis.bar(
                x_position, layer.thickness_m / 1000, bottom=height / 1000,
                width=0.72, color=material_colors[layer.material.value],
                edgecolor="white", linewidth=0.35,
            )
            height += layer.thickness_m

        surface_height = height
        profile = column.vertical_profile
        if profile:
            altitudes = [level["z_m"] for level in profile]
            final_step = altitudes[-1] - altitudes[-2] if len(altitudes) > 1 else 1000.0
            for index, level in enumerate(profile):
                if index + 1 < len(profile):
                    thickness = max(0.0, altitudes[index + 1] - altitudes[index])
                else:
                    thickness = max(0.0, final_step)
                shade = 0.38 + 0.5 * (index / max(1, len(profile) - 1))
                axis.bar(
                    x_position, thickness / 1000, bottom=height / 1000,
                    width=0.72, color=plt.get_cmap("Blues")(shade),
                    edgecolor="white", linewidth=0.2,
                )
                height += thickness

        # The upper atmosphere has a representative bulk height, but no resolved profile.
        representative_atmosphere = next(
            (layer.thickness_m for layer in state.layers
             if layer.material.value == "atmosphere"), 0.0,
        )
        profile_height = height - surface_height
        upper_height = max(0.0, representative_atmosphere - profile_height)
        if upper_height:
            axis.bar(
                x_position, upper_height / 1000, bottom=height / 1000,
                width=0.72, color="#a7cce8", edgecolor="white", linewidth=0.35,
                hatch="//",
            )
            height += upper_height

        maximum_height = max(maximum_height, height)

    handles = [
        Patch(facecolor=color, label=material)
        for material, color in material_colors.items()
    ]
    handles.append(Patch(facecolor="#a7cce8", hatch="//", label="Upper atmosphere (schematic)"))
    legend = None
    if handles:
        legend = axis.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1))
    axis.set_xticks(x_positions, [column.name for column in columns])
    axis.set_xlim(-0.55, max(0.55, len(columns) - 0.45))
    axis.set_ylabel("Height from crust base (km)")
    axis.set_title(f"{planet.name}: crust to atmosphere")
    axis.set_ylim(0, max(1.0, maximum_height / 1000 + 0.6))
    figure.subplots_adjust(left=0.10, right=0.76, bottom=0.10, top=0.90)
    if legend:
        figure.canvas.draw()
        legend_bottom = legend.get_window_extent(figure.canvas.get_renderer())
        legend_bottom = legend_bottom.transformed(figure.transFigure.inverted()).y0
        temperature_lines = []
        for column in columns:
            name = "Midlat" if column.name.casefold() == "midlatitude" else column.name.title()
            temperature_c = column.temperature_k - 273.15
            temperature_lines.append(f"{name} T: {temperature_c:.1f} C")
        figure.text(
            0.78, legend_bottom - 0.035, "\n".join(temperature_lines),
            ha="left", va="top", fontsize=9, linespacing=1.5,
        )
    figure.savefig(output / "layer_cake.png", dpi=160, bbox_inches="tight")
    plt.close(figure)

    moles = {
        species: (state.atm[species] + state.upper_atm[species]) / sim.MM[species]
        for species in sim.SPECIES
    }
    total_moles = sum(moles.values())
    present = [(species, amount) for species, amount in moles.items() if amount > 0]
    figure, axis = plt.subplots(figsize=(8, 6))
    if total_moles:
        wedges, _ = axis.pie(
            [amount for _, amount in present], startangle=90,
            wedgeprops={"edgecolor": "white", "linewidth": 0.8},
        )
        legend_labels = [
            f"{species}: {amount / total_moles * 1_000_000:,.1f} ppm"
            for species, amount in present
        ]
        axis.legend(wedges, legend_labels, title="Atmospheric composition",
                    loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    else:
        axis.text(0.5, 0.5, "No atmospheric gases", ha="center", va="center")
    axis.set_title(f"{planet.name}: atmosphere at {pressure_bar:.3g} bar surface pressure")
    axis.set_aspect("equal")
    figure.tight_layout()
    figure.savefig(output / "atmosphere_composition.png", dpi=160, bbox_inches="tight")
    plt.close(figure)


def run(planets, sim, epochs, output):
    tables = {name: [] for name in (
        "summary", "history", "columns", "layers", "reservoirs", "vertical_profiles",
    )}
    final_states = []
    for planet in planets:
        print(f"\n\nRunning {planet.name} ({epochs} epochs)...", flush=True)
        state, bulk, redox, history = sim.evolve(planet, epochs=epochs)
        row = sim.result_row(planet, state, bulk, redox)
        final_states.append((planet, state, row["P_bar"]))
        row.update({name: getattr(state, name) for name in (
            "formation_candidate_pressure_bar", "initial_target_pressure_bar",
            "formation_retained_fraction", "formation_impact_loss_kg_m2",
        )})
        tables["summary"].append(row)
        tables["columns"].extend(sim.column_rows(planet, state))
        for age, temperature, pressure, material, melt in history:
            tables["history"].append(dict(
                planet=planet.name, age_Gyr=age, T_K=temperature, P_bar=pressure,
                silicate_state=material, melt_fraction=melt,
            ))
        for index, layer in enumerate(state.layers):
            tables["layers"].append(dict(
                planet=planet.name, layer_index=index, material=layer.material.value,
                **{f.name: getattr(layer, f.name) for f in fields(layer) if f.name != "material"},
            ))
        for reservoir in ("atm", "upper_atm", "frost", "clath", "crust_vol", "escaped_kg_m2"):
            for species, amount in getattr(state, reservoir).items():
                tables["reservoirs"].append(dict(
                    planet=planet.name, reservoir=reservoir, species=species, mass_kg_m2=amount,
                ))
        for reservoir, species in (
            ("water_ice", "H2O"), ("ocean", "H2O"), ("hp_ice", "H2O"),
            ("supercrit", "H2O"), ("methane_liquid", "CH4"),
            ("organics", "mixed"), ("refractory_carbon", "C"),
        ):
            tables["reservoirs"].append(dict(
                planet=planet.name, reservoir=reservoir, species=species,
                mass_kg_m2=getattr(state, reservoir),
            ))
        for column in state.climate_columns:
            for index, level in enumerate(column.vertical_profile):
                tables["vertical_profiles"].append(dict(
                    planet=planet.name, column=column.name, level_index=index, **level,
                ))
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        write_csv(output / f"{name}.csv", rows)
    used_directories = set()
    for planet, state, pressure_bar in final_states:
        directory_name = planet_directory_name(planet.name)
        unique_name = directory_name
        suffix = 2
        while unique_name.casefold() in used_directories:
            unique_name = f"{directory_name}_{suffix}"
            suffix += 1
        used_directories.add(unique_name.casefold())
        write_planet_plots(
            planet, state, sim, pressure_bar, output / unique_name,
        )
    print(f"Completed {len(planets)} planet(s). CSV files: {output.resolve()}")


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Evolve JSON-defined planets for 4.5 Gyr using the existing simulator.",
        epilog='Examples: python run_planets.py --list; python run_planets.py --planet "Earth"; '
               'python run_planets.py --all --epochs 128 --output results/full',
    )
    parser.add_argument("--planets", type=Path, default=HERE / "example_planets.json",
                        help="Planet JSON file (default: example_planets.json beside launcher)")
    parser.add_argument("--simulator", type=Path, default=HERE / "atmo-sim.py",
                        help="Simulation library path (default: atmo-sim.py beside launcher)")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--planet", help="Run one planet by case-insensitive exact name")
    selection.add_argument("--all", action="store_true", help="Run all JSON examples")
    selection.add_argument("--list", action="store_true", help="List available planet names")
    parser.add_argument("--epochs", type=positive_int, default=128,
                        help="Evolution steps across 4.5 Gyr (default: 128)")
    parser.add_argument("--output", type=Path, default=Path("results"),
                        help="CSV directory; matching files are overwritten (default: ./results)")
    args = parser.parse_args(argv)
    try:
        sim = load_library(args.simulator.resolve())
        planets = load_planets(args.planets, sim)
        if args.list:
            print("\n".join(p.name for p in planets))
            return 0
        if args.planet:
            selected = [p for p in planets if p.name.casefold() == args.planet.casefold()]
            if not selected:
                raise ValueError(f"Unknown planet {args.planet!r}; use --list to see available names")
            planets = selected
        run(planets, sim, args.epochs, args.output)
    except (OSError, ValueError, TypeError, AssertionError, ImportError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
