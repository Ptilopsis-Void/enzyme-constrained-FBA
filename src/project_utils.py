"""Shared helpers for the luteolin FBA/ecFBA workflow."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "models" / "yeastGEM.xml"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"


try:
    import appdirs

    appdirs.user_cache_dir = lambda *args, **kwargs: str(ROOT / ".cobra_cache")
except Exception:
    pass


if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


import cobra  # noqa: E402
from cobra.io import read_sbml_model  # noqa: E402

from fba_luteolin_burden import configure_environment, detect_growth_reaction  # noqa: E402


def ensure_output_dirs() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_yeast_model(model_path: Path = MODEL_PATH) -> cobra.Model:
    path = Path(model_path)
    try:
        # libSBML can be brittle with non-ASCII absolute Windows paths. When the
        # script is launched from the project root, a relative path is safer.
        path_for_cobra = path.relative_to(Path.cwd())
    except ValueError:
        path_for_cobra = path
    return read_sbml_model(str(path_for_cobra))


def set_default_medium(model: cobra.Model) -> None:
    configure_environment(
        model,
        glucose_reaction="r_1714",
        glucose_uptake=10.0,
        oxygen_reaction="r_1992",
        oxygen_uptake=1000.0,
    )


def set_growth_objective(model: cobra.Model, growth_reaction_id: str | None = None) -> str:
    growth_id = detect_growth_reaction(model, growth_reaction_id)
    model.objective = model.reactions.get_by_id(growth_id)
    model.objective_direction = "max"
    return growth_id


def flux_value(solution: cobra.Solution, reaction_id: str) -> float:
    if solution.status != "optimal" or reaction_id not in solution.fluxes.index:
        return math.nan
    return float(solution.fluxes[reaction_id])


def write_rows(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def get_wt_growth() -> tuple[float, str]:
    model = load_yeast_model()
    set_default_medium(model)
    growth_id = set_growth_objective(model)
    solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError(f"WT growth optimization failed: {solution.status}")
    return float(solution.fluxes[growth_id]), growth_id
