"""Run ordinary FBA baseline scenarios for the luteolin pathway."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from add_luteolin_pathway import add_luteolin_pathway
from project_utils import (
    RESULTS_DIR,
    ensure_output_dirs,
    flux_value,
    get_wt_growth,
    load_yeast_model,
    set_default_medium,
    set_growth_objective,
    write_rows,
)


def parse_float_list(value: str) -> list[float]:
    if value.strip().lower() in {"none", "skip", "off"}:
        return []
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the ordinary FBA baseline scan.")
    parser.add_argument(
        "--luteolin-demands",
        default="0,0.001,0.005,0.01,0.05",
        help="Comma-separated minimum luteolin fluxes (mmol/gDW/h).",
    )
    parser.add_argument(
        "--output",
        default=str(RESULTS_DIR / "fba_baseline.csv"),
        help="Output CSV path.",
    )
    args = parser.parse_args()

    ensure_output_dirs()
    wt_growth, growth_id = get_wt_growth()
    rows = [
        {
            "scenario": "WT",
            "status": "optimal",
            "growth": wt_growth,
            "growth_fraction_of_wt": 1.0,
            "burden_percent_vs_wt": 0.0,
            "luteolin_min": 0.0,
            "luteolin_flux": 0.0,
        }
    ]

    for demand in parse_float_list(args.luteolin_demands):
        model = load_yeast_model()
        set_default_medium(model)
        add_luteolin_pathway(model)
        set_growth_objective(model, growth_id)
        model.reactions.get_by_id("EX_luteolin_e").lower_bound = demand
        solution = model.optimize()
        growth = flux_value(solution, growth_id)
        rows.append(
            {
                "scenario": f"ordinary_FBA_luteolin_min_{demand:g}",
                "status": solution.status,
                "growth": growth,
                "growth_fraction_of_wt": growth / wt_growth if not math.isnan(growth) else math.nan,
                "burden_percent_vs_wt": 100 * (1 - growth / wt_growth)
                if not math.isnan(growth)
                else math.nan,
                "luteolin_min": demand,
                "luteolin_flux": flux_value(solution, "EX_luteolin_e"),
            }
        )

    output = Path(args.output)
    write_rows(rows, output)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
