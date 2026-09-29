"""Run module-by-module growth and endpoint-demand scans."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from add_ec_constraints import (
    add_fixed_expression_atp_cost,
    apply_light_ec_constraints,
    pool_usage_fraction,
)
from add_luteolin_pathway import ENZYME_REACTIONS, add_luteolin_pathway, add_naringenin_feed
from fba_luteolin_burden import apply_capacity_modifiers
from project_utils import (
    RESULTS_DIR,
    ROOT,
    ensure_output_dirs,
    flux_value,
    get_wt_growth,
    load_yeast_model,
    set_default_medium,
    set_growth_objective,
    write_rows,
)


def enabled_reactions_from_row(row) -> set[str]:
    enabled = set()
    for reaction_id in ENZYME_REACTIONS:
        if reaction_id in row and int(row[reaction_id]) == 1:
            enabled.add(reaction_id)
    return enabled


def run_module_row(args, row, wt_growth: float, growth_id: str) -> dict:
    model = load_yeast_model()
    set_default_medium(model)

    enabled = enabled_reactions_from_row(row)
    if int(row.get("apply_capacity_modifiers", 0)) == 1:
        apply_capacity_modifiers(model)

    if enabled or int(row.get("feed_naringenin", 0)) == 1:
        add_luteolin_pathway(model, enabled_reactions=enabled)

    if int(row.get("feed_naringenin", 0)) == 1:
        add_naringenin_feed(model, upper_bound=args.naringenin_feed)

    if enabled:
        apply_light_ec_constraints(
            model,
            args.enzymes,
            heterologous_protein_pool=args.heterologous_pool,
            p450_membrane_pool=args.p450_pool,
            enabled_reactions=enabled,
        )
        fixed_cost = add_fixed_expression_atp_cost(
            model,
            args.enzymes,
            enabled,
            atp_cost_per_relative_expression=args.expression_atp_cost_per_rel,
        )
    else:
        fixed_cost = 0.0

    endpoint = str(row.get("endpoint_reaction", "")).strip()
    endpoint_min = float(row.get("endpoint_min", 0.0) or 0.0)
    if endpoint and endpoint.lower() != "nan" and endpoint_min > 0:
        if endpoint not in model.reactions:
            raise ValueError(f"Endpoint reaction not found for {row['module_id']}: {endpoint}")
        model.reactions.get_by_id(endpoint).lower_bound = endpoint_min

    set_growth_objective(model, growth_id)
    solution = model.optimize()
    growth = flux_value(solution, growth_id)
    endpoint_flux = flux_value(solution, endpoint) if endpoint and endpoint in model.reactions else 0.0

    return {
        "module_id": row["module_id"],
        "description": row.get("description", ""),
        "enabled_reaction_count": len(enabled),
        "feed_naringenin": int(row.get("feed_naringenin", 0)),
        "endpoint_reaction": endpoint,
        "endpoint_min": endpoint_min,
        "status": solution.status,
        "growth": growth,
        "growth_fraction_of_wt": growth / wt_growth if not math.isnan(growth) else math.nan,
        "burden_percent_vs_wt": 100 * (1 - growth / wt_growth)
        if not math.isnan(growth)
        else math.nan,
        "endpoint_flux": endpoint_flux,
        "protein_pool_usage": flux_value(solution, "PROT_POOL_HET_SUPPLY"),
        "protein_pool_usage_fraction": pool_usage_fraction(
            solution, "PROT_POOL_HET_SUPPLY", args.heterologous_pool
        ),
        "p450_pool_usage": flux_value(solution, "P450_POOL_SUPPLY"),
        "p450_pool_usage_fraction": pool_usage_fraction(
            solution, "P450_POOL_SUPPLY", args.p450_pool
        ),
        "fixed_expression_atp_cost": fixed_cost,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare the pathway modules under one ecFBA setting.")
    parser.add_argument(
        "--modules",
        default=str(ROOT / "data" / "module_design.csv"),
        help="CSV describing the reactions enabled in each module.",
    )
    parser.add_argument(
        "--enzymes",
        default=str(ROOT / "data" / "enzymes.csv"),
        help="CSV containing the simplified enzyme parameters.",
    )
    parser.add_argument(
        "--heterologous-pool",
        type=float,
        default=0.02,
        help="Total heterologous protein-pool assumption.",
    )
    parser.add_argument(
        "--p450-pool",
        type=float,
        default=0.005,
        help="P450 protein-pool assumption.",
    )
    parser.add_argument(
        "--naringenin-feed",
        type=float,
        default=1.0,
        help="Diagnostic naringenin supply bound for the late module.",
    )
    parser.add_argument(
        "--expression-atp-cost-per-rel",
        type=float,
        default=0.02,
        help="Fixed ATP cost per relative-expression unit.",
    )
    parser.add_argument(
        "--output",
        default=str(RESULTS_DIR / "module_scan.csv"),
        help="Output CSV path.",
    )
    args = parser.parse_args()

    ensure_output_dirs()
    wt_growth, growth_id = get_wt_growth()
    design = pd.read_csv(args.modules)
    rows = []
    for _, row in design.iterrows():
        try:
            rows.append(run_module_row(args, row, wt_growth, growth_id))
        except Exception as exc:
            rows.append(
                {
                    "module_id": row["module_id"],
                    "description": row.get("description", ""),
                    "enabled_reaction_count": math.nan,
                    "feed_naringenin": row.get("feed_naringenin", 0),
                    "endpoint_reaction": row.get("endpoint_reaction", ""),
                    "endpoint_min": row.get("endpoint_min", 0.0),
                    "status": "error",
                    "growth": math.nan,
                    "growth_fraction_of_wt": math.nan,
                    "burden_percent_vs_wt": math.nan,
                    "endpoint_flux": math.nan,
                    "protein_pool_usage": math.nan,
                    "protein_pool_usage_fraction": math.nan,
                    "p450_pool_usage": math.nan,
                    "p450_pool_usage_fraction": math.nan,
                    "fixed_expression_atp_cost": math.nan,
                    "error": str(exc),
                }
            )

    output = Path(args.output)
    write_rows(rows, output)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
