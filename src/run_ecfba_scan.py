"""Scan lightweight ecFBA protein-pool and P450-pool assumptions."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from add_ec_constraints import (
    add_fixed_expression_atp_cost,
    apply_light_ec_constraints,
    pool_usage_fraction,
)
from add_luteolin_pathway import ENZYME_REACTIONS, add_luteolin_pathway
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


def parse_float_list(value: str) -> list[float]:
    if value.strip().lower() in {"none", "skip", "off"}:
        return []
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def build_model(args, heterologous_pool: float, p450_pool: float):
    model = load_yeast_model()
    set_default_medium(model)
    add_luteolin_pathway(model, fns_type=args.fns_type)
    enabled = set(ENZYME_REACTIONS)
    apply_light_ec_constraints(
        model,
        args.enzymes,
        heterologous_protein_pool=heterologous_pool,
        p450_membrane_pool=p450_pool,
        enabled_reactions=enabled,
    )
    fixed_cost = add_fixed_expression_atp_cost(
        model,
        args.enzymes,
        enabled,
        atp_cost_per_relative_expression=args.expression_atp_cost_per_rel,
    )
    return model, fixed_cost


def run_case(args, wt_growth: float, growth_id: str, heterologous_pool, p450_pool, demand):
    base_model, fixed_cost = build_model(args, heterologous_pool, p450_pool)
    set_growth_objective(base_model, growth_id)
    no_demand_solution = base_model.optimize()
    growth_without_demand = flux_value(no_demand_solution, growth_id)

    base_model.reactions.get_by_id("EX_luteolin_e").lower_bound = demand
    set_growth_objective(base_model, growth_id)
    solution = base_model.optimize()
    growth = flux_value(solution, growth_id)

    return {
        "heterologous_pool": heterologous_pool,
        "p450_pool": p450_pool,
        "luteolin_min": demand,
        "status": solution.status,
        "growth": growth,
        "growth_fraction_of_wt": growth / wt_growth if not math.isnan(growth) else math.nan,
        "burden_percent_vs_wt": 100 * (1 - growth / wt_growth)
        if not math.isnan(growth)
        else math.nan,
        "growth_without_luteolin_demand": growth_without_demand,
        "luteolin_flux": flux_value(solution, "EX_luteolin_e"),
        "protein_pool_usage": flux_value(solution, "PROT_POOL_HET_SUPPLY"),
        "protein_pool_usage_fraction": pool_usage_fraction(
            solution, "PROT_POOL_HET_SUPPLY", heterologous_pool
        ),
        "p450_pool_usage": flux_value(solution, "P450_POOL_SUPPLY"),
        "p450_pool_usage_fraction": pool_usage_fraction(
            solution, "P450_POOL_SUPPLY", p450_pool
        ),
        "fixed_expression_atp_cost": fixed_cost,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Scan the lightweight ecFBA assumptions.")
    parser.add_argument(
        "--enzymes",
        default=str(ROOT / "data" / "enzymes.csv"),
        help="CSV containing MW, kcat, P450 labels and expression weights.",
    )
    parser.add_argument(
        "--fns-type",
        choices=["p450", "2og"],
        default="p450",
        help="Use the P450 FNSII route or the 2-oxoglutarate FNSI route.",
    )
    parser.add_argument(
        "--heterologous-pools",
        default="0.005,0.01,0.02,0.05",
        help="Comma-separated total heterologous protein-pool assumptions.",
    )
    parser.add_argument(
        "--p450-pools",
        default="0.001,0.002,0.005,0.01",
        help="Comma-separated P450 protein-pool assumptions.",
    )
    parser.add_argument(
        "--luteolin-demands",
        default="0,0.001,0.005,0.01,0.05",
        help="Comma-separated minimum luteolin fluxes (mmol/gDW/h).",
    )
    parser.add_argument(
        "--expression-atp-cost-per-rel",
        type=float,
        default=0.02,
        help="Fixed ATP hydrolysis flux per unit relative expression.",
    )
    parser.add_argument(
        "--output",
        default=str(RESULTS_DIR / "ecfba_scan.csv"),
        help="Output CSV path.",
    )
    args = parser.parse_args()

    ensure_output_dirs()
    wt_growth, growth_id = get_wt_growth()
    rows = []
    for heterologous_pool in parse_float_list(args.heterologous_pools):
        for p450_pool in parse_float_list(args.p450_pools):
            for demand in parse_float_list(args.luteolin_demands):
                try:
                    rows.append(
                        run_case(
                            args,
                            wt_growth,
                            growth_id,
                            heterologous_pool,
                            p450_pool,
                            demand,
                        )
                    )
                except Exception as exc:
                    rows.append(
                        {
                            "heterologous_pool": heterologous_pool,
                            "p450_pool": p450_pool,
                            "luteolin_min": demand,
                            "status": "error",
                            "growth": math.nan,
                            "growth_fraction_of_wt": math.nan,
                            "burden_percent_vs_wt": math.nan,
                            "growth_without_luteolin_demand": math.nan,
                            "luteolin_flux": math.nan,
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
