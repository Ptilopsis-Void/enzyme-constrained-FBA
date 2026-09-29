"""Compare ARO7/PHA2 shikimate-end engineering scenarios."""

from __future__ import annotations

import math

from add_luteolin_pathway import add_luteolin_pathway
from fba_luteolin_burden import apply_knockouts
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


def run_scenario(
    name: str,
    add_phea: bool = True,
    knockout: str = "none",
    tyrosine_uptake: float = 0.0,
    aro7_upper_bound=None,
):
    model = load_yeast_model()
    set_default_medium(model)
    if add_phea:
        add_luteolin_pathway(model, enabled_reactions={"LUT_PHEA_fbr"})

    if aro7_upper_bound is not None and "r_0278" in model.reactions:
        reaction = model.reactions.get_by_id("r_0278")
        reaction.upper_bound = min(reaction.upper_bound, float(aro7_upper_bound))

    apply_knockouts(model, knockout)

    if tyrosine_uptake > 0 and "r_1913" in model.reactions:
        reaction = model.reactions.get_by_id("r_1913")
        reaction.lower_bound = -abs(tyrosine_uptake)
        reaction.upper_bound = 1000.0

    growth_id = set_growth_objective(model)
    solution = model.optimize()
    growth = flux_value(solution, growth_id)
    return {
        "scenario": name,
        "status": solution.status,
        "growth": growth,
        "pheA_flux": flux_value(solution, "LUT_PHEA_fbr")
        if "LUT_PHEA_fbr" in model.reactions
        else 0.0,
        "aro7_flux": flux_value(solution, "r_0278"),
        "pha2_flux": flux_value(solution, "r_0938"),
        "tyrosine_exchange_flux": flux_value(solution, "r_1913"),
    }


def main() -> int:
    ensure_output_dirs()
    wt_growth, _ = get_wt_growth()
    wt_row = run_scenario("WT", add_phea=False)
    wt_aro7_flux = wt_row["aro7_flux"]
    wt_row["growth_fraction_of_wt"] = wt_row["growth"] / wt_growth
    wt_row["burden_percent_vs_wt"] = 100 * (1 - wt_row["growth"] / wt_growth)
    rows = [wt_row]
    cases = [
        ("PheA_fbr_only", True, "none", 0.0, None),
        ("PHA2_KO_plus_PheA_fbr", True, "pha2", 0.0, None),
        ("ARO7_KO_plus_PheA_fbr", True, "aro7", 0.0, None),
        # Diagnostic rescue bound, not a measured uptake rate.
        ("ARO7_KO_plus_PheA_fbr_tyrosine_uptake", True, "aro7", 1.0, None),
        (
            "ARO7_weakened_to_10pct_WT_flux_plus_PheA_fbr",
            True,
            "none",
            0.0,
            # Sensitivity scenario: retain 10% of the WT-model ARO7 flux.
            0.1 * wt_aro7_flux,
        ),
    ]

    for name, add_phea, knockout, tyr_uptake, aro7_upper_bound in cases:
        row = run_scenario(name, add_phea, knockout, tyr_uptake, aro7_upper_bound)
        growth = row["growth"]
        row["growth_fraction_of_wt"] = growth / wt_growth if not math.isnan(growth) else math.nan
        row["burden_percent_vs_wt"] = (
            100 * (1 - growth / wt_growth) if not math.isnan(growth) else math.nan
        )
        rows.append(row)

    output = RESULTS_DIR / "aro7_pha2_comparison.csv"
    write_rows(rows, output)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
