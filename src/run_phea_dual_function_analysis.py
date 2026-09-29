"""Corrected FBA analysis for bifunctional PheA_fbr in aro7Delta yeast.

The earlier one-step PheA_fbr representation directly converted chorismate to
phenylpyruvate. Literature/experimental evidence indicates that engineered
PheA_fbr has a functional chorismate mutase domain in yeast and therefore
produces a shared prephenate pool that TYR1 can still use.
"""

from __future__ import annotations

import math

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
from cobra import Reaction
from fba_luteolin_burden import apply_knockouts


def add_old_one_step_phea(model, upper_bound: float = 1000.0) -> None:
    """Old simplified representation: chorismate directly to phenylpyruvate."""

    if "PHEA_fbr_old_one_step" in model.reactions:
        model.reactions.get_by_id("PHEA_fbr_old_one_step").upper_bound = upper_bound
        return

    reaction = Reaction("PHEA_fbr_old_one_step")
    reaction.name = "old simplified PheA_fbr chorismate to phenylpyruvate"
    reaction.lower_bound = 0.0
    reaction.upper_bound = upper_bound
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id("s_0515"): -1.0,  # chorismate
            model.metabolites.get_by_id("s_0794"): -1.0,  # H+
            model.metabolites.get_by_id("s_0951"): 1.0,  # phenylpyruvate
            model.metabolites.get_by_id("s_0456"): 1.0,  # CO2
            model.metabolites.get_by_id("s_0803"): 1.0,  # H2O
        }
    )
    model.add_reactions([reaction])
    reaction.gene_reaction_rule = "pheA_fbr"


def add_phea_dual_function(model, cm_upper_bound: float = 1000.0, pdt_upper_bound: float = 1000.0) -> None:
    """Corrected bifunctional PheA_fbr representation.

    PheA_fbr_CM:
        chorismate -> prephenate

    PheA_fbr_PDT:
        prephenate + H+ -> phenylpyruvate + CO2 + H2O
    """

    if "PHEA_fbr_CM" not in model.reactions:
        cm = Reaction("PHEA_fbr_CM")
        cm.name = "PheA_fbr chorismate mutase domain"
        cm.lower_bound = 0.0
        cm.upper_bound = cm_upper_bound
        cm.add_metabolites(
            {
                model.metabolites.get_by_id("s_0515"): -1.0,  # chorismate
                model.metabolites.get_by_id("s_1377"): 1.0,  # prephenate
            }
        )
        model.add_reactions([cm])
        cm.gene_reaction_rule = "pheA_fbr"
    else:
        model.reactions.get_by_id("PHEA_fbr_CM").upper_bound = cm_upper_bound

    if "PHEA_fbr_PDT" not in model.reactions:
        pdt = Reaction("PHEA_fbr_PDT")
        pdt.name = "PheA_fbr prephenate dehydratase domain"
        pdt.lower_bound = 0.0
        pdt.upper_bound = pdt_upper_bound
        pdt.add_metabolites(
            {
                model.metabolites.get_by_id("s_0794"): -1.0,  # H+
                model.metabolites.get_by_id("s_1377"): -1.0,  # prephenate
                model.metabolites.get_by_id("s_0456"): 1.0,  # CO2
                model.metabolites.get_by_id("s_0803"): 1.0,  # H2O
                model.metabolites.get_by_id("s_0951"): 1.0,  # phenylpyruvate
            }
        )
        model.add_reactions([pdt])
        pdt.gene_reaction_rule = "pheA_fbr"
    else:
        model.reactions.get_by_id("PHEA_fbr_PDT").upper_bound = pdt_upper_bound


def solve_scenario(
    name: str,
    *,
    phea_mode: str = "none",
    knockout: str = "none",
    tyrosine_uptake: float = 0.0,
    tyr1_upper_bound: float | None = None,
) -> dict:
    model = load_yeast_model()
    set_default_medium(model)

    if phea_mode == "old_one_step":
        add_old_one_step_phea(model)
    elif phea_mode == "cm_only":
        add_phea_dual_function(model, pdt_upper_bound=0.0)
    elif phea_mode == "dual":
        add_phea_dual_function(model)
    elif phea_mode != "none":
        raise ValueError(f"Unknown phea_mode: {phea_mode}")

    apply_knockouts(model, knockout)

    if tyrosine_uptake > 0.0 and "r_1913" in model.reactions:
        exchange = model.reactions.get_by_id("r_1913")
        exchange.lower_bound = -abs(tyrosine_uptake)
        exchange.upper_bound = 1000.0

    if tyr1_upper_bound is not None:
        model.reactions.get_by_id("r_0939").upper_bound = tyr1_upper_bound

    growth_id = set_growth_objective(model)
    solution = model.optimize()
    growth = flux_value(solution, growth_id)

    return {
        "scenario": name,
        "phea_mode": phea_mode,
        "knockout": knockout,
        "tyrosine_uptake_bound": tyrosine_uptake,
        "tyr1_upper_bound": tyr1_upper_bound if tyr1_upper_bound is not None else "",
        "status": solution.status,
        "growth": growth,
        "aro7_flux": flux_value(solution, "r_0278"),
        "native_pha2_flux": flux_value(solution, "r_0938"),
        "tyr1_flux": flux_value(solution, "r_0939"),
        "old_one_step_phea_flux": flux_value(solution, "PHEA_fbr_old_one_step")
        if "PHEA_fbr_old_one_step" in model.reactions
        else 0.0,
        "phea_cm_flux": flux_value(solution, "PHEA_fbr_CM")
        if "PHEA_fbr_CM" in model.reactions
        else 0.0,
        "phea_pdt_flux": flux_value(solution, "PHEA_fbr_PDT")
        if "PHEA_fbr_PDT" in model.reactions
        else 0.0,
        "tyrosine_exchange_flux": flux_value(solution, "r_1913"),
    }


def main() -> int:
    ensure_output_dirs()
    wt_growth, _ = get_wt_growth()
    reference = solve_scenario("WT_reference")
    reference_tyr1 = reference["tyr1_flux"]
    # HMA-based sensitivity test; this is not a measured TYR1 capacity.
    tyr1_25pct = 0.25 * reference_tyr1 if not math.isnan(reference_tyr1) else None

    rows = [
        reference,
        solve_scenario("PheA_dual_no_KO", phea_mode="dual"),
        solve_scenario(
            "ARO7_KO_old_one_step_PheA_false_negative",
            phea_mode="old_one_step",
            knockout="aro7",
        ),
        solve_scenario(
            "ARO7_KO_old_one_step_PheA_tyrosine_rescue",
            phea_mode="old_one_step",
            knockout="aro7",
            tyrosine_uptake=1.0,
        ),
        solve_scenario("ARO7_KO_PheA_CM_only", phea_mode="cm_only", knockout="aro7"),
        solve_scenario("ARO7_KO_PheA_dual", phea_mode="dual", knockout="aro7"),
        solve_scenario(
            "ARO7_KO_PheA_dual_TYR1_cap_25pct_WT_flux",
            phea_mode="dual",
            knockout="aro7",
            tyr1_upper_bound=tyr1_25pct,
        ),
        solve_scenario(
            "ARO7_PHA2_KO_PheA_dual",
            phea_mode="dual",
            knockout="both",
        ),
        solve_scenario(
            "PHA2_KO_PheA_dual",
            phea_mode="dual",
            knockout="pha2",
        ),
    ]

    for row in rows:
        growth = row["growth"]
        row["growth_fraction_of_wt"] = (
            growth / wt_growth if not math.isnan(growth) else math.nan
        )
        row["burden_percent_vs_wt"] = (
            100 * (1 - growth / wt_growth) if not math.isnan(growth) else math.nan
        )
        row["tyr1_fraction_of_wt"] = (
            row["tyr1_flux"] / reference_tyr1
            if reference_tyr1 and not math.isnan(row["tyr1_flux"])
            else math.nan
        )
        row["phe_branch_flux_proxy"] = (
            row["native_pha2_flux"]
            + row["phea_pdt_flux"]
            + row["old_one_step_phea_flux"]
        )
        row["ma_hma_branch_proxy"] = (
            row["phe_branch_flux_proxy"] / row["tyr1_flux"]
            if row["tyr1_flux"] and row["tyr1_flux"] > 0
            else math.nan
        )

    output = RESULTS_DIR / "phea_dual_function_analysis.csv"
    write_rows(rows, output)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
