"""Lightweight enzyme-constrained FBA constraints for the luteolin pathway."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
from project_utils import ROOT  # patches COBRApy cache path before cobra import
from cobra import Metabolite, Reaction


REACTION_ALIASES = {
    "PHEA_fbr": "LUT_PHEA_fbr",
    "PAL_Li": "LUT_PAL_Li",
    "C4H_ATR2": "LUT_C4H_ATR2",
    "4CL_Morus": "LUT_4CL_Morus",
    "CHS_Carya": "LUT_CHS_Carya",
    "CHI_virus": "LUT_CHI_virus",
    "FNS_Sm": "LUT_FNS_Sm",
    "F3pH_ATR2_main": "LUT_F3pH_ATR2_main",
    "F3H_ATR2_main": "LUT_F3pH_ATR2_main",
}


def get_or_create_metabolite(model, metabolite_id, name, compartment="c"):
    if metabolite_id in model.metabolites:
        return model.metabolites.get_by_id(metabolite_id)

    metabolite = Metabolite(metabolite_id, name=name, compartment=compartment)
    model.add_metabolites([metabolite])
    return metabolite


def add_supply_reaction(model, metabolite, reaction_id, upper_bound):
    if reaction_id in model.reactions:
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.upper_bound = float(upper_bound)
        return reaction

    reaction = Reaction(reaction_id)
    reaction.name = reaction_id
    reaction.lower_bound = 0.0
    reaction.upper_bound = float(upper_bound)
    reaction.add_metabolites({metabolite: 1.0})
    model.add_reactions([reaction])
    return reaction


def resolve_reaction_id(model, reaction_id: str) -> str | None:
    if reaction_id in model.reactions:
        return reaction_id
    alias = REACTION_ALIASES.get(reaction_id)
    if alias and alias in model.reactions:
        return alias
    return None


def read_enzyme_table(enzyme_csv: str | Path) -> pd.DataFrame:
    enzymes = pd.read_csv(enzyme_csv)
    required = {"reaction_id", "mw_kda", "kcat_s", "is_p450"}
    missing = required.difference(enzymes.columns)
    if missing:
        raise ValueError(f"Missing required enzyme columns: {sorted(missing)}")
    return enzymes


def apply_light_ec_constraints(
    model,
    enzyme_csv,
    heterologous_protein_pool=0.02,
    p450_membrane_pool=0.005,
    enabled_reactions: set[str] | None = None,
):
    """Add a lightweight enzyme-capacity layer.

    Flux is mmol/gDW/h, kcat is s^-1, and MW in kDa is treated as g/mmol.
    The resulting coefficient MW / (kcat * 3600) is an approximate protein
    cost in g protein per mmol reaction flux.
    """

    enzyme_df = read_enzyme_table(enzyme_csv)
    enabled = set(enabled_reactions) if enabled_reactions is not None else None

    protein_pool = get_or_create_metabolite(
        model,
        "prot_pool_heterologous_c",
        "heterologous protein pool",
        "c",
    )
    p450_pool = get_or_create_metabolite(
        model,
        "p450_membrane_pool_c",
        "P450 membrane expression pool",
        "c",
    )

    add_supply_reaction(
        model,
        protein_pool,
        "PROT_POOL_HET_SUPPLY",
        heterologous_protein_pool,
    )
    add_supply_reaction(
        model,
        p450_pool,
        "P450_POOL_SUPPLY",
        p450_membrane_pool,
    )

    applied_rows = []
    for _, row in enzyme_df.iterrows():
        reaction_id = resolve_reaction_id(model, str(row["reaction_id"]))
        if reaction_id is None:
            continue
        if enabled is not None and reaction_id not in enabled:
            continue

        kcat_s = float(row["kcat_s"])
        if kcat_s <= 0:
            raise ValueError(f"kcat must be positive for {row['reaction_id']}")

        mw_kda = float(row["mw_kda"])
        enzyme_cost = mw_kda / (kcat_s * 3600.0)
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.add_metabolites({protein_pool: -enzyme_cost})

        is_p450 = int(row["is_p450"]) == 1
        if is_p450:
            reaction.add_metabolites({p450_pool: -enzyme_cost})

        applied_rows.append(
            {
                "reaction_id": reaction_id,
                "enzyme_cost": enzyme_cost,
                "is_p450": is_p450,
            }
        )

    return pd.DataFrame(applied_rows)


def add_fixed_expression_atp_cost(
    model,
    enzyme_csv,
    enabled_reactions: set[str],
    atp_cost_per_relative_expression: float = 0.0,
) -> float:
    """Approximate expression/folding burden as fixed ATP hydrolysis.

    This complements the flux-coupled enzyme pool. It is intentionally tunable:
    use experimental module growth rates to calibrate the coefficient.
    """

    if atp_cost_per_relative_expression <= 0 or not enabled_reactions:
        return 0.0

    enzyme_df = read_enzyme_table(enzyme_csv)
    relative_total = 0.0
    for _, row in enzyme_df.iterrows():
        reaction_id = resolve_reaction_id(model, str(row["reaction_id"]))
        if reaction_id in enabled_reactions:
            relative_total += float(row.get("relative_expression", 1.0))

    fixed_cost = relative_total * atp_cost_per_relative_expression
    if fixed_cost <= 0:
        return 0.0

    metabolite_ids = {
        "atp": "s_0434",
        "h2o": "s_0803",
        "adp": "s_0394",
        "pi": "s_1322",
        "h": "s_0794",
    }
    if any(mid not in model.metabolites for mid in metabolite_ids.values()):
        return 0.0

    if "FIXED_HET_EXPR_ATP_COST" in model.reactions:
        reaction = model.reactions.get_by_id("FIXED_HET_EXPR_ATP_COST")
        reaction.bounds = (fixed_cost, fixed_cost)
        return fixed_cost

    reaction = Reaction("FIXED_HET_EXPR_ATP_COST")
    reaction.name = "fixed ATP cost for heterologous expression and folding"
    reaction.lower_bound = fixed_cost
    reaction.upper_bound = fixed_cost
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id(metabolite_ids["atp"]): -1.0,
            model.metabolites.get_by_id(metabolite_ids["h2o"]): -1.0,
            model.metabolites.get_by_id(metabolite_ids["adp"]): 1.0,
            model.metabolites.get_by_id(metabolite_ids["pi"]): 1.0,
            model.metabolites.get_by_id(metabolite_ids["h"]): 1.0,
        }
    )
    model.add_reactions([reaction])
    return fixed_cost


def pool_usage_fraction(solution, supply_reaction_id: str, pool_bound: float) -> float:
    if solution.status != "optimal" or supply_reaction_id not in solution.fluxes.index:
        return math.nan
    if pool_bound <= 0:
        return math.nan
    return float(solution.fluxes[supply_reaction_id]) / pool_bound
