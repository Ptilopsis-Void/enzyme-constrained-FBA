"""Add the luteolin pathway and module-specific helper reactions."""

from __future__ import annotations

from project_utils import ROOT  # patches COBRApy cache path before cobra import

from cobra import Metabolite, Reaction

from fba_luteolin_burden import add_luteolin_module


ENZYME_REACTIONS = [
    "LUT_PHEA_fbr",
    "LUT_PAL_Li",
    "LUT_C4H_ATR2",
    "LUT_4CL_Morus",
    "LUT_CHS_Carya",
    "LUT_CHI_virus",
    "LUT_FNS_Sm",
    "LUT_F3pH_ATR2_main",
]

OUTPUT_REACTIONS = ["LUT_LUTEOLINt", "EX_luteolin_e"]

DEMAND_METABOLITES = {
    "DM_p_coumarate_c": "lut_p_coumarate_c",
    "DM_p_coumaroyl_CoA_c": "lut_p_coumaroyl_CoA_c",
    "DM_naringenin_c": "lut_naringenin_c",
    "DM_apigenin_c": "lut_apigenin_c",
    "DM_luteolin_c": "lut_luteolin_c",
}


def add_luteolin_pathway(
    model,
    fns_type: str = "p450",
    include_side_reaction: bool = False,
    enabled_reactions: set[str] | None = None,
):
    """Add the full pathway, then optionally disable enzyme reactions by module."""

    add_luteolin_module(model, fns_type=fns_type, include_side=include_side_reaction)
    add_intermediate_demands(model)

    if enabled_reactions is not None:
        enabled = set(enabled_reactions)
        for reaction_id in ENZYME_REACTIONS:
            if reaction_id in model.reactions and reaction_id not in enabled:
                model.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)

    return model


def add_intermediate_demands(model) -> None:
    for demand_id, metabolite_id in DEMAND_METABOLITES.items():
        if demand_id in model.reactions or metabolite_id not in model.metabolites:
            continue
        metabolite = model.metabolites.get_by_id(metabolite_id)
        reaction = Reaction(demand_id)
        reaction.name = demand_id
        reaction.lower_bound = 0.0
        reaction.upper_bound = 1000.0
        reaction.add_metabolites({metabolite: -1.0})
        model.add_reactions([reaction])


def add_naringenin_feed(model, upper_bound: float = 1.0) -> None:
    """Supply naringenin for validating the late FNS/F3'H module alone."""

    if "lut_naringenin_c" not in model.metabolites:
        model.add_metabolites(
            [
                Metabolite(
                    "lut_naringenin_c",
                    name="naringenin",
                    formula="C15H12O5",
                    charge=0,
                    compartment="c",
                )
            ]
        )

    if "NARINGENIN_FEED" in model.reactions:
        reaction = model.reactions.get_by_id("NARINGENIN_FEED")
        reaction.upper_bound = upper_bound
        return

    reaction = Reaction("NARINGENIN_FEED")
    reaction.name = "naringenin feed for late-module validation"
    reaction.lower_bound = 0.0
    reaction.upper_bound = upper_bound
    reaction.add_metabolites({model.metabolites.get_by_id("lut_naringenin_c"): 1.0})
    model.add_reactions([reaction])
