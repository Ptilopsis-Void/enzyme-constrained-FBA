#!/usr/bin/env python
"""Compare growth burden before/after inserting a luteolin pathway in yeastGEM.

The script loads yeast-GEM SBML, adds a cytosolic phenylpropanoid/flavone
module, and runs FBA with growth as the objective in every burden scenario.

Important modeling note:
    Adding reactions/genes alone does not create a growth burden in classical
    FBA. To quantify burden, constrain either a minimum luteolin secretion flux
    or an optional ATP cost for heterologous expression.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]


# COBRApy 0.31 tries to create a cache under AppData on this machine. Keep all
# generated cache files local to the working directory so the script runs inside
# a workspace sandbox.
try:
    import appdirs

    appdirs.user_cache_dir = lambda *args, **kwargs: str(ROOT / ".cobra_cache")
except Exception:
    pass

try:
    import cobra
    from cobra import Metabolite, Reaction
    from cobra.io import read_sbml_model
except ImportError as exc:  # pragma: no cover - friendly CLI failure
    raise SystemExit(
        "COBRApy is required. Install it with `pip install cobra` and rerun."
    ) from exc


DEFAULT_MODEL = ROOT / "models" / "yeastGEM.xml"
DEFAULT_OUTPUT = ROOT / "results" / "fba_luteolin_burden_results.csv"
EPS = 1e-9


@dataclass(frozen=True)
class ExistingMet:
    id: str
    name: str
    compartment: str


@dataclass(frozen=True)
class NewMet:
    id: str
    name: str
    compartment: str
    formula: str | None = None
    charge: int | None = None


@dataclass(frozen=True)
class HeterologousReaction:
    id: str
    name: str
    stoich: dict[str, float]
    lower_bound: float
    upper_bound: float
    gene_rule: str
    relative_expression: float | None = None


EXISTING_METS: dict[str, ExistingMet] = {
    "chorismate_c": ExistingMet("s_0515", "chorismate", "c"),
    "phenylpyruvate_c": ExistingMet("s_0951", "keto-phenylpyruvate", "c"),
    "phe_l_c": ExistingMet("s_1032", "L-phenylalanine", "c"),
    "nh4_c": ExistingMet("s_0419", "ammonium", "c"),
    "nadph_c": ExistingMet("s_1212", "NADPH", "c"),
    "nadp_c": ExistingMet("s_1207", "NADP(+)", "c"),
    "o2_c": ExistingMet("s_1275", "oxygen", "c"),
    "h_c": ExistingMet("s_0794", "H+", "c"),
    "h2o_c": ExistingMet("s_0803", "H2O", "c"),
    "co2_c": ExistingMet("s_0456", "carbon dioxide", "c"),
    "coa_c": ExistingMet("s_0529", "coenzyme A", "c"),
    "atp_c": ExistingMet("s_0434", "ATP", "c"),
    "adp_c": ExistingMet("s_0394", "ADP", "c"),
    "amp_c": ExistingMet("s_0423", "AMP", "c"),
    "ppi_c": ExistingMet("s_0633", "diphosphate", "c"),
    "pi_c": ExistingMet("s_1322", "phosphate", "c"),
    "malcoa_c": ExistingMet("s_1101", "malonyl-CoA", "c"),
    "akg_c": ExistingMet("s_0180", "2-oxoglutarate", "c"),
    "succ_c": ExistingMet("s_1458", "succinate", "c"),
}


NEW_METS: dict[str, NewMet] = {
    "trans_cinnamate_c": NewMet(
        "lut_trans_cinnamate_c", "trans-cinnamate", "c", "C9H7O2", -1
    ),
    "p_coumarate_c": NewMet(
        "lut_p_coumarate_c", "p-coumarate", "c", "C9H7O3", -1
    ),
    "p_coumaroyl_coa_c": NewMet(
        "lut_p_coumaroyl_CoA_c",
        "p-coumaroyl-CoA",
        "c",
        "C30H37N7O18P3S",
        -4,
    ),
    "naringenin_chalcone_c": NewMet(
        "lut_naringenin_chalcone_c",
        "naringenin chalcone",
        "c",
        "C15H12O5",
        0,
    ),
    "naringenin_c": NewMet(
        "lut_naringenin_c", "naringenin", "c", "C15H12O5", 0
    ),
    "apigenin_c": NewMet("lut_apigenin_c", "apigenin", "c", "C15H10O5", 0),
    "luteolin_c": NewMet("lut_luteolin_c", "luteolin", "c", "C15H10O6", 0),
    "luteolin_e": NewMet("lut_luteolin_e", "luteolin", "e", "C15H10O6", 0),
    "eriodictyol_c": NewMet(
        "lut_eriodictyol_c", "eriodictyol", "c", "C15H12O6", 0
    ),
    "eriodictyol_e": NewMet(
        "lut_eriodictyol_e", "eriodictyol", "e", "C15H12O6", 0
    ),
}


def parse_float_list(value: str) -> list[float]:
    if value.strip().lower() in {"none", "skip", "off", "false"}:
        return []
    if not value:
        return []
    values: list[float] = []
    for part in value.split(","):
        part = part.strip()
        if part:
            values.append(float(part))
    return values


def get_existing_metabolite(model: cobra.Model, key: str) -> Metabolite:
    spec = EXISTING_METS[key]
    if spec.id in model.metabolites:
        return model.metabolites.get_by_id(spec.id)

    hits = [
        met
        for met in model.metabolites
        if met.name.lower() == spec.name.lower() and met.compartment == spec.compartment
    ]
    if len(hits) == 1:
        return hits[0]

    raise KeyError(
        f"Could not find metabolite {spec.name!r} in compartment {spec.compartment!r}"
    )


def ensure_new_metabolite(model: cobra.Model, key: str) -> Metabolite:
    spec = NEW_METS[key]
    if spec.id in model.metabolites:
        return model.metabolites.get_by_id(spec.id)

    met = Metabolite(
        spec.id,
        name=spec.name,
        formula=spec.formula,
        charge=spec.charge,
        compartment=spec.compartment,
    )
    model.add_metabolites([met])
    return met


def build_metabolite_map(model: cobra.Model, include_side: bool) -> dict[str, Metabolite]:
    met_map: dict[str, Metabolite] = {
        key: get_existing_metabolite(model, key) for key in EXISTING_METS
    }

    required_new = [
        "trans_cinnamate_c",
        "p_coumarate_c",
        "p_coumaroyl_coa_c",
        "naringenin_chalcone_c",
        "naringenin_c",
        "apigenin_c",
        "luteolin_c",
        "luteolin_e",
    ]
    if include_side:
        required_new.extend(["eriodictyol_c", "eriodictyol_e"])

    for key in required_new:
        met_map[key] = ensure_new_metabolite(model, key)

    return met_map


def add_reaction(
    model: cobra.Model, met_map: dict[str, Metabolite], reaction: HeterologousReaction
) -> None:
    if reaction.id in model.reactions:
        model.remove_reactions([model.reactions.get_by_id(reaction.id)])

    rxn = Reaction(reaction.id, name=reaction.name)
    rxn.lower_bound = reaction.lower_bound
    rxn.upper_bound = reaction.upper_bound
    rxn.add_metabolites(
        {met_map[met_key]: coefficient for met_key, coefficient in reaction.stoich.items()}
    )
    model.add_reactions([rxn])
    rxn.gene_reaction_rule = reaction.gene_rule


def luteolin_reactions(fns_type: str, include_side: bool) -> list[HeterologousReaction]:
    reactions = [
        HeterologousReaction(
            "LUT_PHEA_fbr",
            "feedback-resistant PheA chorismate to phenylpyruvate",
            {
                "chorismate_c": -1,
                "h_c": -1,
                "phenylpyruvate_c": 1,
                "co2_c": 1,
                "h2o_c": 1,
            },
            0,
            8,
            "pheA_fbr",
            0.80,
        ),
        HeterologousReaction(
            "LUT_PAL_Li",
            "phenylalanine ammonia-lyase",
            {"phe_l_c": -1, "trans_cinnamate_c": 1, "nh4_c": 1},
            0,
            6,
            "PAL",
            0.60,
        ),
        HeterologousReaction(
            "LUT_C4H_ATR2",
            "cinnamate 4-hydroxylase with AtATR2",
            {
                "trans_cinnamate_c": -1,
                "nadph_c": -1,
                "o2_c": -1,
                "h_c": -1,
                "p_coumarate_c": 1,
                "nadp_c": 1,
                "h2o_c": 1,
            },
            0,
            5,
            "C4H_AtATR2",
            0.50,
        ),
        HeterologousReaction(
            "LUT_4CL_Morus",
            "4-coumarate CoA ligase",
            {
                "p_coumarate_c": -1,
                "coa_c": -1,
                "atp_c": -1,
                "p_coumaroyl_coa_c": 1,
                "amp_c": 1,
                "ppi_c": 1,
            },
            0,
            7,
            "4CL",
            0.70,
        ),
        HeterologousReaction(
            "LUT_CHS_Carya",
            "chalcone synthase",
            {
                "p_coumaroyl_coa_c": -1,
                "malcoa_c": -3,
                "naringenin_chalcone_c": 1,
                "coa_c": 4,
                "co2_c": 3,
            },
            0,
            10,
            "CHS",
            1.00,
        ),
        HeterologousReaction(
            "LUT_CHI_virus",
            "chalcone isomerase",
            {"naringenin_chalcone_c": -1, "naringenin_c": 1},
            0,
            8,
            "CHI",
            0.80,
        ),
    ]

    if fns_type == "p450":
        reactions.append(
            HeterologousReaction(
                "LUT_FNS_Sm",
                "P450 flavone synthase II",
                {
                    "naringenin_c": -1,
                    "nadph_c": -1,
                    "o2_c": -1,
                    "h_c": -1,
                    "apigenin_c": 1,
                    "nadp_c": 1,
                    "h2o_c": 2,
                },
                0,
                10,
                "FNS",
                1.00,
            )
        )
    else:
        reactions.append(
            HeterologousReaction(
                "LUT_FNS_Sm",
                "2-oxoglutarate-dependent flavone synthase I",
                {
                    "naringenin_c": -1,
                    "akg_c": -1,
                    "o2_c": -1,
                    "apigenin_c": 1,
                    "succ_c": 1,
                    "co2_c": 1,
                    "h2o_c": 1,
                },
                0,
                10,
                "FNS",
                1.00,
            )
        )

    reactions.extend(
        [
            HeterologousReaction(
                "LUT_F3pH_ATR2_main",
                "flavonoid 3-prime hydroxylase with AtATR2, apigenin to luteolin",
                {
                    "apigenin_c": -1,
                    "nadph_c": -1,
                    "o2_c": -1,
                    "h_c": -1,
                    "luteolin_c": 1,
                    "nadp_c": 1,
                    "h2o_c": 1,
                },
                0,
                4,
                "F3pH_AtATR2",
                0.40,
            ),
            HeterologousReaction(
                "LUT_LUTEOLINt",
                "luteolin transport to extracellular space",
                {"luteolin_c": -1, "luteolin_e": 1},
                0,
                1000,
                "",
                None,
            ),
            HeterologousReaction(
                "EX_luteolin_e",
                "luteolin exchange",
                {"luteolin_e": -1},
                0,
                1000,
                "",
                None,
            ),
        ]
    )

    if include_side:
        reactions.extend(
            [
                HeterologousReaction(
                    "LUT_F3pH_ATR2_side",
                    "flavonoid 3-prime hydroxylase side reaction, naringenin to eriodictyol",
                    {
                        "naringenin_c": -1,
                        "nadph_c": -1,
                        "o2_c": -1,
                        "h_c": -1,
                        "eriodictyol_c": 1,
                        "nadp_c": 1,
                        "h2o_c": 1,
                    },
                    0,
                    0.8,
                    "F3pH_AtATR2",
                    0.08,
                ),
                HeterologousReaction(
                    "LUT_ERIODICTYOLt",
                    "eriodictyol transport to extracellular space",
                    {"eriodictyol_c": -1, "eriodictyol_e": 1},
                    0,
                    1000,
                    "",
                    None,
                ),
                HeterologousReaction(
                    "EX_eriodictyol_e",
                    "eriodictyol exchange",
                    {"eriodictyol_e": -1},
                    0,
                    1000,
                    "",
                    None,
                ),
            ]
        )

    return reactions


def add_luteolin_module(model: cobra.Model, fns_type: str, include_side: bool) -> None:
    met_map = build_metabolite_map(model, include_side)
    for reaction in luteolin_reactions(fns_type, include_side):
        add_reaction(model, met_map, reaction)


def detect_growth_reaction(model: cobra.Model, explicit_id: str | None = None) -> str:
    if explicit_id:
        if explicit_id not in model.reactions:
            raise KeyError(f"Growth reaction {explicit_id!r} not found in model.")
        return explicit_id

    coefficients = model.objective.expression.as_coefficients_dict()
    candidate_ids: list[str] = []
    for var, coefficient in coefficients.items():
        if abs(float(coefficient)) < EPS:
            continue
        name = getattr(var, "name", "")
        if name.endswith("_reverse_"):
            continue
        if name in model.reactions:
            candidate_ids.append(name)

    # yeastGEM v9.0.2 uses r_2111 as the active growth objective.
    if "r_2111" in candidate_ids:
        return "r_2111"
    if len(candidate_ids) == 1:
        return candidate_ids[0]
    if "r_2111" in model.reactions:
        return "r_2111"
    raise RuntimeError(
        "Could not detect a single growth objective. Use --growth-reaction."
    )


def set_growth_objective(model: cobra.Model, growth_reaction_id: str) -> None:
    model.objective = model.reactions.get_by_id(growth_reaction_id)
    model.objective_direction = "max"


def set_exchange_uptake(model: cobra.Model, reaction_id: str, uptake: float) -> None:
    rxn = model.reactions.get_by_id(reaction_id)
    rxn.lower_bound = -abs(uptake)
    if rxn.upper_bound < 0:
        rxn.upper_bound = 1000


def configure_environment(
    model: cobra.Model,
    glucose_reaction: str,
    glucose_uptake: float,
    oxygen_reaction: str,
    oxygen_uptake: float,
) -> None:
    set_exchange_uptake(model, glucose_reaction, glucose_uptake)
    set_exchange_uptake(model, oxygen_reaction, oxygen_uptake)


def scale_upper_bound(model: cobra.Model, reaction_id: str, factor: float) -> None:
    if reaction_id not in model.reactions:
        return
    rxn = model.reactions.get_by_id(reaction_id)
    if math.isfinite(rxn.upper_bound):
        rxn.upper_bound *= factor


def apply_capacity_modifiers(model: cobra.Model) -> None:
    """Apply the capacity modifiers suggested by the project notes.

    In unconstrained yeastGEM, most of these bounds are 1000 and therefore not
    enzyme-capacity estimates. Scaling them is retained as an optional scenario,
    but it should not be over-interpreted unless bounds have been curated.
    """

    for rid in ["r_0042"]:  # DAHP synthase, ARO3/ARO4
        scale_upper_bound(model, rid, 1.5)
    for rid in ["r_1049", "r_1050"]:  # TKL1
        scale_upper_bound(model, rid, 2.0)
    for rid in ["r_1048"]:  # TAL1
        scale_upper_bound(model, rid, 2.0)
    for rid in ["r_0962"]:  # PYK1/PYK2 pyruvate kinase
        scale_upper_bound(model, rid, 0.6)
    for rid in ["r_0884"]:  # PCK1
        scale_upper_bound(model, rid, 1.5)
    for rid in ["r_0109"]:  # cytosolic ACC1 complex
        scale_upper_bound(model, rid, 1.5)


def apply_knockouts(model: cobra.Model, knockout_mode: str) -> None:
    if knockout_mode == "none":
        return
    if knockout_mode in {"pha2", "both"} and "r_0938" in model.reactions:
        model.reactions.get_by_id("r_0938").bounds = (0, 0)
    if knockout_mode in {"aro7", "both"} and "r_0278" in model.reactions:
        model.reactions.get_by_id("r_0278").bounds = (0, 0)


def add_expression_atp_cost(model: cobra.Model, atp_cost: float) -> float:
    heterologous_expression = sum(
        reaction.relative_expression or 0.0 for reaction in luteolin_reactions("p450", True)
    )
    total_cost = heterologous_expression * atp_cost
    if total_cost <= 0:
        return 0.0

    met_map = build_metabolite_map(model, include_side=True)
    reaction = HeterologousReaction(
        "LUT_HET_EXPR_ATP_COST",
        "simplified ATP cost for heterologous expression",
        {"atp_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1, "h_c": 1},
        total_cost,
        total_cost,
        "",
        None,
    )
    add_reaction(model, met_map, reaction)
    return total_cost


def flux_or_nan(solution: cobra.Solution, reaction_id: str) -> float:
    if solution.status != "optimal" or reaction_id not in solution.fluxes.index:
        return math.nan
    return float(solution.fluxes[reaction_id])


def solve_growth_scenario(
    model: cobra.Model,
    scenario: str,
    growth_reaction_id: str,
    wt_growth: float,
    luteolin_min: float | None = None,
) -> dict[str, float | str]:
    set_growth_objective(model, growth_reaction_id)
    if luteolin_min is not None and "EX_luteolin_e" in model.reactions:
        model.reactions.get_by_id("EX_luteolin_e").lower_bound = luteolin_min

    solution = model.optimize()
    growth = flux_or_nan(solution, growth_reaction_id)
    luteolin = flux_or_nan(solution, "EX_luteolin_e")
    chs = flux_or_nan(solution, "LUT_CHS_Carya")
    c4h = flux_or_nan(solution, "LUT_C4H_ATR2")
    fns = flux_or_nan(solution, "LUT_FNS_Sm")
    f3ph = flux_or_nan(solution, "LUT_F3pH_ATR2_main")
    side = flux_or_nan(solution, "LUT_F3pH_ATR2_side")
    if math.isnan(side):
        side = 0.0

    burden = math.nan if math.isnan(growth) else 100 * (1 - growth / wt_growth)
    return {
        "scenario": scenario,
        "status": solution.status,
        "objective": "maximize_growth",
        "growth_reaction": growth_reaction_id,
        "growth": growth,
        "growth_fraction_of_wt": math.nan if math.isnan(growth) else growth / wt_growth,
        "burden_percent_vs_wt": burden,
        "luteolin_min_constraint": luteolin_min if luteolin_min is not None else 0.0,
        "growth_floor_fraction": 0.0,
        "luteolin_exchange_flux": luteolin,
        "malonyl_CoA_demand_from_CHS": math.nan if math.isnan(chs) else 3 * chs,
        "nadph_demand_c4h_fns_f3ph": (
            math.nan
            if any(math.isnan(v) for v in [c4h, fns, f3ph])
            else c4h + fns + f3ph + side
        ),
        "f3ph_side_flux": side,
    }


def solve_luteolin_capacity(
    model: cobra.Model,
    growth_reaction_id: str,
    wt_growth: float,
    growth_floor_fraction: float,
) -> dict[str, float | str]:
    growth_rxn = model.reactions.get_by_id(growth_reaction_id)
    growth_rxn.lower_bound = max(growth_rxn.lower_bound, wt_growth * growth_floor_fraction)
    model.objective = model.reactions.get_by_id("EX_luteolin_e")
    model.objective_direction = "max"

    solution = model.optimize()
    growth = flux_or_nan(solution, growth_reaction_id)
    luteolin = flux_or_nan(solution, "EX_luteolin_e")
    return {
        "scenario": f"engineered_luteolin_capacity_growth_floor_{growth_floor_fraction:g}",
        "status": solution.status,
        "objective": "maximize_luteolin_with_growth_floor",
        "growth_reaction": growth_reaction_id,
        "growth": growth,
        "growth_fraction_of_wt": math.nan if math.isnan(growth) else growth / wt_growth,
        "burden_percent_vs_wt": math.nan if math.isnan(growth) else 100 * (1 - growth / wt_growth),
        "luteolin_min_constraint": 0.0,
        "growth_floor_fraction": growth_floor_fraction,
        "luteolin_exchange_flux": luteolin,
        "malonyl_CoA_demand_from_CHS": math.nan,
        "nadph_demand_c4h_fns_f3ph": math.nan,
        "f3ph_side_flux": math.nan,
    }


def load_prepared_model(args: argparse.Namespace) -> cobra.Model:
    model = read_sbml_model(args.model)
    configure_environment(
        model,
        args.glucose_exchange,
        args.glucose_uptake,
        args.oxygen_exchange,
        args.oxygen_uptake,
    )
    return model


def make_engineered_model(args: argparse.Namespace) -> tuple[cobra.Model, float]:
    model = load_prepared_model(args)
    add_luteolin_module(model, args.fns_type, args.include_side_reaction)
    apply_knockouts(model, args.knockout)
    if args.apply_capacity_modifiers:
        apply_capacity_modifiers(model)
    expression_cost = add_expression_atp_cost(model, args.expression_atp_cost)
    return model, expression_cost


def write_csv(rows: Iterable[dict[str, float | str]], output: Path) -> None:
    rows = list(rows)
    if not rows:
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def print_rows(rows: list[dict[str, float | str]]) -> None:
    cols = [
        "scenario",
        "status",
        "growth",
        "growth_fraction_of_wt",
        "burden_percent_vs_wt",
        "luteolin_exchange_flux",
    ]
    widths = {col: max(len(col), 12) for col in cols}
    for row in rows:
        for col in cols:
            widths[col] = max(widths[col], len(format_value(row[col])))

    header = "  ".join(col.ljust(widths[col]) for col in cols)
    print(header)
    print("-" * len(header))
    for row in rows:
        print("  ".join(format_value(row[col]).ljust(widths[col]) for col in cols))


def format_value(value: float | str) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.6g}"
    return str(value)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run growth-maximizing FBA scenarios to compare WT yeastGEM and a "
            "luteolin-pathway engineered strain."
        )
    )
    parser.add_argument("--model", default=DEFAULT_MODEL, help="Input SBML model.")
    parser.add_argument(
        "--output", default=DEFAULT_OUTPUT, help="CSV file for scenario results."
    )
    parser.add_argument(
        "--growth-reaction",
        default=None,
        help="Growth objective reaction. Default: detect from SBML objective.",
    )
    parser.add_argument(
        "--glucose-exchange",
        default="r_1714",
        help="D-glucose exchange reaction in yeastGEM.",
    )
    parser.add_argument(
        "--oxygen-exchange",
        default="r_1992",
        help="Oxygen exchange reaction in yeastGEM.",
    )
    parser.add_argument(
        "--glucose-uptake",
        type=float,
        default=10.0,
        help="Maximum glucose uptake, mmol/gDW/h. Applied as negative lower bound.",
    )
    parser.add_argument(
        "--oxygen-uptake",
        type=float,
        default=1000.0,
        help="Maximum oxygen uptake, mmol/gDW/h. Applied as negative lower bound.",
    )
    parser.add_argument(
        "--fns-type",
        choices=["p450", "2og"],
        default="p450",
        help="Use P450 FNSII stoichiometry or 2-oxoglutarate FNSI stoichiometry.",
    )
    parser.add_argument(
        "--include-side-reaction",
        action="store_true",
        help="Include low-capacity naringenin to eriodictyol side reaction.",
    )
    parser.add_argument(
        "--knockout",
        choices=["none", "pha2", "aro7", "both"],
        default="none",
        help=(
            "Optional native aromatic-pathway knockout. `pha2` is safer for Phe "
            "rerouting; `aro7` also blocks native tyrosine precursor formation."
        ),
    )
    parser.add_argument(
        "--apply-capacity-modifiers",
        action="store_true",
        help="Apply suggested ARO/TKL/TAL/PYK/PCK/ACC bound-scaling factors.",
    )
    parser.add_argument(
        "--expression-atp-cost",
        type=float,
        default=0.0,
        help=(
            "Fixed ATP hydrolysis flux per unit total relative expression. "
            "Default 0 because no measured expression cost is available."
        ),
    )
    parser.add_argument(
        "--luteolin-demands",
        default="0,0.001,0.01,0.05",
        help=(
            "Comma-separated minimum EX_luteolin_e fluxes. Growth is maximized "
            "under each demand. Use `none` to skip."
        ),
    )
    parser.add_argument(
        "--capacity-growth-floors",
        default="0.1,0.3",
        help=(
            "Comma-separated WT-growth fractions used to estimate maximum "
            "luteolin capacity. Use `none` to skip."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    args.model = str(Path(args.model))

    base = load_prepared_model(args)
    growth_reaction_id = detect_growth_reaction(base, args.growth_reaction)
    set_growth_objective(base, growth_reaction_id)
    wt_solution = base.optimize()
    if wt_solution.status != "optimal":
        print(f"WT optimization failed: {wt_solution.status}", file=sys.stderr)
        return 2

    wt_growth = float(wt_solution.fluxes[growth_reaction_id])
    rows: list[dict[str, float | str]] = [
        {
            "scenario": "wild_type",
            "status": wt_solution.status,
            "objective": "maximize_growth",
            "growth_reaction": growth_reaction_id,
            "growth": wt_growth,
            "growth_fraction_of_wt": 1.0,
            "burden_percent_vs_wt": 0.0,
            "luteolin_min_constraint": 0.0,
            "growth_floor_fraction": 0.0,
            "luteolin_exchange_flux": 0.0,
            "malonyl_CoA_demand_from_CHS": 0.0,
            "nadph_demand_c4h_fns_f3ph": 0.0,
            "f3ph_side_flux": 0.0,
        }
    ]

    expression_cost = math.nan
    for demand in parse_float_list(args.luteolin_demands):
        engineered, expression_cost = make_engineered_model(args)
        scenario = "engineered_inserted_idle" if abs(demand) < EPS else (
            f"engineered_luteolin_min_{demand:g}"
        )
        rows.append(
            solve_growth_scenario(
                engineered, scenario, growth_reaction_id, wt_growth, luteolin_min=demand
            )
        )

    for floor in parse_float_list(args.capacity_growth_floors):
        engineered, expression_cost = make_engineered_model(args)
        rows.append(solve_luteolin_capacity(engineered, growth_reaction_id, wt_growth, floor))

    output = Path(args.output)
    write_csv(rows, output)
    print_rows(rows)
    print()
    print(f"Model: {args.model}")
    print(f"Growth reaction: {growth_reaction_id}")
    print(f"Results CSV: {output}")
    if expression_cost and not math.isnan(expression_cost):
        print(f"Fixed heterologous expression ATP cost: {expression_cost:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
