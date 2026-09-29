"""Create summary figures from FBA/ecFBA result tables."""

from __future__ import annotations

import pandas as pd

from project_utils import FIGURES_DIR, RESULTS_DIR, ensure_output_dirs
from project_utils import ROOT


def plot_growth_vs_luteolin_demand() -> None:
    import matplotlib.pyplot as plt

    path = RESULTS_DIR / "ecfba_scan.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    df = df[df["status"] == "optimal"].copy()
    if df.empty:
        return

    # Keep one heterologous protein pool and compare P450 membrane pools.
    chosen_hp = 0.02 if (df["heterologous_pool"] == 0.02).any() else df["heterologous_pool"].median()
    subset = df[df["heterologous_pool"] == chosen_hp]

    plt.figure(figsize=(7, 4.5))
    for p450_pool, group in subset.groupby("p450_pool"):
        group = group.sort_values("luteolin_min")
        plt.plot(
            group["luteolin_min"],
            group["growth_fraction_of_wt"],
            marker="o",
            label=f"P450 pool={p450_pool:g}",
        )
    plt.xlabel("Minimum luteolin flux (mmol/gDW/h)")
    plt.ylabel("Growth fraction of WT")
    plt.title(f"ecFBA sensitivity at heterologous pool={chosen_hp:g}")
    plt.ylim(bottom=0)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "growth_vs_luteolin_demand.png", dpi=200)
    plt.close()


def plot_module_burden() -> None:
    import matplotlib.pyplot as plt

    path = RESULTS_DIR / "module_scan.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    if df.empty:
        return

    plt.figure(figsize=(8, 4.5))
    colors = ["#4C78A8" if status == "optimal" else "#B54A4A" for status in df["status"]]
    plt.bar(df["module_id"], df["burden_percent_vs_wt"], color=colors)
    plt.ylabel("Predicted burden vs WT (%)")
    plt.title("Module-by-module burden under light ecFBA")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "module_burden_comparison.png", dpi=200)
    plt.close()


def write_model_experiment_gap() -> None:
    experimental_path = ROOT / "data" / "experimental_results.csv"
    module_path = RESULTS_DIR / "module_scan.csv"
    if not experimental_path.exists() or not module_path.exists():
        return

    exp = pd.read_csv(experimental_path)
    if exp.empty or "mu_max" not in exp.columns or "module" not in exp.columns:
        return

    exp = exp.dropna(subset=["mu_max", "module"]).copy()
    if exp.empty:
        return

    wt = exp.loc[exp["module"].astype(str).str.contains("M0|WT", case=False), "mu_max"]
    if wt.empty:
        return

    wt_mu = float(wt.iloc[0])
    exp["experimental_burden_percent_vs_wt"] = 100 * (1 - exp["mu_max"] / wt_mu)

    model = pd.read_csv(module_path)
    merged = exp.merge(
        model[["module_id", "burden_percent_vs_wt", "growth_fraction_of_wt"]],
        left_on="module",
        right_on="module_id",
        how="left",
    )
    merged["model_minus_experiment_burden_percent"] = (
        merged["burden_percent_vs_wt"]
        - merged["experimental_burden_percent_vs_wt"]
    )
    merged.to_csv(RESULTS_DIR / "model_experiment_gap.csv", index=False)


def main() -> int:
    ensure_output_dirs()
    plot_growth_vs_luteolin_demand()
    plot_module_burden()
    write_model_experiment_gap()
    print(f"saved figures in: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
