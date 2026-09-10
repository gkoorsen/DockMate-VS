#!/usr/bin/env python3
"""Regenerate SoftwareX tables and figures from the frozen ESR1 result table."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


SEED = 42
N_BOOTSTRAP = 10_000
ACTIVE = "#0B6E69"
INACTIVE = "#C45131"
INK = "#17324D"
GOLD = "#D49A2A"
GRID = "#D9E0E5"


def bootstrap_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(SEED)
    active_scores = scores[labels == 1]
    inactive_scores = scores[labels == 0]
    auc_values = np.empty(N_BOOTSTRAP)
    ap_values = np.empty(N_BOOTSTRAP)
    for index in range(N_BOOTSTRAP):
        sampled_active = rng.choice(active_scores, size=len(active_scores), replace=True)
        sampled_inactive = rng.choice(inactive_scores, size=len(inactive_scores), replace=True)
        sampled_scores = np.concatenate([sampled_active, sampled_inactive])
        sampled_labels = np.concatenate([
            np.ones(len(sampled_active), dtype=int),
            np.zeros(len(sampled_inactive), dtype=int),
        ])
        auc_values[index] = roc_auc_score(sampled_labels, sampled_scores)
        ap_values[index] = average_precision_score(sampled_labels, sampled_scores)
    return {
        "roc_auc": tuple(np.percentile(auc_values, [2.5, 97.5])),
        "average_precision": tuple(np.percentile(ap_values, [2.5, 97.5])),
    }


def tie_aware_ef(labels: np.ndarray, scores: np.ndarray, percent: float) -> dict[str, float]:
    order = np.argsort(-scores, kind="mergesort")
    ranked_scores = scores[order]
    ranked_labels = labels[order]
    selected = max(1, int(math.ceil(len(scores) * percent / 100.0)))
    cutoff = ranked_scores[selected - 1]
    strict = ranked_labels[ranked_scores > cutoff]
    tied = ranked_labels[ranked_scores == cutoff]
    slots = selected - len(strict)
    strict_actives = int(strict.sum())
    tied_actives = int(tied.sum())
    tied_inactives = len(tied) - tied_actives
    expected = strict_actives + slots * tied_actives / len(tied)
    minimum = strict_actives + max(0, slots - tied_inactives)
    maximum = strict_actives + min(slots, tied_actives)
    baseline = labels.mean()

    def ef(active_count: float) -> float:
        return float((active_count / selected) / baseline)

    return {
        "selected": selected,
        "cutoff": float(cutoff),
        "tie_size": len(tied),
        "expected": ef(expected),
        "minimum": ef(minimum),
        "maximum": ef(maximum),
    }


def metric_row(
    metric: str,
    value: float | int | str,
    ci: tuple[float, float] | None = None,
    note: str = "",
) -> dict[str, object]:
    return {
        "metric": metric,
        "value": value,
        "ci_low": ci[0] if ci else None,
        "ci_high": ci[1] if ci else None,
        "note": note,
    }


def analyze_screening(frame: pd.DataFrame, provenance: dict) -> tuple[pd.DataFrame, dict]:
    complete = frame[frame["docking_completed"].astype(bool)].copy()
    labels = complete["control_label"].astype(int).to_numpy()
    docking_scores = complete["best_score"].astype(float).to_numpy()
    rank_scores = -docking_scores
    confidence = bootstrap_metrics(labels, rank_scores)
    roc_auc = roc_auc_score(labels, rank_scores)
    average_precision = average_precision_score(labels, rank_scores)
    prevalence = float(labels.mean())
    active_scores = docking_scores[labels == 1]
    inactive_scores = docking_scores[labels == 0]
    mann_whitney = mannwhitneyu(
        -active_scores,
        -inactive_scores,
        alternative="greater",
        method="asymptotic",
    )

    ranked = complete.assign(_rank_score=-complete["best_score"]).sort_values(
        ["_rank_score", "dock_name"], ascending=[False, True], kind="mergesort"
    )
    best_active = ranked[ranked["control_label"] == 1].iloc[0]
    best_active_rank = int((ranked["_rank_score"] > best_active["_rank_score"]).sum() + 1)

    ef_metrics = {percent: tie_aware_ef(labels, rank_scores, percent) for percent in (1, 5, 10)}
    rows = [
        metric_row("Planned compounds", len(frame)),
        metric_row("Completed compounds", len(complete)),
        metric_row("Failed compounds", len(frame) - len(complete)),
        metric_row("Assay actives", int(labels.sum())),
        metric_row("Assay inactives", int((labels == 0).sum())),
        metric_row("ROC AUC", roc_auc, confidence["roc_auc"], "Stratified bootstrap, 10,000 replicates"),
        metric_row(
            "Average precision",
            average_precision,
            confidence["average_precision"],
            "Random baseline equals active prevalence",
        ),
        metric_row("Random AP baseline", prevalence),
        metric_row("One-sided Mann-Whitney p-value", float(mann_whitney.pvalue)),
        metric_row("Mean active docking score", float(active_scores.mean())),
        metric_row("Median active docking score", float(np.median(active_scores))),
        metric_row("Mean inactive docking score", float(inactive_scores.mean())),
        metric_row("Median inactive docking score", float(np.median(inactive_scores))),
        metric_row("Best active rank", best_active_rank, note=str(best_active["dock_name"])),
        metric_row("Mean recorded case runtime (s)", float(complete["runtime_sec"].mean())),
        metric_row("Median recorded case runtime (s)", float(complete["runtime_sec"].median())),
        metric_row("Campaign elapsed time (h)", provenance["campaign_timing"]["elapsed_hours"]),
    ]
    for percent, metrics in ef_metrics.items():
        rows.append(
            metric_row(
                f"EF{percent}% tie-aware expected",
                metrics["expected"],
                (metrics["minimum"], metrics["maximum"]),
                f"Cutoff tie: {metrics['tie_size']} compounds; selected: {metrics['selected']}",
            )
        )
    return pd.DataFrame(rows), {
        "labels": labels,
        "docking_scores": docking_scores,
        "rank_scores": rank_scores,
        "ranked": ranked,
        "roc_auc": roc_auc,
        "average_precision": average_precision,
        "prevalence": prevalence,
        "confidence": confidence,
        "ef": ef_metrics,
        "best_active_rank": best_active_rank,
    }


def property_diagnostics(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    properties = ["molecular_weight", "logp", "tpsa", "rotatable_bonds", "ligand_charge"]
    complete = frame[frame["docking_completed"].astype(bool)].copy()
    labels = complete["control_label"].astype(int).to_numpy()
    activity_score = -complete["best_score"].astype(float).to_numpy()
    rows = []
    for column in properties:
        values = complete[column].astype(float).to_numpy()
        correlation, p_value = spearmanr(values, activity_score)
        rows.append({
            "property": column,
            "active_mean": float(values[labels == 1].mean()),
            "inactive_mean": float(values[labels == 0].mean()),
            "score_spearman": float(correlation),
            "score_spearman_p": float(p_value),
        })

    features = complete[properties].astype(float).to_numpy()
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, random_state=SEED),
    )
    folds = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
    auc_values = cross_val_score(model, features, labels, cv=folds, scoring="roc_auc")
    model_summary = {
        "property_only_auc_mean": float(auc_values.mean()),
        "property_only_auc_sd": float(auc_values.std(ddof=1)),
        "fold_auc": auc_values.tolist(),
    }
    return pd.DataFrame(rows), model_summary


def save_screening_figure(output: Path, screening: dict, frame: pd.DataFrame) -> None:
    labels = screening["labels"]
    rank_scores = screening["rank_scores"]
    docking_scores = screening["docking_scores"]
    fpr, tpr, _ = roc_curve(labels, rank_scores)
    precision, recall, _ = precision_recall_curve(labels, rank_scores)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.5), constrained_layout=True)

    axis = axes[0, 0]
    axis.plot(fpr, tpr, color=ACTIVE, linewidth=2.2)
    axis.plot([0, 1], [0, 1], color="#7B8791", linestyle="--", linewidth=1)
    ci = screening["confidence"]["roc_auc"]
    axis.set(title="A  ROC discrimination", xlabel="False-positive rate", ylabel="True-positive rate")
    axis.text(
        0.97,
        0.05,
        f"AUC {screening['roc_auc']:.3f}\n95% CI {ci[0]:.3f}-{ci[1]:.3f}",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        color=INK,
    )

    axis = axes[0, 1]
    axis.plot(recall, precision, color=GOLD, linewidth=2.2)
    axis.axhline(screening["prevalence"], color="#7B8791", linestyle="--", linewidth=1)
    ci = screening["confidence"]["average_precision"]
    axis.set(title="B  Precision-recall", xlabel="Recall", ylabel="Precision")
    axis.text(
        0.97,
        0.95,
        f"AP {screening['average_precision']:.3f}\n95% CI {ci[0]:.3f}-{ci[1]:.3f}\n"
        f"random {screening['prevalence']:.3f}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        color=INK,
    )

    axis = axes[1, 0]
    score_upper = 0.8
    bins = np.linspace(docking_scores.min(), score_upper, 28)
    axis.hist(
        docking_scores[labels == 0], bins=bins, density=True, alpha=0.55,
        color=INACTIVE, label="Assay inactive (n=800)",
    )
    axis.hist(
        docking_scores[labels == 1], bins=bins, density=True, alpha=0.70,
        color=ACTIVE, label="Assay active (n=80)",
    )
    off_scale = int((docking_scores > score_upper).sum())
    axis.set(
        title="C  Docking-score distributions",
        xlabel="Vina docking score (kcal/mol; lower is better)",
        ylabel="Density",
        xlim=(docking_scores.min() - 0.5, score_upper),
    )
    if off_scale:
        axis.text(
            0.97, 0.92, f"{off_scale} inactive score > {score_upper:g} off-scale",
            transform=axis.transAxes, ha="right", va="top", color=INK, fontsize=8,
        )
    axis.legend(frameon=False, fontsize=8)

    axis = axes[1, 1]
    active = frame[frame["control_label"] == 1]
    inactive = frame[frame["control_label"] == 0]
    axis.scatter(
        inactive["molecular_weight"], inactive["best_score"], s=12,
        alpha=0.35, color=INACTIVE, edgecolors="none", label="Inactive",
    )
    axis.scatter(
        active["molecular_weight"], active["best_score"], s=20,
        alpha=0.75, color=ACTIVE, edgecolors="none", label="Active",
    )
    correlation = spearmanr(frame["molecular_weight"], -frame["best_score"]).statistic
    axis.set(
        title="D  Molecular-weight diagnostic",
        xlabel="Molecular weight (Da)",
        ylabel="Vina docking score (kcal/mol)",
        ylim=(docking_scores.min() - 0.7, score_upper),
    )
    axis.text(
        0.97, 0.05, f"Spearman with -score: {correlation:.2f}",
        transform=axis.transAxes, ha="right", va="bottom", color=INK,
    )
    axis.legend(frameon=False, fontsize=8, loc="upper right")

    for axis in axes.flat:
        axis.grid(color=GRID, linewidth=0.6, alpha=0.65)
        axis.set_axisbelow(True)
    figure.suptitle(
        "LIT-PCBA ESR1 antagonist screening: deterministic 80/800 subset",
        color=INK,
        fontsize=14,
        fontweight="bold",
    )
    figure.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def save_protocol_figure(output: Path, seeds: pd.DataFrame) -> None:
    seeds = seeds.sort_values("seed")
    positions = np.arange(len(seeds))
    width = 0.34
    figure, axis = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    axis.bar(
        positions - width / 2,
        seeds["best_rmsd"],
        width,
        color=ACTIVE,
        label="Lowest generated RMSD",
    )
    axis.bar(
        positions + width / 2,
        seeds["top1_rmsd"],
        width,
        color=INACTIVE,
        label="Top-scoring pose RMSD",
    )
    axis.axhline(2.0, color=GOLD, linestyle="--", linewidth=1.4, label="2 A threshold")
    for position, (_, row) in zip(positions, seeds.iterrows()):
        axis.text(
            position - width / 2,
            row["best_rmsd"] + 0.18,
            f"rank {int(row['best_rmsd_rank'])}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=INK,
        )
    axis.set(
        title="ESR1 1XP1/AIH pose recovery across Vina seeds",
        xlabel="Vina random seed",
        ylabel="Heavy-atom RMSD (A)",
        xticks=positions,
        xticklabels=[str(int(seed)) for seed in seeds["seed"]],
    )
    axis.grid(axis="y", color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.15))
    figure.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def write_markdown_tables(
    path: Path,
    metrics: pd.DataFrame,
    properties: pd.DataFrame,
    property_model: dict,
    seeds: pd.DataFrame,
) -> None:
    lookup = metrics.set_index("metric")

    def value(metric: str, digits: int = 3) -> str:
        row = lookup.loc[metric]
        result = f"{float(row['value']):.{digits}f}"
        if pd.notna(row["ci_low"]):
            result += f" ({float(row['ci_low']):.{digits}f}-{float(row['ci_high']):.{digits}f})"
        return result

    lines = [
        "# Frozen ESR1 case-study tables",
        "",
        "## Screening performance",
        "",
        "| Measure | Result |",
        "| --- | ---: |",
        f"| Completed | {int(lookup.loc['Completed compounds', 'value'])}/"
        f"{int(lookup.loc['Planned compounds', 'value'])} |",
        f"| ROC AUC (95% bootstrap CI) | {value('ROC AUC')} |",
        f"| Average precision (95% bootstrap CI) | {value('Average precision')} |",
        f"| Random AP baseline | {value('Random AP baseline')} |",
        f"| Tie-aware EF1% (min-max) | {value('EF1% tie-aware expected', 2)} |",
        f"| Tie-aware EF5% (min-max) | {value('EF5% tie-aware expected', 2)} |",
        f"| Tie-aware EF10% (min-max) | {value('EF10% tie-aware expected', 2)} |",
        f"| Best assay-active rank | {int(lookup.loc['Best active rank', 'value'])} |",
        f"| End-to-end elapsed time | {value('Campaign elapsed time (h)', 1)} h |",
        "",
        "## Pose recovery across seeds",
        "",
        "| Seed | Best RMSD (A) | Top-1 RMSD (A) | Top-5 RMSD (A) | Best-pose rank | Poses |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, row in seeds.sort_values("seed").iterrows():
        lines.append(
            f"| {int(row['seed'])} | {row['best_rmsd']:.2f} | {row['top1_rmsd']:.2f} | "
            f"{row['top5_rmsd']:.2f} | {int(row['best_rmsd_rank'])} | {int(row['pose_count'])} |"
        )
    lines.extend([
        "",
        "All three seeds generated a pose below 2 A, but the Vina top-scoring pose was "
        "not the native-like pose. The native-like pose occurred within the top five.",
        "",
        "## Physicochemical diagnostics",
        "",
        "| Property | Active mean | Inactive mean | Spearman with -score |",
        "| --- | ---: | ---: | ---: |",
    ])
    for _, row in properties.iterrows():
        lines.append(
            f"| {row['property']} | {row['active_mean']:.2f} | "
            f"{row['inactive_mean']:.2f} | {row['score_spearman']:.2f} |"
        )
    lines.extend([
        "",
        f"A property-only logistic model produced a mean 10-fold cross-validated ROC AUC "
        f"of {property_model['property_only_auc_mean']:.3f} "
        f"(SD {property_model['property_only_auc_sd']:.3f}). This diagnostic shows that "
        "the case study demonstrates workflow execution and transparent evaluation, not "
        "receptor-specific scoring superiority.",
        "",
    ])
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frozen",
        type=Path,
        help="Frozen evidence directory",
    )
    parser.add_argument(
        "--manuscript-root",
        type=Path,
        help="Directory containing figures/ and tables/",
    )
    args = parser.parse_args()

    benchmark_root = Path(__file__).resolve().parent
    manuscript_root = args.manuscript_root or benchmark_root.parents[1]
    frozen = args.frozen or benchmark_root / "frozen_esr1_seed42"
    figures = manuscript_root / "figures"
    tables = manuscript_root / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(frozen / "results/redock_results.csv")
    seeds = pd.read_csv(frozen / "protocol_development/seed_results.csv")
    provenance = json.loads((frozen / "provenance/provenance.json").read_text())
    metrics, screening = analyze_screening(frame, provenance)
    properties, property_model = property_diagnostics(frame)

    metrics.to_csv(tables / "esr1_seed42_screening_metrics.csv", index=False)
    properties.to_csv(tables / "esr1_seed42_property_diagnostics.csv", index=False)
    seeds.to_csv(tables / "esr1_seed42_protocol_seeds.csv", index=False)
    screening["ranked"].head(20).drop(columns="_rank_score").to_csv(
        tables / "esr1_seed42_top20.csv", index=False
    )
    (tables / "esr1_seed42_property_model.json").write_text(
        json.dumps(property_model, indent=2) + "\n"
    )
    write_markdown_tables(
        tables / "esr1_case_study_tables.md",
        metrics,
        properties,
        property_model,
        seeds,
    )
    save_screening_figure(figures / "esr1_screening_performance", screening, frame)
    save_protocol_figure(figures / "esr1_protocol_seed_recovery", seeds)
    print(f"Wrote tables to {tables}")
    print(f"Wrote figures to {figures}")


if __name__ == "__main__":
    main()
