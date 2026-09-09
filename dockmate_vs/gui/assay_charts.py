"""Raw-score diagnostics and the four-panel assay screening figure."""

from typing import Iterable, Mapping, Optional

import numpy as np
from matplotlib.figure import Figure
from scipy.stats import spearmanr


ACTIVE = "#177E89"
INACTIVE = "#A6ADB4"
GOLD = "#D97732"


def _number(value) -> Optional[float]:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def docking_diagnostics(records: Iterable[Mapping]) -> dict:
    """Keep raw docking scores separate from any selected rescoring function."""
    points = []
    engines = set()
    for record in records:
        label = _number(record.get("control_label"))
        score = _number(record.get("best_score"))
        completed = record.get("docking_completed")
        if (label not in (0, 1) or score is None
                or str(completed).lower() in ("false", "0", "0.0")):
            continue
        weight = _number(record.get("molecular_weight"))
        points.append([score, int(label), weight if weight and weight > 0 else None])
        engine = record.get("engine")
        if engine:
            engines.add(str(engine).lower())
    paired = [point for point in points if point[2] is not None]
    correlation = None
    if (len(paired) > 1 and len({p[0] for p in paired}) > 1
            and len({p[2] for p in paired}) > 1):
        correlation = float(spearmanr(
            [p[2] for p in paired], [-p[0] for p in paired]
        ).statistic)
    score_label = "Docking score"
    if engines and engines <= {"vina", "smina"}:
        engine_label = next(iter(engines)).capitalize() if len(engines) == 1 else "Vina/Smina"
        score_label = f"{engine_label} docking score (kcal/mol)"
    return {"points": points, "score_label": score_label,
            "spearman_negative_score": correlation}


def screening_figure(summary: dict, score_upper: Optional[float] = None) -> Figure:
    """Render assay curves and raw-score diagnostics without global pyplot state."""
    data = summary["assay_benchmark_charts"]
    diagnostics = data.get("docking_diagnostics") or {}
    figure = Figure(figsize=(10.2, 7.5), dpi=100, constrained_layout=True)
    axes = figure.subplots(2, 2)
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=8)

    axis = axes[0, 0]
    roc = np.asarray(data["roc_curve"])
    axis.plot(roc[:, 0], roc[:, 1], color=ACTIVE, linewidth=2.2)
    axis.plot([0, 1], [0, 1], color="#7B8791", linestyle="--", linewidth=1)
    axis.set(title="A  ROC discrimination", xlabel="False-positive rate",
             ylabel="True-positive rate", xlim=(0, 1), ylim=(0, 1.02))
    auc = summary.get("roc_auc")
    if auc is not None:
        axis.text(.97, .05, f"AUC {auc:.3f}", transform=axis.transAxes,
                  ha="right", va="bottom", fontsize=9)

    axis = axes[0, 1]
    pr = np.asarray(data["precision_recall_curve"])
    axis.plot(pr[:, 0], pr[:, 1], color=GOLD, linewidth=2.2)
    prevalence = data["prevalence"]
    axis.axhline(prevalence, color="#7B8791", linestyle="--", linewidth=1)
    axis.set(title="B  Precision-recall", xlabel="Recall", ylabel="Precision",
             xlim=(0, 1), ylim=(0, 1.02))
    ap = summary.get("average_precision")
    annotation = f"AP {ap:.3f}\n" if ap is not None else ""
    axis.text(.97, .95, annotation + f"random {prevalence:.3f}",
              transform=axis.transAxes, ha="right", va="top", fontsize=9)

    points = diagnostics.get("points") or []
    score_label = diagnostics.get("score_label", "Docking score")
    hist, scatter = axes[1]
    hist.set(title="C  Docking-score distributions",
             xlabel=score_label + "; lower is better", ylabel="Density")
    scatter.set(title="D  Molecular-weight diagnostic",
                xlabel="Molecular weight (Da)", ylabel=score_label)
    if not points:
        for axis in (hist, scatter):
            axis.text(.5, .5, "Raw docking scores unavailable", ha="center",
                      transform=axis.transAxes)
        return figure

    scores = np.asarray([p[0] for p in points])
    lower, upper = float(scores.min()), float(scores.max())
    if score_upper is not None:
        if not np.isfinite(score_upper) or score_upper <= lower:
            raise ValueError(f"Score upper limit must be greater than {lower:g}.")
        upper = score_upper
    if lower == upper:
        lower -= .5
        upper += .5
    edges = np.linspace(lower, upper, 28)
    off_scale = []
    paired_count = 0
    for label, name, color, alpha in (
        (0, "Inactive", INACTIVE, .55), (1, "Active", ACTIVE, .70)
    ):
        group = [p for p in points if p[1] == label]
        visible = [p[0] for p in group if p[0] <= upper]
        if visible:
            hist.hist(visible, bins=edges, density=True, alpha=alpha,
                      color=color, label=f"{name} (n={len(group)})")
        excluded = len(group) - len(visible)
        if excluded:
            off_scale.append(f"{excluded} {name.lower()}")
        paired = [p for p in group if p[2] is not None]
        paired_count += len(paired)
        if paired:
            scatter.scatter([p[2] for p in paired], [p[0] for p in paired],
                            s=20 if label else 12, alpha=.75 if label else .35,
                            color=color, edgecolors="none", label=name)
    hist.set_xlim(lower - .5, upper)
    scatter.set_ylim(lower - .7, upper)
    if hist.get_legend_handles_labels()[0]:
        hist.legend(frameon=False, fontsize=8, loc="upper left")
    if paired_count:
        scatter.legend(frameon=False, fontsize=8, loc="upper right")
    else:
        scatter.text(.5, .5, "Molecular weights unavailable", ha="center",
                     transform=scatter.transAxes)
    if off_scale:
        hist.text(.97, .92, ", ".join(off_scale) + f" > {upper:g} off-scale",
                  transform=hist.transAxes, ha="right", va="top", fontsize=8)
        hidden_paired = sum(p[2] is not None and p[0] > upper for p in points)
        if hidden_paired:
            scatter.text(.03, .97, f"{hidden_paired} point(s) > {upper:g} off-scale",
                         transform=scatter.transAxes, va="top", fontsize=8)
    correlation = diagnostics.get("spearman_negative_score")
    if correlation is not None:
        scatter.text(.97, .05, f"Spearman with -score: {correlation:.2f}",
                     transform=scatter.transAxes, ha="right", va="bottom", fontsize=9)
    if paired_count < len(points):
        scatter.text(.03, .02, f"MW available: {paired_count}/{len(points)}",
                     transform=scatter.transAxes, fontsize=8)
    return figure


def populate_assay_charts(parent, summary: dict) -> None:
    """Embed the publication panels with a score-range control and export toolbar."""
    import tkinter as tk
    from tkinter import ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    parent.grid_rowconfigure(1, weight=1)
    parent.grid_rowconfigure(2, weight=0)
    parent.grid_columnconfigure(0, weight=1)
    parent.grid_columnconfigure(1, weight=0)
    controls = ttk.Frame(parent)
    controls.grid(row=0, column=0, sticky="ew", padx=10, pady=5)
    ttk.Label(controls, text="Score upper limit:").pack(side="left")
    upper_var = tk.StringVar(master=parent)
    entry = ttk.Entry(controls, textvariable=upper_var, width=9)
    entry.pack(side="left", padx=5)
    error = ttk.Label(controls, foreground="#A32929")
    figure = screening_figure(summary)
    # Retina canvas resize requests must not enlarge the surrounding notebook.
    plot_frame = ttk.Frame(parent, width=640, height=480)
    plot_frame.grid(row=1, column=0, sticky="nsew")
    plot_frame.grid_propagate(False)
    plot_frame.grid_rowconfigure(0, weight=1)
    plot_frame.grid_columnconfigure(0, weight=1)
    canvas = FigureCanvasTkAgg(figure, master=plot_frame)
    widget = canvas.get_tk_widget()
    widget.configure(width=640, height=480)
    widget.grid(row=0, column=0, sticky="nsew")
    toolbar_frame = ttk.Frame(parent)
    toolbar_frame.grid(row=2, column=0, sticky="ew")
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame, pack_toolbar=False)
    toolbar.pack(fill="x")

    def redraw(_event=None):
        try:
            upper = float(upper_var.get()) if upper_var.get().strip() else None
            updated = screening_figure(summary, upper)
        except ValueError as exc:
            error.configure(text=str(exc))
            return
        error.configure(text="")
        old = canvas.figure
        updated.set_size_inches(old.get_size_inches(), forward=False)
        updated.set_dpi(old.dpi)
        old.set_canvas(None)
        canvas.figure = updated
        updated.set_canvas(canvas)
        toolbar.update()
        canvas.draw_idle()

    def full_range():
        upper_var.set("")
        redraw()

    ttk.Button(controls, text="Apply", command=redraw).pack(side="left")
    ttk.Button(controls, text="Full range", command=full_range).pack(side="left", padx=5)
    error.pack(side="left", padx=5)
    entry.bind("<Return>", redraw)
    canvas.draw_idle()
