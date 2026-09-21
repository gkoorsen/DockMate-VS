"""Selected docking scores for unknown screening compounds, grouped by structure."""

import math
import textwrap
from collections import defaultdict

import numpy as np
from matplotlib.figure import Figure


COLORS = ("#1876A3", "#C56B23", "#348447", "#C53D4B", "#8660A8", "#428C8C", "#BD599B")


def unknown_docking_data(records) -> list:
    groups = {}
    for record in records:
        label = record.get("control_label")
        unlabelled = label is None or str(label).strip().lower() in ("", "nan")
        if record.get("mode") != "screening" or not unlabelled:
            continue
        key = tuple(str(record.get(field) or "") for field in (
            "target_name", "pdb_id", "ligand_resname", "engine",
        ))
        group = groups.setdefault(key, {
            "target": key[0], "pdb_id": key[1], "ligand": key[2],
            "engine": key[3], "cases": 0, "points": [], "score_source": None,
        })
        group["cases"] += 1
        score_value = None
        score_source = None
        for field, label in (
            ("rescore_cnn_affinity", "GNINA CNN affinity"),
            ("rescore_cnn_score", "GNINA CNN score"),
            ("rescore_score", "Vinardo score-only"),
            ("best_score", "docking score"),
        ):
            if record.get(field) is not None:
                score_value, score_source = record.get(field), label
                break
        try:
            score = float(score_value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(score) or str(record.get("docking_completed")).lower() in ("false", "0", "0.0"):
            continue
        group["points"].append({
            "compound": str(record.get("dock_name") or record.get("case_id") or "Unnamed"),
            "score": score, "score_source": score_source,
            "case_id": record.get("case_id"),
        })
        group["score_source"] = score_source
    return [groups[key] for key in sorted(groups)]


def attach_pose_quality(groups, quality_rows) -> list:
    """Attach top-hit PoseBusters and PLIP results to matching chart points."""
    quality_rows = list(quality_rows or [])

    def _text(value) -> str:
        return str(value or "").strip()

    def _score(value):
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        return round(number, 8) if math.isfinite(number) else None

    by_case = {
        _text(row.get("case_id")): row
        for row in quality_rows
        if _text(row.get("case_id"))
    }
    by_identity = {}
    for row in quality_rows:
        key = (
            _text(row.get("target_name")),
            _text(row.get("pdb_id")),
            _text(row.get("ligand")),
            _text(row.get("compound")),
            _score(row.get("score")),
        )
        by_identity.setdefault(key, []).append(row)

    enriched = []
    for group in groups or []:
        updated_group = dict(group)
        updated_points = []
        for point in group.get("points") or []:
            updated = dict(point)
            case_id = _text(point.get("case_id"))
            quality = by_case.get(case_id) if case_id else None
            if quality is None:
                key = (
                    _text(group.get("target")),
                    _text(group.get("pdb_id")),
                    _text(group.get("ligand")),
                    _text(point.get("compound")),
                    _score(point.get("score")),
                )
                matches = by_identity.get(key) or []
                compatible = [
                    row
                    for row in matches
                    if not case_id
                    or not _text(row.get("case_id"))
                    or _text(row.get("case_id")) == case_id
                ]
                quality = compatible[0] if compatible else None
            updated["pose_quality_analyzed"] = quality is not None
            updated["posebusters_available"] = (
                quality.get("posebusters_available") if quality is not None else None
            )
            updated["posebusters_pass"] = (
                quality.get("posebusters_pass") if quality is not None else None
            )
            updated["native_contact_recovery"] = (
                quality.get("native_contact_recovery") if quality is not None else None
            )
            updated["plip_available"] = (
                quality.get("plip_available") if quality is not None else None
            )
            updated_points.append(updated)
        updated_group["points"] = updated_points
        enriched.append(updated_group)
    return enriched


def structure_label(group):
    return f"{group['target'] or group['pdb_id']} | {group['pdb_id']}/{group['ligand']} | {group['engine'] or 'unspecified engine'}"


def unknown_score_figure(groups):
    """Build an exportable strip plot and a point index for hover hit testing."""
    figure = Figure(figsize=(10.5, max(4.2, 1.5 + .75 * len(groups))), dpi=100)
    axis = figure.add_subplot()
    figure.subplots_adjust(left=.20, right=.70, top=.89, bottom=.16)
    axis.spines[["top", "right"]].set_visible(False)
    sources = {g["score_source"] for g in groups if g.get("score_source")}
    source_label = next(iter(sources), "Selected score") if len(sources) == 1 else "Selected score"
    axis.set_title(f"Unknown-compound {source_label.lower()}s", loc="left", fontsize=12, pad=18)
    engines = {g["engine"].lower() for g in groups}
    units = " (kcal/mol)" if engines and engines <= {"smina", "vina"} else ""
    axis.set_xlabel(f"{source_label}{units}; lower is better", fontsize=9)
    axis.tick_params(labelsize=8)
    axis.grid(axis="x", color="#E7E7E7", linewidth=.7)
    axis.set_axisbelow(True)
    points = []
    labels = []
    for y, group in enumerate(groups):
        label = f"{group['pdb_id']}/{group['ligand']}"
        if group["target"] and group["target"] != group["pdb_id"]:
            label = textwrap.shorten(group["target"], width=23, placeholder="...") + "\n" + label
        if len(engines) > 1:
            label += f" ({group['engine'] or 'unspecified'})"
        labels.append(label)
        ordered = sorted(group["points"], key=lambda p: (p["score"], p["compound"]))
        color = COLORS[y % len(COLORS)]
        tied = defaultdict(list)
        for point in ordered:
            tied[point["score"]].append(point)
        for score, members in tied.items():
            offsets = np.linspace(-.16, .16, len(members)) if len(members) > 1 else [0.0]
            for point, offset in zip(members, offsets):
                points.append(dict(point, x=score, y=y + float(offset), group=group,
                                   best=score == ordered[0]["score"]))
        group_points = points[-len(ordered):] if ordered else []
        axis.scatter([p["x"] for p in group_points], [p["y"] for p in group_points],
                     s=25, color="#90979B", alpha=.48, edgecolors="none")
        if ordered:
            best = ordered[0]
            winners = [p for p in group_points if p["score"] == best["score"]]
            axis.scatter([p["x"] for p in winners], [p["y"] for p in winners],
                         s=85, color=color, edgecolors="white", linewidth=.8, zorder=3)
            name = textwrap.shorten(best["compound"], width=48, placeholder="...")
            text = "\n".join(textwrap.wrap(name, width=27))
            text += f"\n{best['score']:.3f}"
            if len(winners) > 1:
                text += f" ({len(winners)} tied)"
        else:
            text = "No scored unknowns"
        text += f"\nScored {len(ordered)}/{group['cases']}"
        axis.text(1.04, y, text, transform=axis.get_yaxis_transform(),
                  va="center", ha="left", fontsize=8, color=color, clip_on=False)
    axis.set_yticks(range(len(groups)), labels)
    axis.set_ylim(-.65, max(len(groups) - .35, .65))
    if points:
        low, high = min(p["x"] for p in points), max(p["x"] for p in points)
        pad = max((high - low) * .08, .5)
        axis.set_xlim(high + pad, low - pad)
    else:
        axis.set_xlim(1, -1)
    return figure, axis, points


class ScoreHover:
    """Find the nearest visible point in screen pixels, including tied scores."""

    def __init__(self, canvas, axis, points):
        self.canvas, self.axis, self.points = canvas, axis, points
        self.annotation = axis.annotate(
            "", xy=(0, 0), xytext=(12, 12), textcoords="offset points",
            bbox={"boxstyle": "round,pad=.5", "fc": "white", "ec": "#666666"},
            fontsize=9, zorder=10, annotation_clip=False,
        )
        self.annotation.set_visible(False)
        self.connection = canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.leave_connection = canvas.mpl_connect("figure_leave_event", self.on_motion)

    def disconnect(self):
        self.canvas.mpl_disconnect(self.connection)
        self.canvas.mpl_disconnect(self.leave_connection)

    def on_motion(self, event):
        selected = None
        if event.inaxes is self.axis and self.points and event.name != "figure_leave_event":
            pixels = self.axis.transData.transform([(p["x"], p["y"]) for p in self.points])
            distance = np.hypot(pixels[:, 0] - event.x, pixels[:, 1] - event.y)
            index = int(np.argmin(distance))
            scale = self.canvas.figure.dpi / 72
            # Highlighted markers sit above nearby gray points, so hit-test them first.
            highlighted = [i for i, p in enumerate(self.points) if p["best"] and distance[i] <= 5 * scale]
            if highlighted:
                index = min(highlighted, key=lambda i: distance[i])
            if distance[index] <= 6 * scale:
                selected = self.points[index]
        if selected is None:
            if self.annotation.get_visible():
                self.annotation.set_visible(False)
                self.canvas.draw_idle()
            return
        group = selected["group"]
        name = "\n".join(textwrap.wrap(selected["compound"], width=40))
        analyzed = bool(selected.get("pose_quality_analyzed"))
        pass_value = selected.get("posebusters_pass")
        if not analyzed:
            posebusters = "Not analyzed"
        elif pass_value is True or str(pass_value).strip().lower() == "true":
            posebusters = "Pass"
        elif pass_value is False or str(pass_value).strip().lower() == "false":
            posebusters = "Fail"
        else:
            posebusters = "Unavailable"
        recovery_value = selected.get("native_contact_recovery")
        try:
            recovery = float(recovery_value)
        except (TypeError, ValueError):
            recovery = None
        if recovery is not None and math.isfinite(recovery):
            contacts = f"{100.0 * recovery:.1f}%"
        else:
            contacts = "Unavailable" if analyzed else "Not analyzed"
        text = (f"{name}\n{group['pdb_id']}/{group['ligand']}"
                f"\nTarget: {group['target'] or 'unspecified'}"
                f"\n{selected['score_source']}: {selected['score']:.4f}"
                f"\nEngine: {group['engine'] or 'unspecified'}"
                f"\nPoseBusters: {posebusters}"
                f"\nNative contacts recovered: {contacts}")
        self.annotation.xy = (selected["x"], selected["y"])
        right = event.x > self.axis.bbox.x0 + self.axis.bbox.width / 2
        above = event.y < self.axis.bbox.y0 + self.axis.bbox.height / 2
        self.annotation.set_position((-12 if right else 12, 12 if above else -12))
        self.annotation.set_ha("right" if right else "left")
        self.annotation.set_va("bottom" if above else "top")
        self.annotation.set_text(text)
        self.annotation.set_visible(True)
        self.canvas.draw_idle()


def populate_unknown_charts(parent, groups):
    import tkinter as tk
    from tkinter import ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    parent.grid_columnconfigure(0, weight=1)
    parent.grid_rowconfigure(1, weight=1)
    controls = ttk.Frame(parent)
    controls.grid(row=0, column=0, sticky="ew", padx=10, pady=6)
    ttk.Label(controls, text="Structure:").pack(side="left", padx=(0, 6))
    choices = ["All structures"] + [structure_label(g) for g in groups]
    selection = ttk.Combobox(controls, values=choices, state="readonly", width=48)
    selection.current(0)
    selection.pack(side="left", fill="x", expand=True)

    viewport = tk.Canvas(parent, highlightthickness=0, width=640, height=480)
    viewport.grid(row=1, column=0, sticky="nsew")
    scrollbar = ttk.Scrollbar(parent, orient="vertical", command=viewport.yview)
    scrollbar.grid(row=1, column=1, sticky="ns")
    viewport.configure(yscrollcommand=scrollbar.set)
    figure, axis, points = unknown_score_figure(groups)
    canvas = FigureCanvasTkAgg(figure, master=viewport)
    widget = canvas.get_tk_widget()
    window = viewport.create_window(0, 0, anchor="nw", window=widget)
    toolbar_frame = ttk.Frame(parent)
    toolbar_frame.grid(row=2, column=0, sticky="ew")
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame, pack_toolbar=False)
    toolbar.pack(fill="x")
    hover = ScoreHover(canvas, axis, points)
    visible = groups

    def resize(_event=None):
        width = max(viewport.winfo_width(), 400)
        height = max(viewport.winfo_height(), 140 + 50 * len(visible), 420)
        viewport.itemconfigure(window, width=width, height=height)
        viewport.configure(scrollregion=(0, 0, width, height))

    def redraw(_event=None):
        nonlocal hover, visible
        index = selection.current()
        visible = groups if index <= 0 else [groups[index - 1]]
        updated, axis, points = unknown_score_figure(visible)
        old = canvas.figure
        updated.set_size_inches(old.get_size_inches(), forward=False)
        updated.set_dpi(old.dpi)
        hover.disconnect()
        old.set_canvas(None)
        canvas.figure = updated
        updated.set_canvas(canvas)
        hover = ScoreHover(canvas, axis, points)
        toolbar.update()
        viewport.yview_moveto(0)
        resize()
        canvas.draw_idle()

    def destroy(event):
        if event.widget is parent:
            hover.disconnect()

    selection.bind("<<ComboboxSelected>>", redraw)
    viewport.bind("<Configure>", resize)
    parent.bind("<Destroy>", destroy, add="+")
    parent.after_idle(resize)
    canvas.draw_idle()
