import json

import numpy as np
from matplotlib.backend_bases import MouseEvent
from matplotlib.backends.backend_agg import FigureCanvasAgg

from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.gui.unknown_charts import (
    ScoreHover,
    attach_pose_quality,
    unknown_docking_data,
    unknown_score_figure,
)


def record(name="sample", **overrides):
    row = dict(
        mode="screening", control_label=None, pdb_id="1ABC", ligand_resname="LIG",
        target_name="Target", engine="smina", dock_name=name, best_score=-7.0,
        docking_completed=True, rescore_score=-12.0,
    )
    row.update(overrides)
    return row


def test_selected_unknown_scores_exclude_controls_and_failures_without_truncation():
    rows = [record(f"sample {i}", best_score=-i, rescore_score=-i) for i in range(77)]
    rows += [
        record("active", control_label=1), record("decoy", control_label=0),
        record("native", mode="single"), record("failed", docking_completed=False),
        record("missing", best_score=None, rescore_score=None),
        record("nan", best_score=float("nan"), rescore_score=float("nan")),
        record("infinite", best_score=float("inf"), rescore_score=float("inf")),
    ]
    groups = unknown_docking_data(rows)
    assert len(groups) == 1
    assert groups[0]["cases"] == 81
    assert len(groups[0]["points"]) == 77
    assert groups[0]["points"][0]["score"] == 0
    assert groups[0]["points"][-1]["score"] == -76
    json.dumps(groups, allow_nan=False)


def test_groups_preserve_structures_engines_and_missing_score_rows():
    groups = unknown_docking_data([
        record(), record("other site", ligand_resname="ABC"),
        record("other receptor", pdb_id="2ABC"), record("rdock", engine="rdock"),
        record("unscored receptor", pdb_id="3ABC", best_score=None, rescore_score=None),
        record("csv unknown", control_label=float("nan")),
    ])
    assert len(groups) == 5
    assert sum(len(g["points"]) for g in groups) == 5
    assert groups[-1]["cases"] == 1
    assert groups[-1]["points"] == []


def test_figure_plots_every_score_highlights_minimum_and_keeps_ties_hoverable():
    groups = unknown_docking_data([
        record("worse", best_score=-3, rescore_score=-3),
        record("best", best_score=-9, rescore_score=-9),
        record("tied best", best_score=-9, rescore_score=-9),
        record("missing", pdb_id="2ABC", best_score=None, rescore_score=None),
    ])
    figure, axis, points = unknown_score_figure(groups)
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    assert len(points) == 3
    assert len({(p["x"], p["y"]) for p in points}) == 3
    assert axis.xaxis_inverted()
    assert np.asarray(axis.collections[1].get_offsets())[:, 0].tolist() == [-9, -9]
    assert any("No scored unknowns" in t.get_text() for t in axis.texts)
    assert any("Scored 3/3" in t.get_text() for t in axis.texts)
    assert np.asarray(canvas.buffer_rgba())[:, :, :3].std() > 5

    hover = ScoreHover(canvas, axis, points)
    for point in points:
        x, y = axis.transData.transform((point["x"], point["y"]))
        event = MouseEvent("motion_notify_event", canvas, x, y)
        canvas.callbacks.process("motion_notify_event", event)
        assert hover.annotation.get_visible()
        assert point["compound"] in hover.annotation.get_text()
        assert f"{point['score']:.4f}" in hover.annotation.get_text()
        assert "1ABC/LIG" in hover.annotation.get_text()
    canvas.callbacks.process("motion_notify_event", MouseEvent("motion_notify_event", canvas, 1, 1))
    assert not hover.annotation.get_visible()
    hover.disconnect()


def test_csv_only_saved_results_receive_unknown_chart_data(tmp_path):
    import pandas as pd

    path = tmp_path / "redock_results.csv"
    pd.DataFrame([record()]).to_csv(path, index=False)
    summary_path = tmp_path / "redock_summary.json"
    summary_path.write_text(json.dumps({"n_samples": 1}))
    app = object.__new__(DockMateVSApp)
    summary = app._summary_for_display(path, summary_path)
    assert summary["unknown_docking_scores"][0]["points"][0]["score"] == -12


def test_hover_prefers_visible_highlight_over_nearly_overlapping_gray_point():
    groups = unknown_docking_data([
        record("highlighted best", best_score=-8.744, rescore_score=-8.744),
        record("almost best", best_score=-8.7439, rescore_score=-8.7439),
        record("other", best_score=-3, rescore_score=-3),
    ])
    figure, axis, points = unknown_score_figure(groups)
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    hover = ScoreHover(canvas, axis, points)
    x, y = axis.transData.transform((-8.744, 0))
    hover.on_motion(MouseEvent("motion_notify_event", canvas, x, y))
    assert "highlighted best" in hover.annotation.get_text()
    hover.disconnect()


def test_unknown_hover_shows_posebusters_and_native_contact_recovery():
    groups = unknown_docking_data([
        record("analyzed", case_id="case-1", rescore_score=-9.0),
        record("not analyzed", case_id="case-2", rescore_score=-7.0),
    ])
    groups = attach_pose_quality(groups, [{
        "target_name": "Target",
        "pdb_id": "1ABC",
        "ligand": "LIG",
        "compound": "analyzed",
        "case_id": "case-1",
        "score": -9.0,
        "posebusters_available": True,
        "posebusters_pass": False,
        "plip_available": True,
        "native_contact_recovery": 0.75,
    }])

    analyzed = next(point for point in groups[0]["points"] if point["compound"] == "analyzed")
    missing = next(point for point in groups[0]["points"] if point["compound"] == "not analyzed")
    assert analyzed["posebusters_pass"] is False
    assert analyzed["native_contact_recovery"] == 0.75
    assert missing["pose_quality_analyzed"] is False

    figure, axis, points = unknown_score_figure(groups)
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    hover = ScoreHover(canvas, axis, points)

    for point, expected in (
        (next(point for point in points if point["compound"] == "analyzed"),
         ("PoseBusters: Fail", "Native contacts recovered: 75.0%")),
        (next(point for point in points if point["compound"] == "not analyzed"),
         ("PoseBusters: Not analyzed", "Native contacts recovered: Not analyzed")),
    ):
        x, y = axis.transData.transform((point["x"], point["y"]))
        hover.on_motion(MouseEvent("motion_notify_event", canvas, x, y))
        assert expected[0] in hover.annotation.get_text()
        assert expected[1] in hover.annotation.get_text()
    hover.disconnect()
