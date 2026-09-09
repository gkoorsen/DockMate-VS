"""Scientific contracts for the GUI's publication-style assay panels."""

from dataclasses import asdict
import json

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg

from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.gui.assay_charts import docking_diagnostics, screening_figure
from dockmate_vs.tests.test_screening_pipeline import _app_without_tk, _result


def _results():
    return [
        _result(pdb_id="1XP1", ligand_resname="AIH", engine="vina",
                dock_name=name, control_label=label, best_score=score,
                molecular_weight=weight, docking_completed=True,
                rescore_score=rescore)
        for name, label, score, weight, rescore in (
            ("active_1", 1, -10., 500., -1.),
            ("active_2", 1, -8., 400., -2.),
            ("inactive_1", 0, -9., 450., -3.),
            ("inactive_2", 0, 20., 200., -4.),
        )
    ]


def test_diagnostics_use_raw_scores_and_tie_aware_negative_score_correlation():
    data = docking_diagnostics([
        {"best_score": -10, "control_label": 1, "molecular_weight": 400,
         "rescore_score": 100, "engine": "vina"},
        {"best_score": -8, "control_label": 0, "molecular_weight": 300},
        {"best_score": -8, "control_label": 1, "molecular_weight": 300},
    ])
    assert [p[0] for p in data["points"]] == [-10, -8, -8]
    assert data["spearman_negative_score"] == pytest.approx(1.)
    assert data["score_label"] == "Vina docking score (kcal/mol)"


def test_diagnostics_filter_invalid_values_without_losing_scores_with_missing_mw():
    records = [
        {"best_score": -8, "control_label": 1, "molecular_weight": None},
        {"best_score": "-7", "control_label": "0", "molecular_weight": "nan"},
        {"best_score": -6, "control_label": 0, "molecular_weight": -1},
        {"best_score": np.inf, "control_label": 0},
        {"best_score": "bad", "control_label": 1},
        {"best_score": -9, "control_label": None},
        {"best_score": -9, "control_label": 1, "docking_completed": False},
        {"best_score": -9, "control_label": 1, "docking_completed": "False"},
    ]
    data = docking_diagnostics(records)
    assert data["points"] == [[-8., 1, None], [-7., 0, None], [-6., 0, None]]
    assert data["spearman_negative_score"] is None
    json.dumps(data, allow_nan=False)


def test_summary_rebuild_adds_diagnostics_to_older_json_and_csv(tmp_path):
    app = _app_without_tk()
    results_path = tmp_path / "redock_results.json"
    summary_path = tmp_path / "redock_summary.json"
    results = _results()
    results_path.write_text(json.dumps({"results": [asdict(r) for r in results]}))
    summary = app._summary_for_display(results_path, summary_path)
    expected = summary["assay_benchmark_charts"]["docking_diagnostics"]
    assert len(expected["points"]) == 4
    assert expected["points"][0][0] == -10.
    assert summary["roc_auc"] == pytest.approx(0.)
    assert summary["assay_benchmark_charts"]["roc_curve"] == (
        DockMateVSApp._assay_benchmark_chart_data(
            [(-r.rescore_score, r.control_label) for r in results]
        )["roc_curve"]
    )
    # Legacy CSV plus summary remains loadable without a sibling result JSON.
    results_path.unlink()
    import pandas as pd
    csv_path = tmp_path / "redock_results.csv"
    pd.DataFrame([asdict(r) for r in results]).to_csv(csv_path, index=False)
    del summary["assay_benchmark_charts"]["docking_diagnostics"]
    summary_path.write_text(json.dumps(summary))
    loaded = app._summary_for_display(csv_path, summary_path)
    assert loaded["assay_benchmark_charts"]["docking_diagnostics"] == expected


def test_figure_density_scatter_and_axis_limit_preserve_underlying_data():
    summary = _app_without_tk()._build_summary(_results(), 2.)
    before = json.dumps(summary)
    full = screening_figure(summary)
    assert full.axes[3].get_ylim()[1] == 20.
    figure = screening_figure(summary, score_upper=.8)
    FigureCanvasAgg(figure).draw()
    roc, pr, hist, scatter = figure.axes
    assert [ax.get_title() for ax in figure.axes] == [
        "A  ROC discrimination", "B  Precision-recall",
        "C  Docking-score distributions", "D  Molecular-weight diagnostic",
    ]
    assert hist.get_ylabel() == "Density"
    assert "lower is better" in hist.get_xlabel()
    # Each overlaid class histogram integrates to one over the visible range.
    for container in hist.containers:
        assert sum(p.get_width() * p.get_height() for p in container) == pytest.approx(1.)
    assert sum(len(c.get_offsets()) for c in scatter.collections) == 4
    assert any("1 inactive > 0.8 off-scale" in t.get_text() for t in hist.texts)
    assert any("1 point(s)" in t.get_text() for t in scatter.texts)
    assert scatter.get_ylim()[1] == .8
    assert json.dumps(summary) == before
    assert np.array_equal(roc.lines[0].get_xydata(),
                          np.array(summary["assay_benchmark_charts"]["roc_curve"]))
    assert np.array_equal(pr.lines[0].get_xydata(),
                          np.array(summary["assay_benchmark_charts"]["precision_recall_curve"]))


@pytest.mark.parametrize("limit", [-10., -20., float("nan"), float("inf")])
def test_invalid_axis_limits_are_rejected(limit):
    summary = _app_without_tk()._build_summary(_results(), 2.)
    with pytest.raises(ValueError, match="Score upper limit"):
        screening_figure(summary, limit)


def test_missing_weights_and_constant_scores_render_without_false_correlation():
    results = _results()
    for result in results:
        result.best_score = -8.
        result.molecular_weight = None
    summary = _app_without_tk()._build_summary(results, 2.)
    figure = screening_figure(summary)
    FigureCanvasAgg(figure).draw()
    assert any(t.get_text() == "Molecular weights unavailable"
               for t in figure.axes[3].texts)
    assert not any("Spearman" in t.get_text() for t in figure.axes[3].texts)


def test_multiple_receptors_do_not_get_pooled_assay_panels():
    results = _results()
    results[-1].pdb_id = "2AAA"
    summary = _app_without_tk()._build_summary(results, 2.)
    assert not summary.get("assay_benchmark_charts")
