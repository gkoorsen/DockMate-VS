import json
from types import SimpleNamespace
from unittest.mock import Mock

from dockmate_vs.gui.ui_scale import FONT_ROLES, font_spec, load_ui_scale, save_ui_scale


def test_font_scale_preserves_readable_role_order_without_tk():
    for scale in (1.0, 1.25, 1.5):
        sizes = [font_spec(role, scale)[1] for role in FONT_ROLES]
        assert sizes == sorted(set(sizes))
    for role in FONT_ROLES:
        assert font_spec(role, 1.5)[1] > font_spec(role)[1]
        assert font_spec(role, 1.5)[1] == round(font_spec(role)[1] * 1.5)


def test_preference_roundtrip_preserves_other_settings_and_recovers_invalid_json(tmp_path):
    path = tmp_path / "preferences.json"
    assert load_ui_scale(path) == 1.0
    path.write_text('{"other_setting": true}')
    save_ui_scale(1.25, path)
    assert load_ui_scale(path) == 1.25
    assert json.loads(path.read_text())["other_setting"] is True
    path.write_text("broken")
    assert load_ui_scale(path) == 1.0
    save_ui_scale(1.5, path)
    assert load_ui_scale(path) == 1.5


def test_chart_control_redraw_keeps_scale_without_display():
    from dockmate_vs.gui.app import DockMateVSApp

    pending = []
    draws = []
    canvas = SimpleNamespace(
        _root=lambda: SimpleNamespace(ui_scale=1.5),
        cget=lambda name: 300,
        configure=lambda **kwargs: None,
        bind=lambda *args, **kwargs: None,
        after_idle=lambda callback: pending.append(callback),
        winfo_exists=lambda: True,
        scale=lambda tag, x, y, sx, sy: draws.append((sx, sy)),
        find_all=lambda: [],
    )
    redraw = DockMateVSApp._install_chart(canvas, lambda: draws.append("draw"))
    pending.pop()()
    redraw()
    pending.pop()()
    assert draws == ["draw", (1.5, 1.5), "draw", (1.5, 1.5)]


def test_protocol_charts_install_and_switch_all_cutoffs_without_display(monkeypatch, tmp_path):
    import dockmate_vs.gui.app as app_module

    class Variable:
        def __init__(self, master=None, value=None):
            self.value = value

        def get(self):
            return self.value

        def set(self, value):
            self.value = value

    for name in ("Frame", "Label", "Canvas"):
        monkeypatch.setattr(app_module.tk, name, lambda *args, **kwargs: Mock())
    for name in ("StringVar", "IntVar"):
        monkeypatch.setattr(app_module.tk, name, Variable)
    monkeypatch.setattr(app_module.ttk, "Notebook", lambda *args, **kwargs: Mock())
    controls = []

    def radio(*args, **kwargs):
        controls.append(kwargs)
        return Mock()

    monkeypatch.setattr(app_module.ttk, "Radiobutton", radio)
    path = tmp_path / "protocol_development_results.csv"
    path.write_text(
        "pdb_id,ligand_resname,engine,status,best_rmsd,top1_rmsd,top5_rmsd,top10_rmsd\n"
        "1ABC,LIG,smina,complete,0.5,1.5,0.7,0.5\n"
    )
    app = object.__new__(app_module.DockMateVSApp)
    for name in ("_draw_protocol_scatter", "_draw_factor_effect_chart", "_draw_grouped_horizontal_chart"):
        monkeypatch.setattr(app, name, Mock())
    redraws = []

    def install(canvas, draw):
        redraw = Mock(side_effect=draw)
        redraws.append(redraw)
        redraw()
        return redraw

    monkeypatch.setattr(app, "_install_chart", install)
    parent = Mock()
    app._populate_protocol_charts(parent, path)

    assert len(redraws) == 6
    assert [control["value"] for control in controls] == [1, 5, 10]
    assert parent._protocol_top_n_var.get() == 1
    for control in controls[1:]:
        control["variable"].set(control["value"])
        control["command"]()
        assert f"Top-{control['value']}" in parent._protocol_pose_note_var.get()
    assert [redraw.call_count for redraw in redraws] == [3] * 6
