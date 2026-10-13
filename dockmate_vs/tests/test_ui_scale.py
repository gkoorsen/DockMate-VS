import json
from types import SimpleNamespace

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
