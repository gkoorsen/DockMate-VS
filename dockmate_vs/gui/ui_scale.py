"""UI preferences and typography; importing this module does not create Tk."""

import json
import math
from pathlib import Path


SCALE_CHOICES = (1.0, 1.25, 1.5)
PREFERENCES_PATH = Path.home() / ".dockmate-vs" / "preferences.json"
FONT_ROLES = {"small": (9, False), "body": (10, False), "heading": (12, True), "title": (16, True)}


def font_spec(role, scale=1.0, family="TkDefaultFont"):
    """Return a font specification without consulting a display."""
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("UI scale must be a positive finite number")
    size, bold = FONT_ROLES[role]
    return (family, max(1, round(size * scale)), "bold" if bold else "normal")


def ui_font(role):
    """Named fonts update existing labels and canvas text when scale changes."""
    if role not in FONT_ROLES:
        raise ValueError(f"Unknown font role: {role}")
    return f"DockMate{role.title()}Font"


def load_ui_scale(path=None):
    try:
        value = float(json.loads(Path(path or PREFERENCES_PATH).read_text())["ui_scale"])
        return value if value in SCALE_CHOICES else 1.0
    except (OSError, ValueError, TypeError, KeyError):
        return 1.0


def save_ui_scale(scale, path=None):
    if scale not in SCALE_CHOICES:
        raise ValueError("Unsupported UI scale")
    path = Path(path or PREFERENCES_PATH)
    try:
        data = json.loads(path.read_text())
        if not isinstance(data, dict):
            data = {}
    except (OSError, ValueError):
        data = {}
    data["ui_scale"] = scale
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def widget_scale(widget):
    return widget._root().__dict__.get("ui_scale", 1.0)


def canvas_size(canvas, minimum=(300, 200)):
    """Logical chart dimensions; chart geometry is scaled after each redraw."""
    scale = widget_scale(canvas)
    return (
        max(int(canvas.winfo_width() / scale), minimum[0]),
        max(int(canvas.winfo_height() / scale), minimum[1]),
    )


def scale_figure(figure, scale):
    """Apply UI zoom once to a newly constructed Matplotlib figure."""
    from matplotlib.text import Text

    for text in figure.findobj(Text):
        text.set_fontsize(text.get_fontsize() * scale)


def scale_table(table):
    from tkinter import font

    scale = widget_scale(table)
    columns = ("#0", *table["columns"])
    widths = table.__dict__.setdefault("_unscaled_widths", {
        column: int(table.column(column, "width")) for column in columns
    })
    heading_font = font.Font(root=table, name="TkHeadingFont", exists=True)
    for column, width in widths.items():
        heading_width = heading_font.measure(table.heading(column, "text")) + round(20 * scale)
        table.column(column, width=max(round(width * scale), heading_width))


class UIScale:
    def __init__(self, root):
        from tkinter import font, ttk

        self.root = root
        self.fonts = {
            name: font.Font(root=root, name=name, exists=True)
            for name in root.tk.splitlist(root.tk.call("font", "names")) if name.startswith("Tk")
        }
        self.sizes = {name: abs(item.cget("size")) for name, item in self.fonts.items()}
        family = self.fonts["TkDefaultFont"].cget("family")
        for role in FONT_ROLES:
            name = ui_font(role)
            self.fonts[name] = font.Font(root=root, name=name, font=font_spec(role, family=family))
            self.sizes[name] = font_spec(role)[1]
        self.style = ttk.Style(root)
        self.apply(load_ui_scale())

    def apply(self, scale):
        # Tk already applies the display DPI to point sizes. Apply user zoom
        # here only, leaving native `tk scaling` untouched to avoid double zoom.
        self.root.ui_scale = scale
        for name, item in self.fonts.items():
            item.configure(size=max(1, round(self.sizes[name] * scale)))
        self.style.configure(".", font="TkDefaultFont")
        self.style.configure("Treeview", font="TkDefaultFont", rowheight=self.fonts["TkDefaultFont"].metrics("linespace") + round(8 * scale))
        self.style.configure("Treeview.Heading", font="TkHeadingFont")
        self.style.configure("TNotebook.Tab", padding=(round(8 * scale), round(4 * scale)))
        self._notify(self.root, scale)
        self.root.update_idletasks()

    def _notify(self, widget, scale):
        if "font" in widget.keys():
            # Some Tk 8.6 Aqua widgets repaint a changed named font without
            # recalculating their requested size until the option is reapplied.
            widget.configure(font=widget.cget("font"))
        if widget.winfo_class() == "Treeview":
            scale_table(widget)
        widget.event_generate("<<UIScaleChanged>>")
        for child in widget.winfo_children():
            self._notify(child, scale)
