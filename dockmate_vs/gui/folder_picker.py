"""A directory picker with separate navigation and selection actions."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import tkinter as tk
from tkinter import ttk


def subdirectories(directory: Path) -> list[Path]:
    """Return the directory's immediate child folders in display order."""
    with os.scandir(directory) as entries:
        children = [
            Path(entry.path)
            for entry in entries
            if entry.is_dir(follow_symlinks=True)
        ]
    return sorted(children, key=lambda path: (path.name.casefold(), path.name))


class FolderPickerDialog:
    """Modal folder browser that does not conflate opening with accepting."""

    def __init__(self, parent: tk.Misc, title: str, initial_dir: Path) -> None:
        self.parent = parent
        self.result: Optional[Path] = None
        self.current_dir = self._usable_initial_directory(initial_dir)
        self._paths: dict[str, Path] = {}

        self.window = tk.Toplevel(parent)
        self.window.title(title)
        self.window.geometry("760x500")
        self.window.minsize(540, 360)
        self.window.transient(parent)
        self.window.protocol("WM_DELETE_WINDOW", self._cancel)

        self.path_var = tk.StringVar(master=self.window, value=str(self.current_dir))
        self.status_var = tk.StringVar(master=self.window)
        self._build()
        self._show_directory(self.current_dir)

    @staticmethod
    def _usable_initial_directory(initial_dir: Path) -> Path:
        candidate = Path(initial_dir).expanduser()
        if candidate.is_file():
            candidate = candidate.parent
        while not candidate.is_dir() and candidate != candidate.parent:
            candidate = candidate.parent
        return candidate if candidate.is_dir() else Path.home()

    def _build(self) -> None:
        self.window.grid_columnconfigure(0, weight=1)
        self.window.grid_rowconfigure(1, weight=1)

        navigation = tk.Frame(self.window)
        navigation.grid(row=0, column=0, sticky="ew", padx=10, pady=(10, 6))
        navigation.grid_columnconfigure(2, weight=1)

        tk.Button(navigation, text="Home", command=self._go_home, width=8).grid(
            row=0, column=0, padx=(0, 6)
        )
        tk.Button(navigation, text="Up", command=self._go_up, width=8).grid(
            row=0, column=1, padx=(0, 6)
        )
        path_entry = tk.Entry(navigation, textvariable=self.path_var)
        path_entry.grid(row=0, column=2, sticky="ew")
        path_entry.bind("<Return>", self._go_to_typed_path)
        tk.Button(navigation, text="Go", command=self._go_to_typed_path, width=8).grid(
            row=0, column=3, padx=(6, 0)
        )

        browser = tk.Frame(self.window)
        browser.grid(row=1, column=0, sticky="nsew", padx=10)
        browser.grid_columnconfigure(0, weight=1)
        browser.grid_rowconfigure(0, weight=1)

        self.tree = ttk.Treeview(browser, show="tree", selectmode="browse")
        self.tree.heading("#0", text="Subfolders", anchor="w")
        self.tree.column("#0", anchor="w", stretch=True)
        self.tree.grid(row=0, column=0, sticky="nsew")
        self.tree.bind("<Double-1>", self._open_selected)
        self.tree.bind("<Return>", self._open_selected)

        scrollbar = ttk.Scrollbar(browser, orient="vertical", command=self.tree.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=scrollbar.set)

        tk.Label(
            self.window,
            textvariable=self.status_var,
            anchor="w",
            fg="#a12b2b",
        ).grid(row=2, column=0, sticky="ew", padx=10, pady=(5, 0))

        actions = tk.Frame(self.window)
        actions.grid(row=3, column=0, sticky="ew", padx=10, pady=10)
        tk.Button(actions, text="Cancel", command=self._cancel, width=12).pack(
            side="right"
        )
        tk.Button(
            actions,
            text="Select Folder",
            command=self._select_folder,
            width=14,
        ).pack(side="right", padx=(0, 6))
        tk.Button(
            actions,
            text="Open Selected",
            command=self._open_selected,
            width=14,
        ).pack(side="right", padx=(0, 6))

    def _show_directory(self, directory: Path) -> None:
        directory = Path(directory).expanduser()
        try:
            children = subdirectories(directory)
        except OSError as exc:
            self.status_var.set(f"Cannot open {directory}: {exc}")
            return

        self.current_dir = directory
        self.path_var.set(str(directory))
        self.status_var.set("")
        self._paths.clear()
        for item in self.tree.get_children():
            self.tree.delete(item)
        for index, child in enumerate(children):
            item_id = f"folder-{index}"
            self._paths[item_id] = child
            self.tree.insert("", "end", iid=item_id, text=child.name)

    def _selected_path(self) -> Optional[Path]:
        selected = self.tree.selection()
        return self._paths.get(selected[0]) if selected else None

    def _open_selected(self, _event: object = None) -> None:
        selected = self._selected_path()
        if selected is not None:
            self._show_directory(selected)

    def _select_folder(self) -> None:
        self.result = self._selected_path() or self.current_dir
        self.window.destroy()

    def _go_home(self) -> None:
        self._show_directory(Path.home())

    def _go_up(self) -> None:
        self._show_directory(self.current_dir.parent)

    def _go_to_typed_path(self, _event: object = None) -> None:
        candidate = Path(self.path_var.get().strip()).expanduser()
        if candidate.is_dir():
            self._show_directory(candidate)
        else:
            self.status_var.set(f"Folder does not exist: {candidate}")

    def _cancel(self) -> None:
        self.result = None
        self.window.destroy()

    def show(self) -> Optional[Path]:
        self.window.grab_set()
        self.tree.focus_set()
        self.parent.wait_window(self.window)
        return self.result


def choose_directory(
    parent: tk.Misc,
    *,
    title: str,
    initial_dir: Path,
) -> Optional[Path]:
    """Show the navigable directory picker and return the accepted folder."""
    return FolderPickerDialog(parent, title, initial_dir).show()
