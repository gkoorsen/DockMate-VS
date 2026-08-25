"""
GUI for redock analysis using single or adaptive docking.
"""

import base64
import copy
import importlib.metadata
import json
import os
import math
import random
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.request
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
import platform
from typing import Callable, Dict, List, Optional, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import pandas as pd
import numpy as np
from loguru import logger
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolAlign, rdFMCS, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D
from docking_platform_gui.utils.redock_results import RedockAnalyzer
from Bio import PDB

from docking_platform_gui.adaptive_docking import AdaptiveDockingPipeline
from docking_platform_gui.binding_site.cocrystal import BindingSite, BindingSiteDefinition
from docking_platform_gui.docking.smina import SminaDockingEngine
from docking_platform_gui.gui.utils import download_pdb_structure
from docking_platform_gui.gui.widgets.progress_dialog import ProgressDialog
from docking_platform_gui.utils.rmsd import calculate_rmsd


COFACTORS = {
    "HEM", "FAD", "NAD", "NAP", "ADP", "ATP", "GDP", "GTP",
    "FMN", "NDP", "ANP", "AMP", "CMP", "UMP",
}

KNOWN_ADDITIVES = {
    "NAG", "MAN", "FUC", "GAL", "GLC", "BMA", "NDG", "FUL",
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PG", "P3G", "PG6",
    "MPD", "DMS", "ACT", "FMT", "MES", "TRS", "EPE",
    "CIT", "SO4", "PO4", "ACE", "NH4", "CA", "MG", "ZN",
    "FE", "CL", "NA", "K",
    "HEM", "FAD", "NAD", "NAP", "ADP", "ATP", "GDP", "GTP",
    "FMN", "NDP", "ANP", "AMP", "CMP", "UMP",
    "MSE", "SEP", "TPO", "PTR", "MLY", "ALY", "CSO", "CSX",
    "HYP", "MLZ", "TYS", "PCA",
    "OLC", "OLA", "PLM", "MYR", "STE", "PC",
    "DAO", "FTT", "TAR", "BOG", "LBN",
    "FLC", "CLR", "IPH", "DTT", "3PE", "P1L",
    "LLP", "BCN", "PEE", "OLB"
}

ADDITIVES_ONLY = KNOWN_ADDITIVES - COFACTORS
FILTERS_PATH = Path.home() / ".docking_platform_gui" / "redock_filters.json"


@dataclass
class RedockResult:
    pdb_id: str
    ligand_resname: str
    ligand_chain: str
    mode: str
    engine: str
    protocol: str
    best_rmsd: float
    success: bool
    runtime_sec: float
    output_file: Optional[str] = None
    pose_count: Optional[int] = None
    best_score: Optional[float] = None
    dock_name: Optional[str] = None
    top1_rmsd: Optional[float] = None
    top5_rmsd: Optional[float] = None
    top10_rmsd: Optional[float] = None
    best_rmsd_rank: Optional[int] = None
    rmsd_best_score: Optional[float] = None
    rmsd_mean: Optional[float] = None
    rmsd_median: Optional[float] = None
    rmsd_std: Optional[float] = None
    near_native_fraction: Optional[float] = None
    score_rmsd_pearson: Optional[float] = None
    score_rmsd_spearman: Optional[float] = None
    protocols_tried: Optional[int] = None
    protocols_success: Optional[int] = None
    protocol_attempts: Optional[List[dict]] = None
    error_message: Optional[str] = None
    control_label: Optional[int] = None
    ligand_charge: Optional[int] = None
    rescore_method: Optional[str] = None
    rescore_score: Optional[float] = None
    rescore_pose_count: Optional[int] = None
    rescore_top1_rmsd: Optional[float] = None
    rescore_top5_rmsd: Optional[float] = None
    rescore_top10_rmsd: Optional[float] = None
    rescore_best_rmsd_rank: Optional[int] = None
    rescore_rmsd_best_score: Optional[float] = None
    rescore_score_rmsd_pearson: Optional[float] = None
    rescore_score_rmsd_spearman: Optional[float] = None
    rescore_cnn_score: Optional[float] = None
    rescore_cnn_affinity: Optional[float] = None
    rescore_error: Optional[str] = None
    site_method: Optional[str] = None
    docking_completed: Optional[bool] = None
    case_id: Optional[str] = None
    variants_prepared: Optional[int] = None
    variants_docked: Optional[int] = None
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    tpsa: Optional[float] = None
    rotatable_bonds: Optional[int] = None
    target_name: Optional[str] = None


class RedockAnalysisApp(tk.Tk):
    """Standalone GUI for redock analysis."""

    def __init__(
        self,
        use_vina_default: bool = False,
        vina_binary_default: Optional[str] = None,
        smina_binary_default: Optional[str] = None
    ):
        super().__init__()

        self.title("Docking Analysis")
        self.geometry("1000x780")
        self.resizable(True, True)

        self.file_var = tk.StringVar()
        self.output_var = tk.StringVar(value=str(Path("output/redock_analysis").resolve()))
        self.mode_var = tk.StringVar(value="protocol_development")
        self.exclude_additives_var = tk.BooleanVar(value=False)
        self.exclude_cofactors_var = tk.BooleanVar(value=False)
        self.sample_enable_var = tk.BooleanVar(value=False)
        self.sample_size_var = tk.StringVar(value="")
        self.sample_seed_var = tk.StringVar(value="")
        self.sample_include_controls_var = tk.BooleanVar(value=True)
        self.use_smiles_var = tk.BooleanVar(value=False)
        self.variant_mode_var = tk.StringVar(value="adaptive") 
        self.max_tautomers_var = tk.StringVar(value="8")
        self.max_conformers_var = tk.StringVar(value="10")
        self.enable_rescore_var = tk.BooleanVar(value=False)
        self.rescore_scoring_var = tk.StringVar(value="vina")

        self.engine_var = tk.StringVar(value="vina")
        self.box_margin_var = tk.StringVar(value="4.0")
        self.apo_site_mode_var = tk.StringVar(value="auto")
        self.site_definition_mode_var = tk.StringVar(value="auto")
        self.site_residues_var = tk.StringVar()
        self.size_x_var = tk.StringVar()
        self.size_y_var = tk.StringVar()
        self.size_z_var = tk.StringVar()
        self.water_handling_var = tk.StringVar(value="remove_all")
        self.exhaustiveness_var = tk.StringVar(value="16")
        self.num_modes_var = tk.StringVar(value="20")
        self.energy_range_var = tk.StringVar(value="3.0")
        self.cpu_var = tk.StringVar(value="4")
        self.seed_var = tk.StringVar(value="42")
        self.timeout_var = tk.StringVar(value="1200")
        self.scoring_var = tk.StringVar(value="vina")
        self.smina_bin_var = tk.StringVar(
            value=smina_binary_default or str(Path.home() / "Documents/apps/smina/smina")
        )
        self.vina_bin_var = tk.StringVar(
            value=vina_binary_default or "/usr/local/bin/vina"
        )

        self.enable_rdock_var = tk.BooleanVar(value=True)
        self.use_vina_var = tk.BooleanVar(value=use_vina_default)
        self.rdock_root_var = tk.StringVar(value=str(Path.home() / "Documents/apps/rdock"))
        self.rdock_runs_var = tk.StringVar(value="20")
        self.rdock_seed_var = tk.StringVar(value="42")
        self.rdock_radius_var = tk.StringVar()

        self.threshold_var = tk.StringVar(value="2.0")
        self.protocol_water_vars = {
            "remove_all": tk.BooleanVar(value=True),
            "retain_all": tk.BooleanVar(value=True),
            "selective": tk.BooleanVar(value=True),
        }
        self.protocol_exhaustiveness_var = tk.StringVar(value="8, 16, 32")
        self.protocol_seeds_var = tk.StringVar(value="42")
        self.protocol_engine_vars = {
            "smina": tk.BooleanVar(value=True),
            "vina": tk.BooleanVar(value=False),
            "rdock": tk.BooleanVar(value=False),
        }
        self.protocol_box_definitions_var = tk.StringVar(value="margin:4.0")
        self.protocol_rescore_methods_var = tk.StringVar(value="none")

        self.pairs_label_var = tk.StringVar(value="Loaded pairs: 0")

        self.progress_dialog = None
        self._queue = queue.Queue()
        self._ui_queue = queue.Queue()
        self._ui_poll_job = None
        self._worker = None
        self.last_results_path: Optional[Path] = None
        self._run_button: Optional[tk.Widget] = None
        self._busy_widgets: List[tk.Widget] = []
        self._busy_widget_states: Dict[tk.Widget, Optional[str]] = {}
        self._pair_count_job = None
        self._pair_count_request_id = 0
        self._variant_all_rmsd_btn: Optional[tk.Radiobutton] = None
        self._rmsd_variant_available = True
        self._network_phase_complete = False
        self._protocol_swept_widgets: List[tk.Widget] = []

        self._load_filter_config()
        self._build_ui()
        self._update_mode()
        self._update_engine()
        self._start_ui_queue()
        self.after(200, self._bring_to_front)

    def _build_ui(self) -> None:
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=0)

        canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)

        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        canvas.configure(yscrollcommand=scrollbar.set)

        container = tk.Frame(canvas, padx=20, pady=20)
        container_id = canvas.create_window((0, 0), window=container, anchor="nw")

        def _on_container_configure(event: tk.Event) -> None:
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_configure(event: tk.Event) -> None:
            canvas.itemconfigure(container_id, width=event.width)

        container.bind("<Configure>", _on_container_configure)
        canvas.bind("<Configure>", _on_canvas_configure)

        def _on_mousewheel(event: tk.Event) -> None:
            delta = int(-1 * (event.delta / 120)) if event.delta else 0
            if delta:
                canvas.yview_scroll(delta, "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        canvas.bind_all("<Button-4>", lambda e: canvas.yview_scroll(-1, "units"))
        canvas.bind_all("<Button-5>", lambda e: canvas.yview_scroll(1, "units"))
        container.grid_columnconfigure(1, weight=1)

        row = 0
        tk.Label(container, text="Input Excel file:").grid(row=row, column=0, sticky="w")
        file_entry = tk.Entry(container, textvariable=self.file_var)
        file_entry.grid(row=row, column=1, sticky="ew", padx=5)
        self._register_busy_widget(file_entry)
        file_btn = tk.Button(container, text="Browse", command=self._safe_call(self._browse_excel))
        file_btn.grid(row=row, column=2, padx=5)
        self._register_busy_widget(file_btn)

        row += 1
        tk.Label(container, textvariable=self.pairs_label_var).grid(row=row, column=1, sticky="w", pady=(5, 10))

        row += 1
        filter_frame = tk.LabelFrame(container, text="Filters", padx=8, pady=8)
        filter_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        add_chk = tk.Checkbutton(
            filter_frame,
            text="Exclude known additives",
            variable=self.exclude_additives_var,
            command=self._update_pair_count
        )
        add_chk.grid(row=0, column=0, sticky="w", padx=(0, 10))
        self._register_busy_widget(add_chk)
        cof_chk = tk.Checkbutton(
            filter_frame,
            text="Exclude cofactors",
            variable=self.exclude_cofactors_var,
            command=self._update_pair_count
        )
        cof_chk.grid(row=0, column=1, sticky="w")
        self._register_busy_widget(cof_chk)
        edit_btn = tk.Button(
            filter_frame,
            text="Edit lists...",
            command=self._safe_call(self._edit_filters)
        )
        edit_btn.grid(row=0, column=2, sticky="e", padx=(20, 0))
        self._register_busy_widget(edit_btn)
        sample_chk = tk.Checkbutton(
            filter_frame,
            text="Random sample",
            variable=self.sample_enable_var,
            command=self._update_pair_count
        )
        sample_chk.grid(row=1, column=0, sticky="w", pady=(8, 0))
        self._register_busy_widget(sample_chk)
        tk.Label(filter_frame, text="Sample size:").grid(row=1, column=1, sticky="e", pady=(8, 0))
        sample_entry = tk.Entry(filter_frame, textvariable=self.sample_size_var, width=8)
        sample_entry.grid(row=1, column=2, sticky="w", pady=(8, 0))
        self._register_busy_widget(sample_entry)
        tk.Label(filter_frame, text="Seed:").grid(row=1, column=3, sticky="e", padx=(20, 0), pady=(8, 0))
        seed_entry = tk.Entry(filter_frame, textvariable=self.sample_seed_var, width=10)
        seed_entry.grid(row=1, column=4, sticky="w", pady=(8, 0))
        self._register_busy_widget(seed_entry)
        include_controls_chk = tk.Checkbutton(
            filter_frame,
            text="Always include controls",
            variable=self.sample_include_controls_var,
            command=self._update_pair_count
        )
        include_controls_chk.grid(row=2, column=0, sticky="w", pady=(8, 0))
        self._register_busy_widget(include_controls_chk)
        use_smiles_chk = tk.Checkbutton(
            filter_frame,
            text="Use SMILES column",
            variable=self.use_smiles_var,
            command=self._update_pair_count
        )
        use_smiles_chk.grid(row=2, column=1, sticky="w", pady=(8, 0))
        self._register_busy_widget(use_smiles_chk)

        row += 1
        tk.Label(container, text="Output directory:").grid(row=row, column=0, sticky="w")
        output_entry = tk.Entry(container, textvariable=self.output_var)
        output_entry.grid(row=row, column=1, sticky="ew", padx=5)
        self._register_busy_widget(output_entry)
        output_btn = tk.Button(container, text="Browse", command=self._safe_call(self._browse_output))
        output_btn.grid(row=row, column=2, padx=5)
        self._register_busy_widget(output_btn)

        row += 1
        self.workflow_notebook = ttk.Notebook(container)
        self.workflow_notebook.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(10, 5))
        self.protocol_tab = tk.Frame(self.workflow_notebook, padx=10, pady=10)
        self.screening_tab = tk.Frame(self.workflow_notebook, padx=10, pady=10)
        self.workflow_notebook.add(self.protocol_tab, text="Protocol Development")
        self.workflow_notebook.add(self.screening_tab, text="Screening")
        self.workflow_notebook.bind("<<NotebookTabChanged>>", self._on_workflow_changed)
        tk.Label(
            self.protocol_tab,
            text="Benchmark protocol combinations on control actives only; completed conditions resume automatically.",
            anchor="w"
        ).grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 8))
        tk.Label(self.protocol_tab, text="Water handling:").grid(row=1, column=0, sticky="w")
        for column, (value, label) in enumerate((
            ("remove_all", "Remove all"),
            ("retain_all", "Retain all"),
            ("selective", "Selective"),
        ), start=1):
            widget = tk.Checkbutton(
                self.protocol_tab, text=label, variable=self.protocol_water_vars[value]
            )
            widget.grid(row=1, column=column, sticky="w", padx=(5, 10))
            self._register_busy_widget(widget)
        tk.Label(self.protocol_tab, text="Exhaustiveness values:").grid(row=2, column=0, sticky="w", pady=(8, 0))
        protocol_exhaustiveness = tk.Entry(
            self.protocol_tab, textvariable=self.protocol_exhaustiveness_var, width=18
        )
        protocol_exhaustiveness.grid(row=2, column=1, columnspan=2, sticky="w", pady=(8, 0))
        self._register_busy_widget(protocol_exhaustiveness)
        tk.Label(self.protocol_tab, text="Seeds:").grid(row=2, column=3, sticky="e", pady=(8, 0))
        protocol_seeds = tk.Entry(self.protocol_tab, textvariable=self.protocol_seeds_var, width=14)
        protocol_seeds.grid(row=2, column=4, sticky="w", padx=5, pady=(8, 0))
        self._register_busy_widget(protocol_seeds)
        tk.Label(self.protocol_tab, text="Docking engines:").grid(
            row=3, column=0, sticky="w", pady=(8, 0)
        )
        for column, engine in enumerate(("smina", "vina", "rdock"), start=1):
            widget = tk.Checkbutton(
                self.protocol_tab, text=engine, variable=self.protocol_engine_vars[engine]
            )
            widget.grid(row=3, column=column, sticky="w", padx=(5, 10), pady=(8, 0))
            self._register_busy_widget(widget)
        tk.Label(self.protocol_tab, text="Box definitions:").grid(
            row=4, column=0, sticky="w", pady=(8, 0)
        )
        protocol_boxes = tk.Entry(
            self.protocol_tab, textvariable=self.protocol_box_definitions_var, width=42
        )
        protocol_boxes.grid(row=4, column=1, columnspan=4, sticky="w", pady=(8, 0))
        self._register_busy_widget(protocol_boxes)
        tk.Label(
            self.protocol_tab,
            text="Separate with ';', e.g. margin:4; margin:6; 20x20x20",
            fg="#666666",
        ).grid(row=5, column=1, columnspan=5, sticky="w")
        tk.Label(self.protocol_tab, text="Rescoring methods:").grid(
            row=6, column=0, sticky="w", pady=(8, 0)
        )
        protocol_rescore = tk.Entry(
            self.protocol_tab, textvariable=self.protocol_rescore_methods_var, width=42
        )
        protocol_rescore.grid(row=6, column=1, columnspan=4, sticky="w", pady=(8, 0))
        self._register_busy_widget(protocol_rescore)
        tk.Label(
            self.protocol_tab,
            text=("Comma-separated: none, vina, vinardo, ad4_scoring, dkoes_fast, dkoes_scoring. "
                  "Every method is compared with the original docking ranking."),
            fg="#666666",
        ).grid(row=7, column=1, columnspan=5, sticky="w")

        tk.Label(
            self.screening_tab,
            text=("Dock controls and screening samples with one validated protocol. "
                  "Controls are used for enrichment; unlabelled samples are excluded from AUC."),
            anchor="w", justify="left", wraplength=850
        ).pack(anchor="w")

        row += 1
        tk.Label(container, text="RMSD threshold (A):").grid(row=row, column=0, sticky="w")
        threshold_entry = tk.Entry(container, textvariable=self.threshold_var, width=10)
        threshold_entry.grid(row=row, column=1, sticky="w", padx=5)
        self._register_busy_widget(threshold_entry)

        row += 1
        variant_frame = tk.LabelFrame(container, text="Ligand variants", padx=8, pady=8)
        variant_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(5, 5))
        variant_frame.grid_columnconfigure(0, weight=1)

        # Adaptive mode (NEW - recommended default)
        variant_adaptive = tk.Radiobutton(
            variant_frame,
            text="⭐ Adaptive (recommended) - Smart variant selection",
            variable=self.variant_mode_var,
            value="adaptive",
            font=("Arial", 10, "bold")
        )
        variant_adaptive.grid(row=0, column=0, sticky="w", pady=(0, 2))
        self._register_busy_widget(variant_adaptive)

        # Add help text for adaptive mode
        adaptive_help = tk.Label(
            variant_frame,
            text="Automatically selects 1-10 variants based on molecular flexibility",
            font=("Arial", 8, "italic"),
            foreground="gray"
        )
        adaptive_help.grid(row=1, column=0, sticky="w", padx=20, pady=(0, 8))

        # Best variant only (renamed for clarity)
        variant_best = tk.Radiobutton(
            variant_frame,
            text="Best variant only (fast) - Dock lowest energy variant",
            variable=self.variant_mode_var,
            value="best"
        )
        variant_best.grid(row=2, column=0, sticky="w", pady=2)
        self._register_busy_widget(variant_best)

        # Thorough sampling (NEW)
        variant_thorough = tk.Radiobutton(
            variant_frame,
            text="Thorough sampling (slow) - Dock 10-15 diverse variants",
            variable=self.variant_mode_var,
            value="thorough"
        )
        variant_thorough.grid(row=3, column=0, sticky="w", pady=2)
        self._register_busy_widget(variant_thorough)

        # Keep best RMSD (for backward compatibility)
        variant_all_rmsd = tk.Radiobutton(
            variant_frame,
            text="Dock all variants, keep best RMSD (very slow)",
            variable=self.variant_mode_var,
            value="all_rmsd"
        )
        variant_all_rmsd.grid(row=4, column=0, sticky="w", pady=2)
        self._variant_all_rmsd_btn = variant_all_rmsd
        self._register_busy_widget(variant_all_rmsd)

        # Keep best score (for backward compatibility)
        variant_all_score = tk.Radiobutton(
            variant_frame,
            text="Dock all variants, keep best score (very slow)",
            variable=self.variant_mode_var,
            value="all_score"
        )
        variant_all_score.grid(row=5, column=0, sticky="w", pady=2)
        self._register_busy_widget(variant_all_score)

        # Configuration inputs (keep these)
        tk.Label(variant_frame, text="Max tautomers:").grid(row=2, column=1, sticky="e", padx=(20, 5))
        max_taut_entry = tk.Entry(variant_frame, textvariable=self.max_tautomers_var, width=6)
        max_taut_entry.grid(row=2, column=2, sticky="w")
        self._register_busy_widget(max_taut_entry)

        tk.Label(variant_frame, text="Max conformers:").grid(row=3, column=1, sticky="e", padx=(20, 5))
        max_conf_entry = tk.Entry(variant_frame, textvariable=self.max_conformers_var, width=6)
        max_conf_entry.grid(row=3, column=2, sticky="w")
        self._register_busy_widget(max_conf_entry)

        row += 1
        self.adaptive_frame = tk.LabelFrame(container, text="Adaptive search settings", padx=10, pady=10)
        self.adaptive_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(10, 5))
        self.adaptive_frame.grid_columnconfigure(1, weight=1)

        adaptive_chk = tk.Checkbutton(
            self.adaptive_frame,
            text="Enable rDock (if installed)",
            variable=self.enable_rdock_var
        )
        adaptive_chk.grid(row=0, column=0, sticky="w")
        self._register_busy_widget(adaptive_chk)
        vina_chk = tk.Checkbutton(
            self.adaptive_frame,
            text="Use Vina (no Smina)",
            variable=self.use_vina_var
        )
        vina_chk.grid(row=0, column=1, sticky="w", padx=(10, 0))
        self._register_busy_widget(vina_chk)
        tk.Label(self.adaptive_frame, text="rDock root:").grid(row=1, column=0, sticky="w", pady=(5, 0))
        adaptive_root = tk.Entry(self.adaptive_frame, textvariable=self.rdock_root_var)
        adaptive_root.grid(row=1, column=1, sticky="ew", padx=5)
        self._register_busy_widget(adaptive_root)
        adaptive_btn = tk.Button(self.adaptive_frame, text="Browse", command=self._safe_call(self._browse_rdock))
        adaptive_btn.grid(row=1, column=2, padx=5)
        self._register_busy_widget(adaptive_btn)

        tk.Label(self.adaptive_frame, text="CPU cores:").grid(
            row=2, column=0, sticky="w", pady=(8, 0)
        )
        adaptive_cpu_entry = tk.Entry(
            self.adaptive_frame, textvariable=self.cpu_var, width=8
        )
        adaptive_cpu_entry.grid(row=2, column=1, sticky="w", padx=5, pady=(8, 0))
        self._register_busy_widget(adaptive_cpu_entry)
        tk.Label(
            self.adaptive_frame,
            text="Used by ligand preparation and each Vina/Smina docking step",
            fg="#666666"
        ).grid(row=2, column=2, sticky="w", pady=(8, 0))

        row += 1
        self.single_frame = tk.LabelFrame(
            container, text="Docking protocol (single protocol / screening)", padx=10, pady=10
        )
        self.single_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(10, 5))
        self.single_frame.grid_columnconfigure(3, weight=1)

        engine_label = tk.Label(self.single_frame, text="Engine:")
        engine_label.grid(row=0, column=0, sticky="w")
        engine_menu = ttk.Combobox(
            self.single_frame,
            textvariable=self.engine_var,
            values=["smina", "vina", "rdock"],
            state="readonly",
            width=12
        )
        engine_menu.grid(row=0, column=1, sticky="w", padx=5)
        engine_menu.bind("<<ComboboxSelected>>", lambda _e: self._update_engine())
        self._register_busy_widget(engine_menu)

        box_margin_label = tk.Label(self.single_frame, text="Box margin (A):")
        box_margin_label.grid(row=0, column=2, sticky="w")
        box_entry = tk.Entry(self.single_frame, textvariable=self.box_margin_var, width=8)
        box_entry.grid(row=0, column=3, sticky="w")
        self._register_busy_widget(box_entry)

        box_size_label = tk.Label(self.single_frame, text="Box size override (x,y,z):")
        box_size_label.grid(row=1, column=0, sticky="w", pady=(5, 0))
        size_x_entry = tk.Entry(self.single_frame, textvariable=self.size_x_var, width=6)
        size_x_entry.grid(row=1, column=1, sticky="w")
        self._register_busy_widget(size_x_entry)
        size_y_entry = tk.Entry(self.single_frame, textvariable=self.size_y_var, width=6)
        size_y_entry.grid(row=1, column=2, sticky="w")
        self._register_busy_widget(size_y_entry)
        size_z_entry = tk.Entry(self.single_frame, textvariable=self.size_z_var, width=6)
        size_z_entry.grid(row=1, column=3, sticky="w")
        self._register_busy_widget(size_z_entry)

        water_label = tk.Label(self.single_frame, text="Water handling:")
        water_label.grid(row=2, column=0, sticky="w", pady=(8, 0))
        water_menu = ttk.Combobox(
            self.single_frame,
            textvariable=self.water_handling_var,
            values=["remove_all", "retain_all", "selective"],
            state="readonly",
            width=12
        )
        water_menu.grid(row=2, column=1, sticky="w")
        self._register_busy_widget(water_menu)

        exhaust_label = tk.Label(self.single_frame, text="Exhaustiveness:")
        exhaust_label.grid(row=3, column=0, sticky="w", pady=(8, 0))
        exhaust_entry = tk.Entry(self.single_frame, textvariable=self.exhaustiveness_var, width=8)
        exhaust_entry.grid(row=3, column=1, sticky="w")
        self._register_busy_widget(exhaust_entry)
        tk.Label(self.single_frame, text="Num modes:").grid(row=3, column=2, sticky="w")
        modes_entry = tk.Entry(self.single_frame, textvariable=self.num_modes_var, width=8)
        modes_entry.grid(row=3, column=3, sticky="w")
        self._register_busy_widget(modes_entry)

        tk.Label(self.single_frame, text="Energy range:").grid(row=4, column=0, sticky="w", pady=(5, 0))
        energy_entry = tk.Entry(self.single_frame, textvariable=self.energy_range_var, width=8)
        energy_entry.grid(row=4, column=1, sticky="w")
        self._register_busy_widget(energy_entry)
        tk.Label(self.single_frame, text="CPU:").grid(row=4, column=2, sticky="w")
        cpu_entry = tk.Entry(self.single_frame, textvariable=self.cpu_var, width=8)
        cpu_entry.grid(row=4, column=3, sticky="w")
        self._register_busy_widget(cpu_entry)

        docking_seed_label = tk.Label(self.single_frame, text="Seed:")
        docking_seed_label.grid(row=5, column=0, sticky="w", pady=(5, 0))
        seed_entry = tk.Entry(self.single_frame, textvariable=self.seed_var, width=8)
        seed_entry.grid(row=5, column=1, sticky="w")
        self._register_busy_widget(seed_entry)
        self._protocol_swept_widgets = [
            engine_label, engine_menu,
            box_margin_label, box_entry,
            box_size_label, size_x_entry, size_y_entry, size_z_entry,
            water_label, water_menu,
            exhaust_label, exhaust_entry,
            docking_seed_label, seed_entry,
        ]
        tk.Label(self.single_frame, text="Timeout (s):").grid(row=5, column=2, sticky="w")
        timeout_entry = tk.Entry(self.single_frame, textvariable=self.timeout_var, width=8)
        timeout_entry.grid(row=5, column=3, sticky="w")
        self._register_busy_widget(timeout_entry)

        tk.Label(self.single_frame, text="Scoring (smina only):").grid(row=6, column=0, sticky="w", pady=(5, 0))
        scoring_menu = ttk.Combobox(
            self.single_frame,
            textvariable=self.scoring_var,
            values=["vina", "vinardo", "ad4_scoring", "dkoes_fast", "dkoes_scoring"],
            state="readonly",
            width=12
        )
        scoring_menu.grid(row=6, column=1, sticky="w")
        self._register_busy_widget(scoring_menu)

        rescore_enable_label = tk.Label(self.single_frame, text="Rescore (smina --score_only):")
        rescore_enable_label.grid(row=6, column=2, sticky="w", pady=(5, 0))
        rescore_chk = tk.Checkbutton(
            self.single_frame,
            text="Enable",
            variable=self.enable_rescore_var
        )
        rescore_chk.grid(row=6, column=3, sticky="w")
        self._register_busy_widget(rescore_chk)

        rescore_method_label = tk.Label(self.single_frame, text="Rescore scoring:")
        rescore_method_label.grid(row=7, column=0, sticky="w", pady=(5, 0))
        rescore_menu = ttk.Combobox(
            self.single_frame,
            textvariable=self.rescore_scoring_var,
            values=["vina", "vinardo", "ad4_scoring", "dkoes_fast", "dkoes_scoring"],
            state="readonly",
            width=12
        )
        rescore_menu.grid(row=7, column=1, sticky="w", padx=5, pady=(5, 0))
        self._register_busy_widget(rescore_menu)
        self._protocol_swept_widgets.extend([
            rescore_enable_label, rescore_chk, rescore_method_label, rescore_menu,
        ])

        tk.Label(self.single_frame, text="Smina binary:").grid(row=8, column=0, sticky="w", pady=(8, 0))
        smina_entry = tk.Entry(self.single_frame, textvariable=self.smina_bin_var)
        smina_entry.grid(row=8, column=1, columnspan=2, sticky="ew", padx=5)
        self._register_busy_widget(smina_entry)
        smina_btn = tk.Button(self.single_frame, text="Browse", command=self._safe_call(self._browse_smina))
        smina_btn.grid(row=8, column=3, padx=5)
        self._register_busy_widget(smina_btn)

        tk.Label(self.single_frame, text="Vina binary:").grid(row=9, column=0, sticky="w", pady=(5, 0))
        vina_entry = tk.Entry(self.single_frame, textvariable=self.vina_bin_var)
        vina_entry.grid(row=9, column=1, columnspan=2, sticky="ew", padx=5)
        self._register_busy_widget(vina_entry)
        vina_btn = tk.Button(self.single_frame, text="Browse", command=self._safe_call(self._browse_vina))
        vina_btn.grid(row=9, column=3, padx=5)
        self._register_busy_widget(vina_btn)

        self.rdock_single_frame = tk.LabelFrame(self.single_frame, text="rDock settings", padx=8, pady=8)
        self.rdock_single_frame.grid(row=10, column=0, columnspan=4, sticky="ew", pady=(10, 0))
        self.rdock_single_frame.grid_columnconfigure(1, weight=1)

        tk.Label(self.rdock_single_frame, text="rDock root:").grid(row=0, column=0, sticky="w")
        rdock_entry = tk.Entry(self.rdock_single_frame, textvariable=self.rdock_root_var)
        rdock_entry.grid(row=0, column=1, sticky="ew", padx=5)
        self._register_busy_widget(rdock_entry)
        rdock_btn = tk.Button(self.rdock_single_frame, text="Browse", command=self._safe_call(self._browse_rdock))
        rdock_btn.grid(row=0, column=2, padx=5)
        self._register_busy_widget(rdock_btn)

        tk.Label(self.rdock_single_frame, text="Runs:").grid(row=1, column=0, sticky="w", pady=(5, 0))
        runs_entry = tk.Entry(self.rdock_single_frame, textvariable=self.rdock_runs_var, width=10)
        runs_entry.grid(row=1, column=1, sticky="w")
        self._register_busy_widget(runs_entry)
        tk.Label(self.rdock_single_frame, text="Seed:").grid(row=1, column=2, sticky="w")
        rdock_seed_entry = tk.Entry(self.rdock_single_frame, textvariable=self.rdock_seed_var, width=10)
        rdock_seed_entry.grid(row=1, column=3, sticky="w")
        self._register_busy_widget(rdock_seed_entry)

        tk.Label(self.rdock_single_frame, text="Radius override (A):").grid(row=2, column=0, sticky="w", pady=(5, 0))
        rdock_radius_entry = tk.Entry(self.rdock_single_frame, textvariable=self.rdock_radius_var, width=10)
        rdock_radius_entry.grid(row=2, column=1, sticky="w")
        self._register_busy_widget(rdock_radius_entry)

        tk.Label(self.single_frame, text="Site definition:").grid(row=11, column=0, sticky="w", pady=(8, 0))
        site_definition_menu = ttk.Combobox(
            self.single_frame,
            textvariable=self.site_definition_mode_var,
            values=["auto", "cocrystal_ligand", "detected_pocket", "specified_residues", "protein_centroid"],
            state="readonly",
            width=18
        )
        site_definition_menu.grid(row=11, column=1, sticky="w", pady=(8, 0))
        self._register_busy_widget(site_definition_menu)
        tk.Label(
            self.single_frame,
            text="Auto uses Ligand column, then residues, then pocket detection",
            fg="#666666"
        ).grid(row=11, column=2, columnspan=2, sticky="w", pady=(8, 0))

        tk.Label(self.single_frame, text="Site residues:").grid(row=12, column=0, sticky="w", pady=(5, 0))
        site_residues_entry = tk.Entry(self.single_frame, textvariable=self.site_residues_var)
        site_residues_entry.grid(row=12, column=1, columnspan=3, sticky="ew", pady=(5, 0))
        self._register_busy_widget(site_residues_entry)

        row += 1
        button_frame = tk.Frame(container)
        button_frame.grid(row=row, column=0, columnspan=3, sticky="e", pady=(15, 0))
        run_btn = tk.Button(button_frame, text="Run Analysis", command=self._safe_call(self._start_run), width=16)
        run_btn.pack(side="right")
        self._register_busy_widget(run_btn)
        self._run_button = run_btn

        row += 1
        self.status_var = tk.StringVar(value="Ready")
        status_label = tk.Label(container, textvariable=self.status_var, anchor="w", fg="#555555")
        status_label.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(10, 0))

        row += 1
        self.results_frame = tk.LabelFrame(container, text="Results", padx=10, pady=10)
        self.results_frame.grid(row=row, column=0, columnspan=3, sticky="nsew", pady=(10, 0))
        self.results_frame.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(row, weight=1)

        results_actions = tk.Frame(self.results_frame)
        results_actions.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        load_results_btn = tk.Button(
            results_actions,
            text="Load Run Folder...",
            command=self._safe_call(self._browse_results_folder)
        )
        load_results_btn.pack(side="left")
        self._register_busy_widget(load_results_btn)
        load_results_file_btn = tk.Button(
            results_actions,
            text="Load Results File...",
            command=self._safe_call(self._browse_results_file)
        )
        load_results_file_btn.pack(side="left", padx=6)
        self._register_busy_widget(load_results_file_btn)
        pose_viewer_btn = tk.Button(
            results_actions,
            text="Pose Viewer",
            command=self._safe_call(self._open_pose_viewer_from_last)
        )
        pose_viewer_btn.pack(side="left")
        self._register_busy_widget(pose_viewer_btn)
        tk.Label(
            results_actions,
            text="Choose one completed campaign folder, not the parent output folder.",
            fg="#666666",
        ).pack(side="left", padx=10)

        self.results_notebook = ttk.Notebook(self.results_frame)
        self.results_notebook.grid(row=1, column=0, sticky="nsew")
        self.results_summary_tab = tk.Frame(self.results_notebook)
        self.results_charts_tab = tk.Frame(self.results_notebook)
        self.results_notebook.add(self.results_summary_tab, text="Summary")
        self.results_notebook.add(self.results_charts_tab, text="Charts")
        self.results_frame.grid_rowconfigure(1, weight=1)
        self._populate_empty_results()

    def _browse_excel(self) -> None:
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            self.file_var.set(path)
            self._update_pair_count()
            self._set_status(f"Loaded: {Path(path).name}")

    def _browse_output(self) -> None:
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            self.output_var.set(path)
            self._set_status("Output directory updated")

    def _browse_rdock(self) -> None:
        path = filedialog.askdirectory(title="Select rDock root directory")
        if path:
            self.rdock_root_var.set(path)
            self._set_status("rDock root updated")

    def _browse_smina(self) -> None:
        path = filedialog.askopenfilename(title="Select smina binary")
        if path:
            self.smina_bin_var.set(path)
            self._set_status("Smina binary updated")

    def _browse_vina(self) -> None:
        path = filedialog.askopenfilename(title="Select vina binary")
        if path:
            self.vina_bin_var.set(path)
            self._set_status("Vina binary updated")

    def _on_workflow_changed(self, _event: Optional[tk.Event] = None) -> None:
        selected = self.workflow_notebook.select()
        if selected == str(self.screening_tab):
            self.mode_var.set("screening")
        else:
            self.mode_var.set("protocol_development")
        self._update_mode()
        self._update_pair_count()

    def _update_mode(self) -> None:
        mode = self.mode_var.get()
        if mode == "adaptive":
            self.adaptive_frame.grid()
            self.single_frame.grid_remove()
        else:
            self.adaptive_frame.grid_remove()
            self.single_frame.grid()
        protocol_mode = mode == "protocol_development"
        for widget in self._protocol_swept_widgets:
            if protocol_mode:
                widget.grid_remove()
            else:
                widget.grid()
        if protocol_mode:
            self.single_frame.configure(text="Base docking protocol (shared by all sweep conditions)")
        else:
            self.single_frame.configure(text="Docking protocol (single protocol / screening)")
        if self._run_button:
            labels = {
                "adaptive": "Run Redock",
                "single": "Run Redock",
                "protocol_development": "Run Protocol Sweep",
                "screening": "Run Screening",
            }
            self._run_button.configure(text=labels.get(mode, "Run Analysis"))

    def _update_engine(self) -> None:
        engine = self.engine_var.get()
        if engine == "rdock":
            self.rdock_single_frame.grid()
        else:
            self.rdock_single_frame.grid_remove()

        if engine == "vina":
            self.scoring_var.set("vina")

    def _update_pair_count(self) -> None:
        self._schedule_pair_count_update()

    def _schedule_pair_count_update(self) -> None:
        if self._pair_count_job:
            self.after_cancel(self._pair_count_job)
        logger.debug("Scheduling pair count update")
        self._pair_count_job = self.after(300, self._run_pair_count_update)

    def _run_pair_count_update(self) -> None:
        self._pair_count_request_id += 1
        request_id = self._pair_count_request_id

        excel_path = Path(self.file_var.get())
        exclude_additives = self.exclude_additives_var.get()
        exclude_cofactors = self.exclude_cofactors_var.get()
        sample_enabled = self.mode_var.get() == "screening" and self.sample_enable_var.get()
        sample_size_raw = self.sample_size_var.get()
        include_controls = self.sample_include_controls_var.get()
        use_smiles = self.use_smiles_var.get() or self.mode_var.get() in (
            "screening", "protocol_development"
        )

        if not excel_path.exists():
            self.pairs_label_var.set("Loaded pairs: 0")
            self._set_rmsd_variant_available(True)
            return

        self.pairs_label_var.set("Loaded pairs: ...")

        def _worker():
            try:
                pairs, _ = self._load_pairs_from_excel(
                    excel_path,
                    exclude_additives=exclude_additives,
                    exclude_cofactors=exclude_cofactors,
                    use_smiles=use_smiles,
                    include_controls=include_controls
                )
                total = len(pairs)
                rmsd_variant_available = any(p.get("site_ligand") for p in pairs)
                controls_count = sum(1 for p in pairs if p.get("control_label") is not None)
                non_controls_count = total - controls_count
                if sample_enabled:
                    sample_size = None
                    sample_valid = True
                    raw = sample_size_raw.strip()
                    if raw:
                        try:
                            sample_size = int(raw)
                            if sample_size <= 0:
                                sample_valid = False
                        except ValueError:
                            sample_valid = False

                    if sample_size and sample_valid:
                        if include_controls:
                            sampled_non_controls = min(sample_size, non_controls_count)
                            label = (
                                f"Loaded pairs: {total} "
                                f"(sample {sampled_non_controls} + controls {controls_count})"
                            )
                        else:
                            sample_size = min(sample_size, total)
                            label = f"Loaded pairs: {total} (sample {sample_size})"
                    elif raw:
                        label = f"Loaded pairs: {total} (sample size invalid)"
                    else:
                        label = f"Loaded pairs: {total}"
                else:
                    label = f"Loaded pairs: {total}"
            except Exception as exc:
                label = "Loaded pairs: 0"
                rmsd_variant_available = True
                logger.debug("Failed to parse Excel: {}", exc)

            def _apply():
                if request_id != self._pair_count_request_id:
                    return
                self.pairs_label_var.set(label)
                self._set_rmsd_variant_available(rmsd_variant_available)

            self._run_on_ui(_apply)

        threading.Thread(target=_worker, daemon=True).start()

    def _start_run(self) -> None:
        if self._worker and self._worker.is_alive():
            self._set_status("Run already active")
            return
        self._set_status("Validating inputs...")
        self.update_idletasks()
        logger.info("Run analysis requested")

        excel_path = Path(self.file_var.get())
        if not excel_path.exists():
            messagebox.showerror("Missing file", "Please choose an Excel file.")
            self._set_status("Missing Excel file")
            return

        output_dir = Path(self.output_var.get()).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        sample_size = None
        sample_seed = None
        include_controls = self.sample_include_controls_var.get()
        sampling_enabled = self.mode_var.get() == "screening" and self.sample_enable_var.get()
        if sampling_enabled:
            sample_size = self._parse_sample_size(silent=False)
            if sample_size is None:
                self._set_status("Invalid sample size")
                return
            sample_seed = self._parse_sample_seed(silent=False)
            if sample_seed is None and self.sample_seed_var.get().strip():
                self._set_status("Invalid sample seed")
                return

        threshold = self._parse_float(self.threshold_var.get(), "RMSD threshold")
        if threshold is None:
            self._set_status("Invalid threshold")
            return

        config = {
            "mode": self.mode_var.get(),
            "threshold": threshold,
            "output_dir": output_dir,
            "input_file": str(excel_path.resolve()),
            "filters": {
                "exclude_additives": self.exclude_additives_var.get(),
                "exclude_cofactors": self.exclude_cofactors_var.get(),
                "use_smiles": self.use_smiles_var.get() or self.mode_var.get() in (
                    "screening", "protocol_development"
                ),
            },
            "sampling": {
                "enabled": sampling_enabled,
                "size": sample_size,
                "seed": sample_seed,
                "include_all_controls": include_controls,
                "strategy": "stratified_by_structure" if include_controls else "global",
            },
            "single": self._collect_single_config(),
            "adaptive": self._collect_adaptive_config(),
            "rescore": {
                "enable": self.enable_rescore_var.get(),
                "scoring": self.rescore_scoring_var.get(),
                "smina_binary": self.smina_bin_var.get()
            }
        }
        if config["single"] is None or config["adaptive"] is None:
            logger.error(
                "Configuration invalid | mode={} single={} adaptive={}",
                config["mode"],
                config["single"],
                config["adaptive"]
            )
            self._set_status("Configuration invalid - check log for details")
            return
        if config["mode"] == "protocol_development":
            try:
                config["protocol_sweep"] = self._collect_protocol_sweep_config()
            except ValueError as exc:
                messagebox.showerror("Protocol settings", str(exc))
                self._set_status("Invalid protocol sweep settings")
                return

        self._set_busy(True, "Parsing Excel...")
        self.update_idletasks()

        def _load_pairs_worker():
            try:
                pairs, col_info = self._load_pairs_from_excel(
                    excel_path,
                    exclude_additives=self.exclude_additives_var.get(),
                    exclude_cofactors=self.exclude_cofactors_var.get(),
                    use_smiles=self.use_smiles_var.get() or config["mode"] in (
                        "screening", "protocol_development"
                    ),
                    include_controls=include_controls
                )
            except Exception as exc:
                message = str(exc)
                self._run_on_ui(lambda msg=message: self._start_run_failed("Excel error", msg))
                return

            if not pairs:
                self._run_on_ui(lambda: self._start_run_failed("No pairs", "No valid PDB/ligand pairs found."))
                return

            if config["mode"] == "protocol_development":
                actives = self._protocol_active_pairs(pairs)
                if not actives:
                    msg = (
                        "Protocol Development requires at least one control active. "
                        "In the current template format, a row with a decoy SMILES "
                        "creates an active/decoy control pair."
                    )
                    self._run_on_ui(lambda m=msg: self._start_run_failed("No control actives", m))
                    return
                self._run_on_ui(lambda: self._begin_protocol_run(actives, config))
                return

            if sample_size and sample_size > 0:
                pairs = self._apply_random_sample(
                    pairs,
                    sample_size,
                    sample_seed,
                    include_controls=include_controls
                )

            config["planned_cases"] = {
                "total": len(pairs),
                "actives": sum(p.get("control_label") == 1 for p in pairs),
                "decoys": sum(p.get("control_label") == 0 for p in pairs),
                "samples": sum(p.get("control_label") is None for p in pairs),
            }

            rmsd_variant_available = any(p.get("site_ligand") for p in pairs)
            if not rmsd_variant_available and self.variant_mode_var.get() == "all_rmsd":
                msg = (
                    "Dock all variants, keep best RMSD requires at least one row with a co-crystal Ligand. "
                    "Use a score-based variant option for blank-ligand screening sheets."
                )
                self._run_on_ui(lambda m=msg: self._start_run_failed("Mode mismatch", m))
                return

            if config["mode"] == "adaptive":
                non_cocrystal_rows = sum(1 for p in pairs if p.get("site_mode") != "cocrystal")
                if non_cocrystal_rows:
                    msg = (
                        f"Found {non_cocrystal_rows} rows without a co-crystal site ligand. "
                        "Redock (adaptive search) requires a co-crystal ligand. "
                        "Use Redock (single protocol) or Screen compounds."
                    )
                    self._run_on_ui(lambda m=msg: self._start_run_failed("Mode mismatch", m))
                    return

                # Decoy/control rows have no reference pose, so the adaptive
                # cascade cannot terminate early and would run the full ladder
                # per compound for a score-only result. Caught here, before any
                # structure is downloaded.
                control_rows = sum(1 for p in pairs if p.get("control_label") is not None)
                if control_rows:
                    msg = (
                        f"This sheet contains {control_rows} decoy/control rows, which have "
                        "no reference pose to reproduce. The adaptive cascade would run every "
                        "protocol for each one and still report no RMSD.\n\n"
                        "Use 'Screen compounds' for compound sheets, or 'Redock' modes on a "
                        "sheet containing only each structure's own crystal ligand."
                    )
                    self._run_on_ui(lambda m=msg: self._start_run_failed("Mode mismatch", m))
                    return

            self._run_on_ui(lambda: self._begin_run(pairs, col_info, config))

        threading.Thread(target=_load_pairs_worker, daemon=True).start()

    def _start_run_failed(self, title: str, message: str) -> None:
        if title == "Excel error":
            messagebox.showerror(title, message)
        else:
            messagebox.showwarning(title, message)
        self._set_status(message)
        self._set_busy(False)

    def _begin_run(self, pairs: List[Dict[str, str]], col_info: dict, config: dict) -> None:
        config["columns"] = col_info

        self.progress_dialog = ProgressDialog(self, total_ligands=len(pairs))
        self._queue = queue.Queue()

        self._worker = threading.Thread(
            target=self._run_worker,
            args=(pairs, config),
            daemon=True
        )
        self._worker.start()
        self.after(200, self._poll_queue)
        self._set_status("Run started")

    @staticmethod
    def _parse_positive_int_list(value: str, field: str) -> List[int]:
        parts = [part.strip() for part in value.split(",") if part.strip()]
        if not parts:
            raise ValueError(f"{field} must contain at least one integer.")
        try:
            values = list(dict.fromkeys(int(part) for part in parts))
        except ValueError as exc:
            raise ValueError(f"{field} must be a comma-separated list of integers.") from exc
        if any(item <= 0 for item in values):
            raise ValueError(f"{field} values must be greater than zero.")
        return values

    @staticmethod
    def _parse_protocol_box_definitions(value: str) -> List[dict]:
        """Parse adaptive margins or fixed XYZ boxes separated by semicolons."""
        parts = [part.strip() for part in value.split(";") if part.strip()]
        if not parts:
            raise ValueError("Box definitions must contain at least one value.")
        definitions = []
        seen = set()
        for part in parts:
            lowered = part.lower().replace(" ", "")
            if lowered.startswith("margin:"):
                try:
                    margin = float(lowered.split(":", 1)[1])
                except ValueError as exc:
                    raise ValueError(f"Invalid box margin definition: {part}") from exc
                if margin <= 0:
                    raise ValueError("Box margins must be greater than zero.")
                definition = {
                    "label": f"margin:{margin:g}", "box_margin": margin,
                    "size_override": None,
                }
            else:
                dimensions = lowered.replace(",", "x").split("x")
                if len(dimensions) != 3:
                    raise ValueError(
                        f"Invalid box definition '{part}'. Use margin:4 or 20x20x20."
                    )
                try:
                    size = tuple(float(item) for item in dimensions)
                except ValueError as exc:
                    raise ValueError(f"Invalid box size definition: {part}") from exc
                if any(item <= 0 for item in size):
                    raise ValueError("Box dimensions must be greater than zero.")
                definition = {
                    "label": "x".join(f"{item:g}" for item in size),
                    "box_margin": None, "size_override": size,
                }
            if definition["label"] not in seen:
                seen.add(definition["label"])
                definitions.append(definition)
        return definitions

    @staticmethod
    def _parse_protocol_rescore_methods(value: str) -> List[str]:
        allowed = {
            "none", "vina", "vinardo", "ad4_scoring", "dkoes_fast", "dkoes_scoring"
        }
        methods = list(dict.fromkeys(
            part.strip().lower() for part in value.split(",") if part.strip()
        ))
        if not methods:
            raise ValueError("Rescoring methods must contain at least one value.")
        unsupported = [method for method in methods if method not in allowed]
        if unsupported:
            raise ValueError(
                "Unsupported rescoring method(s): " + ", ".join(unsupported)
            )
        return methods

    @staticmethod
    def _protocol_active_pairs(pairs: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """Return one positive control for each structure/crystal-ligand pair."""
        active_pairs = []
        seen = set()
        for pair in pairs:
            if pair.get("control_label") != 1 or not pair.get("site_ligand"):
                continue
            key = (pair.get("pdb_id"), pair.get("site_ligand"), pair.get("chain"))
            if key not in seen:
                seen.add(key)
                active_pairs.append(pair)
        return active_pairs

    def _collect_protocol_sweep_config(self) -> dict:
        water_modes = [
            mode for mode, variable in self.protocol_water_vars.items() if variable.get()
        ]
        if not water_modes:
            raise ValueError("Select at least one water-handling method.")
        engines = [
            engine for engine, variable in self.protocol_engine_vars.items() if variable.get()
        ]
        if not engines:
            raise ValueError("Select at least one docking engine.")
        return {
            "water_modes": water_modes,
            "engines": engines,
            "box_definitions": self._parse_protocol_box_definitions(
                self.protocol_box_definitions_var.get()
            ),
            "rescore_methods": self._parse_protocol_rescore_methods(
                self.protocol_rescore_methods_var.get()
            ),
            "exhaustiveness": self._parse_positive_int_list(
                self.protocol_exhaustiveness_var.get(), "Exhaustiveness"
            ),
            "seeds": self._parse_positive_int_list(self.protocol_seeds_var.get(), "Seeds"),
        }

    @staticmethod
    def _expand_protocol_conditions(
        actives: List[Dict[str, str]], sweep: dict
    ) -> List[tuple]:
        """Expand meaningful sweep conditions without rDock-only duplicates."""
        conditions = []
        for pair in actives:
            for engine in sweep["engines"]:
                # Exhaustiveness is a Vina/Smina parameter. rDock sampling is
                # controlled independently by rdock_runs, so sweep it once.
                exhaustiveness_values = (
                    [None] if engine == "rdock" else sweep["exhaustiveness"]
                )
                for box in sweep["box_definitions"]:
                    for water in sweep["water_modes"]:
                        for exhaustiveness in exhaustiveness_values:
                            for seed in sweep["seeds"]:
                                for rescore_method in sweep["rescore_methods"]:
                                    conditions.append((
                                        pair, engine, box, rescore_method, water,
                                        exhaustiveness, seed,
                                    ))
        return conditions

    def _begin_protocol_run(self, actives: List[Dict[str, str]], config: dict) -> None:
        sweep = config["protocol_sweep"]
        total = len(self._expand_protocol_conditions(actives, sweep))
        self.progress_dialog = ProgressDialog(self, total_ligands=total)
        self._queue = queue.Queue()
        self._worker = threading.Thread(
            target=self._run_protocol_worker, args=(actives, config), daemon=True
        )
        self._worker.start()
        self.after(200, self._poll_queue)
        self._set_status(f"Protocol sweep started ({total} conditions)")

    def _run_protocol_worker(self, actives: List[Dict[str, str]], config: dict) -> None:
        output_dir = config["output_dir"] / "protocol_development"
        output_dir.mkdir(parents=True, exist_ok=True)
        results_path = output_dir / "protocol_development_results.csv"
        manifest_path = output_dir / "protocol_development_manifest.json"
        rows = pd.read_csv(results_path).to_dict("records") if results_path.exists() else []

        sweep = config["protocol_sweep"]
        default_box = sweep["box_definitions"][0]["label"]

        def condition_identity(row: dict) -> Tuple[str, str, str, str, str, str, int, int]:
            engine = str(row.get("engine", config["single"].get("engine", "")))
            raw_exhaustiveness = row.get("exhaustiveness", -1)
            try:
                exhaustiveness = int(float(raw_exhaustiveness))
            except (TypeError, ValueError):
                exhaustiveness = -1
            if engine == "rdock":
                exhaustiveness = -1
            return (
                str(row.get("pdb_id", "")), str(row.get("ligand_resname", row.get("site_ligand", ""))),
                engine,
                str(row.get("box_definition", default_box)),
                str(row.get("rescore_method") or "none"),
                str(row.get("water_handling", "")), exhaustiveness,
                int(row.get("seed", -1)),
            )

        def condition_key(row: dict) -> tuple:
            identity = condition_identity(row)
            engine = identity[2]
            rescore_method = identity[4]
            cavity_version = row.get("rdock_cavity_version")
            if not cavity_version:
                cavity_version = "legacy" if engine == "rdock" else "not_applicable"
            pose_version = row.get("rescore_pose_version")
            if not pose_version:
                pose_version = "legacy"
            return (*identity, str(cavity_version), str(pose_version))

        completed = {
            condition_key(row) for row in rows
            if row.get("status") in ("complete", "unsupported")
        }
        conditions = self._expand_protocol_conditions(actives, sweep)
        pending_actives = []
        for pair in actives:
            has_pending_condition = any(
                (pair["pdb_id"], pair["site_ligand"], engine, box["label"], rescore_method,
                 water, -1 if engine == "rdock" else exhaustiveness, seed,
                 "crystal_or_center_v2" if engine == "rdock" else "not_applicable",
                 "all_poses_v2" if rescore_method != "none" else "baseline_v2")
                not in completed
                for candidate_pair, engine, box, rescore_method, water, exhaustiveness, seed
                in self._expand_protocol_conditions([pair], sweep)
            )
            if has_pending_condition:
                pending_actives.append(pair)
        if pending_actives and not (self.progress_dialog and self.progress_dialog.cancelled):
            try:
                self._prefetch_remote_inputs(pending_actives, output_dir)
            except Exception as exc:
                self._queue.put(("preflight_failed", str(exc)))
                return
        self._write_json_atomic(manifest_path, {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "software": self._software_provenance(config),
            "input_file": config["input_file"],
            "protocol_sweep": sweep,
            "single_protocol": config["single"],
            "active_count": len(actives),
        })

        needs_rescore = any(method != "none" for method in sweep["rescore_methods"])
        rescore_binary = (
            self._resolve_smina_binary(config.get("rescore", {}).get("smina_binary"))
            if needs_rescore else None
        )
        docking_cache = {}

        for index, (pair, engine_name, box, rescore_method, water, exhaustiveness, seed) in enumerate(conditions, 1):
            if self.progress_dialog and self.progress_dialog.cancelled:
                self._queue.put(("cancelled", results_path))
                return
            pdb_id = pair["pdb_id"]
            site_ligand = pair["site_ligand"]
            key = (
                pdb_id, site_ligand, engine_name, box["label"], rescore_method,
                water, -1 if engine_name == "rdock" else exhaustiveness, seed,
                "crystal_or_center_v2" if engine_name == "rdock" else "not_applicable",
                "all_poses_v2" if rescore_method != "none" else "baseline_v2",
            )
            label = (
                f"{pdb_id} {site_ligand}: {engine_name}, {box['label']}, "
                f"rescore={rescore_method}, "
                f"{water}, "
                f"{'rDock runs=' + str(config['single'].get('rdock_runs', 20)) if engine_name == 'rdock' else 'e' + str(exhaustiveness)}, "
                f"seed {seed}"
            )
            self._queue.put(("progress", index - 1, len(conditions), label))
            if key in completed:
                self._queue.put(("log", f"Skipping completed condition: {label}"))
                self._queue.put(("progress", index, len(conditions), f"{label} (resumed)"))
                continue

            row = {
                "pdb_id": pdb_id, "ligand_resname": site_ligand,
                "engine": engine_name, "box_definition": box["label"],
                "rescore_method": rescore_method,
                "rdock_cavity_version": (
                    "crystal_or_center_v2" if engine_name == "rdock" else "not_applicable"
                ),
                "rescore_pose_version": (
                    "all_poses_v2" if rescore_method != "none" else "baseline_v2"
                ),
                "water_handling": water, "exhaustiveness": exhaustiveness,
                "seed": seed, "status": "failed",
            }
            try:
                pdb_file = self._download_pdb(pdb_id, output_dir / "pdbs")
                ligand_chain = pair.get("chain") or self._detect_ligand_chain(pdb_file, site_ligand)
                if not ligand_chain:
                    raise ValueError(f"Ligand chain not found for {pdb_id}/{site_ligand}")
                if self._has_covalent_ligand_link(pdb_file, site_ligand, ligand_chain):
                    row["status"] = "unsupported"
                    row["error_message"] = (
                        "Covalent receptor-ligand LINK: standard Vina/Smina docking is invalid"
                    )
                    self._queue.put(("log", f"{label}: skipped covalent complex"))
                    raise StopIteration
                smiles = pair.get("smiles") or self._get_ligand_smiles(
                    pdb_file, site_ligand, ligand_chain, output_dir
                )
                if not smiles:
                    raise ValueError(f"Could not resolve SMILES for {pdb_id}/{site_ligand}")
                single_cfg = dict(config["single"])
                single_cfg.update({
                    "engine": engine_name,
                    "box_margin": (
                        box["box_margin"] if box["box_margin"] is not None
                        else single_cfg["box_margin"]
                    ),
                    "size_override": box["size_override"],
                    "water_handling": water,
                    "seed": seed,
                    "variant_select_by": "rmsd",
                })
                if exhaustiveness is not None:
                    single_cfg["exhaustiveness"] = exhaustiveness
                case_id = self._safe_case_id(
                    f"{pdb_id}_{site_ligand}_{engine_name}_{box['label']}_"
                    f"{water}_"
                    f"{'runs' + str(single_cfg.get('rdock_runs', 20)) if engine_name == 'rdock' else 'e' + str(exhaustiveness)}_s{seed}"
                )
                docking_key = (
                    pdb_id, site_ligand, engine_name, box["label"], water,
                    -1 if engine_name == "rdock" else exhaustiveness, seed,
                )
                if docking_key in docking_cache:
                    result = copy.deepcopy(docking_cache[docking_key])
                    self._queue.put(("log", f"Reusing docking poses for {label}"))
                else:
                    result = self._run_single_case(
                        pdb_file=pdb_file,
                        ligand_name=pair.get("dock_name") or site_ligand,
                        ligand_chain=ligand_chain,
                        smiles=smiles,
                        case_dir=output_dir / case_id,
                        threshold=config["threshold"],
                        single_cfg=single_cfg,
                        ligand_resname=site_ligand,
                        site_mode="cocrystal",
                        run_mode="protocol_development",
                        control_label=1,
                    )
                    docking_cache[docking_key] = copy.deepcopy(result)
                result.target_name = pair.get("target_name")
                row.update(asdict(result))
                if rescore_method == "none":
                    result.rescore_score = result.best_score
                    result.rescore_pose_count = result.pose_count
                    result.rescore_top1_rmsd = result.top1_rmsd
                    result.rescore_top5_rmsd = result.top5_rmsd
                    result.rescore_top10_rmsd = result.top10_rmsd
                    result.rescore_best_rmsd_rank = result.best_rmsd_rank
                    result.rescore_rmsd_best_score = result.rmsd_best_score
                    result.rescore_score_rmsd_pearson = result.score_rmsd_pearson
                    result.rescore_score_rmsd_spearman = result.score_rmsd_spearman
                    row.update(asdict(result))
                else:
                    if not rescore_binary:
                        result.rescore_error = "Smina binary not found"
                    elif result.output_file:
                        rescored = self._rescore_with_smina(
                            output_file=Path(result.output_file),
                            case_dir=output_dir / case_id,
                            smina_binary=rescore_binary,
                            scoring=rescore_method,
                        )
                        if rescored:
                            result.rescore_score = rescored.get("score")
                            result.rescore_error = rescored.get("error")
                            if rescored.get("scores"):
                                rescored_metrics = self._compute_rescored_pose_metrics(
                                    crystal_ligand_pdb=(output_dir / case_id / "crystal_ligand.pdb"),
                                    docked_file=Path(result.output_file),
                                    scores=rescored["scores"],
                                    threshold=config["threshold"],
                                    has_reference_pose=True,
                                )
                                for field, value in rescored_metrics.items():
                                    setattr(result, field, value)
                    row.update(asdict(result))
                row.update({
                    "engine": engine_name, "box_definition": box["label"],
                    "rescore_method": rescore_method,
                    "water_handling": water, "exhaustiveness": exhaustiveness,
                    "seed": seed, "status": "complete",
                })
                self._queue.put(("log", f"{label}: best RMSD {result.best_rmsd:.2f} A"))
            except StopIteration:
                pass
            except Exception as exc:
                row["error_message"] = str(exc)
                self._queue.put(("log", f"{label} failed: {exc}"))

            rows = [old for old in rows if condition_identity(old) != key[:-2]]
            rows.append(row)
            temporary = results_path.with_suffix(".tmp")
            pd.DataFrame(rows).to_csv(temporary, index=False)
            temporary.replace(results_path)
            self._queue.put(("progress", index, len(conditions), label))

        report_path = self._write_protocol_report(rows, output_dir, config["threshold"])
        self._queue.put(("protocol_done", results_path, report_path))

    @staticmethod
    def _write_protocol_report(
        rows: List[dict], output_dir: Path, threshold: float = 2.0
    ) -> Path:
        report_path = output_dir / "protocol_development_summary.md"
        frame = pd.DataFrame(rows)
        complete = frame[frame.get("status") == "complete"].copy() if not frame.empty else frame
        lines = ["# Protocol Development Summary", ""]
        if not frame.empty and "status" in frame:
            unsupported = int((frame["status"] == "unsupported").sum())
            if unsupported:
                lines.extend([
                    f"- Unsupported covalent conditions skipped: {unsupported}",
                    "- Covalent complexes require a dedicated covalent-docking method.",
                    "",
                ])
        if complete.empty:
            lines.append("No protocol conditions completed successfully.")
        else:
            if "engine" not in complete:
                complete["engine"] = "unknown"
            if "box_definition" not in complete:
                complete["box_definition"] = "legacy"
            if "rescore_method" not in complete:
                complete["rescore_method"] = "none"
            complete["rescore_method"] = complete["rescore_method"].fillna("none")
            if "rescore_error" not in complete:
                complete["rescore_error"] = None
            for column in (
                "top1_rmsd", "top5_rmsd", "top10_rmsd", "best_rmsd_rank",
                "best_score", "score_rmsd_pearson", "score_rmsd_spearman",
                "rescore_score", "rescore_top1_rmsd", "rescore_top5_rmsd",
                "rescore_top10_rmsd", "rescore_best_rmsd_rank",
                "rescore_score_rmsd_pearson", "rescore_score_rmsd_spearman",
            ):
                if column not in complete:
                    complete[column] = np.nan
            # Missing/unmappable RMSDs use 999.9 internally; never average the
            # sentinel into protocol-quality statistics.
            rmsd_columns = [
                "best_rmsd", "top1_rmsd", "top5_rmsd", "top10_rmsd",
                "rescore_top1_rmsd", "rescore_top5_rmsd", "rescore_top10_rmsd",
            ]
            for column in rmsd_columns:
                complete.loc[complete[column] >= 900, column] = np.nan

            # Older runs expanded rDock once for every Vina exhaustiveness
            # value even though rDock ignores that parameter. Normalize and
            # collapse those legacy duplicates before calculating statistics.
            complete.loc[complete["engine"] == "rdock", "exhaustiveness"] = np.nan
            legacy_identity = [
                column for column in (
                    "pdb_id", "ligand_resname", "ligand_chain", "target_name", "engine",
                    "box_definition", "rescore_method", "water_handling", "seed",
                    "exhaustiveness",
                )
                if column in complete.columns
            ]
            if legacy_identity:
                complete = complete.drop_duplicates(subset=legacy_identity, keep="last")

            complex_identity = [
                column for column in ("pdb_id", "ligand_resname", "ligand_chain")
                if column in complete.columns
            ]
            unique_complexes = (
                complete[complex_identity].drop_duplicates().shape[0]
                if complex_identity else 0
            )
            grouped = complete.groupby(
                ["engine", "box_definition", "rescore_method", "water_handling", "exhaustiveness"],
                dropna=False,
            ).agg(
                cases=("status", "size"),
                rmsd_cases=("best_rmsd", "count"),
                success_rate=("success", "mean"),
                mean_best_rmsd=("best_rmsd", "mean"),
                median_best_rmsd=("best_rmsd", "median"),
                mean_runtime_sec=("runtime_sec", "mean"),
                mean_baseline_score=("best_score", "mean"),
                mean_baseline_top1_rmsd=("top1_rmsd", "mean"),
                mean_baseline_top5_rmsd=("top5_rmsd", "mean"),
                mean_baseline_top10_rmsd=("top10_rmsd", "mean"),
                mean_baseline_best_rmsd_rank=("best_rmsd_rank", "mean"),
                mean_baseline_pearson=("score_rmsd_pearson", "mean"),
                mean_baseline_spearman=("score_rmsd_spearman", "mean"),
                rescore_failures=("rescore_error", lambda values: values.notna().sum()),
                mean_rescore_score=("rescore_score", "mean"),
                mean_rescore_top1_rmsd=("rescore_top1_rmsd", "mean"),
                mean_rescore_top5_rmsd=("rescore_top5_rmsd", "mean"),
                mean_rescore_top10_rmsd=("rescore_top10_rmsd", "mean"),
                mean_rescore_best_rmsd_rank=("rescore_best_rmsd_rank", "mean"),
                mean_rescore_pearson=("rescore_score_rmsd_pearson", "mean"),
                mean_rescore_spearman=("rescore_score_rmsd_spearman", "mean"),
            ).reset_index().sort_values(["success_rate", "mean_best_rmsd"], ascending=[False, True])
            lines.extend([
                f"- Unique crystal complexes evaluated: {unique_complexes}",
                f"- Completed protocol-ranking conditions: {len(complete)}",
                "",
                "A protocol condition is one parameter/ranking combination evaluated against a crystal complex. Multiple conditions for the same complex are not independent validation cases.",
                "",
                "The baseline columns use the docking engine's original ranking. When Smina/Vina docks with Vina scoring, these are the Vina results; the rescored columns show the selected Smina scoring function on the same poses.",
                "",
                "| Engine | Box | Rescoring | Water | Exhaust. | Conditions (with RMSD) | Mean best RMSD | Baseline Top-1/5/10 | Rescored Top-1/5/10 | Baseline / rescored best-RMSD rank | Baseline / rescored score | Baseline / rescored Spearman | Failures | Runtime (s) |",
                "| --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | ---: | ---: |",
            ])
            for _, item in grouped.iterrows():
                def _value(name: str) -> str:
                    value = getattr(item, name)
                    return "N/A" if pd.isna(value) else f"{value:.2f}"
                exhaustiveness_label = (
                    "N/A" if pd.isna(item.exhaustiveness)
                    else str(int(item.exhaustiveness))
                )
                lines.append(
                    f"| {item.engine} | {item.box_definition} | {item.rescore_method} | "
                    f"{item.water_handling} | "
                    f"{exhaustiveness_label} | {int(item.cases)} ({int(item.rmsd_cases)}) | "
                    f"{_value('mean_best_rmsd')} | "
                    f"{_value('mean_baseline_top1_rmsd')}/{_value('mean_baseline_top5_rmsd')}/{_value('mean_baseline_top10_rmsd')} | "
                    f"{_value('mean_rescore_top1_rmsd')}/{_value('mean_rescore_top5_rmsd')}/{_value('mean_rescore_top10_rmsd')} | "
                    f"{_value('mean_baseline_best_rmsd_rank')} / {_value('mean_rescore_best_rmsd_rank')} | "
                    f"{_value('mean_baseline_score')} / {_value('mean_rescore_score')} | "
                    f"{_value('mean_baseline_spearman')} / {_value('mean_rescore_spearman')} | "
                    f"{int(item.rescore_failures)} | "
                    f"{item.mean_runtime_sec:.1f} |"
                )

            lines.extend([
                "",
                "## Pose recovery across conditions",
                "",
                f"Success rates use an RMSD threshold of {threshold:g} A.",
                "",
                "| Water | Rescoring | Conditions | Best-pose success | Baseline Top-1 | Rescored Top-1 | Baseline Top-5 | Rescored Top-5 |",
                "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ])
            for (water, method), group in complete.groupby(
                ["water_handling", "rescore_method"], dropna=False
            ):
                valid_best = group["best_rmsd"].dropna()

                def _success_rate(column: str, limit: Optional[int] = None) -> str:
                    values = group[column].dropna()
                    if limit is not None:
                        # Top-N columns already contain the best RMSD within N.
                        values = values
                    if values.empty:
                        return "N/A"
                    return f"{(values < threshold).mean():.1%}"

                best_success = (
                    f"{(valid_best < threshold).mean():.1%}" if not valid_best.empty else "N/A"
                )
                lines.append(
                    f"| {water} | {method} | {len(group)} | {best_success} | "
                    f"{_success_rate('top1_rmsd')} | {_success_rate('rescore_top1_rmsd')} | "
                    f"{_success_rate('top5_rmsd', 5)} | {_success_rate('rescore_top5_rmsd', 5)} |"
                )

            recommendation_frame = complete.copy()
            if "target_name" not in recommendation_frame:
                recommendation_frame["target_name"] = None
            if "pdb_id" not in recommendation_frame:
                recommendation_frame["pdb_id"] = "Unknown"
            recommendation_frame["report_target"] = recommendation_frame["target_name"].where(
                recommendation_frame["target_name"].notna(), recommendation_frame["pdb_id"]
            )
            protocol_columns = [
                "report_target", "engine", "box_definition", "rescore_method",
                "water_handling", "exhaustiveness",
            ]
            recommendations = []
            for protocol_key, group in recommendation_frame.groupby(protocol_columns, dropna=False):
                for ranking, top1_column, top5_column, rank_column in (
                    ("baseline", "top1_rmsd", "top5_rmsd", "best_rmsd_rank"),
                    ("rescored", "rescore_top1_rmsd", "rescore_top5_rmsd", "rescore_best_rmsd_rank"),
                ):
                    top1 = group[top1_column].dropna()
                    top5 = group[top5_column].dropna()
                    if top1.empty:
                        continue
                    recommendations.append({
                        "target": protocol_key[0], "engine": protocol_key[1],
                        "box": protocol_key[2], "rescore": protocol_key[3],
                        "water": protocol_key[4], "exhaustiveness": protocol_key[5],
                        "ranking": ranking, "cases": len(top1),
                        "top1_success": float((top1 < threshold).mean()),
                        "top5_success": float((top5 < threshold).mean()) if not top5.empty else 0.0,
                        "mean_top1": float(top1.mean()),
                        "mean_best": float(group["best_rmsd"].mean()),
                        "mean_rank": float(group[rank_column].mean()),
                    })

            best_by_target = {}
            for candidate in recommendations:
                key = (
                    -candidate["top1_success"], -candidate["top5_success"],
                    candidate["mean_rank"], candidate["mean_best"], candidate["mean_top1"],
                )
                current = best_by_target.get(candidate["target"])
                if current is None or key < current[0]:
                    best_by_target[candidate["target"]] = (key, candidate)

            lines.extend([
                "",
                "## Recommended protocol per target",
                "",
                "Recommendations prioritize Top-1 success, Top-5 success, best-pose rank, mean best RMSD, and then mean Top-1 RMSD. This prevents tiny differences between failed Top-1 poses from outranking a protocol that places the native-like pose earlier. Recommendations validate pose recovery only; confirm enrichment with matched decoys before screening.",
                "",
                "| Target | Engine | Box | Water | Exhaust. | Ranking | Rescorer | Conditions | Top-1 success | Top-5 success | Mean Top-1 RMSD | Mean best RMSD | Mean best-pose rank |",
                "| --- | --- | --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ])
            for target in sorted(best_by_target, key=lambda value: str(value)):
                item = best_by_target[target][1]
                rescorer = item["rescore"] if item["ranking"] == "rescored" else "original score"
                exhaustiveness_label = (
                    "N/A" if pd.isna(item["exhaustiveness"])
                    else str(int(item["exhaustiveness"]))
                )
                lines.append(
                    f"| {target} | {item['engine']} | {item['box']} | {item['water']} | "
                    f"{exhaustiveness_label} | {item['ranking']} | {rescorer} | "
                    f"{item['cases']} | {item['top1_success']:.1%} | {item['top5_success']:.1%} | "
                    f"{item['mean_top1']:.2f} | {item['mean_best']:.2f} | {item['mean_rank']:.2f} |"
                )
        report_path.write_text("\n".join(lines) + "\n")
        return report_path

    def _poll_queue(self) -> None:
        if not self.progress_dialog:
            return

        while not self._queue.empty():
            msg = self._queue.get()
            msg_type = msg[0]
            if msg_type == "progress":
                _, current, total, text = msg
                self.progress_dialog.update_progress(current, total, text)
            elif msg_type == "log":
                _, text = msg
                self.progress_dialog.log(text)
            elif msg_type == "done":
                _, results_path = msg
                self.progress_dialog.log(f"Results saved to {results_path}")
                self.progress_dialog.destroy()
                self.progress_dialog = None
                messagebox.showinfo("Completed", f"Docking analysis completed.\n{results_path}")
                self._set_status("Run completed")
                self.last_results_path = Path(results_path)
                self._safe_call(self._render_results_from_path)(results_path)
                self._safe_call(self._show_results)(results_path)
                self._set_busy(False)
                return
            elif msg_type == "protocol_done":
                _, results_path, report_path = msg
                self.progress_dialog.log(f"Protocol results saved to {results_path}")
                self.progress_dialog.destroy()
                self.progress_dialog = None
                messagebox.showinfo(
                    "Completed", f"Protocol development completed.\n{report_path}"
                )
                self._set_status("Protocol sweep completed")
                self.last_results_path = Path(results_path)
                self._safe_call(self._render_protocol_results)(results_path, report_path)
                self._safe_call(self._show_protocol_results)(results_path, report_path)
                self._set_busy(False)
                return
            elif msg_type == "preflight_failed":
                _, error = msg
                self.progress_dialog.destroy()
                self.progress_dialog = None
                messagebox.showerror(
                    "Input download failed",
                    "The required online inputs could not be downloaded before docking.\n\n"
                    f"{error}\n\nNo docking calculations were started."
                )
                self._set_status("Input download failed; docking not started")
                self._set_busy(False)
                return
            elif msg_type == "cancelled":
                _, results_path = msg
                self.progress_dialog.log(f"Partial results saved to {results_path}")
                self.progress_dialog.destroy()
                self.progress_dialog = None
                messagebox.showinfo(
                    "Cancelled",
                    "Docking campaign cancelled. Completed cases were saved and "
                    "will be skipped when the same run is restarted."
                )
                self._set_status("Run cancelled; progress saved")
                self.last_results_path = Path(results_path)
                self._set_busy(False)
                return

        if self._worker and self._worker.is_alive():
            self.after(200, self._poll_queue)
        else:
            if self.progress_dialog:
                self.progress_dialog.destroy()
                self.progress_dialog = None
                self._set_status("Run finished")
            self._set_busy(False)

    def _start_ui_queue(self) -> None:
        def _poll() -> None:
            while True:
                try:
                    func = self._ui_queue.get_nowait()
                except queue.Empty:
                    break
                try:
                    func()
                except Exception as exc:
                    logger.error("UI task failed: {}", exc)
            self._ui_poll_job = self.after(100, _poll)

        if self._ui_poll_job is None:
            self._ui_poll_job = self.after(100, _poll)

    def _run_on_ui(self, func: Callable[[], None]) -> None:
        if threading.current_thread() is threading.main_thread():
            func()
        else:
            self._ui_queue.put(func)

    def _collect_single_config(self) -> Optional[dict]:
        if self.mode_var.get() not in ("single", "screening", "protocol_development"):
            return {}

        variant_mode, variant_select_by = self._variant_config()
        requested_selection = variant_select_by
        variant_select_by = self._variant_selection_for_mode(
            self.mode_var.get(), variant_select_by
        )
        if variant_select_by != requested_selection:
            logger.warning(
                "Screening mode requires score-based ligand-variant selection; "
                "overriding '{}' with 'score'",
                requested_selection
            )
        max_tautomers = self._parse_int(self.max_tautomers_var.get(), "Max tautomers")
        max_conformers = self._parse_int(self.max_conformers_var.get(), "Max conformers")
        if None in (max_tautomers, max_conformers):
            logger.error("Single config invalid: ligand variant counts")
            return None
        if not (1 <= max_tautomers <= 20):
            messagebox.showerror("Input error", "Max tautomers must be between 1 and 20.")
            return None
        if not (1 <= max_conformers <= 30):
            messagebox.showerror("Input error", "Max conformers must be between 1 and 30.")
            return None
        engine = self.engine_var.get()
        logger.info(
            "Collecting single config | engine={} box_margin={} size_override={},{},{} "
            "exhaustiveness={} num_modes={} energy_range={} cpu={} seed={} timeout={} "
            "scoring={} water={}",
            engine,
            self.box_margin_var.get(),
            self.size_x_var.get(),
            self.size_y_var.get(),
            self.size_z_var.get(),
            self.exhaustiveness_var.get(),
            self.num_modes_var.get(),
            self.energy_range_var.get(),
            self.cpu_var.get(),
            self.seed_var.get(),
            self.timeout_var.get(),
            self.scoring_var.get(),
            self.water_handling_var.get()
        )
        box_margin = self._parse_float(self.box_margin_var.get(), "Box margin")
        if box_margin is None:
            logger.error("Single config invalid: box margin")
            return None

        size_override, size_override_ok = self._parse_size_override()
        if not size_override_ok:
            logger.error("Single config invalid: box size override")
            return None

        exhaustiveness = self._parse_int(self.exhaustiveness_var.get(), "Exhaustiveness")
        num_modes = self._parse_int(self.num_modes_var.get(), "Num modes")
        energy_range = self._parse_float(self.energy_range_var.get(), "Energy range")
        cpu = self._parse_int(self.cpu_var.get(), "CPU")
        seed = self._parse_int(self.seed_var.get(), "Seed")
        timeout = self._parse_int(self.timeout_var.get(), "Timeout")
        if None in (exhaustiveness, num_modes, energy_range, cpu, seed, timeout):
            logger.error("Single config invalid: numeric fields")
            return None

        rdock_runs = None
        rdock_seed = None
        rdock_radius = None
        if engine == "rdock":
            logger.info(
                "Collecting rDock config | root={} runs={} seed={} radius={}",
                self.rdock_root_var.get(),
                self.rdock_runs_var.get(),
                self.rdock_seed_var.get(),
                self.rdock_radius_var.get()
            )
            rdock_runs = self._parse_int(self.rdock_runs_var.get(), "rDock runs")
            rdock_seed = self._parse_int(self.rdock_seed_var.get(), "rDock seed")
            if None in (rdock_runs, rdock_seed):
                logger.error("Single config invalid: rDock runs/seed")
                return None

            rdock_radius = self._parse_float(self.rdock_radius_var.get(), "rDock radius", allow_blank=True)
            if rdock_radius is None and self.rdock_radius_var.get().strip():
                logger.error("Single config invalid: rDock radius")
                return None

        return {
            "engine": engine,
            "box_margin": box_margin,
            "apo_site_mode": self.apo_site_mode_var.get(),
            "site_definition_mode": self.site_definition_mode_var.get(),
            "site_residues": self.site_residues_var.get().strip(),
            "size_override": size_override,
            "water_handling": self.water_handling_var.get(),
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "energy_range": energy_range,
            "cpu": cpu,
            "seed": seed,
            "timeout": timeout,
            "scoring": self.scoring_var.get(),
            "smina_binary": self.smina_bin_var.get(),
            "vina_binary": self.vina_bin_var.get(),
            "rdock_root": self.rdock_root_var.get(),
            "rdock_runs": rdock_runs or 20,
            "rdock_seed": rdock_seed or 42,
            "rdock_radius": rdock_radius,
            "ligand_variant_mode": variant_mode,
            "variant_select_by": variant_select_by,
            "max_tautomers": max_tautomers,
            "max_conformers": max_conformers,
            "n_cpus": cpu
        }

    def _collect_adaptive_config(self) -> Optional[dict]:
        if self.mode_var.get() != "adaptive":
            return {}

        variant_mode, variant_select_by = self._variant_config()
        max_tautomers = self._parse_int(self.max_tautomers_var.get(), "Max tautomers")
        max_conformers = self._parse_int(self.max_conformers_var.get(), "Max conformers")
        if None in (max_tautomers, max_conformers):
            logger.error("Adaptive config invalid: ligand variant counts")
            return None
        if not (1 <= max_tautomers <= 20):
            messagebox.showerror("Input error", "Max tautomers must be between 1 and 20.")
            return None
        if not (1 <= max_conformers <= 30):
            messagebox.showerror("Input error", "Max conformers must be between 1 and 30.")
            return None
        rdock_root = Path(self.rdock_root_var.get())
        enable_rdock = self.enable_rdock_var.get() and rdock_root.exists()
        if self.enable_rdock_var.get() and not rdock_root.exists():
            messagebox.showwarning(
                "rDock not found",
                f"rDock root not found: {rdock_root}\nAdaptive docking will run without rDock."
            )
            logger.warning("Adaptive config: rDock root not found: {}", rdock_root)
            enable_rdock = False

        return {
            "enable_rdock": enable_rdock,
            "rdock_root": rdock_root,
            "use_vina": self.use_vina_var.get(),
            "smina_binary": self.smina_bin_var.get(),
            "vina_binary": self.vina_bin_var.get(),
            "ligand_variant_mode": variant_mode,
            "variant_select_by": variant_select_by,
            "max_tautomers": max_tautomers,
            "max_conformers": max_conformers,
            "n_cpus": self._parse_int(self.cpu_var.get(), "CPU") or 4
        }

    def _run_worker(self, pairs: List[Dict[str, str]], config: dict) -> None:
        output_dir = config["output_dir"]
        output_dir.mkdir(parents=True, exist_ok=True)
        progress_path = output_dir / "redock_progress.json"
        results_path = output_dir / "redock_results.json"
        results_csv = output_dir / "redock_results.csv"
        manifest_path = output_dir / "run_manifest.json"
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "software": self._software_provenance(config),
            "config": config,
            "cases": [
                {
                    "pdb_id": p.get("pdb_id"),
                    "site_ligand": p.get("site_ligand"),
                    "dock_name": p.get("dock_name"),
                    "control_label": p.get("control_label"),
                    "case_id": p.get("case_id"),
                    "reference_pose_policy": "molecular_identity_v1",
                    **({"target_name": p.get("target_name")} if p.get("target_name") else {}),
                }
                for p in pairs
            ],
        }
        resumed = self._load_resumable_results(manifest_path, progress_path, manifest)
        results_by_case = {result.case_id: result for result in resumed if result.case_id}
        if resumed:
            self._queue.put((
                "log",
                f"Resuming compatible run: {len(resumed)} completed cases will be skipped"
            ))
        self._write_json_atomic(manifest_path, manifest)

        pending_pairs = [
            pair for pair in pairs
            if (pair.get("case_id") or f"{pair['pdb_id']}_{pair.get('dock_name') or pair['ligand']}")
            not in results_by_case
        ]
        if pending_pairs and not (self.progress_dialog and self.progress_dialog.cancelled):
            try:
                self._prefetch_remote_inputs(pending_pairs, output_dir)
            except Exception as exc:
                self._queue.put(("preflight_failed", str(exc)))
                return

        rescore_cfg = config.get("rescore", {})
        rescore_enabled = bool(rescore_cfg.get("enable"))
        rescore_binary = None
        if rescore_enabled:
            rescore_binary = self._resolve_smina_binary(rescore_cfg.get("smina_binary"))
            if not rescore_binary:
                self._queue.put(("log", "Smina binary not found; rescoring disabled"))
                rescore_enabled = False

        cancelled = False
        for idx, item in enumerate(pairs, 1):
            if self.progress_dialog and self.progress_dialog.cancelled:
                self._queue.put(("log", "Run cancelled by user"))
                cancelled = True
                break

            pdb_id = item["pdb_id"]
            ligand = item["ligand"]
            site_ligand = item.get("site_ligand")
            dock_name = item.get("dock_name") or ligand
            site_mode = item.get("site_mode") or ("cocrystal" if site_ligand else "prediction")
            pocket_center = item.get("pocket_center")
            site_residues = item.get("site_residues")
            single_cfg = dict(config.get("single", {}))
            single_site_mode = single_cfg.get("site_definition_mode", "auto")
            if single_site_mode == "cocrystal_ligand":
                site_mode = "cocrystal"
            elif single_site_mode == "detected_pocket":
                site_mode = "prediction"
                single_cfg["apo_site_mode"] = "fpocket"
            elif single_site_mode == "specified_residues":
                site_mode = "residues"
                site_residues = site_residues or single_cfg.get("site_residues")
            elif single_site_mode == "protein_centroid":
                site_mode = "prediction"
                single_cfg["apo_site_mode"] = "protein_centroid"
            chain = item.get("chain")
            control_label = item.get("control_label")
            target_name = item.get("target_name")
            # Report the case as STARTING (idx-1 complete), not as finished.
            # Emitting idx here drove the bar to 100% while the final case was
            # still docking.
            self._queue.put(("progress", idx - 1, len(pairs), f"{pdb_id} {ligand}"))

            case_id = item.get("case_id") or f"{pdb_id}_{dock_name}"
            if case_id in results_by_case:
                self._queue.put(("log", f"Skipping completed case {case_id}"))
                self._queue.put(("progress", idx, len(pairs), f"{pdb_id} {ligand} resumed"))
                continue
            case_dir = output_dir / self._safe_case_id(case_id)
            case_dir.mkdir(parents=True, exist_ok=True)

            try:
                pdb_file = self._download_pdb(pdb_id, output_dir / "pdbs")
                if site_mode == "cocrystal":
                    if not site_ligand:
                        raise ValueError("Co-crystal mode requires a ligand residue name")
                    ligand_chain = chain or self._detect_ligand_chain(pdb_file, site_ligand)
                    if not ligand_chain:
                        raise ValueError("Ligand chain not found")
                else:
                    ligand_chain = chain or ""

                smiles = item.get("smiles")
                if smiles:
                    smiles = str(smiles).strip()
                    if not smiles or smiles.lower() == "nan":
                        smiles = None
                if not smiles and site_ligand:
                    smiles = self._get_ligand_smiles(pdb_file, site_ligand, ligand_chain, output_dir)
                if not smiles:
                    raise ValueError(
                        "Could not resolve ligand SMILES. Provide a SMILES column for apo/pocket-detection rows."
                    )
                ligand_charge = self._get_ligand_charge(smiles)
                ligand_properties = self._get_ligand_properties(smiles)
                is_reference_case = (
                    config["mode"] != "screening"
                    or item.get("is_reference_ligand") is True
                )

                if config["mode"] == "adaptive":
                    if site_mode != "cocrystal":
                        raise ValueError(
                            "Redock (adaptive search) requires a co-crystal ligand. "
                            "Use Redock (single protocol) for apo pocket-detection workflows."
                        )
                    # The adaptive ladder escalates until RMSD clears the threshold.
                    # If the docked ligand is not the crystal ligand there is no
                    # reference pose, RMSD is undefined, the early-exit can never
                    # fire, and every compound burns the full cascade to produce a
                    # score the search was not optimising for. Refuse up front.
                    is_ref = self._is_reference_ligand(
                        pdb_file, site_ligand, ligand_chain, smiles, output_dir
                    )
                    if control_label is not None or is_ref is False:
                        raise ValueError(
                            f"'{dock_name}' is not the co-crystal ligand of {pdb_id} "
                            f"({site_ligand}), so RMSD is undefined and the adaptive "
                            "cascade cannot terminate early — every compound would run "
                            "the full protocol ladder for a score-only result. "
                            "Use Screen compounds for novel ligands and decoys; use "
                            "Redock modes only on the structure's own crystal ligand."
                        )
                    if is_ref is None:
                        logger.warning(
                            "Could not confirm %s is the co-crystal ligand of %s; "
                            "proceeding, but check the RMSD is meaningful.",
                            dock_name, pdb_id
                        )
                    result = self._run_adaptive_case(
                        pdb_file,
                        dock_name,
                        ligand_chain,
                        smiles,
                        case_dir,
                        config["threshold"],
                        config["adaptive"],
                        ligand_resname=site_ligand,
                        control_label=control_label
                    )
                else:
                    result = self._run_single_case(
                        pdb_file,
                        dock_name,
                        ligand_chain,
                        smiles,
                        case_dir,
                        config["threshold"],
                        single_cfg,
                        ligand_resname=site_ligand,
                        site_mode=site_mode,
                        pocket_center=pocket_center,
                        site_residues=site_residues,
                        run_mode=config["mode"],
                        control_label=control_label,
                        is_reference_ligand=is_reference_case,
                    )
                result.control_label = control_label
                result.target_name = target_name
                result.ligand_charge = ligand_charge
                result.molecular_weight = ligand_properties.get("molecular_weight")
                result.logp = ligand_properties.get("logp")
                result.tpsa = ligand_properties.get("tpsa")
                result.rotatable_bonds = ligand_properties.get("rotatable_bonds")
                result.dock_name = dock_name
                result.docking_completed = bool(result.output_file)
                result.case_id = case_id

                if rescore_enabled and rescore_binary and result.output_file:
                    out_path = Path(result.output_file)
                    if out_path.exists():
                        rescore = self._rescore_with_smina(
                            output_file=out_path,
                            case_dir=case_dir,
                            smina_binary=rescore_binary,
                            scoring=rescore_cfg.get("scoring", "vina")
                        )
                        if rescore:
                            result.rescore_method = rescore.get("method")
                            result.rescore_score = rescore.get("score")
                            result.rescore_error = rescore.get("error")
                            if rescore.get("scores"):
                                rescored_metrics = self._compute_rescored_pose_metrics(
                                    crystal_ligand_pdb=case_dir / "crystal_ligand.pdb",
                                    docked_file=out_path,
                                    scores=rescore["scores"],
                                    threshold=config["threshold"],
                                    has_reference_pose=is_reference_case,
                                )
                                for field, value in rescored_metrics.items():
                                    setattr(result, field, value)
                results_by_case[case_id] = result
                if result.best_rmsd >= 900 and result.best_score is not None:
                    self._queue.put((
                        "log",
                        f"{pdb_id} {dock_name} score={result.best_score:.2f} site={result.site_method or site_mode}"
                    ))
                else:
                    self._queue.put(("log", f"{pdb_id} {ligand} RMSD={result.best_rmsd:.2f}"))

            except Exception as exc:
                results_by_case[case_id] = RedockResult(
                    pdb_id=pdb_id,
                    ligand_resname=site_ligand or ligand,
                    ligand_chain=chain or "",
                    mode=config["mode"],
                    engine=config["single"].get("engine", "adaptive"),
                    protocol="N/A",
                    best_rmsd=999.9,
                    success=False,
                    runtime_sec=0.0,
                    error_message=str(exc),
                    control_label=control_label,
                    dock_name=dock_name,
                    site_method=site_mode,
                    docking_completed=False,
                    case_id=case_id,
                    target_name=target_name
                )
                self._queue.put(("log", f"{pdb_id} {ligand} failed: {exc}"))

            ordered_results = self._ordered_results(pairs, results_by_case)
            self._write_progress(progress_path, ordered_results)
            # Case finished — advance the bar only now.
            self._queue.put(("progress", idx, len(pairs), f"{pdb_id} {ligand} done"))

        results = self._ordered_results(pairs, results_by_case)
        self._write_results(results_path, results_csv, results, config["threshold"])
        self._queue.put(("cancelled" if cancelled else "done", results_path))

    def _run_adaptive_case(
        self,
        pdb_file: Path,
        ligand_name: str,
        ligand_chain: str,
        smiles: str,
        case_dir: Path,
        threshold: float,
        adaptive_cfg: dict,
        ligand_resname: Optional[str] = None,
        control_label: Optional[int] = None
    ) -> RedockResult:
        if self._has_covalent_ligand_link(pdb_file, ligand_resname or ligand_name, ligand_chain):
            raise ValueError(
                f"{pdb_file.stem.upper()}/{ligand_resname or ligand_name} is covalently "
                "linked to the receptor. Standard Vina/Smina redocking is not valid for "
                "this complex; use a covalent-docking protocol or another validation structure."
            )
        pipeline = AdaptiveDockingPipeline(
            output_dir=case_dir,
            rmsd_threshold=threshold,
            rdock_available=adaptive_cfg["enable_rdock"],
            rdock_root=adaptive_cfg["rdock_root"],
            use_vina=adaptive_cfg.get("use_vina", False),
            smina_binary=adaptive_cfg.get("smina_binary", "smina"),
            vina_binary=adaptive_cfg.get("vina_binary", "vina"),
            ligand_variant_mode=adaptive_cfg.get("ligand_variant_mode", "first"),
            variant_select_by=adaptive_cfg.get("variant_select_by", "rmsd"),
            max_tautomers=adaptive_cfg.get("max_tautomers", 8),
            max_conformers=adaptive_cfg.get("max_conformers", 10),
            n_cpus=adaptive_cfg.get("n_cpus")
        )

        start = time.time()
        best_result, all_results = pipeline.run_adaptive_docking(
            pdb_file=pdb_file,
            ligand_smiles=smiles,
            ligand_name=ligand_name,
            ligand_resname=ligand_resname or ligand_name,
            ligand_chain=ligand_chain
        )
        runtime = time.time() - start

        if best_result is None:
            raise ValueError("Adaptive docking returned no result")

        metrics = self._compute_pose_metrics(
            crystal_ligand_pdb=case_dir / "crystal_ligand.pdb",
            docked_file=best_result.output_file,
            threshold=threshold,
            # The adaptive path rejects any labelled row upstream, so reaching
            # here means this IS the structure's own crystal ligand.
            has_reference_pose=True,
        )
        engine = "rdock" if best_result.output_file.suffix in (".sd", ".sdf") else "smina"
        protocol_attempts = [
            {
                "protocol": attempt.protocol_name,
                "rmsd": attempt.rmsd,
                "success": attempt.success
            }
            for attempt in all_results
        ]

        best_rmsd = metrics.get("best_rmsd")
        if best_rmsd is None:
            best_rmsd = best_result.rmsd if best_result.rmsd is not None else 999.9
        best_score_value = metrics.get("best_score")
        if best_score_value is None:
            best_score_value = best_result.score if hasattr(best_result, "score") else None

        return RedockResult(
            pdb_id=pdb_file.stem.upper(),
            ligand_resname=ligand_resname or ligand_name,
            ligand_chain=ligand_chain,
            mode="adaptive",
            engine=engine,
            protocol=best_result.protocol_name,
            best_rmsd=best_rmsd,
            success=best_rmsd < threshold,
            runtime_sec=runtime,
            output_file=str(best_result.output_file),
            pose_count=metrics.get("pose_count"),
            best_score=best_score_value,
            top1_rmsd=metrics.get("top1_rmsd"),
            top5_rmsd=metrics.get("top5_rmsd"),
            top10_rmsd=metrics.get("top10_rmsd"),
            best_rmsd_rank=metrics.get("best_rmsd_rank"),
            rmsd_best_score=metrics.get("rmsd_best_score"),
            rmsd_mean=metrics.get("rmsd_mean"),
            rmsd_median=metrics.get("rmsd_median"),
            rmsd_std=metrics.get("rmsd_std"),
            near_native_fraction=metrics.get("near_native_fraction"),
            score_rmsd_pearson=metrics.get("score_rmsd_pearson"),
            score_rmsd_spearman=metrics.get("score_rmsd_spearman"),
            protocols_tried=len(all_results),
            protocols_success=sum(1 for attempt in all_results if attempt.success),
            protocol_attempts=protocol_attempts,
            error_message=best_result.error_message,
            dock_name=ligand_name,
            site_method="cocrystal"
        )

    def _run_single_case(
        self,
        pdb_file: Path,
        ligand_name: str,
        ligand_chain: Optional[str],
        smiles: str,
        case_dir: Path,
        threshold: float,
        single_cfg: dict,
        ligand_resname: Optional[str] = None,
        site_mode: str = "cocrystal",
        pocket_center: Optional[Tuple[float, float, float]] = None,
        site_residues: Optional[str] = None,
        run_mode: str = "single",
        control_label: Optional[int] = None,
        is_reference_ligand: Optional[bool] = None,
    ) -> RedockResult:
        if site_mode == "cocrystal" and self._has_covalent_ligand_link(
            pdb_file, ligand_resname or ligand_name, ligand_chain
        ):
            raise ValueError(
                f"{pdb_file.stem.upper()}/{ligand_resname or ligand_name} is covalently "
                "linked to the receptor. Standard Vina/Smina redocking is not valid for "
                "this complex; use a covalent-docking protocol or another validation structure."
            )
        engine_name = single_cfg["engine"]
        pipeline = AdaptiveDockingPipeline(
            output_dir=case_dir,
            rdock_available=True,
            rdock_root=Path(single_cfg["rdock_root"]),
            smina_binary=single_cfg["smina_binary"],
            vina_binary=single_cfg["vina_binary"],
            ligand_variant_mode=single_cfg.get("ligand_variant_mode", "first"),
            variant_select_by=single_cfg.get("variant_select_by", "rmsd"),
            max_tautomers=single_cfg.get("max_tautomers", 8),
            max_conformers=single_cfg.get("max_conformers", 10),
            n_cpus=single_cfg.get("n_cpus")
        )

        receptor_pdbqt, receptor_pdb = pipeline._prepare_receptor(
            pdb_file,
            water_handling=single_cfg["water_handling"],
            site_ligand_resname=(
                ligand_resname or ligand_name if site_mode == "cocrystal" else None
            )
        )

        enumerate_states = not pipeline._contains_metal(smiles)
        variants = pipeline._prepare_ligand_variants(
            ligand_smiles=smiles,
            ligand_name=ligand_name,
            enumerate_states=enumerate_states
        )
        variant_mode = single_cfg.get("ligand_variant_mode", "first")
        variants_prepared = len(variants)
        if variant_mode == "best":
            variants = [pipeline._select_best_variant(variants)]
        elif variant_mode == "first":
            variants = variants[:1]
        elif variant_mode == "adaptive":
            variants = pipeline._adaptive_variant_selection(
                variants, ligand_smiles=smiles, ligand_name=ligand_name
            )
        elif variant_mode == "thorough":
            variants = pipeline._select_diverse_variants(
                variants, min(15, len(variants))
            )
        variants_docked = len(variants)

        binding_site: BindingSite
        size_override = single_cfg["size_override"]
        if site_mode == "cocrystal":
            crystal_ligand_pdb = pipeline._extract_crystal_ligand(
                pdb_file,
                ligand_resname or ligand_name,
                ligand_chain
            )
            binding_site = BindingSiteDefinition(margin=single_cfg["box_margin"]).from_cocrystal(
                pdb_file,
                ligand_resname=ligand_resname or ligand_name,
                ligand_chain=ligand_chain
            )
            if size_override is not None:
                binding_site.size = size_override
            site_method = "cocrystal"
        elif site_mode == "residues":
            crystal_ligand_pdb = case_dir / "crystal_ligand.pdb"
            binding_site, site_method = self._binding_site_from_residues(
                pdb_path=pdb_file,
                residues_text=site_residues or "",
                box_margin=single_cfg["box_margin"],
                size_override=size_override,
                default_chain=ligand_chain
            )
            self._queue.put(
                ("log", f"{pdb_file.stem.upper()} {ligand_name} site={site_method} "
                        f"center=({binding_site.center[0]:.2f},{binding_site.center[1]:.2f},{binding_site.center[2]:.2f})")
            )
        else:
            crystal_ligand_pdb = case_dir / "crystal_ligand.pdb"
            binding_site, site_method = self._predict_binding_site(
                pdb_path=pdb_file,
                case_dir=case_dir,
                box_margin=single_cfg["box_margin"],
                size_override=size_override,
                apo_site_mode=single_cfg.get("apo_site_mode", "auto"),
                manual_center=pocket_center
            )
            self._queue.put(
                ("log", f"{pdb_file.stem.upper()} {ligand_name} site={site_method} "
                        f"center=({binding_site.center[0]:.2f},{binding_site.center[1]:.2f},{binding_site.center[2]:.2f})")
            )

        start = time.time()
        variant_results = []
        for variant in variants:
            variant_label = variant["label"]
            ligand_pdbqt = variant["pdbqt"]
            variant_smiles = variant.get("smiles") or smiles
            variant_dir = case_dir / "variants" / variant_label
            variant_dir.mkdir(parents=True, exist_ok=True)
            output_file = None
            docking_result = None

            if engine_name in ("smina", "vina"):
                binary = single_cfg["smina_binary"] if engine_name == "smina" else single_cfg["vina_binary"]
                scoring = single_cfg["scoring"] if engine_name == "smina" else "vina"

                engine = SminaDockingEngine(
                    smina_binary=binary,
                    scoring_function=scoring,
                    exhaustiveness=single_cfg["exhaustiveness"],
                    num_modes=single_cfg["num_modes"],
                    energy_range=single_cfg["energy_range"],
                    cpu=single_cfg["cpu"],
                    seed=single_cfg["seed"],
                    timeout_sec=single_cfg["timeout"]
                )
                output_path = variant_dir / "docked.pdbqt"
                docking_result = engine.dock(
                    receptor_path=str(receptor_pdbqt),
                    ligand_path=str(ligand_pdbqt),
                    center=binding_site.center,
                    size=binding_site.size,
                    output_path=str(output_path)
                )
                if docking_result.success:
                    output_file = Path(output_path)
            else:
                docking_result = pipeline._dock_with_rdock(
                    receptor_pdb=pdb_file,
                    ligand_smiles=variant_smiles,
                    ligand_name=variant_label,
                    binding_site=binding_site,
                    output_dir=variant_dir,
                    prepared_receptor_pdb=receptor_pdb,
                    reference_ligand_pdb=crystal_ligand_pdb,
                    radius_override=single_cfg["rdock_radius"],
                    runs=single_cfg["rdock_runs"],
                    seed=single_cfg["rdock_seed"]
                )
                if docking_result.get("success"):
                    output_file = Path(docking_result["output_file"])

            metrics = {}
            rmsd = 999.9
            best_score = None
            if output_file and output_file.exists():
                metrics = self._compute_pose_metrics(
                    crystal_ligand_pdb=crystal_ligand_pdb,
                    docked_file=output_file,
                    threshold=threshold,
                    has_reference_pose=(
                        crystal_ligand_pdb.exists()
                        and (run_mode != "screening" or is_reference_ligand is True)
                    ),
                )
                rmsd = metrics.get("best_rmsd", 999.9) if metrics else 999.9
                best_score = metrics.get("best_score")
            if best_score is None:
                if engine_name in ("smina", "vina") and docking_result:
                    best_score = docking_result.poses[0].score if docking_result.poses else None
                elif docking_result:
                    best_score = docking_result.get("score")

            variant_results.append({
                "label": variant_label,
                "output_file": output_file,
                "metrics": metrics,
                "rmsd": rmsd,
                "best_score": best_score
            })

        runtime = time.time() - start
        valid_variants = [v for v in variant_results if v["output_file"] is not None]
        if not valid_variants:
            raise ValueError("Docking failed for all ligand variants")

        select_by = single_cfg.get("variant_select_by", "rmsd")
        if select_by == "score":
            best_variant = min(valid_variants, key=lambda v: v["best_score"] if v["best_score"] is not None else float("inf"))
        else:
            best_variant = min(valid_variants, key=lambda v: v["rmsd"])

        output_file = best_variant["output_file"]
        metrics = best_variant["metrics"] or {}
        rmsd = best_variant["rmsd"]
        best_score = best_variant["best_score"]

        best_score_value = metrics.get("best_score")
        if best_score_value is None:
            best_score_value = best_score

        return RedockResult(
            pdb_id=pdb_file.stem.upper(),
            ligand_resname=ligand_resname or ligand_name,
            ligand_chain=ligand_chain or "",
            mode=run_mode,
            engine=engine_name,
            protocol="single",
            best_rmsd=rmsd,
            success=rmsd < threshold,
            runtime_sec=runtime,
            output_file=str(output_file),
            pose_count=metrics.get("pose_count"),
            best_score=best_score_value,
            dock_name=ligand_name,
            top1_rmsd=metrics.get("top1_rmsd"),
            top5_rmsd=metrics.get("top5_rmsd"),
            top10_rmsd=metrics.get("top10_rmsd"),
            best_rmsd_rank=metrics.get("best_rmsd_rank"),
            rmsd_best_score=metrics.get("rmsd_best_score"),
            rmsd_mean=metrics.get("rmsd_mean"),
            rmsd_median=metrics.get("rmsd_median"),
            rmsd_std=metrics.get("rmsd_std"),
            near_native_fraction=metrics.get("near_native_fraction"),
            score_rmsd_pearson=metrics.get("score_rmsd_pearson"),
            score_rmsd_spearman=metrics.get("score_rmsd_spearman"),
            error_message=None,
            site_method=site_method,
            docking_completed=True,
            variants_prepared=variants_prepared,
            variants_docked=variants_docked
        )

    @staticmethod
    def _ordered_results(
        pairs: List[Dict[str, str]], results_by_case: Dict[str, RedockResult]
    ) -> List[RedockResult]:
        ordered = []
        for pair in pairs:
            case_id = pair.get("case_id") or f"{pair['pdb_id']}_{pair.get('dock_name') or pair['ligand']}"
            if case_id in results_by_case:
                ordered.append(results_by_case[case_id])
        return ordered

    def _load_resumable_results(
        self, manifest_path: Path, progress_path: Path, current_manifest: dict
    ) -> List[RedockResult]:
        if not manifest_path.exists() or not progress_path.exists():
            return []
        try:
            previous_manifest = json.loads(manifest_path.read_text())
            expected = self._json_normalize({
                "config": current_manifest.get("config"),
                "cases": current_manifest.get("cases"),
            })
            previous = {
                "config": previous_manifest.get("config"),
                "cases": previous_manifest.get("cases"),
            }
            if previous != expected:
                logger.info("Existing progress is incompatible with this run; starting fresh")
                return []

            payload = json.loads(progress_path.read_text())
            cases = current_manifest.get("cases", [])
            resumed = []
            for item in payload.get("results", []):
                result = RedockResult(**item)
                if not result.case_id:
                    matches = [
                        case for case in cases
                        if case.get("pdb_id") == result.pdb_id
                        and case.get("site_ligand") == result.ligand_resname
                        and case.get("dock_name") == result.dock_name
                        and case.get("control_label") == result.control_label
                    ]
                    if len(matches) == 1:
                        result.case_id = matches[0].get("case_id")
                output_exists = bool(result.output_file and Path(result.output_file).exists())
                completed = result.docking_completed is True or (
                    result.docking_completed is None
                    and output_exists
                    and not result.error_message
                )
                if result.case_id and completed and output_exists:
                    result.docking_completed = True
                    resumed.append(result)
            return resumed
        except Exception as exc:
            logger.warning("Could not resume existing progress: {}", exc)
            return []

    @classmethod
    def _json_normalize(cls, value: object):
        return json.loads(json.dumps(value, default=cls._json_default))

    @classmethod
    def _write_json_atomic(cls, path: Path, payload: object) -> None:
        temporary = path.with_name(f".{path.name}.tmp")
        temporary.write_text(json.dumps(payload, indent=2, default=cls._json_default))
        temporary.replace(path)

    @staticmethod
    def _write_text_atomic(path: Path, text: str) -> None:
        temporary = path.with_name(f".{path.name}.tmp")
        temporary.write_text(text)
        temporary.replace(path)

    @staticmethod
    def _command_version(command: List[str]) -> Optional[str]:
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=20,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        output = "\n".join(
            part.strip()
            for part in (completed.stdout, completed.stderr)
            if part and part.strip()
        )
        return output or None

    @classmethod
    def _software_provenance(cls, config: dict) -> dict:
        """Capture code, dependency, and selected binary versions for a run."""
        try:
            package_version = importlib.metadata.version("docking_platform_gui")
        except importlib.metadata.PackageNotFoundError:
            package_version = "development"

        repository_root = Path(__file__).resolve().parents[2]
        git_commit = cls._command_version(
            ["git", "-C", str(repository_root), "rev-parse", "HEAD"]
        )
        git_status = cls._command_version(
            ["git", "-C", str(repository_root), "status", "--porcelain"]
        )
        dependencies = {}
        for distribution in (
            "numpy", "pandas", "openpyxl", "rdkit", "loguru",
            "biopython", "pydantic", "MDAnalysis", "gemmi", "scipy",
            "scikit-learn", "matplotlib", "seaborn", "openmm",
            "pdbfixer", "openbabel",
        ):
            try:
                dependencies[distribution] = importlib.metadata.version(distribution)
            except importlib.metadata.PackageNotFoundError:
                dependencies[distribution] = None

        single = config.get("single", {})
        binaries = {}

        def short_version(command: List[str]) -> Optional[str]:
            output = cls._command_version(command)
            return "\n".join(output.splitlines()[:5]) if output else None

        for name, key in (("vina", "vina_binary"), ("smina", "smina_binary")):
            raw_path = single.get(key)
            if not raw_path:
                continue
            binary_path = str(Path(raw_path).expanduser())
            binaries[name] = {
                "path": binary_path,
                "version": short_version([binary_path, "--version"]),
            }
        rdock_root = single.get("rdock_root")
        if rdock_root:
            rdock_binary = Path(rdock_root).expanduser() / "bin/rbdock"
            if rdock_binary.exists():
                binaries["rdock"] = {
                    "path": str(rdock_binary),
                    "version": short_version([str(rdock_binary), "-h"]),
                }

        return {
            "package": "docking_platform_gui",
            "version": package_version,
            "git_commit": git_commit.splitlines()[0] if git_commit else None,
            "git_dirty": bool(git_status),
            "python": sys.version,
            "platform": platform.platform(),
            "dependencies": dependencies,
            "external_binaries": binaries,
        }

    def _write_summary_files(self, results_path: Path, summary: dict) -> None:
        """Keep the machine-readable and Markdown summaries in sync."""
        summary_path = results_path.with_name("redock_summary.json")
        self._write_json_atomic(summary_path, summary)
        self._write_text_atomic(
            results_path.with_name("redock_summary.md"),
            self._summary_to_markdown(summary),
        )

    def _write_progress(self, progress_path: Path, results: List[RedockResult]) -> None:
        payload = {"results": [asdict(r) for r in results]}
        self._write_json_atomic(progress_path, payload)

    @staticmethod
    def _json_default(value: object):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")

    def _write_results(
        self,
        json_path: Path,
        csv_path: Path,
        results: List[RedockResult],
        threshold: float
    ) -> None:
        payload = {"results": [asdict(r) for r in results]}
        self._write_json_atomic(json_path, payload)

        if results:
            df = pd.DataFrame([asdict(r) for r in results])
            df.to_csv(csv_path, index=False)

        summary = self._build_summary(results, threshold)
        self._write_summary_files(json_path, summary)

    def _resolve_results_path(self, allow_csv: bool = False) -> Optional[Path]:
        if self.last_results_path and self.last_results_path.exists():
            return self.last_results_path

        raw = self.output_var.get().strip()
        if not raw:
            return None
        output_path = Path(raw).expanduser()
        if output_path.is_file():
            if output_path.name == "redock_results.json":
                return output_path
            if allow_csv and output_path.name == "redock_results.csv":
                return output_path
            return None

        candidate = output_path / "redock_results.json"
        if candidate.exists():
            return candidate
        if allow_csv:
            csv_candidate = output_path / "redock_results.csv"
            if csv_candidate.exists():
                return csv_candidate
            direct_protocol_candidate = output_path / "protocol_development_results.csv"
            if direct_protocol_candidate.exists():
                return direct_protocol_candidate
            protocol_candidate = (
                output_path / "protocol_development" / "protocol_development_results.csv"
            )
            if protocol_candidate.exists():
                return protocol_candidate
        return None

    def _load_last_results(self) -> None:
        results_path = self._resolve_results_path(allow_csv=True)
        if not results_path or not results_path.exists():
            messagebox.showwarning(
                "Results missing",
                "No results file found in the output directory."
            )
            return
        self.last_results_path = results_path
        if results_path.name == "protocol_development_results.csv":
            report_path = results_path.with_name("protocol_development_summary.md")
            self._render_protocol_results(results_path, report_path)
        else:
            self._render_results_from_path(results_path)
        self._set_status(f"Loaded results from {results_path.parent}")

    @staticmethod
    def _result_file_for_selection(
        selected_path: Path,
        preferred_mode: Optional[str] = None,
    ) -> Optional[Path]:
        """Resolve a selected run directory or exact results file."""
        selected_path = Path(selected_path).expanduser()
        valid_names = {
            "redock_results.json",
            "redock_results.csv",
            "protocol_development_results.csv",
        }
        if selected_path.is_file():
            return selected_path if selected_path.name in valid_names else None
        if not selected_path.is_dir():
            return None

        screening_candidates = [
            selected_path / "redock_results.json",
            selected_path / "redock_results.csv",
        ]
        protocol_candidates = [
            selected_path / "protocol_development_results.csv",
            selected_path / "protocol_development" / "protocol_development_results.csv",
        ]
        candidates = (
            protocol_candidates + screening_candidates
            if preferred_mode == "protocol_development"
            else screening_candidates + protocol_candidates
        )
        return next((candidate for candidate in candidates if candidate.is_file()), None)

    def _load_results_selection(self, selected_path: Path) -> None:
        """Render results discovered from a user-selected folder or file."""
        preferred_mode = self.mode_var.get() if hasattr(self, "mode_var") else None
        results_path = self._result_file_for_selection(selected_path, preferred_mode)
        if results_path is None:
            messagebox.showerror(
                "Results not found",
                "No supported results were found in that selection.\n\n"
                "For Screening, select the completed run folder containing "
                "redock_results.json (or redock_results.csv).\n\n"
                "For Protocol Development, select either the completed run folder "
                "containing protocol_development/, or the protocol_development "
                "folder itself.\n\n"
                "Do not select the parent output folder when it contains multiple runs.",
            )
            return

        self.last_results_path = results_path
        if results_path.name == "protocol_development_results.csv":
            self._render_protocol_results(
                results_path, results_path.with_name("protocol_development_summary.md")
            )
            result_type = "Protocol Development"
        else:
            self._render_results_from_path(results_path)
            result_type = "Screening"
        self._set_status(f"Loaded {result_type} results from {results_path.parent}")

    def _browse_results_folder(self) -> None:
        """Load a campaign by selecting its run folder."""
        initial = Path(self.output_var.get()).expanduser()
        initial_dir = initial if initial.is_dir() else initial.parent
        selected = filedialog.askdirectory(
            title="Select one completed docking run folder",
            initialdir=str(initial_dir) if initial_dir.exists() else None,
            mustexist=True,
        )
        if selected:
            self._load_results_selection(Path(selected))

    def _browse_results_file(self) -> None:
        """Load a campaign by selecting its exact results CSV or JSON file."""
        initial = Path(self.output_var.get()).expanduser()
        initial_dir = initial if initial.is_dir() else initial.parent
        selected = filedialog.askopenfilename(
            title="Select docking results",
            initialdir=str(initial_dir) if initial_dir.exists() else None,
            filetypes=[
                ("Docking results", ("*.json", "*.csv")),
                ("JSON files", "*.json"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        if not selected:
            return
        self._load_results_selection(Path(selected))

    def _open_pose_viewer_from_last(self) -> None:
        results_path = self._resolve_results_path(allow_csv=True)
        if not results_path or not results_path.exists():
            messagebox.showwarning(
                "Pose viewer",
                "No results file found in the output directory."
            )
            return
        self.last_results_path = results_path
        self._show_pose_viewer(results_path)

    def _show_results(self, results_path: Path) -> None:
        summary_path = Path(results_path).with_name("redock_summary.json")
        csv_path = Path(results_path).with_name("redock_results.csv")
        summary = self._summary_for_display(Path(results_path), summary_path)
        if not summary:
            messagebox.showwarning("Results missing", "Results could not be summarized.")
            return
        rmsd_values = []
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            if "best_rmsd" in df.columns:
                rmsd_values = [v for v in df["best_rmsd"].tolist() if isinstance(v, (int, float)) and v < 900]

        dialog = tk.Toplevel(self)
        dialog.title("Docking Analysis Results")
        dialog.geometry("900x700")
        dialog.transient(self)

        container = tk.Frame(dialog)
        container.pack(fill="both", expand=True)
        canvas = tk.Canvas(container, borderwidth=0)
        vscroll = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        hscroll = ttk.Scrollbar(container, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vscroll.set, xscrollcommand=hscroll.set)

        vscroll.pack(side="right", fill="y")
        hscroll.pack(side="bottom", fill="x")
        canvas.pack(side="left", fill="both", expand=True)

        content = tk.Frame(canvas)
        window_id = canvas.create_window((0, 0), window=content, anchor="nw")

        def _on_frame_configure(_event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_configure(event):
            canvas.itemconfigure(window_id, width=event.width)

        content.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_configure)

        notebook = ttk.Notebook(content)
        notebook.pack(fill="both", expand=True, padx=10, pady=10)

        summary_frame = tk.Frame(notebook)
        charts_frame = tk.Frame(notebook)
        notebook.add(summary_frame, text="Summary")
        notebook.add(charts_frame, text="Charts")

        self._populate_summary_tab(summary_frame, summary)
        self._populate_charts_tab(charts_frame, summary, rmsd_values)

    def _render_protocol_results(self, results_path: Path, report_path: Path) -> None:
        """Render a protocol-development report in the main Results card."""
        self._clear_frame(self.results_summary_tab)
        self._clear_frame(self.results_charts_tab)
        self._populate_protocol_report(
            self.results_summary_tab, Path(results_path), Path(report_path)
        )
        self._populate_protocol_charts(self.results_charts_tab, Path(results_path))
        self.results_notebook.select(self.results_summary_tab)

    def _show_protocol_results(self, results_path: Path, report_path: Path) -> None:
        """Open the completion summary window for a protocol-development run."""
        dialog = tk.Toplevel(self)
        dialog.title("Protocol Development Results")
        dialog.geometry("1100x720")
        dialog.transient(self)

        content = tk.Frame(dialog, padx=10, pady=10)
        content.pack(fill="both", expand=True)
        notebook = ttk.Notebook(content)
        notebook.pack(fill="both", expand=True)
        summary_tab = tk.Frame(notebook, padx=8, pady=8)
        charts_tab = tk.Frame(notebook, padx=8, pady=8)
        notebook.add(summary_tab, text="Summary")
        notebook.add(charts_tab, text="Charts")
        self._populate_protocol_report(summary_tab, Path(results_path), Path(report_path))
        self._populate_protocol_charts(charts_tab, Path(results_path))

        actions = tk.Frame(content)
        actions.pack(fill="x", pady=(8, 0))
        tk.Button(
            actions,
            text="Pose Viewer",
            command=self._safe_call(lambda: self._show_pose_viewer(Path(results_path))),
        ).pack(side="left")
        tk.Button(actions, text="Close", command=dialog.destroy).pack(side="right")
        dialog.lift()
        dialog.focus_force()

    @staticmethod
    def _parse_protocol_markdown(report: str) -> Tuple[List[str], List[str], List[List[str]]]:
        """Extract prose and the comparison table from a generated report."""
        prose: List[str] = []
        headers: List[str] = []
        rows: List[List[str]] = []
        table_started = False
        for raw_line in report.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("|") and line.endswith("|"):
                cells = [cell.strip() for cell in line.strip("|").split("|")]
                if not headers:
                    headers = cells
                    table_started = True
                    continue
                if all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells):
                    continue
                rows.append(cells)
                continue
            if not table_started:
                prose.append(line.lstrip("# "))
        return prose, headers, rows

    @staticmethod
    def _parse_protocol_markdown_sections(
        report: str,
    ) -> Tuple[List[str], List[Tuple[str, List[str], List[List[str]]]]]:
        """Parse all Markdown tables so the GUI can place each in its own tab."""
        prose: List[str] = []
        tables: List[Tuple[str, List[str], List[List[str]]]] = []
        section = "Protocol comparison"
        headers: List[str] = []
        rows: List[List[str]] = []

        def _finish_table() -> None:
            nonlocal headers, rows
            if headers:
                tables.append((section, headers, rows))
                headers, rows = [], []

        for raw_line in report.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("## "):
                _finish_table()
                section = line[3:].strip()
                continue
            if line.startswith("|") and line.endswith("|"):
                cells = [cell.strip() for cell in line.strip("|").split("|")]
                if not headers:
                    headers = cells
                elif not all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells):
                    rows.append(cells)
                continue
            if not headers:
                cleaned = line.lstrip("# ")
                if cleaned and cleaned not in prose:
                    prose.append(cleaned)
        _finish_table()
        return prose, tables

    @staticmethod
    def _populate_protocol_report(parent: tk.Widget, results_path: Path, report_path: Path) -> None:
        if report_path.exists():
            report = report_path.read_text()
        else:
            report = f"Protocol summary was not found:\n{report_path}"
        RedockAnalysisApp._populate_markdown_report(
            parent, report, source_text=f"Results: {results_path}"
        )

    @staticmethod
    def _populate_markdown_report(
        parent: tk.Widget, report: str, source_text: Optional[str] = None
    ) -> None:
        """Render generated Markdown as labels and sortable-looking tables."""
        if source_text:
            tk.Label(
                parent,
                text=source_text,
                anchor="w",
                justify="left",
                fg="#555555",
            ).pack(fill="x", pady=(0, 6))
        prose, tables = RedockAnalysisApp._parse_protocol_markdown_sections(report)
        for index, line in enumerate(prose):
            if index == 0:
                tk.Label(
                    parent, text=line, anchor="w", font=("TkDefaultFont", 13, "bold")
                ).pack(fill="x", pady=(0, 5))
            else:
                is_bullet = line.startswith("- ")
                rendered = line[2:].strip() if is_bullet else line
                tk.Label(
                    parent, text=rendered, anchor="w", justify="left", wraplength=1050
                ).pack(fill="x", padx=(14 if is_bullet else 0, 0), pady=(0, 3))

        if not tables:
            if not prose:
                tk.Label(parent, text="No summary data available.", anchor="nw").pack(
                    fill="both", expand=True
                )
            return

        notebook = ttk.Notebook(parent)
        notebook.pack(fill="both", expand=True, pady=(8, 0))
        numeric_headers = {
            "Exhaust.", "Exhaustiveness", "Mean best RMSD", "Mean Top-1 RMSD",
            "Mean best-pose rank", "Failures", "Runtime (s)", "Cases", "Conditions",
            "Conditions (with RMSD)",
        }
        for table_index, (title, headers, rows) in enumerate(tables):
            frame = tk.Frame(notebook)
            notebook.add(frame, text=title[:32])
            column_ids = [f"table_{table_index}_column_{index}" for index in range(len(headers))]
            table = ttk.Treeview(
                frame, columns=column_ids, show="headings", height=max(6, min(25, len(rows)))
            )
            vertical = ttk.Scrollbar(frame, orient="vertical", command=table.yview)
            horizontal = ttk.Scrollbar(frame, orient="horizontal", command=table.xview)
            table.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
            for column_id, heading in zip(column_ids, headers):
                table.heading(column_id, text=heading)
                width = max(90, min(230, 8 * len(heading) + 24))
                table.column(
                    column_id, width=width, minwidth=75,
                    anchor="e" if heading in numeric_headers else "center", stretch=False,
                )
            for row in rows:
                padded = row + [""] * (len(headers) - len(row))
                table.insert("", "end", values=padded[:len(headers)])
            vertical.pack(side="right", fill="y")
            horizontal.pack(side="bottom", fill="x")
            table.pack(side="left", fill="both", expand=True)

    def _render_results_from_path(self, results_path: Path) -> None:
        summary_path = Path(results_path).with_name("redock_summary.json")
        csv_path = Path(results_path).with_name("redock_results.csv")
        summary = self._summary_for_display(Path(results_path), summary_path)
        if not summary:
            messagebox.showwarning("Results missing", "Results could not be summarized.")
            return
        rmsd_values = []
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            if "best_rmsd" in df.columns:
                rmsd_values = [
                    v for v in df["best_rmsd"].tolist()
                    if isinstance(v, (int, float)) and v < 900
                ]

        self._render_results(summary, rmsd_values)

    def _summary_for_display(self, results_path: Path, summary_path: Path) -> dict:
        """Rebuild and persist metrics so older runs use current reporting."""
        saved_summary = {}
        if summary_path.exists():
            try:
                saved_summary = json.loads(summary_path.read_text())
            except Exception as exc:
                logger.warning("Could not read saved summary: {}", exc)
        json_path = (
            results_path.with_name("redock_results.json")
            if results_path.suffix.lower() == ".csv" else results_path
        )
        if not json_path.exists():
            return saved_summary
        try:
            payload = json.loads(json_path.read_text())
            results = [RedockResult(**item) for item in payload.get("results", [])]
            if not results:
                return saved_summary
            rebuilt_summary = self._build_summary(
                results, float(saved_summary.get("threshold", 2.0))
            )
            try:
                self._write_summary_files(json_path, rebuilt_summary)
            except Exception as exc:
                logger.warning("Could not update saved summary files: {}", exc)
            return rebuilt_summary
        except Exception as exc:
            logger.warning("Could not rebuild summary for display: {}", exc)
            return saved_summary

    @staticmethod
    def _protocol_chart_data(
        frame: pd.DataFrame, threshold: float = 2.0
    ) -> dict:
        """Build continuous-metric chart data for a factorial protocol sweep."""
        empty = {
            "pose_ranking_points": [], "rescore_points": [],
            "runtime_points": [], "factor_effects": [],
            "source_rows": len(frame), "condition_count": 0,
            "complex_count": 0, "top1_success": 0, "top5_success": 0,
            "best_pose_success": 0, "rescore_improved": 0,
            "rescore_unchanged": 0, "rescore_worse": 0,
            "median_rescore_delta": None,
        }
        if frame.empty or "status" not in frame:
            return empty

        working = frame.copy()
        condition_columns = (
            "engine", "box_definition", "rescore_method", "water_handling",
            "exhaustiveness",
        )
        defaults = {
            "engine": "unknown", "box_definition": "unknown",
            "rescore_method": "none", "water_handling": "unknown",
            "exhaustiveness": "N/A", "seed": "N/A",
        }
        for column in condition_columns:
            if column not in working:
                working[column] = defaults[column]
            working[column] = working[column].fillna(defaults[column])
        if "seed" not in working:
            working["seed"] = defaults["seed"]
        working["seed"] = working["seed"].fillna(defaults["seed"])

        metric_columns = (
            "best_rmsd", "top1_rmsd", "top5_rmsd", "rescore_top1_rmsd",
            "rescore_top5_rmsd", "runtime_sec",
        )
        for column in metric_columns:
            if column not in working:
                working[column] = np.nan
            working[column] = pd.to_numeric(working[column], errors="coerce")
            if "rmsd" in column:
                working.loc[working[column] >= 900, column] = np.nan

        working["_status"] = working["status"].astype(str).str.lower()
        working = working[
            ~working["_status"].str.contains("unsupported|skipped", na=False)
        ]
        if working.empty:
            return empty

        # Older runs repeated identical rDock work for each Vina
        # exhaustiveness value even though rDock does not use that parameter.
        working.loc[
            working["engine"].astype(str).str.lower() == "rdock", "exhaustiveness"
        ] = "N/A"
        identity_columns = [
            column for column in (
                "pdb_id", "ligand_resname", "ligand_chain", "target_name",
                "engine", "box_definition", "rescore_method", "water_handling",
                "exhaustiveness", "seed",
            )
            if column in working.columns
        ]
        if identity_columns:
            rdock_mask = working["engine"].astype(str).str.lower() == "rdock"
            working = pd.concat(
                [
                    working[~rdock_mask],
                    working[rdock_mask].drop_duplicates(
                        subset=identity_columns, keep="last"
                    ),
                ],
                axis=0,
            ).sort_index()

        def _compact(value: object) -> str:
            if isinstance(value, (int, np.integer)):
                return str(value)
            if isinstance(value, (float, np.floating)) and np.isfinite(value):
                return str(int(value)) if float(value).is_integer() else f"{value:g}"
            return str(value).replace("_", " ")

        def _condition_label(key: tuple, seed_label: str) -> str:
            engine, box, rescore, water, exhaustiveness = key
            return " | ".join((
                _compact(engine), _compact(box), _compact(water),
                f"e{_compact(exhaustiveness)}", seed_label,
                f"rescore:{_compact(rescore)}",
            ))

        records = []
        grouped = working.groupby(list(condition_columns), dropna=False, sort=False)
        for raw_key, all_rows in grouped:
            key = raw_key if isinstance(raw_key, tuple) else (raw_key,)
            complete = all_rows[all_rows["_status"] == "complete"]
            method = str(key[2]).strip().lower()
            has_rescore = method not in {"", "none", "nan", "n/a"}
            selected_top1 = "rescore_top1_rmsd" if has_rescore else "top1_rmsd"
            selected_top5 = "rescore_top5_rmsd" if has_rescore else "top5_rmsd"
            runtime_values = complete["runtime_sec"].dropna()
            seeds = sorted({_compact(value) for value in all_rows["seed"].dropna().tolist()})
            seed_label = f"s{seeds[0]}" if len(seeds) == 1 else f"{len(seeds)} seeds"
            selected_top1_values = complete[selected_top1].dropna()
            best_values = complete["best_rmsd"].dropna()
            baseline_values = complete["top1_rmsd"].dropna()
            rescored_values = complete["rescore_top1_rmsd"].dropna()
            selected_top5_values = complete[selected_top5].dropna()
            label = _condition_label(key, seed_label)
            records.append({
                "key": key,
                "label": label,
                "engine": str(key[0]),
                "has_rescore": has_rescore,
                "rows": complete,
                "completion": 100.0 * len(complete) / len(all_rows),
                "top1_rmsd": (
                    float(selected_top1_values.median())
                    if not selected_top1_values.empty else None
                ),
                "top5_rmsd": (
                    float(selected_top5_values.median())
                    if not selected_top5_values.empty else None
                ),
                "best_rmsd": float(best_values.median()) if not best_values.empty else None,
                "baseline_top1_rmsd": (
                    float(baseline_values.median()) if not baseline_values.empty else None
                ),
                "rescored_top1_rmsd": (
                    float(rescored_values.median())
                    if has_rescore and not rescored_values.empty else None
                ),
                "runtime": float(runtime_values.median()) if not runtime_values.empty else None,
            })

        pose_ranking_points = []
        rescore_points = []
        runtime_points = []
        for row in records:
            common = {
                "engine": row["engine"], "rescore": row["has_rescore"],
                "label": row["label"],
            }
            if row["best_rmsd"] is not None and row["top1_rmsd"] is not None:
                pose_ranking_points.append({
                    **common, "x": row["best_rmsd"], "y": row["top1_rmsd"],
                    "tooltip": (
                        f"{row['label']}\nBest pose: {row['best_rmsd']:.2f} A; "
                        f"selected Top-1: {row['top1_rmsd']:.2f} A"
                    ),
                })
            if (
                row["has_rescore"]
                and row["baseline_top1_rmsd"] is not None
                and row["rescored_top1_rmsd"] is not None
            ):
                delta = row["rescored_top1_rmsd"] - row["baseline_top1_rmsd"]
                rescore_points.append({
                    **common,
                    "x": row["baseline_top1_rmsd"],
                    "y": row["rescored_top1_rmsd"],
                    "delta": delta,
                    "tooltip": (
                        f"{row['label']}\nDocking: {row['baseline_top1_rmsd']:.2f} A; "
                        f"rescored: {row['rescored_top1_rmsd']:.2f} A; "
                        f"change: {delta:+.2f} A"
                    ),
                })
            if row["runtime"] is not None and row["top1_rmsd"] is not None:
                runtime_points.append({
                    **common, "x": row["runtime"], "y": row["top1_rmsd"],
                    "tooltip": (
                        f"{row['label']}\nRuntime: {row['runtime']:.1f} s; "
                        f"selected Top-1: {row['top1_rmsd']:.2f} A"
                    ),
                })

        for point in runtime_points:
            point["pareto"] = not any(
                other["x"] <= point["x"]
                and other["y"] <= point["y"]
                and (other["x"] < point["x"] or other["y"] < point["y"])
                for other in runtime_points
            )

        factor_specs = (
            ("Engine", 0), ("Box", 1), ("Water", 3),
            ("Exhaustiveness", 4), ("Ranking", 2),
        )
        factor_values: Dict[Tuple[str, str], List[float]] = {}
        for row in records:
            if row["top1_rmsd"] is None:
                continue
            for factor, key_index in factor_specs:
                raw_level = row["key"][key_index]
                if factor == "Exhaustiveness" and str(raw_level).lower() in {"n/a", "nan"}:
                    continue
                if factor == "Ranking":
                    level = "rescored" if row["has_rescore"] else "docking"
                else:
                    level = _compact(raw_level)
                factor_values.setdefault((factor, level), []).append(row["top1_rmsd"])

        factor_effects = []
        preferred_order = {
            "Engine": {"smina": 0, "vina": 1, "rdock": 2},
            "Water": {"remove all": 0, "selective": 1, "retain all": 2},
            "Ranking": {"docking": 0, "rescored": 1},
        }

        def _factor_level_key(factor: str, level: str) -> Tuple[int, object]:
            if factor == "Exhaustiveness":
                try:
                    return 0, float(level)
                except (TypeError, ValueError):
                    return 1, level
            order = preferred_order.get(factor, {})
            return (0, order[level]) if level in order else (1, str(level))

        for factor, _ in factor_specs:
            levels = sorted(
                [
                    (level, values)
                    for (record_factor, level), values in factor_values.items()
                    if record_factor == factor
                ],
                key=lambda item, current_factor=factor:
                _factor_level_key(current_factor, item[0]),
            )
            for level, values in levels:
                array = np.asarray(values, dtype=float)
                factor_effects.append({
                    "factor": factor, "label": f"{factor}: {level}",
                    "minimum": float(np.min(array)),
                    "q1": float(np.percentile(array, 25)),
                    "median": float(np.median(array)),
                    "q3": float(np.percentile(array, 75)),
                    "maximum": float(np.max(array)), "n": len(array),
                })

        complex_columns = [
            column for column in ("pdb_id", "ligand_resname", "ligand_chain")
            if column in working.columns
        ]
        complex_count = (
            working[complex_columns].drop_duplicates().shape[0]
            if complex_columns else 0
        )
        deltas = np.asarray([point["delta"] for point in rescore_points], dtype=float)
        return {
            "pose_ranking_points": pose_ranking_points,
            "rescore_points": rescore_points,
            "runtime_points": runtime_points,
            "factor_effects": factor_effects,
            "source_rows": len(frame),
            "condition_count": len(records),
            "complex_count": complex_count,
            "top1_success": sum(
                row["top1_rmsd"] is not None and row["top1_rmsd"] < threshold
                for row in records
            ),
            "top5_success": sum(
                row["top5_rmsd"] is not None and row["top5_rmsd"] < threshold
                for row in records
            ),
            "best_pose_success": sum(
                row["best_rmsd"] is not None and row["best_rmsd"] < threshold
                for row in records
            ),
            "rescore_improved": int((deltas < -1e-6).sum()) if deltas.size else 0,
            "rescore_unchanged": int((np.abs(deltas) <= 1e-6).sum()) if deltas.size else 0,
            "rescore_worse": int((deltas > 1e-6).sum()) if deltas.size else 0,
            "median_rescore_delta": float(np.median(deltas)) if deltas.size else None,
        }

    def _populate_protocol_charts(self, parent: tk.Frame, results_path: Path) -> None:
        """Render protocol charts from the selected CSV, not stale summary state."""
        if not results_path.exists():
            tk.Label(parent, text="Protocol results CSV not found.").pack(anchor="w", padx=10, pady=10)
            return
        try:
            frame = pd.read_csv(results_path)
        except Exception as exc:
            tk.Label(parent, text=f"Could not load chart data: {exc}").pack(anchor="w", padx=10, pady=10)
            return
        data = self._protocol_chart_data(frame)
        condition_count = data["condition_count"]
        note = (
            f"{data['source_rows']} saved rows represent {condition_count} distinct "
            f"protocol conditions across {data['complex_count']} crystal complex(es). "
            f"At 2 A: Top-1 {data['top1_success']}/{condition_count}, "
            f"Top-5 {data['top5_success']}/{condition_count}, and a generated pose "
            f"{data['best_pose_success']}/{condition_count}."
        )
        if data["rescore_points"]:
            note += (
                f" Rescoring improved {data['rescore_improved']}, left "
                f"{data['rescore_unchanged']} unchanged, and worsened "
                f"{data['rescore_worse']} conditions; median Top-1 change "
                f"{data['median_rescore_delta']:+.2f} A."
            )
        tk.Label(
            parent, text=note, anchor="w", justify="left", wraplength=1050,
            fg="#555555",
        ).grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=(10, 0))
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(1, weight=1)
        parent.grid_rowconfigure(2, weight=1)

        pose_canvas = tk.Canvas(parent, height=310, bg="white", highlightthickness=1)
        pose_canvas.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)
        self._install_chart(
            pose_canvas,
            lambda: self._draw_protocol_scatter(
                pose_canvas,
                "Pose generation vs selected ranking",
                data["pose_ranking_points"],
                "Best generated-pose RMSD (A)",
                "Selected Top-1 RMSD (A)",
                x_threshold=2.0,
                y_threshold=2.0,
                aggregate_points=True,
                quadrant_labels=True,
            ),
        )

        rescore_canvas = tk.Canvas(parent, height=310, bg="white", highlightthickness=1)
        rescore_canvas.grid(row=1, column=1, sticky="nsew", padx=10, pady=10)
        self._install_chart(
            rescore_canvas,
            lambda: self._draw_protocol_scatter(
                rescore_canvas,
                "Paired Top-1 RMSD: docking vs rescoring",
                data["rescore_points"],
                "Docking-ranked Top-1 RMSD (A)",
                "Rescored Top-1 RMSD (A)",
                x_threshold=2.0,
                y_threshold=2.0,
                diagonal=True,
                equal_axes=True,
                aggregate_points=True,
            ),
        )

        runtime_canvas = tk.Canvas(parent, height=310, bg="white", highlightthickness=1)
        runtime_canvas.grid(row=2, column=0, sticky="nsew", padx=10, pady=10)
        self._install_chart(
            runtime_canvas,
            lambda: self._draw_protocol_scatter(
                runtime_canvas,
                "Accuracy-runtime trade-off",
                data["runtime_points"],
                "Recorded docking runtime (s, log scale)",
                "Selected Top-1 RMSD (A)",
                y_threshold=2.0,
                pareto=True,
                log_x=True,
            ),
        )

        effects_canvas = tk.Canvas(parent, height=310, bg="white", highlightthickness=1)
        effects_canvas.grid(row=2, column=1, sticky="nsew", padx=10, pady=10)
        self._install_chart(
            effects_canvas,
            lambda: self._draw_factor_effect_chart(
                effects_canvas,
                "Top-1 RMSD by protocol factor",
                data["factor_effects"],
                threshold=2.0,
            ),
        )

    def _render_results(self, summary: dict, rmsd_values: List[float]) -> None:
        self._clear_frame(self.results_summary_tab)
        self._clear_frame(self.results_charts_tab)
        self._populate_summary_tab(self.results_summary_tab, summary)
        self._populate_charts_tab(self.results_charts_tab, summary, rmsd_values)
        self.results_notebook.select(self.results_summary_tab)

    def _populate_empty_results(self) -> None:
        self._clear_frame(self.results_summary_tab)
        self._clear_frame(self.results_charts_tab)
        tk.Label(
            self.results_summary_tab,
            text="No results yet. Run an analysis or load existing results.",
            fg="#555555"
        ).pack(anchor="w", padx=10, pady=10)

    def _clear_frame(self, frame: tk.Widget) -> None:
        for child in frame.winfo_children():
            child.destroy()

    def _populate_summary_tab(self, parent: tk.Frame, summary: dict) -> None:
        self._populate_markdown_report(parent, self._summary_to_markdown(summary))

    @staticmethod
    def _result_is_completed(result: RedockResult) -> bool:
        output_exists = bool(result.output_file and Path(result.output_file).exists())
        return bool(
            result.docking_completed is True
            or (
                result.docking_completed is None
                and output_exists
                and not result.error_message
            )
        )

    @staticmethod
    def _result_is_failed(result: RedockResult) -> bool:
        output_exists = bool(result.output_file and Path(result.output_file).exists())
        return bool(
            result.docking_completed is False
            or (
                result.docking_completed is None
                and not output_exists
            )
        )

    @staticmethod
    def _screening_chart_data(summary: dict) -> dict:
        screening_rows = summary.get("per_structure_screening") or []
        enrichment_rows = summary.get("per_structure_enrichment") or []

        def _structure_label(row: dict) -> str:
            structure = f"{row.get('pdb_id', 'N/A')}/{row.get('ligand', 'N/A')}"
            target = str(row.get("target_name") or "").strip()
            return f"{target} | {structure}" if target and target != row.get("pdb_id") else structure

        structure_auc = []
        score_margin = []
        active_rank_percentile = []
        for row in enrichment_rows:
            rank = row.get("active_rank")
            actives = int(row.get("actives") or 1)
            decoys = int(row.get("decoys") or 0)
            total = actives + decoys
            rank_text = f"rank {rank}/{total}" if rank is not None else "rank N/A"
            label = f"{_structure_label(row)} ({rank_text})"
            if row.get("roc_auc") is not None:
                structure_auc.append((label, float(row["roc_auc"])))
            if row.get("score_margin") is not None:
                score_margin.append((label, float(row["score_margin"])))
            if rank is not None:
                denominator = max(total - 1, 1)
                percentile = 100.0 * max(0.0, total - float(rank)) / denominator
                active_rank_percentile.append((label, percentile))
        structure_auc.sort(key=lambda item: item[1], reverse=True)
        score_margin.sort(key=lambda item: item[1], reverse=True)
        active_rank_percentile.sort(key=lambda item: item[1], reverse=True)

        target_auc = []
        for row in summary.get("per_target_enrichment") or []:
            structures = int(row.get("structures") or 0)
            structure_text = "structure" if structures == 1 else "structures"
            target_auc.append((
                f"{row.get('target_name') or 'N/A'} ({structures} {structure_text})",
                [row.get("macro_roc_auc"), row.get("pooled_roc_auc")],
            ))
        target_auc.sort(
            key=lambda item: -(item[1][0] if item[1][0] is not None else -1.0)
        )

        hit_stats: Dict[str, Dict[str, float]] = {}
        for hit in summary.get("screening_top_hits") or []:
            compound = str(hit.get("compound") or "Unnamed")
            rank = int(hit.get("rank") or 0)
            if rank < 1:
                continue
            stats = hit_stats.setdefault(compound, {"rank1": 0.0, "top5": 0.0, "rr": 0.0})
            stats["top5"] += 1.0
            stats["rr"] += 1.0 / rank
            if rank == 1:
                stats["rank1"] += 1.0
        recurrence_rows = sorted(
            (
                (compound, stats["rank1"], stats["top5"], stats["rr"])
                for compound, stats in hit_stats.items() if stats["rank1"] > 0
            ),
            key=lambda item: (-item[1], -item[3], item[0].lower()),
        )
        top_hit_recurrence = [
            (f"{compound} (top-5 hits {int(top5)})", rank1)
            for compound, rank1, top5, _ in recurrence_rows[:12]
        ]

        direction = summary.get("screening_score_direction")
        score_advantage = []
        for row in screening_rows:
            best = row.get("best_score")
            median = row.get("median_score")
            if best is None or median is None:
                continue
            advantage = float(best) - float(median)
            if direction == "lower":
                advantage = -advantage
            score_advantage.append((_structure_label(row), advantage))
        score_advantage.sort(key=lambda item: item[1], reverse=True)

        early_enrichment = []
        for label, key in (("EF 1%", "ef_1_percent"), ("EF 5%", "ef_5_percent"), ("EF 10%", "ef_10_percent")):
            if summary.get(key) is not None:
                early_enrichment.append((label, float(summary[key])))

        campaign_coverage = []
        total_cases = int(summary.get("total_cases") or 0)
        if total_cases:
            campaign_coverage.append((
                "Docking completed",
                100.0 * float(summary.get("docking_completed") or 0) / total_cases,
            ))
        control_total = int(summary.get("n_actives") or 0) + int(summary.get("n_decoys") or 0)
        if control_total:
            control_scored = int(summary.get("control_actives") or 0) + int(
                summary.get("control_decoys") or 0
            )
            campaign_coverage.append((
                "Controls with ranking scores", 100.0 * control_scored / control_total,
            ))
        sample_total = int(summary.get("n_samples") or 0)
        if sample_total:
            campaign_coverage.append((
                "Unknowns with ranking scores",
                100.0 * float(summary.get("screening_score_count") or 0) / sample_total,
            ))

        completion = [
            (_structure_label(row), float(row["completion_rate"]))
            for row in screening_rows
        ]
        score_coverage = [
            (
                _structure_label(row),
                100.0 * float(row["scored"]) / float(row["cases"])
                if row.get("cases") else 0.0,
            )
            for row in screening_rows
        ]
        failures = [
            (_structure_label(row), float(row["cases"] - row["completed"]))
            for row in screening_rows
        ]
        return {
            "assay_benchmark": summary.get("assay_benchmark_charts") or {},
            "structure_auc": structure_auc[:15],
            "score_margin": score_margin[:15],
            "active_rank_percentile": active_rank_percentile[:15],
            "target_auc": target_auc[:15],
            "top_hit_recurrence": top_hit_recurrence,
            "score_advantage": score_advantage[:15],
            "early_enrichment": early_enrichment,
            "campaign_coverage": campaign_coverage,
            "has_coverage_gap": any(value < 99.999 for _, value in campaign_coverage),
            "completion": completion,
            "score_coverage": score_coverage,
            "failures": failures,
        }

    @staticmethod
    def _is_screening_summary(summary: dict) -> bool:
        """Recognize screening reports even when control ligands provide RMSDs."""
        return bool(
            summary.get("per_structure_screening")
            or summary.get("per_structure_enrichment")
            or summary.get("n_samples")
        )

    def _populate_charts_tab(self, parent: tk.Frame, summary: dict, rmsd_values: List[float]) -> None:
        if self._is_screening_summary(summary):
            chart_data = self._screening_chart_data(summary)
            parent.grid_rowconfigure(0, weight=0)
            parent.grid_rowconfigure(1, weight=1)
            parent.grid_rowconfigure(2, weight=1)
            parent.grid_columnconfigure(0, weight=1)
            parent.grid_columnconfigure(1, weight=1)

            assay_data = chart_data["assay_benchmark"]
            if (
                summary.get("enrichment_dataset_type") == "assay_benchmark"
                and assay_data.get("roc_curve")
            ):
                tk.Label(
                    parent,
                    text=(
                        f"Assay benchmark: {assay_data['actives']} actives and "
                        f"{assay_data['inactives']} inactives ranked at one receptor. "
                        "Curves and distributions use the same ranking score; higher "
                        "values rank better. Histogram classes are normalized separately."
                    ),
                    anchor="w",
                    justify="left",
                    wraplength=1050,
                    fg="#555555",
                ).grid(
                    row=0, column=0, columnspan=2, sticky="ew",
                    padx=10, pady=(10, 0),
                )

                auc = summary.get("roc_auc")
                ap = summary.get("average_precision")
                prevalence = assay_data.get("prevalence")
                ef_values = [
                    summary.get("ef_1_percent"),
                    summary.get("ef_5_percent"),
                    summary.get("ef_10_percent"),
                ]
                chart_specs = [
                    (
                        "curve",
                        f"ROC curve (AUC = {auc:.3f})" if auc is not None else "ROC curve",
                        assay_data["roc_curve"],
                        "False-positive rate", "True-positive rate",
                        True, None, (), "#177E89",
                    ),
                    (
                        "curve",
                        (
                            f"Precision-recall curve (AP = {ap:.3f})"
                            if ap is not None else "Precision-recall curve"
                        ),
                        assay_data["precision_recall_curve"],
                        "Recall", "Precision",
                        False, prevalence, (), "#D97732",
                    ),
                    (
                        "histogram", "Ranking-score distributions",
                        assay_data["score_histogram"],
                        "Ranking score (higher is better)", "Within-class frequency (%)",
                        False, None, (), "",
                    ),
                    (
                        "curve",
                        "Cumulative active recovery (EF 1/5/10% = "
                        + "/".join(
                            f"{value:.2f}" if value is not None else "N/A"
                            for value in ef_values
                        )
                        + ")",
                        assay_data["cumulative_recovery_curve"],
                        "Fraction of library screened", "Fraction of actives recovered",
                        True, None, (0.01, 0.05, 0.10), "#5B8E3E",
                    ),
                ]
                for index, spec in enumerate(chart_specs):
                    canvas = tk.Canvas(
                        parent, height=300, bg="white", highlightthickness=1
                    )
                    canvas.grid(
                        row=1 + index // 2, column=index % 2,
                        sticky="nsew", padx=10, pady=10,
                    )
                    if spec[0] == "histogram":
                        self._install_chart(
                            canvas,
                            lambda c=canvas, s=spec: self._draw_score_distribution_chart(
                                c, s[1], s[2], s[3], s[4]
                            ),
                        )
                    else:
                        self._install_chart(
                            canvas,
                            lambda c=canvas, s=spec: self._draw_xy_line_chart(
                                c, s[1], s[2], s[3], s[4],
                                reference_diagonal=s[5],
                                horizontal_reference=s[6],
                                vertical_guides=s[7],
                                color=s[8],
                            ),
                        )
                return

            tk.Label(
                parent,
                text=(
                    "Enrichment is assessed within each receptor structure. Positive score "
                    "margins mean the active beats the best decoy; target-pooled AUC is "
                    "diagnostic when several receptor structures are combined."
                    if chart_data["structure_auc"] else
                    "Unknown-compound scores are compared within each receptor structure; "
                    "raw scores should not be pooled across structures."
                ),
                anchor="w",
                justify="left",
                wraplength=1050,
                fg="#555555",
            ).grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=(10, 0))

            if chart_data["structure_auc"]:
                if chart_data["has_coverage_gap"]:
                    fourth = (
                        "Campaign completion and score coverage (%)",
                        "horizontal", chart_data["campaign_coverage"], (), 100.0, None,
                    )
                elif chart_data["top_hit_recurrence"]:
                    fourth = (
                        "Rank-1 compound recurrence across structures",
                        "horizontal", chart_data["top_hit_recurrence"], (), None, None,
                    )
                elif (
                    summary.get("enrichment_dataset_type") == "assay_benchmark"
                    and chart_data["early_enrichment"]
                ):
                    fourth = (
                        "Early enrichment factor",
                        "horizontal", chart_data["early_enrichment"], (), None, None,
                    )
                else:
                    fourth = (
                        "Active rank percentile (100% = rank 1)",
                        "horizontal", chart_data["active_rank_percentile"], (), 100.0, None,
                    )
                margin_title = (
                    "Best active vs best inactive score margin"
                    if summary.get("enrichment_dataset_type") == "assay_benchmark"
                    else "Reference active vs best-decoy score margin"
                )
                screening_charts = [
                    (
                        "ROC AUC by receptor structure",
                        "horizontal", chart_data["structure_auc"], (), 1.0,
                        [(0.5, "random"), (0.7, "review")],
                    ),
                    (
                        margin_title,
                        "diverging", chart_data["score_margin"], (), None, None,
                    ),
                    (
                        "Target enrichment: macro vs pooled AUC",
                        "grouped", chart_data["target_auc"],
                        ("Macro AUC", "Target-pooled AUC"), 1.0,
                        [(0.5, "random"), (0.7, "review")],
                    ),
                    fourth,
                ]
            else:
                screening_charts = [
                    (
                        "Best-score advantage over structure median",
                        "horizontal", chart_data["score_advantage"], (), None, None,
                    ),
                    (
                        "Rank-1 compound recurrence across structures",
                        "horizontal", chart_data["top_hit_recurrence"], (), None, None,
                    ),
                    (
                        "Completion rate by receptor structure (%)",
                        "horizontal", chart_data["completion"], (), 100.0, None,
                    ),
                    (
                        "Score coverage by receptor structure (%)",
                        "horizontal", chart_data["score_coverage"], (), 100.0, None,
                    ),
                ]

            for index, (title, chart_type, values, series, maximum, thresholds) in enumerate(screening_charts):
                canvas = tk.Canvas(parent, height=300, bg="white", highlightthickness=1)
                canvas.grid(
                    row=1 + index // 2, column=index % 2, sticky="nsew", padx=10, pady=10
                )
                if chart_type == "grouped":
                    self._install_chart(
                        canvas,
                        lambda c=canvas, t=title, v=values, s=series, m=maximum, h=thresholds:
                        self._draw_grouped_horizontal_chart(
                            c, t, v, s, max_value=m, thresholds=h
                        ),
                    )
                elif chart_type == "diverging":
                    self._install_chart(
                        canvas,
                        lambda c=canvas, t=title, v=values:
                        self._draw_diverging_horizontal_chart(c, t, v),
                    )
                else:
                    self._install_chart(
                        canvas,
                        lambda c=canvas, t=title, v=values, m=maximum, h=thresholds:
                        self._draw_horizontal_chart(
                            c, t, v, max_value=m, thresholds=h
                        ),
                    )
            return

        parent.grid_rowconfigure(1, weight=1)
        parent.grid_rowconfigure(2, weight=0)
        parent.grid_columnconfigure(0, weight=1)

        protocol_rates = []
        for proto, stats in summary.get("by_protocol", {}).items():
            count = stats.get("count") or 0
            success = stats.get("success") or 0
            rate = 100.0 * success / count if count else 0.0
            protocol_rates.append((proto, rate))

        chart_frame = tk.Frame(parent)
        chart_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        chart_frame.grid_columnconfigure(0, weight=1)

        canvas1 = tk.Canvas(chart_frame, height=240, bg="white")
        canvas1.grid(row=0, column=0, sticky="ew")
        self._draw_bar_chart(canvas1, "Success Rate by Protocol (%)", protocol_rates)

        hist_frame = tk.Frame(parent)
        hist_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)
        hist_frame.grid_columnconfigure(0, weight=1)
        hist_frame.grid_rowconfigure(0, weight=1)
        canvas2 = tk.Canvas(hist_frame, height=260, bg="white")
        canvas2.grid(row=0, column=0, sticky="nsew")
        self._draw_histogram(canvas2, "Best RMSD Distribution", rmsd_values)

        if summary.get("charge_count"):
            charge_frame = tk.Frame(parent)
            charge_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=(0, 10))
            charge_frame.grid_columnconfigure(0, weight=1)
            charge_frame.grid_columnconfigure(1, weight=1)

            overall = [
                ("Pos", summary.get("charge_frac_positive") or 0.0),
                ("Neutral", summary.get("charge_frac_neutral") or 0.0),
                ("Neg", summary.get("charge_frac_negative") or 0.0)
            ]
            top10 = [
                ("Pos", summary.get("charge_top10_frac_positive") or 0.0),
                ("Neutral", summary.get("charge_top10_frac_neutral") or 0.0),
                ("Neg", summary.get("charge_top10_frac_negative") or 0.0)
            ]

            canvas3 = tk.Canvas(charge_frame, height=220, bg="white")
            canvas3.grid(row=0, column=0, sticky="ew", padx=(0, 10))
            self._draw_bar_chart(canvas3, "Charge Distribution (%)", overall)

            canvas4 = tk.Canvas(charge_frame, height=220, bg="white")
            canvas4.grid(row=0, column=1, sticky="ew")
            self._draw_bar_chart(canvas4, "Charge Distribution Top-10% (%)", top10)

    def _build_table(self, parent: tk.Frame, columns: List[str], rows: List[Tuple]) -> ttk.Treeview:
        table = ttk.Treeview(parent, columns=columns, show="headings", height=6)
        for col in columns:
            table.heading(col, text=col)
            table.column(col, width=120, anchor="center")
        for row in rows:
            table.insert("", "end", values=row)
        return table

    @staticmethod
    def _install_chart(canvas: tk.Canvas, draw: Callable[[], None]) -> None:
        """Draw now and redraw when a result window is resized."""
        pending = {"job": None}

        def _draw() -> None:
            pending["job"] = None
            draw()

        def _redraw(_event=None) -> None:
            if pending["job"] is not None:
                canvas.after_cancel(pending["job"])
            pending["job"] = canvas.after_idle(_draw)

        canvas.bind("<Configure>", _redraw)
        _redraw()

    @staticmethod
    def _chart_dimensions(canvas: tk.Canvas) -> Tuple[int, int, int, int, int, int]:
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 420)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 260)
        label_width = min(245, max(145, int(width * 0.40)))
        left = label_width
        right = 58
        top = 54
        bottom = 28
        return width, height, left, right, top, bottom

    @staticmethod
    def _short_chart_label(label: str, limit: int = 48) -> str:
        cleaned = str(label).replace("\n", " ")
        return cleaned if len(cleaned) <= limit else f"{cleaned[:limit - 3]}..."

    @staticmethod
    def _protocol_engine_color(engine: str) -> str:
        return {
            "smina": "#177E89",
            "vina": "#D97732",
            "rdock": "#5B8E3E",
        }.get(str(engine).lower(), "#60758A")

    @staticmethod
    def _draw_chart_tooltip(
        canvas: tk.Canvas, event: tk.Event, text: str, width: int, height: int
    ) -> None:
        canvas.delete("chart_tooltip")
        x = min(max(event.x + 12, 8), max(width - 310, 8))
        y = min(max(event.y + 12, 8), max(height - 66, 8))
        text_id = canvas.create_text(
            x + 7, y + 6, text=text, anchor="nw", width=285,
            fill="#202020", font=("TkDefaultFont", 8), tags="chart_tooltip",
        )
        bounds = canvas.bbox(text_id)
        if bounds:
            rectangle = canvas.create_rectangle(
                bounds[0] - 5, bounds[1] - 4, bounds[2] + 5, bounds[3] + 4,
                fill="#FFFBEA", outline="#5D5A4F", tags="chart_tooltip",
            )
            canvas.tag_lower(rectangle, text_id)

    def _draw_protocol_scatter(
        self,
        canvas: tk.Canvas,
        title: str,
        points: List[dict],
        x_label: str,
        y_label: str,
        x_threshold: Optional[float] = None,
        y_threshold: Optional[float] = None,
        diagonal: bool = False,
        equal_axes: bool = False,
        pareto: bool = False,
        aggregate_points: bool = False,
        quadrant_labels: bool = False,
        log_x: bool = False,
    ) -> None:
        """Render a responsive condition-level scatter plot with hover details."""
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 460)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 280)
        left, right, top, bottom = 64, 24, 67, 48
        plot_width = max(100, width - left - right)
        plot_height = max(90, height - top - bottom)
        canvas.create_text(
            14, 12, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold")
        )
        if not points:
            canvas.create_text(
                left, height // 2, text="No applicable data", anchor="w", fill="#666666"
            )
            return

        raw_valid_points = [
            dict(point) for point in points
            if np.isfinite(point.get("x", np.nan))
            and np.isfinite(point.get("y", np.nan))
            and (not log_x or float(point.get("x", 0.0)) > 0)
        ]
        if not raw_valid_points:
            canvas.create_text(
                left, height // 2, text="No finite values", anchor="w", fill="#666666"
            )
            return

        if aggregate_points:
            grouped_points: Dict[Tuple[str, bool, float, float], List[dict]] = {}
            for point in raw_valid_points:
                key = (
                    str(point.get("engine") or "unknown"),
                    bool(point.get("rescore")),
                    round(float(point["x"]), 2),
                    round(float(point["y"]), 2),
                )
                grouped_points.setdefault(key, []).append(point)
            valid_points = []
            for group in grouped_points.values():
                point = dict(group[0])
                point["count"] = len(group)
                if len(group) > 1:
                    point["tooltip"] = (
                        f"{len(group)} overlapping conditions at this coordinate\n"
                        f"Example: {point.get('tooltip') or point.get('label', '')}"
                    )
                valid_points.append(point)
        else:
            valid_points = raw_valid_points

        x_values = [float(point["x"]) for point in valid_points]
        y_values = [float(point["y"]) for point in valid_points]
        raw_x_max = max(
            x_values + ([float(x_threshold)] if x_threshold is not None else [])
        )
        y_max = max(y_values + ([float(y_threshold)] if y_threshold is not None else []))
        if equal_axes:
            shared_max = max(raw_x_max, y_max) * 1.06
            x_axis_min, x_axis_max = 0.0, shared_max
            y_max = shared_max
        elif log_x:
            log_values = [math.log10(value) for value in x_values]
            x_axis_min, x_axis_max = min(log_values), max(log_values)
            padding = max((x_axis_max - x_axis_min) * 0.04, 0.03)
            x_axis_min -= padding
            x_axis_max += padding
            y_max = max(y_max * 1.06, 1e-6)
        else:
            x_axis_min = 0.0
            x_axis_max = max(raw_x_max * 1.06, 1e-6)
            y_max = max(y_max * 1.06, 1e-6)

        def _x(value: float) -> float:
            transformed = math.log10(float(value)) if log_x else float(value)
            return left + (
                (transformed - x_axis_min) / (x_axis_max - x_axis_min) * plot_width
            )

        def _y(value: float) -> float:
            return top + plot_height - float(value) / y_max * plot_height

        def _tick_text(value: float, maximum: float) -> str:
            return f"{value:.0f}" if maximum >= 20 else f"{value:.1f}"

        for index in range(6):
            fraction = index / 5
            x_axis_value = x_axis_min + (x_axis_max - x_axis_min) * fraction
            x_value = 10 ** x_axis_value if log_x else x_axis_value
            y_value = y_max * fraction
            x_coord = left + plot_width * fraction
            y_coord = top + plot_height * (1 - fraction)
            canvas.create_line(
                x_coord, top, x_coord, top + plot_height, fill="#ECECE8"
            )
            canvas.create_line(
                left, y_coord, left + plot_width, y_coord, fill="#ECECE8"
            )
            canvas.create_text(
                x_coord, top + plot_height + 5,
                text=_tick_text(x_value, max(x_values)),
                anchor="n", fill="#555555", font=("TkDefaultFont", 8),
            )
            canvas.create_text(
                left - 7, y_coord, text=_tick_text(y_value, y_max),
                anchor="e", fill="#555555", font=("TkDefaultFont", 8),
            )

        if x_threshold is not None:
            x_coord = _x(float(x_threshold))
            canvas.create_line(
                x_coord, top, x_coord, top + plot_height,
                fill="#B5452D", dash=(5, 3), width=1,
            )
            canvas.create_text(
                x_coord + 3, top + 2, text=f"{x_threshold:g} A", anchor="nw",
                fill="#9C3D24", font=("TkDefaultFont", 8),
            )
        if y_threshold is not None:
            y_coord = _y(float(y_threshold))
            canvas.create_line(
                left, y_coord, left + plot_width, y_coord,
                fill="#B5452D", dash=(5, 3), width=1,
            )
            canvas.create_text(
                left + 3, y_coord - 2, text=f"{y_threshold:g} A", anchor="sw",
                fill="#9C3D24", font=("TkDefaultFont", 8),
            )
        if diagonal:
            diagonal_max = min(x_axis_max, y_max)
            canvas.create_line(
                _x(0.0), _y(0.0), _x(diagonal_max), _y(diagonal_max),
                fill="#595959", dash=(3, 3), width=1,
            )
            canvas.create_text(
                _x(diagonal_max * 0.74), _y(diagonal_max * 0.74) - 5,
                text="no change", anchor="sw", fill="#555555",
                font=("TkDefaultFont", 8),
            )
            canvas.create_text(
                _x(diagonal_max * 0.22), _y(diagonal_max * 0.70),
                text="rescoring worsens", anchor="center", fill="#8A4A3B",
                font=("TkDefaultFont", 8),
            )
            canvas.create_text(
                _x(diagonal_max * 0.72), _y(diagonal_max * 0.23),
                text="rescoring improves", anchor="center", fill="#356B39",
                font=("TkDefaultFont", 8),
            )
        if quadrant_labels and x_threshold is not None and y_threshold is not None:
            canvas.create_text(
                left + 8, top + 8, text="pose generated; ranking failed",
                anchor="nw", fill="#666666", font=("TkDefaultFont", 8),
            )
            canvas.create_text(
                left + 8, _y(y_threshold) + 8, text="Top-1 recovered",
                anchor="nw", fill="#356B39", font=("TkDefaultFont", 8),
            )
            canvas.create_text(
                _x(x_threshold) + 8, top + 8, text="search and ranking failed",
                anchor="nw", fill="#8A4A3B", font=("TkDefaultFont", 8),
            )

        if pareto:
            frontier = sorted(
                (point for point in valid_points if point.get("pareto")),
                key=lambda point: float(point["x"]),
            )
            if len(frontier) > 1:
                coordinates = []
                for point in frontier:
                    coordinates.extend((_x(point["x"]), _y(point["y"])))
                canvas.create_line(
                    *coordinates, fill="#313131", width=2, dash=(5, 2)
                )

        engines = sorted({str(point.get("engine") or "unknown") for point in valid_points})
        legend_x = left
        for engine in engines:
            color = self._protocol_engine_color(engine)
            canvas.create_oval(
                legend_x, 39, legend_x + 9, 48, fill=color, outline=color
            )
            canvas.create_text(
                legend_x + 13, 43.5, text=engine, anchor="w",
                font=("TkDefaultFont", 8),
            )
            legend_x += max(62, len(engine) * 7 + 28)
        if len({bool(point.get("rescore")) for point in valid_points}) > 1:
            canvas.create_text(
                legend_x, 43.5, text="circle=docking  diamond=rescored",
                anchor="w", fill="#555555", font=("TkDefaultFont", 8),
            )
            legend_x += 175
        if any(int(point.get("count", 1)) > 1 for point in valid_points):
            canvas.create_text(
                legend_x, 43.5, text="size=overlap count",
                anchor="w", fill="#555555", font=("TkDefaultFont", 8),
            )
        if pareto:
            canvas.create_text(
                width - right, 43.5, text="outlined = Pareto-efficient",
                anchor="e", fill="#333333", font=("TkDefaultFont", 8),
            )

        for index, point in enumerate(valid_points):
            x_coord, y_coord = _x(point["x"]), _y(point["y"])
            color = self._protocol_engine_color(point.get("engine", "unknown"))
            tag = f"protocol_point_{index}"
            outline = "#171717" if point.get("pareto") else "#FFFFFF"
            line_width = 2 if point.get("pareto") else 1
            radius = 4.0 + min(5.0, max(0.0, math.sqrt(point.get("count", 1)) - 1.0))
            if point.get("rescore"):
                marker = canvas.create_polygon(
                    x_coord, y_coord - radius, x_coord + radius, y_coord,
                    x_coord, y_coord + radius, x_coord - radius, y_coord,
                    fill=color, outline=outline, width=line_width, tags=(tag,),
                )
            else:
                marker = canvas.create_oval(
                    x_coord - radius, y_coord - radius,
                    x_coord + radius, y_coord + radius,
                    fill=color, outline=outline, width=line_width, tags=(tag,),
                )
            tooltip = point.get("tooltip") or point.get("label") or "Protocol condition"
            canvas.tag_bind(
                marker, "<Enter>",
                lambda event, text=tooltip, w=width, h=height:
                self._draw_chart_tooltip(canvas, event, text, w, h),
            )
            canvas.tag_bind(
                marker, "<Leave>", lambda _event: canvas.delete("chart_tooltip")
            )

        canvas.create_line(
            left, top, left, top + plot_height, fill="#4F4F4F", width=1
        )
        canvas.create_line(
            left, top + plot_height, left + plot_width, top + plot_height,
            fill="#4F4F4F", width=1,
        )
        canvas.create_text(
            left + plot_width / 2, height - 9, text=x_label, anchor="s",
            font=("TkDefaultFont", 8),
        )
        canvas.create_text(
            13, top + plot_height / 2, text=y_label, anchor="s", angle=90,
            font=("TkDefaultFont", 8),
        )

    def _draw_factor_effect_chart(
        self,
        canvas: tk.Canvas,
        title: str,
        rows: List[dict],
        threshold: float = 2.0,
    ) -> None:
        """Show median, IQR, and range for each protocol factor level."""
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 460)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 280)
        left, right, top, bottom = 142, 42, 48, 34
        plot_width = max(100, width - left - right)
        plot_height = max(100, height - top - bottom)
        canvas.create_text(
            14, 12, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold")
        )
        canvas.create_text(
            14, 29,
            text="Each level pools other tested factors. Line=range; box=IQR; tick=median.",
            anchor="nw", fill="#666666", font=("TkDefaultFont", 8),
        )
        if not rows:
            canvas.create_text(
                left, height // 2, text="No applicable data", anchor="w", fill="#666666"
            )
            return

        maximum = max(max(float(row["maximum"]), threshold) for row in rows) * 1.05
        row_height = plot_height / len(rows)

        def _x(value: float) -> float:
            return left + float(value) / maximum * plot_width

        for index in range(6):
            fraction = index / 5
            x_coord = left + fraction * plot_width
            value = fraction * maximum
            canvas.create_line(
                x_coord, top, x_coord, top + plot_height, fill="#ECECE8"
            )
            canvas.create_text(
                x_coord, top + plot_height + 4, text=f"{value:.1f}",
                anchor="n", fill="#555555", font=("TkDefaultFont", 8),
            )

        threshold_x = _x(threshold)
        canvas.create_line(
            threshold_x, top, threshold_x, top + plot_height,
            fill="#B5452D", dash=(5, 3), width=1,
        )
        canvas.create_text(
            threshold_x + 3, top + 2, text=f"{threshold:g} A", anchor="nw",
            fill="#9C3D24", font=("TkDefaultFont", 8),
        )

        factor_colors = {
            "Engine": "#177E89", "Box": "#D97732", "Water": "#5B8E3E",
            "Exhaustiveness": "#60758A", "Ranking": "#B45B73",
        }
        previous_factor = None
        for index, row in enumerate(rows):
            y_coord = top + (index + 0.5) * row_height
            factor = str(row["factor"])
            if previous_factor is not None and factor != previous_factor:
                canvas.create_line(
                    8, y_coord - row_height / 2, width - right,
                    y_coord - row_height / 2, fill="#D6D6D0",
                )
            previous_factor = factor
            color = factor_colors.get(factor, "#60758A")
            canvas.create_text(
                left - 7, y_coord, text=self._short_chart_label(row["label"], 24),
                anchor="e", fill="#252525", font=("TkDefaultFont", 8),
            )
            canvas.create_line(
                _x(row["minimum"]), y_coord, _x(row["maximum"]), y_coord,
                fill=color, width=2,
            )
            half_height = max(2.0, min(6.0, row_height * 0.30))
            canvas.create_rectangle(
                _x(row["q1"]), y_coord - half_height,
                _x(row["q3"]), y_coord + half_height,
                fill=color, outline=color,
            )
            canvas.create_line(
                _x(row["median"]), y_coord - half_height - 2,
                _x(row["median"]), y_coord + half_height + 2,
                fill="#151515", width=2,
            )
            value_x = _x(row["maximum"])
            value_anchor = "e" if value_x > width - right - 70 else "w"
            value_x += -4 if value_anchor == "e" else 4
            canvas.create_text(
                value_x, y_coord,
                text=f"{row['median']:.2f} (n={row['n']})", anchor=value_anchor,
                fill="#333333", font=("TkDefaultFont", 7),
            )
        canvas.create_text(
            left + plot_width / 2, height - 4, text="Selected Top-1 RMSD (A)",
            anchor="s", font=("TkDefaultFont", 8),
        )

    def _draw_xy_line_chart(
        self,
        canvas: tk.Canvas,
        title: str,
        points: List[List[float]],
        x_label: str,
        y_label: str,
        reference_diagonal: bool = False,
        horizontal_reference: Optional[float] = None,
        vertical_guides: Tuple[float, ...] = (),
        color: str = "#177E89",
    ) -> None:
        """Render a unit-scale screening-performance curve."""
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 460)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 280)
        left, right, top, bottom = 62, 24, 54, 48
        plot_width = max(100, width - left - right)
        plot_height = max(90, height - top - bottom)
        canvas.create_text(
            14, 12, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold")
        )

        valid = [
            (float(point[0]), float(point[1])) for point in points
            if len(point) >= 2
            and np.isfinite(point[0])
            and np.isfinite(point[1])
        ]
        if not valid:
            canvas.create_text(
                left, height // 2, text="No applicable data", anchor="w",
                fill="#666666",
            )
            return

        def _x(value: float) -> float:
            return left + min(max(value, 0.0), 1.0) * plot_width

        def _y(value: float) -> float:
            return top + (1.0 - min(max(value, 0.0), 1.0)) * plot_height

        for tick in np.linspace(0.0, 1.0, 6):
            x_coord = _x(float(tick))
            y_coord = _y(float(tick))
            canvas.create_line(
                x_coord, top, x_coord, top + plot_height, fill="#E7E7E2"
            )
            canvas.create_line(
                left, y_coord, left + plot_width, y_coord, fill="#E7E7E2"
            )
            canvas.create_text(
                x_coord, top + plot_height + 5, text=f"{tick:.1f}",
                anchor="n", fill="#555555", font=("TkDefaultFont", 7),
            )
            canvas.create_text(
                left - 7, y_coord, text=f"{tick:.1f}", anchor="e",
                fill="#555555", font=("TkDefaultFont", 7),
            )

        if reference_diagonal:
            canvas.create_line(
                _x(0.0), _y(0.0), _x(1.0), _y(1.0),
                fill="#8B8B83", dash=(4, 3),
            )
            canvas.create_text(
                _x(0.76), _y(0.70), text="random", anchor="nw",
                fill="#777777", font=("TkDefaultFont", 7),
            )
        if horizontal_reference is not None and np.isfinite(horizontal_reference):
            baseline = min(max(float(horizontal_reference), 0.0), 1.0)
            canvas.create_line(
                _x(0.0), _y(baseline), _x(1.0), _y(baseline),
                fill="#8B8B83", dash=(4, 3),
            )
            canvas.create_text(
                _x(1.0) - 3, _y(baseline) - 3,
                text=f"prevalence {baseline:.3f}", anchor="se",
                fill="#666666", font=("TkDefaultFont", 7),
            )
        for guide in vertical_guides:
            if not np.isfinite(guide) or guide <= 0 or guide >= 1:
                continue
            x_coord = _x(float(guide))
            canvas.create_line(
                x_coord, top, x_coord, top + plot_height,
                fill="#A5A59D", dash=(2, 3),
            )
            canvas.create_text(
                x_coord + 2, top + 2, text=f"{100 * guide:g}%",
                anchor="nw", fill="#666666", font=("TkDefaultFont", 7),
            )

        coordinates = []
        for x_value, y_value in valid:
            coordinates.extend((_x(x_value), _y(y_value)))
        if len(coordinates) >= 4:
            canvas.create_line(*coordinates, fill=color, width=2.5, smooth=False)
        else:
            canvas.create_oval(
                coordinates[0] - 2, coordinates[1] - 2,
                coordinates[0] + 2, coordinates[1] + 2,
                fill=color, outline=color,
            )
        canvas.create_line(
            left, top, left, top + plot_height, fill="#555555", width=1
        )
        canvas.create_line(
            left, top + plot_height, left + plot_width, top + plot_height,
            fill="#555555", width=1,
        )
        canvas.create_text(
            left + plot_width / 2, height - 4, text=x_label,
            anchor="s", font=("TkDefaultFont", 8),
        )
        canvas.create_text(
            12, top + plot_height / 2, text=y_label, angle=90,
            anchor="center", font=("TkDefaultFont", 8),
        )

    def _draw_score_distribution_chart(
        self,
        canvas: tk.Canvas,
        title: str,
        histogram: dict,
        x_label: str,
        y_label: str,
    ) -> None:
        """Compare active and inactive score distributions on equal class scales."""
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 460)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 280)
        left, right, top, bottom = 62, 24, 62, 48
        plot_width = max(100, width - left - right)
        plot_height = max(90, height - top - bottom)
        canvas.create_text(
            14, 12, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold")
        )

        edges = [float(value) for value in histogram.get("edges") or []]
        active = [float(value) for value in histogram.get("active_percent") or []]
        inactive = [float(value) for value in histogram.get("inactive_percent") or []]
        if len(edges) < 2 or len(active) != len(edges) - 1 or len(inactive) != len(active):
            canvas.create_text(
                left, height // 2, text="No applicable data", anchor="w",
                fill="#666666",
            )
            return

        score_min, score_max = edges[0], edges[-1]
        score_span = max(score_max - score_min, 1e-9)
        observed_max = max(active + inactive + [1.0])
        y_max = max(5.0, math.ceil(observed_max / 5.0) * 5.0)

        def _x(value: float) -> float:
            return left + (value - score_min) / score_span * plot_width

        def _y(value: float) -> float:
            return top + (1.0 - value / y_max) * plot_height

        for tick in np.linspace(0.0, y_max, 5):
            y_coord = _y(float(tick))
            canvas.create_line(
                left, y_coord, left + plot_width, y_coord, fill="#E7E7E2"
            )
            canvas.create_text(
                left - 7, y_coord, text=f"{tick:.0f}", anchor="e",
                fill="#555555", font=("TkDefaultFont", 7),
            )
        for tick in np.linspace(score_min, score_max, 5):
            x_coord = _x(float(tick))
            canvas.create_line(
                x_coord, top, x_coord, top + plot_height, fill="#F0F0EC"
            )
            canvas.create_text(
                x_coord, top + plot_height + 5, text=f"{tick:.1f}",
                anchor="n", fill="#555555", font=("TkDefaultFont", 7),
            )

        colors = (("Actives", "#177E89"), ("Inactives", "#D97732"))
        legend_x = left
        for label, fill in colors:
            canvas.create_rectangle(
                legend_x, 37, legend_x + 10, 47, fill=fill, outline=fill
            )
            canvas.create_text(
                legend_x + 14, 42, text=label, anchor="w",
                font=("TkDefaultFont", 8),
            )
            legend_x += 82

        for index, (active_value, inactive_value) in enumerate(zip(active, inactive)):
            x0 = _x(edges[index])
            x1 = _x(edges[index + 1])
            group_width = max(x1 - x0, 1.0)
            padding = min(1.0, group_width * 0.08)
            half_width = max((group_width - 2 * padding) / 2.0, 0.5)
            baseline = _y(0.0)
            canvas.create_rectangle(
                x0 + padding, _y(active_value),
                x0 + padding + half_width, baseline,
                fill="#177E89", outline="#0E5961",
            )
            canvas.create_rectangle(
                x0 + padding + half_width, _y(inactive_value),
                x1 - padding, baseline,
                fill="#D97732", outline="#9C4D18",
            )

        canvas.create_line(
            left, top, left, top + plot_height, fill="#555555", width=1
        )
        canvas.create_line(
            left, top + plot_height, left + plot_width, top + plot_height,
            fill="#555555", width=1,
        )
        canvas.create_text(
            left + plot_width / 2, height - 4, text=x_label,
            anchor="s", font=("TkDefaultFont", 8),
        )
        canvas.create_text(
            12, top + plot_height / 2, text=y_label, angle=90,
            anchor="center", font=("TkDefaultFont", 8),
        )

    def _draw_horizontal_chart(
        self,
        canvas: tk.Canvas,
        title: str,
        data: List[Tuple[str, Optional[float]]],
        max_value: Optional[float] = None,
        thresholds: Optional[List[Tuple[float, str]]] = None,
    ) -> None:
        canvas.delete("all")
        width, height, left, right, top, bottom = self._chart_dimensions(canvas)
        canvas.create_text(14, 14, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold"))
        if not data:
            canvas.create_text(left, height // 2, text="No applicable data", anchor="w", fill="#666666")
            return

        valid = [float(value) for _, value in data if value is not None and np.isfinite(value)]
        axis_max = float(max_value) if max_value is not None else (max(valid) if valid else 1.0)
        if thresholds:
            axis_max = max(axis_max, max(float(value) for value, _ in thresholds))
        axis_max = max(axis_max, 1e-9)
        plot_width = max(80, width - left - right)
        plot_height = max(70, height - top - bottom)
        row_height = plot_height / max(len(data), 1)

        for threshold, label in thresholds or []:
            x = left + min(max(float(threshold) / axis_max, 0.0), 1.0) * plot_width
            canvas.create_line(x, top - 7, x, height - bottom, fill="#8B8B83", dash=(4, 3))
            canvas.create_text(x + 3, top - 9, text=label, anchor="sw", fill="#666666", font=("TkDefaultFont", 8))

        for index, (label, raw_value) in enumerate(data):
            y = top + (index + 0.5) * row_height
            canvas.create_text(
                left - 8, y, text=self._short_chart_label(label), anchor="e",
                fill="#252525", font=("TkDefaultFont", 8),
            )
            if raw_value is None or not np.isfinite(raw_value):
                canvas.create_text(left + 4, y, text="N/A", anchor="w", fill="#888888")
                continue
            value = float(raw_value)
            x1 = left + min(max(value / axis_max, 0.0), 1.0) * plot_width
            half_height = max(2.0, min(9.0, row_height * 0.30))
            canvas.create_rectangle(
                left, y - half_height, x1, y + half_height,
                fill="#177E89", outline="#0E5961",
            )
            canvas.create_text(
                min(x1 + 5, width - right + 3), y, text=f"{value:.2f}",
                anchor="w", fill="#252525", font=("TkDefaultFont", 8),
            )
        canvas.create_line(left, height - bottom, width - right, height - bottom, fill="#777777")

    def _draw_grouped_horizontal_chart(
        self,
        canvas: tk.Canvas,
        title: str,
        data: List[Tuple[str, List[Optional[float]]]],
        series_names: Tuple[str, ...],
        max_value: Optional[float] = None,
        thresholds: Optional[List[Tuple[float, str]]] = None,
    ) -> None:
        canvas.delete("all")
        width, height, left, right, top, bottom = self._chart_dimensions(canvas)
        canvas.create_text(14, 14, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold"))
        if not data:
            canvas.create_text(left, height // 2, text="No applicable data", anchor="w", fill="#666666")
            return

        colors = (("#177E89", "#0E5961"), ("#D97732", "#9C4D18"), ("#5B8E3E", "#365E23"))
        legend_x = left
        for index, name in enumerate(series_names):
            fill, _ = colors[index % len(colors)]
            canvas.create_rectangle(legend_x, 34, legend_x + 10, 44, fill=fill, outline=fill)
            canvas.create_text(legend_x + 14, 39, text=name, anchor="w", font=("TkDefaultFont", 8))
            legend_x += max(95, len(name) * 7 + 30)

        valid = [
            float(value) for _, values in data for value in values
            if value is not None and np.isfinite(value)
        ]
        axis_max = float(max_value) if max_value is not None else (max(valid) if valid else 1.0)
        if thresholds:
            axis_max = max(axis_max, max(float(value) for value, _ in thresholds))
        axis_max = max(axis_max, 1e-9)
        plot_width = max(80, width - left - right)
        plot_height = max(70, height - top - bottom)
        group_height = plot_height / max(len(data), 1)

        for threshold, label in thresholds or []:
            x = left + min(max(float(threshold) / axis_max, 0.0), 1.0) * plot_width
            canvas.create_line(x, top - 2, x, height - bottom, fill="#8B8B83", dash=(4, 3))
            canvas.create_text(x + 3, top, text=label, anchor="nw", fill="#666666", font=("TkDefaultFont", 8))

        series_count = max(len(series_names), 1)
        for row_index, (label, values) in enumerate(data):
            center_y = top + (row_index + 0.5) * group_height
            canvas.create_text(
                left - 8, center_y, text=self._short_chart_label(label), anchor="e",
                fill="#252525", font=("TkDefaultFont", 8),
            )
            bar_height = max(2.0, min(7.0, group_height * 0.28))
            separation = bar_height * 2 + 2
            first_y = center_y - separation * (series_count - 1) / 2
            for series_index in range(series_count):
                value = values[series_index] if series_index < len(values) else None
                y = first_y + series_index * separation
                if value is None or not np.isfinite(value):
                    canvas.create_text(left + 4, y, text="N/A", anchor="w", fill="#888888", font=("TkDefaultFont", 8))
                    continue
                numeric = float(value)
                x1 = left + min(max(numeric / axis_max, 0.0), 1.0) * plot_width
                fill, outline = colors[series_index % len(colors)]
                canvas.create_rectangle(left, y - bar_height, x1, y + bar_height, fill=fill, outline=outline)
                canvas.create_text(
                    min(x1 + 4, width - right + 3), y, text=f"{numeric:.2f}",
                    anchor="w", fill="#252525", font=("TkDefaultFont", 8),
                )
        canvas.create_line(left, height - bottom, width - right, height - bottom, fill="#777777")

    def _draw_diverging_horizontal_chart(
        self, canvas: tk.Canvas, title: str, data: List[Tuple[str, Optional[float]]]
    ) -> None:
        canvas.delete("all")
        width, height, left, right, top, bottom = self._chart_dimensions(canvas)
        canvas.create_text(14, 14, text=title, anchor="nw", font=("TkDefaultFont", 10, "bold"))
        if not data:
            canvas.create_text(left, height // 2, text="No applicable data", anchor="w", fill="#666666")
            return

        valid = [abs(float(value)) for _, value in data if value is not None and np.isfinite(value)]
        extent = max(valid) if valid else 1.0
        extent = max(extent, 1e-9)
        plot_width = max(80, width - left - right)
        center = left + plot_width / 2
        row_height = max(70, height - top - bottom) / max(len(data), 1)
        canvas.create_line(center, top - 8, center, height - bottom, fill="#555555", width=2)
        canvas.create_text(center - 5, top - 10, text="decoy wins", anchor="se", fill="#9C3D24", font=("TkDefaultFont", 8))
        canvas.create_text(center + 5, top - 10, text="active wins", anchor="sw", fill="#0E5961", font=("TkDefaultFont", 8))

        for index, (label, raw_value) in enumerate(data):
            y = top + (index + 0.5) * row_height
            canvas.create_text(
                left - 8, y, text=self._short_chart_label(label), anchor="e",
                fill="#252525", font=("TkDefaultFont", 8),
            )
            if raw_value is None or not np.isfinite(raw_value):
                continue
            value = float(raw_value)
            delta = value / extent * (plot_width / 2)
            x0, x1 = sorted((center, center + delta))
            half_height = max(2.0, min(9.0, row_height * 0.30))
            fill, outline = ("#177E89", "#0E5961") if value >= 0 else ("#D95D39", "#9C3D24")
            canvas.create_rectangle(x0, y - half_height, x1, y + half_height, fill=fill, outline=outline)
            anchor = "w" if value >= 0 else "e"
            offset = 4 if value >= 0 else -4
            canvas.create_text(
                center + delta + offset, y, text=f"{value:+.2f}", anchor=anchor,
                fill="#252525", font=("TkDefaultFont", 8),
            )

    def _draw_bar_chart(self, canvas: tk.Canvas, title: str, data: List[Tuple[str, float]]) -> None:
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 300)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 200)
        margin = 40

        canvas.create_text(margin, 10, text=title, anchor="w")

        if not data:
            canvas.create_text(margin, height // 2, text="No data", anchor="w")
            return

        max_val = max(v for _, v in data) if data else 1.0
        bar_width = max(20, int((width - 2 * margin) / max(1, len(data)) - 10))
        scale = (height - 2 * margin) / max(max_val, 1.0)

        for idx, (label, value) in enumerate(data):
            x0 = margin + idx * (bar_width + 10)
            x1 = x0 + bar_width
            y1 = height - margin
            y0 = y1 - value * scale
            canvas.create_rectangle(x0, y0, x1, y1, fill="#4C89FF", outline="#2C5FA8")
            canvas.create_text((x0 + x1) / 2, y0 - 5, text=f"{value:.1f}", anchor="s")
            canvas.create_text((x0 + x1) / 2, y1 + 5, text=label, anchor="n", angle=45)

    def _draw_histogram(self, canvas: tk.Canvas, title: str, values: List[float]) -> None:
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), int(canvas.winfo_reqwidth()), 300)
        height = max(int(canvas.winfo_height()), int(canvas.winfo_reqheight()), 200)
        margin = 40

        canvas.create_text(margin, 10, text=title, anchor="w")

        if not values:
            canvas.create_text(margin, height // 2, text="No RMSD data", anchor="w")
            return

        bins = [0, 1, 2, 3, 4, 5, 10]
        counts = [0] * (len(bins) - 1)
        for v in values:
            for i in range(len(bins) - 1):
                if bins[i] <= v < bins[i + 1]:
                    counts[i] += 1
                    break
            else:
                if v >= bins[-1]:
                    counts[-1] += 1

        max_count = max(counts) if counts else 1
        bar_width = max(20, int((width - 2 * margin) / len(counts) - 10))
        scale = (height - 2 * margin) / max(max_count, 1)

        for idx, count in enumerate(counts):
            x0 = margin + idx * (bar_width + 10)
            x1 = x0 + bar_width
            y1 = height - margin
            y0 = y1 - count * scale
            canvas.create_rectangle(x0, y0, x1, y1, fill="#7BC96F", outline="#3C7F35")
            canvas.create_text((x0 + x1) / 2, y0 - 5, text=str(count), anchor="s")
            label = f"{bins[idx]}-{bins[idx+1]}"
            canvas.create_text((x0 + x1) / 2, y1 + 5, text=label, anchor="n")

    def _download_pdb(self, pdb_id: str, pdb_dir: Path) -> Path:
        pdb_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = pdb_dir / f"{pdb_id.upper()}.pdb"
        if pdb_path.exists():
            return pdb_path
        success, message, file_path = download_pdb_structure(pdb_id, str(pdb_dir))
        if not success or not file_path:
            raise ValueError(message)
        return Path(file_path)

    @staticmethod
    def _has_covalent_ligand_link(
        pdb_path: Path, ligand_resname: str, ligand_chain: Optional[str]
    ) -> bool:
        """Detect PDB LINK records joining the site ligand to an amino acid."""
        amino_acids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        }
        ligand_resname = ligand_resname.upper()
        with open(pdb_path) as handle:
            for line in handle:
                if not line.startswith("LINK"):
                    continue
                first = (line[17:20].strip().upper(), line[21:22].strip())
                second = (line[47:50].strip().upper(), line[51:52].strip())
                for ligand, partner in ((first, second), (second, first)):
                    chain_matches = not ligand_chain or ligand[1] == ligand_chain
                    if ligand[0] == ligand_resname and chain_matches and partner[0] in amino_acids:
                        return True
        return False

    def _prefetch_remote_inputs(self, pairs: List[Dict[str, str]], output_dir: Path) -> None:
        """Download every campaign input before any docking calculation starts."""
        self._network_phase_complete = False
        pdb_dir = output_dir / "pdbs"
        unique_pdb_ids = list(dict.fromkeys(str(pair["pdb_id"]).upper() for pair in pairs))
        self._queue.put((
            "log",
            f"Online preflight: acquiring {len(unique_pdb_ids)} structure(s) and ligand records"
        ))
        pdb_files = {}
        for pdb_id in unique_pdb_ids:
            self._queue.put(("log", f"Online preflight: {pdb_id}"))
            pdb_files[pdb_id] = self._download_pdb(pdb_id, pdb_dir)

        resolved_crystal_ligands = {}
        for pair in pairs:
            pdb_id = str(pair["pdb_id"]).upper()
            site_ligand = pair.get("site_ligand")
            if not site_ligand:
                smiles = pair.get("smiles")
                if smiles is None or not str(smiles).strip() or str(smiles).lower() == "nan":
                    raise ValueError(
                        f"No SMILES was supplied for an apo/predicted-site case at {pdb_id}."
                    )
                pair["is_reference_ligand"] = False
                continue
            chain = pair.get("chain") or self._detect_ligand_chain(
                pdb_files[pdb_id], site_ligand
            )
            if not chain:
                raise ValueError(f"Ligand chain not found for {pdb_id}/{site_ligand}")
            pair["chain"] = chain
            ligand_key = (pdb_id, site_ligand, chain)
            if ligand_key not in resolved_crystal_ligands:
                resolved_crystal_ligands[ligand_key] = self._get_ligand_smiles(
                    pdb_files[pdb_id], site_ligand, chain, output_dir
                )

            smiles = pair.get("smiles")
            if smiles is None or not str(smiles).strip() or str(smiles).lower() == "nan":
                smiles = resolved_crystal_ligands[ligand_key]
                if not smiles:
                    raise ValueError(
                        f"Could not resolve ligand SMILES for {pdb_id}/{site_ligand}. "
                        "Provide it in the template before running offline."
                    )
                pair["smiles"] = smiles
            pair["is_reference_ligand"] = self._same_ligand_smiles(
                str(pair["smiles"]), resolved_crystal_ligands[ligand_key]
            )

        self._network_phase_complete = True
        self._queue.put((
            "log",
            "Online preflight complete. All remaining preparation, docking, and analysis are local."
        ))

    def _detect_ligand_chain(self, pdb_path: Path, ligand_resname: str) -> Optional[str]:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname == ligand_resname:
                        return chain.id.strip()
        return None

    def _is_reference_ligand(
        self,
        pdb_path: Path,
        ligand_resname: Optional[str],
        ligand_chain: str,
        smiles: str,
        output_dir: Path
    ) -> Optional[bool]:
        """Is `smiles` the structure's own co-crystal ligand?

        Returns True (same molecule), False (different), or None when it
        cannot be determined — callers should warn rather than block on None.

        Used to keep the adaptive cascade off ligands that have no reference
        pose, where RMSD is undefined and the escalation cannot terminate.
        """
        if not ligand_resname or not smiles:
            return None
        try:
            crystal_smiles = self._get_ligand_smiles(
                pdb_path, ligand_resname, ligand_chain, output_dir
            )
            return self._same_ligand_smiles(smiles, crystal_smiles)
        except Exception as exc:
            logger.warning("Reference-ligand check failed for %s: %s", ligand_resname, exc)
            return None

    @staticmethod
    def _same_ligand_smiles(smiles: Optional[str], crystal_smiles: Optional[str]) -> Optional[bool]:
        """Compare molecular identity without conflating it with an activity label."""
        if not smiles or not crystal_smiles:
            return None
        docked_mol = Chem.MolFromSmiles(smiles)
        crystal_mol = Chem.MolFromSmiles(crystal_smiles)
        if docked_mol is None or crystal_mol is None:
            return None
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize

            def identity_key(molecule: Chem.Mol) -> str:
                parent = rdMolStandardize.FragmentParent(Chem.Mol(molecule))
                parent = rdMolStandardize.Uncharger().uncharge(parent)
                parent = rdMolStandardize.TautomerEnumerator().Canonicalize(parent)
                return Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)

            return identity_key(docked_mol) == identity_key(crystal_mol)
        except Exception as exc:
            logger.warning("Ligand identity standardization failed: %s", exc)
            return None

    def _get_ligand_smiles(
        self,
        pdb_path: Path,
        ligand_resname: str,
        ligand_chain: str,
        output_dir: Path
    ) -> Optional[str]:
        known = {
            "BUM": "CNC(=O)[C@@H](N)C(C)(C)C",
            "STI": "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
        }
        if ligand_resname in known:
            return known[ligand_resname]

        sdf_file = output_dir / f"{ligand_resname}_ideal.sdf"
        try:
            if not sdf_file.exists() and not getattr(self, "_network_phase_complete", False):
                url = f"https://files.rcsb.org/ligands/view/{ligand_resname}_ideal.sdf"
                urllib.request.urlretrieve(url, sdf_file)
            if not sdf_file.exists():
                raise FileNotFoundError(sdf_file)
            suppl = Chem.SDMolSupplier(str(sdf_file))
            mol = next(suppl)
            if mol is not None:
                return Chem.MolToSmiles(mol)
        except Exception:
            pass

        return self._extract_smiles_from_pdb(pdb_path, ligand_resname, ligand_chain)

    def _extract_smiles_from_pdb(
        self,
        pdb_file: Path,
        ligand_resname: str,
        ligand_chain: str
    ) -> Optional[str]:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_file))

        ligand_residue = None
        for model in structure:
            for chain in model:
                if ligand_chain and chain.id != ligand_chain:
                    continue
                for residue in chain:
                    if residue.resname == ligand_resname:
                        ligand_residue = residue
                        break
                if ligand_residue:
                    break
            if ligand_residue:
                break

        if not ligand_residue:
            return None

        temp_pdb = pdb_file.parent / f"{ligand_resname}_tmp.pdb"
        io = PDB.PDBIO()
        io.set_structure(ligand_residue)
        io.save(str(temp_pdb))
        try:
            mol = Chem.MolFromPDBFile(str(temp_pdb), removeHs=True)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)
        finally:
            if temp_pdb.exists():
                temp_pdb.unlink()

    def _get_ligand_charge(self, smiles: str) -> Optional[int]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return int(Chem.GetFormalCharge(mol))
        except Exception:
            return None

    def _get_ligand_properties(self, smiles: str) -> dict:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {}
            return {
                "molecular_weight": float(Descriptors.MolWt(mol)),
                "logp": float(Crippen.MolLogP(mol)),
                "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
                "rotatable_bonds": int(Descriptors.NumRotatableBonds(mol)),
            }
        except Exception:
            return {}

    def _parse_pdb_coords(
        self,
        pdb_path: Path,
        record_types: Tuple[str, ...] = ("ATOM", "HETATM")
    ) -> np.ndarray:
        coords: List[List[float]] = []
        with open(pdb_path) as handle:
            for line in handle:
                if not line.startswith(record_types):
                    continue
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except ValueError:
                    continue
                coords.append([x, y, z])
        return np.asarray(coords, dtype=float) if coords else np.empty((0, 3), dtype=float)

    def _protein_centroid(self, pdb_path: Path) -> np.ndarray:
        coords = self._parse_pdb_coords(pdb_path, record_types=("ATOM",))
        if coords.size == 0:
            raise ValueError(f"No ATOM coordinates found in {pdb_path.name}")
        return coords.mean(axis=0)

    def _parse_site_residue_specs(
        self,
        residues_text: str,
        default_chain: Optional[str] = None
    ) -> List[Tuple[Optional[str], int]]:
        specs: List[Tuple[Optional[str], int]] = []
        for raw in re.split(r"[,;\s]+", residues_text.strip()):
            token = raw.strip()
            if not token:
                continue
            if token.isdigit():
                specs.append((default_chain or None, int(token)))
                continue
            match = re.match(r"^([A-Za-z0-9])[:_-](\d+)$", token)
            if match:
                specs.append((match.group(1), int(match.group(2))))
                continue
            match = re.match(r"^([A-Za-z])(\d+)$", token)
            if match:
                specs.append((match.group(1), int(match.group(2))))
                continue
            match = re.match(r"^(\d+)([A-Za-z0-9])$", token)
            if match:
                specs.append((match.group(2), int(match.group(1))))
                continue
            raise ValueError(f"Invalid residue spec '{token}'. Use A:118, A118, 118A, or 118.")
        if not specs:
            raise ValueError("No site residues provided")
        return specs

    def _binding_site_from_residues(
        self,
        pdb_path: Path,
        residues_text: str,
        box_margin: float,
        size_override: Optional[np.ndarray],
        default_chain: Optional[str] = None
    ) -> Tuple[BindingSite, str]:
        specs = self._parse_site_residue_specs(residues_text, default_chain=default_chain)
        spec_set = set(specs)
        coords: List[List[float]] = []

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        for model in structure:
            for chain in model:
                chain_id = chain.id.strip() or None
                for residue in chain:
                    resseq = residue.id[1]
                    if (chain_id, resseq) not in spec_set and (None, resseq) not in spec_set:
                        continue
                    for atom in residue.get_atoms():
                        if atom.element == "H":
                            continue
                        coord = atom.get_coord()
                        coords.append([float(coord[0]), float(coord[1]), float(coord[2])])
            break

        if not coords:
            formatted = ",".join(f"{chain or '*'}:{resnum}" for chain, resnum in specs)
            raise ValueError(f"No atoms found for site residues: {formatted}")

        arr = np.asarray(coords, dtype=float)
        center = arr.mean(axis=0)
        extent = arr.max(axis=0) - arr.min(axis=0)
        size = size_override if size_override is not None else np.clip(extent + (2.0 * float(box_margin)), 18.0, 40.0)
        method = "residues:" + ",".join(f"{chain or '*'}:{resnum}" for chain, resnum in specs)

        return (
            BindingSite(
                center=np.asarray(center, dtype=float),
                size=np.asarray(size, dtype=float),
                ligand_resname="RES",
                ligand_chain=default_chain or "",
                ligand_atoms=len(coords),
                source_pdb=str(pdb_path)
            ),
            method
        )

    def _detect_fpocket_center(
        self,
        pdb_path: Path,
        case_dir: Path,
        box_margin: float
    ) -> Optional[Tuple[np.ndarray, np.ndarray, str]]:
        fpocket_bin = shutil.which("fpocket")
        if not fpocket_bin:
            return None

        fpocket_input = case_dir / "fpocket_input.pdb"
        shutil.copy2(pdb_path, fpocket_input)

        fpocket_out = case_dir / f"{fpocket_input.stem}_out"
        if fpocket_out.exists():
            shutil.rmtree(fpocket_out, ignore_errors=True)

        try:
            result = subprocess.run(
                [fpocket_bin, "-f", str(fpocket_input)],
                cwd=case_dir,
                capture_output=True,
                text=True,
                timeout=600
            )
        except Exception as exc:
            logger.warning("fpocket failed to launch: {}", exc)
            return None

        if result.returncode != 0:
            stderr = (result.stderr or "").strip().splitlines()
            msg = stderr[-1] if stderr else "unknown error"
            logger.warning("fpocket failed: {}", msg)
            return None

        pockets_dir = fpocket_out / "pockets"
        if not pockets_dir.exists():
            logger.warning("fpocket output missing pockets directory: {}", pockets_dir)
            return None

        pocket_files = sorted(pockets_dir.glob("pocket*_atm.pdb"), key=lambda p: p.name)
        if not pocket_files:
            logger.warning("fpocket produced no pocket*_atm.pdb files")
            return None

        def _pocket_rank(path: Path) -> int:
            m = re.search(r"pocket(\d+)_atm\.pdb", path.name)
            return int(m.group(1)) if m else 10_000

        pocket_file = sorted(pocket_files, key=_pocket_rank)[0]
        pocket_coords = self._parse_pdb_coords(pocket_file, record_types=("ATOM", "HETATM"))
        if pocket_coords.size == 0:
            logger.warning("fpocket top pocket had no coordinates: {}", pocket_file.name)
            return None

        center = pocket_coords.mean(axis=0)
        extent = pocket_coords.max(axis=0) - pocket_coords.min(axis=0)
        size = np.clip(extent + (2.0 * float(box_margin)), 18.0, 40.0)
        return center, size, f"fpocket:{pocket_file.name}"

    def _predict_binding_site(
        self,
        pdb_path: Path,
        case_dir: Path,
        box_margin: float,
        size_override: Optional[np.ndarray],
        apo_site_mode: str,
        manual_center: Optional[Tuple[float, float, float]]
    ) -> Tuple[BindingSite, str]:
        if manual_center is not None:
            center = np.asarray(manual_center, dtype=float)
            size = size_override if size_override is not None else np.array([22.0, 22.0, 22.0], dtype=float)
            method = "manual_center"
        else:
            fpocket_result = None
            if apo_site_mode in {"auto", "fpocket"}:
                fpocket_result = self._detect_fpocket_center(
                    pdb_path=pdb_path,
                    case_dir=case_dir,
                    box_margin=box_margin
                )

            if fpocket_result is not None:
                center, fpocket_size, method = fpocket_result
                size = size_override if size_override is not None else fpocket_size
            else:
                if apo_site_mode == "fpocket":
                    raise ValueError(
                        "Apo site mode is 'fpocket' but fpocket failed or is unavailable. "
                        "Install fpocket or switch to 'auto'/'protein_centroid'."
                    )
                center = self._protein_centroid(pdb_path)
                size = size_override if size_override is not None else np.array([24.0, 24.0, 24.0], dtype=float)
                method = "protein_centroid"

        site = BindingSite(
            center=np.asarray(center, dtype=float),
            size=np.asarray(size, dtype=float),
            ligand_resname="POC",
            ligand_chain="A",
            ligand_atoms=0,
            source_pdb=str(pdb_path)
        )
        return site, method

    def _resolve_smina_binary(self, candidate: Optional[str]) -> Optional[str]:
        if candidate:
            path = Path(candidate)
            if path.exists():
                return str(path)
        resolved = shutil.which("smina")
        return resolved

    def _find_receptor_pdbqt(self, case_dir: Path, output_file: Path) -> Optional[Path]:
        candidates = [
            case_dir / "receptor_prepared.pdbqt",
            output_file.parent / "receptor.pdbqt",
            output_file.parent / "receptor_prepared.pdbqt",
            output_file.parent.parent / "receptor_prepared.pdbqt"
        ]
        for path in candidates:
            if path.exists():
                return path
        return None

    def _convert_to_pdbqt(self, input_file: Path, output_file: Path) -> bool:
        try:
            result = subprocess.run(
                ["obabel", str(input_file), "-O", str(output_file)],
                capture_output=True,
                text=True
            )
            return result.returncode == 0 and output_file.exists()
        except Exception:
            return False

    def _parse_smina_score(self, stdout: str) -> Optional[dict]:
        scores: List[float] = []
        for line in stdout.splitlines():
            if "Affinity:" not in line:
                continue
            match = re.search(r"Affinity:\s*([-0-9.]+)", line)
            if not match:
                continue
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                continue
        if not scores:
            return None
        return {"score": min(scores), "scores": scores}

    def _single_model_pdbqt(self, source: Path, destination: Path) -> Path:
        """Return a score-only compatible PDBQT containing the first pose."""
        lines = source.read_text().splitlines(keepends=True)
        if not any(line.lstrip().startswith("MODEL") for line in lines):
            return source

        pose_lines = []
        in_first_model = False
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith("MODEL"):
                if in_first_model:
                    break
                in_first_model = True
                continue
            if in_first_model and stripped.startswith("ENDMDL"):
                break
            if in_first_model:
                pose_lines.append(line)

        if not pose_lines:
            return source
        destination.write_text("".join(pose_lines))
        return destination

    def _rescore_with_smina(
        self,
        output_file: Path,
        case_dir: Path,
        smina_binary: str,
        scoring: str
    ) -> Optional[dict]:
        receptor_pdbqt = self._find_receptor_pdbqt(case_dir, output_file)
        if receptor_pdbqt is None:
            return {"error": "Receptor PDBQT not found"}

        ligand_pdbqt = output_file
        if output_file.suffix.lower() not in (".pdbqt",):
            candidate = case_dir / f"rescore_{output_file.stem}.pdbqt"
            if self._convert_to_pdbqt(output_file, candidate):
                ligand_pdbqt = candidate
            else:
                return {"error": f"Cannot rescore non-PDBQT output: {output_file.name}"}

        pose_files = self._split_pdbqt_models(
            ligand_pdbqt, case_dir / f"rescore_{scoring}_poses"
        )
        scores = []
        for pose_index, pose_file in enumerate(pose_files, 1):
            try:
                result = subprocess.run(
                    [
                        smina_binary,
                        "--receptor", str(receptor_pdbqt),
                        "--ligand", str(pose_file),
                        "--score_only",
                        "--scoring", scoring,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
            except Exception as exc:
                return {"error": f"Smina rescoring failed for pose {pose_index}: {exc}"}
            if result.returncode != 0:
                message = result.stderr.strip() or "Smina rescoring failed"
                return {"error": f"Pose {pose_index}: {message}"}
            parsed = self._parse_smina_score(result.stdout)
            if not parsed:
                return {"error": f"Pose {pose_index}: Smina output did not contain scores"}
            scores.append(float(parsed["score"]))

        if not scores:
            return {"error": "No poses available for rescoring"}
        return {
            "method": f"smina_score_only:{scoring}",
            "score": min(scores),
            "scores": scores,
            "pose_count": len(scores),
        }

    @staticmethod
    def _split_pdbqt_models(source: Path, output_dir: Path) -> List[Path]:
        """Split a Vina/Smina multi-model PDBQT into score-only inputs."""
        lines = source.read_text().splitlines(keepends=True)
        models: List[List[str]] = []
        current: Optional[List[str]] = None
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith("MODEL"):
                if current:
                    models.append(current)
                current = []
                continue
            if stripped.startswith("ENDMDL"):
                if current is not None:
                    models.append(current)
                    current = None
                continue
            if current is not None:
                current.append(line)
        if current:
            models.append(current)
        if not models:
            models = [lines]

        output_dir.mkdir(parents=True, exist_ok=True)
        outputs = []
        for index, model in enumerate(models, 1):
            destination = output_dir / f"pose_{index:03d}.pdbqt"
            destination.write_text("".join(model))
            outputs.append(destination)
        return outputs

    @staticmethod
    def _aggregate_pose_metrics(results: list, threshold: float) -> dict:
        """Aggregate per-case pose metrics into summary-level statistics.

        Only cases carrying a real RMSD contribute: decoys and failed cases
        have no reference pose, and the 999.9 sentinel would otherwise poison
        every mean. Returns None for a field when no case supplies it, so the
        summary still renders N/A rather than a misleading zero.
        """
        def _vals(attr):
            out = []
            for r in results:
                v = getattr(r, attr, None)
                if v is None:
                    continue
                if attr.endswith("rmsd") and v >= 900:
                    continue  # sentinel for "no RMSD", not a 999 Å pose
                out.append(float(v))
            return out

        def _mean(vals):
            return float(sum(vals) / len(vals)) if vals else None

        def _median(vals):
            if not vals:
                return None
            s = sorted(vals)
            mid = len(s) // 2
            return float(s[mid]) if len(s) % 2 else float((s[mid - 1] + s[mid]) / 2)

        def _rate(attr):
            vals = _vals(attr)
            if not vals:
                return None
            return 100.0 * sum(1 for v in vals if v < threshold) / len(vals)

        top1 = _vals("top1_rmsd")
        return {
            "success_rate_top1": _rate("top1_rmsd"),
            "success_rate_top5": _rate("top5_rmsd"),
            "success_rate_top10": _rate("top10_rmsd"),
            "mean_top1_rmsd": _mean(top1),
            "median_top1_rmsd": _median(top1),
            "mean_rmsd_best_score": _mean(_vals("rmsd_best_score")),
            "mean_near_native_fraction": _mean(_vals("near_native_fraction")),
            "mean_pose_count": _mean(_vals("pose_count")),
            "mean_score_rmsd_pearson": _mean(_vals("score_rmsd_pearson")),
            "mean_score_rmsd_spearman": _mean(_vals("score_rmsd_spearman")),
            "median_runtime_sec": _median(_vals("runtime_sec")),
        }

    def _compute_pose_metrics(
        self,
        crystal_ligand_pdb: Path,
        docked_file: Path,
        threshold: float,
        has_reference_pose: bool = True,
    ) -> dict:
        """Pose-level metrics for one docked case.

        RMSD is only defined when the docked molecule IS the reference molecule.
        Decoys and most screening actives are different molecules with no known pose, so any
        "RMSD" for them is meaningless. It is not merely unknown: computing it
        anyway previously let the MCS fallback match a fragment of a small
        compound onto a large crystal ligand and report sub-angstrom agreement
        (an 11-atom compound scored 0.06 A against 68-atom imatinib), which then
        propagated into the summary's success rate and mean RMSD.

        Activity labels are deliberately not accepted here: they classify
        compounds for enrichment but do not establish a reference pose.
        has_reference_pose=False returns score-only metrics. Redock callers set
        it because identity is guaranteed by construction; screening callers
        set it only after a molecular-identity check.
        """
        if not docked_file.exists():
            return {}

        poses, scores = self._load_poses_and_scores(docked_file)
        if not poses:
            return {}

        pose_count = len(poses)
        best_score = None
        if scores:
            valid_scores = [s for s in scores if s is not None]
            if valid_scores:
                best_score = min(valid_scores)

        if not has_reference_pose:
            return {
                "pose_count": pose_count,
                "best_score": best_score
            }

        ref_mol = self._load_reference_mol(crystal_ligand_pdb)
        if ref_mol is None:
            return {
                "pose_count": pose_count,
                "best_score": best_score
            }

        rmsds: List[Optional[float]] = []
        for mol in poses:
            rmsd = self._pose_rmsd(ref_mol, mol)
            rmsds.append(rmsd)

        valid_rmsds = [r for r in rmsds if r is not None]
        if not valid_rmsds:
            return {
                "pose_count": pose_count,
                "best_score": best_score
            }
        top1_rmsd = rmsds[0] if rmsds else None
        top5_candidates = [r for r in rmsds[:min(5, pose_count)] if r is not None]
        top10_candidates = [r for r in rmsds[:min(10, pose_count)] if r is not None]
        top5_rmsd = min(top5_candidates) if top5_candidates else None
        top10_rmsd = min(top10_candidates) if top10_candidates else None
        best_rmsd = min(valid_rmsds)

        best_rmsd_rank = None
        rmsd_best_score = None
        rmsd_mean = float(np.mean(valid_rmsds)) if valid_rmsds else None
        rmsd_median = float(np.median(valid_rmsds)) if valid_rmsds else None
        rmsd_std = float(np.std(valid_rmsds)) if valid_rmsds else None
        near_native_fraction = None
        if pose_count and valid_rmsds:
            near_native_fraction = sum(r < threshold for r in valid_rmsds) / len(valid_rmsds)

        if valid_rmsds:
            best_rmsd_index = min(
                (i for i, r in enumerate(rmsds) if r is not None),
                key=lambda i: rmsds[i]
            )
            best_rmsd_rank = best_rmsd_index + 1

        score_rmsd_pearson = None
        score_rmsd_spearman = None
        if scores:
            score_rmsd_pearson, score_rmsd_spearman = self._score_rmsd_corr(scores, rmsds)
            valid_score_indices = [i for i, s in enumerate(scores) if s is not None]
            if valid_score_indices:
                best_score_index = min(valid_score_indices, key=lambda i: scores[i])
                rmsd_best_score = rmsds[best_score_index]

        return {
            "pose_count": pose_count,
            "best_score": best_score,
            "best_rmsd": best_rmsd,
            "top1_rmsd": top1_rmsd,
            "top5_rmsd": top5_rmsd,
            "top10_rmsd": top10_rmsd,
            "best_rmsd_rank": best_rmsd_rank,
            "rmsd_best_score": rmsd_best_score,
            "rmsd_mean": rmsd_mean,
            "rmsd_median": rmsd_median,
            "rmsd_std": rmsd_std,
            "near_native_fraction": near_native_fraction,
            "score_rmsd_pearson": score_rmsd_pearson,
            "score_rmsd_spearman": score_rmsd_spearman
        }

    def _compute_rescored_pose_metrics(
        self,
        crystal_ligand_pdb: Path,
        docked_file: Path,
        scores: List[float],
        threshold: float,
        has_reference_pose: bool,
    ) -> dict:
        """Calculate pose-recovery metrics after ranking poses by rescored values."""
        poses, _ = self._load_poses_and_scores(docked_file)
        count = min(len(poses), len(scores))
        if not count:
            return {}
        scores = [float(value) for value in scores[:count]]
        metrics = {
            "rescore_score": min(scores),
            "rescore_pose_count": count,
        }
        if not has_reference_pose:
            return metrics
        reference = self._load_reference_mol(crystal_ligand_pdb)
        if reference is None:
            return metrics
        rmsds = [self._pose_rmsd(reference, pose) for pose in poses[:count]]
        ranked_indices = sorted(range(count), key=lambda index: scores[index])
        ranked_rmsds = [rmsds[index] for index in ranked_indices]
        valid = [value for value in ranked_rmsds if value is not None]
        if not valid:
            return metrics

        def _top_best(limit: int) -> Optional[float]:
            values = [value for value in ranked_rmsds[:limit] if value is not None]
            return min(values) if values else None

        best_rmsd_index = min(
            (index for index, value in enumerate(rmsds) if value is not None),
            key=lambda index: rmsds[index],
        )
        metrics.update({
            "rescore_top1_rmsd": ranked_rmsds[0],
            "rescore_top5_rmsd": _top_best(5),
            "rescore_top10_rmsd": _top_best(10),
            "rescore_best_rmsd_rank": ranked_indices.index(best_rmsd_index) + 1,
            "rescore_rmsd_best_score": ranked_rmsds[0],
        })
        pearson, spearman = self._score_rmsd_corr(scores, rmsds)
        metrics["rescore_score_rmsd_pearson"] = pearson
        metrics["rescore_score_rmsd_spearman"] = spearman
        return metrics

    def _load_reference_mol(self, pdb_file: Path) -> Optional[Chem.Mol]:
        if not pdb_file.exists() or pdb_file.stat().st_size == 0:
            return None
        try:
            mol = Chem.MolFromPDBFile(str(pdb_file), removeHs=False)
            if mol is None:
                mol = Chem.MolFromPDBFile(str(pdb_file), removeHs=False, sanitize=False)
            return mol
        except Exception:
            return None

    def _load_poses_and_scores(self, docked_file: Path) -> Tuple[List[Chem.Mol], List[Optional[float]]]:
        suffix = docked_file.suffix.lower()
        poses = []
        scores: List[Optional[float]] = []

        if suffix in (".pdbqt",):
            scores = self._parse_pdbqt_scores(docked_file)
            from docking_platform_gui.utils.rmsd import _parse_pdbqt_models  # type: ignore
            for block in _parse_pdbqt_models(docked_file):
                mol = Chem.MolFromPDBBlock(block, removeHs=False)
                if mol is None:
                    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
                if mol is not None:
                    poses.append(mol)
            if len(scores) < len(poses):
                scores.extend([None] * (len(poses) - len(scores)))
        elif suffix in (".sd", ".sdf"):
            suppl = Chem.SDMolSupplier(str(docked_file), removeHs=False, sanitize=False)
            for mol in suppl:
                if mol is None:
                    continue
                poses.append(mol)
                props = mol.GetPropsAsDict()
                score = props.get("SCORE", props.get("SCORE.INTER"))
                scores.append(float(score) if score is not None else None)

        return poses, scores

    def _parse_pdbqt_scores(self, docked_file: Path) -> List[Optional[float]]:
        scores: List[Optional[float]] = []
        with open(docked_file, "r") as handle:
            for line in handle:
                # AutoDock Vina writes "REMARK VINA RESULT: <affinity> ..."
                if line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            scores.append(float(parts[3]))
                        except ValueError:
                            scores.append(None)
                # smina writes "REMARK minimizedAffinity <affinity>". Without
                # this branch every smina run yields no scores, which silently
                # blanks rmsd_best_score and both score-RMSD correlations.
                elif line.startswith("REMARK minimizedAffinity"):
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            scores.append(float(parts[2]))
                        except ValueError:
                            scores.append(None)
        return scores

    def _case_dir_from_output_file(self, output_file: Path) -> Path:
        if output_file.parent.name.startswith("Protocol"):
            candidate = output_file.parent.parent
        else:
            candidate = output_file.parent
        if candidate.parent and candidate.parent.name == "variants":
            candidate = candidate.parent.parent
        if (candidate / "crystal_ligand.pdb").exists():
            return candidate
        if candidate.parent and (candidate.parent / "crystal_ligand.pdb").exists():
            return candidate.parent
        return candidate

    def _ensure_receptor_pdb(self, case_dir: Path) -> Optional[Path]:
        receptor_pdb = case_dir / "receptor_prepared.pdb"
        if receptor_pdb.exists():
            return receptor_pdb
        receptor_pdbqt = case_dir / "receptor_prepared.pdbqt"
        if receptor_pdbqt.exists():
            lines = []
            with open(receptor_pdbqt, "r") as handle:
                for line in handle:
                    if line.startswith(("ATOM", "HETATM")):
                        lines.append(line[:66] + "\n")
                    elif line.startswith(("TER", "END")):
                        lines.append(line)
            if not lines or not lines[-1].startswith("END"):
                lines.append("END\n")
            receptor_pdb.write_text("".join(lines))
            return receptor_pdb
        return None

    def _write_ligand_pdb(
        self,
        mol: Chem.Mol,
        ligand_resname: str,
        output_pdb: Path,
        chain: str = "L",
        resnum: int = 1
    ) -> None:
        block = Chem.MolToPDBBlock(mol)
        resname = (ligand_resname or "LIG")[:3].ljust(3)
        chain_id = (chain or "L")[0]
        resnum_str = f"{resnum:>4}"
        out_lines = []
        element_counts: Dict[str, int] = {}
        atom_idx = 0
        for line in block.splitlines():
            if line.startswith("CONECT"):
                out_lines.append(line)
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            atom_idx += 1
            element = atom.GetSymbol().upper()
            element_counts[element] = element_counts.get(element, 0) + 1
            atom_name = f"{element}{element_counts[element]}"[-4:]
            line = list(line)
            while len(line) < 54:
                line.append(" ")
            line[12:16] = list(f"{atom_name:>4}")
            line[17:20] = list(resname)
            line[21] = chain_id
            line[22:26] = list(resnum_str)
            line[0:6] = list("HETATM")
            out_lines.append("".join(line))
        out_lines.append("END")
        output_pdb.write_text("\n".join(out_lines) + "\n")

    def _mol_to_photoimage(self, mol: Chem.Mol, size: Tuple[int, int]) -> Optional[tk.PhotoImage]:
        try:
            draw_mol = Chem.RemoveHs(Chem.Mol(mol))
            Chem.rdDepictor.Compute2DCoords(draw_mol)
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, draw_mol)
            drawer.FinishDrawing()
            png = drawer.GetDrawingText()
            if isinstance(png, str):
                png = png.encode("utf-8")
            encoded = base64.b64encode(png).decode("ascii")
            return tk.PhotoImage(data=encoded)
        except Exception as exc:
            logger.warning("Failed to render molecule image: {}", exc)
            return None

    def _mol_to_interaction_photoimage(
        self,
        mol: Chem.Mol,
        contacts: List[Dict[str, object]],
        size: Tuple[int, int],
        max_labels: int = 12
    ) -> Optional[tk.PhotoImage]:
        try:
            draw_mol = Chem.RemoveHs(Chem.Mol(mol))
            Chem.rdDepictor.Compute2DCoords(draw_mol)
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, draw_mol)

            atom_coords = [drawer.GetDrawCoords(idx) for idx in range(draw_mol.GetNumAtoms())]
            if atom_coords:
                centroid_x = sum(pt.x for pt in atom_coords) / len(atom_coords)
                centroid_y = sum(pt.y for pt in atom_coords) / len(atom_coords)
            else:
                centroid_x = 0.0
                centroid_y = 0.0

            for contact in contacts[:max_labels]:
                atom_idx = int(contact["atom_idx"])
                if atom_idx >= draw_mol.GetNumAtoms():
                    continue
                atom_pt = drawer.GetDrawCoords(atom_idx)
                dx = atom_pt.x - centroid_x
                dy = atom_pt.y - centroid_y
                norm = (dx * dx + dy * dy) ** 0.5
                if norm == 0:
                    dx, dy = 1.0, 0.0
                    norm = 1.0
                scale = 24.0
                label_pt = Point2D(
                    atom_pt.x + (dx / norm) * scale,
                    atom_pt.y + (dy / norm) * scale
                )
                drawer.DrawLine(atom_pt, label_pt)
                label = f"{contact['resname']}{contact['resnum']}"
                drawer.DrawString(label, label_pt)

            drawer.FinishDrawing()
            png = drawer.GetDrawingText()
            if isinstance(png, str):
                png = png.encode("utf-8")
            encoded = base64.b64encode(png).decode("ascii")
            return tk.PhotoImage(data=encoded)
        except Exception as exc:
            logger.warning("Failed to render interaction image: {}", exc)
            return None

    def _compute_contact_summary(
        self,
        receptor_pdb: Path,
        ligand_mol: Chem.Mol,
        cutoff: float = 4.0
    ) -> List[Dict[str, object]]:
        parser = PDB.PDBParser(QUIET=True)
        try:
            receptor = parser.get_structure("receptor", str(receptor_pdb))
        except Exception as exc:
            logger.warning("Failed to parse PDB for contacts: {}", exc)
            return []

        if ligand_mol is None or ligand_mol.GetNumAtoms() == 0:
            return []

        try:
            conf = ligand_mol.GetConformer()
        except Exception:
            return []

        heavy_indices = [
            idx for idx, atom in enumerate(ligand_mol.GetAtoms())
            if atom.GetSymbol() != "H"
        ]
        if not heavy_indices:
            return []

        lig_coords = np.array(
            [list(conf.GetAtomPosition(i)) for i in heavy_indices],
            dtype=float
        )

        contacts: List[Dict[str, object]] = []
        for residue in receptor.get_residues():
            if residue.id[0] != " ":
                continue
            min_dist = None
            min_idx = None
            for atom in residue.get_atoms():
                name = atom.get_name()
                element = (atom.element or "").strip()
                if element.upper() == "H" or name.startswith("H"):
                    continue
                diff = lig_coords - atom.coord
                distances = np.sqrt((diff * diff).sum(axis=1))
                local_min = float(distances.min())
                if min_dist is None or local_min < min_dist:
                    min_dist = local_min
                    min_idx = int(distances.argmin())
            if min_dist is not None and min_dist <= cutoff and min_idx is not None:
                chain_id = residue.get_parent().id
                contacts.append(
                    {
                        "chain": chain_id,
                        "resname": residue.resname,
                        "resnum": residue.id[1],
                        "dist": min_dist,
                        "atom_idx": min_idx
                    }
                )

        contacts.sort(key=lambda item: float(item["dist"]))
        return contacts

    def _resolve_pymol_bin(self) -> Optional[str]:
        pymol_bin = shutil.which("pymol")
        if pymol_bin:
            return pymol_bin
        fallback = Path("/Applications/PyMOL.app/Contents/MacOS/PyMOL")
        if fallback.exists():
            return str(fallback)
        return None

    def _resolve_ligplot_bin(self) -> Optional[str]:
        candidates: List[Path] = []
        configured = os.environ.get("LIGPLOT_BIN")
        if configured:
            candidates.append(Path(configured).expanduser())

        on_path = shutil.which("ligplot")
        if on_path:
            candidates.append(Path(on_path))

        roots = [
            os.environ.get("LIGPLUS_ROOT"),
            os.environ.get("LIGPLOT_HOME"),
            str(Path.home() / "LigPlus"),
            "/Applications/LigPlus",
        ]
        for raw_root in roots:
            if not raw_root:
                continue
            root = Path(raw_root).expanduser()
            candidates.extend(sorted((root / "lib").glob("exe_*/ligplot")))

        for candidate in candidates:
            if candidate.is_file():
                return str(candidate.resolve())
        return None

    def _resolve_ligplus_tool(self, ligplot_bin: str, tool_name: str) -> Optional[str]:
        ligplot_path = Path(ligplot_bin)
        candidate = ligplot_path.with_name(tool_name)
        if candidate.exists():
            return str(candidate)
        return None

    def _resolve_components_cif(self, ligplot_bin: Optional[str] = None) -> Optional[Path]:
        env_path = os.environ.get("HET_GROUP_DICTIONARY")
        if env_path:
            path = Path(env_path).expanduser()
            if path.is_file():
                return path.resolve()

        roots: List[Path] = []
        if ligplot_bin:
            ligplot_path = Path(ligplot_bin).expanduser()
            if len(ligplot_path.parents) >= 3:
                roots.append(ligplot_path.parents[2])
        for variable in ("LIGPLUS_ROOT", "LIGPLOT_HOME"):
            value = os.environ.get(variable)
            if value:
                roots.append(Path(value).expanduser())
        roots.extend((Path.home() / "LigPlus", Path("/Applications/LigPlus")))

        for root in roots:
            candidate = root / "lib" / "data" / "components.cif"
            if candidate.is_file():
                return candidate.resolve()
        return None

    def _combine_complex(self, receptor_pdb: Path, ligand_pdb: Path, output_pdb: Path) -> None:
        # LigPlot/HBPLUS identify atoms by PDB serial number. RDKit ligand PDBs
        # restart numbering at 1, so concatenating files directly makes ligand
        # CONECT records point at receptor atoms and fragments the ligand.
        combined = []
        next_serial = 1
        with open(receptor_pdb) as f_in:
            for line in f_in:
                if line.startswith(("ATOM", "HETATM")):
                    combined.append(f"{line[:6]}{next_serial:5d}{line[11:]}".rstrip("\n"))
                    next_serial += 1
                elif line.startswith("TER"):
                    combined.append("TER")

        ligand_atoms = []
        ligand_conect = []
        with open(ligand_pdb) as f_in:
            for line in f_in:
                if line.startswith(("ATOM", "HETATM")):
                    ligand_atoms.append(line.rstrip("\n"))
                elif line.startswith("CONECT"):
                    ligand_conect.append(line.rstrip("\n"))

        serial_map = {}
        for line in ligand_atoms:
            try:
                old_serial = int(line[6:11])
            except ValueError:
                continue
            serial_map[old_serial] = next_serial
            combined.append(f"{line[:6]}{next_serial:5d}{line[11:]}")
            next_serial += 1

        for line in ligand_conect:
            try:
                old_serials = [int(value) for value in line[6:].split()]
            except ValueError:
                continue
            new_serials = [serial_map[value] for value in old_serials if value in serial_map]
            if len(new_serials) >= 2:
                combined.append("CONECT" + "".join(f"{value:5d}" for value in new_serials))

        combined.append("END")
        output_pdb.write_text("\n".join(combined) + "\n")

    def _ligplot_ps_to_png(self, ps_file: Path, png_file: Path) -> bool:
        gs_bin = shutil.which("gs")
        if not gs_bin:
            return False
        cmd = [
            gs_bin,
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-sDEVICE=pngalpha",
            "-r200",
            f"-sOutputFile={png_file}",
            str(ps_file)
        ]
        subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return png_file.exists()

    def _run_ligplot(
        self,
        ligplot_bin: str,
        complex_pdb: Path,
        ligand_resname: str,
        resnum: int,
        chain: str,
        out_dir: Path
    ) -> Optional[Path]:
        out_dir.mkdir(parents=True, exist_ok=True)
        # Do not mistake files from an earlier successful render for the
        # current run when LigPlot or Ghostscript fails.
        for stale_name in ("ligplot.ps", "ligplot.png", "ligplot.sum"):
            stale_path = out_dir / stale_name
            if stale_path.exists():
                stale_path.unlink()
        ligplot_path = Path(ligplot_bin)
        ligplus_root = ligplot_path.parents[2] if len(ligplot_path.parents) >= 3 else None
        hbadd_bin = self._resolve_ligplus_tool(ligplot_bin, "hbadd")
        hbplus_bin = self._resolve_ligplus_tool(ligplot_bin, "hbplus")
        if ligplus_root:
            prm = ligplus_root / "lib" / "params" / "ligplot.prm"
            if prm.exists():
                target_prm = out_dir / "ligplot.prm"
                if not target_prm.exists():
                    target_prm.write_text(prm.read_text())

        env = os.environ.copy()
        env["PATH"] = f"{ligplot_path.parent}:{env.get('PATH', '')}"
        if ligplus_root:
            env["LIGPLUS_DIR"] = str(ligplus_root)
            env["LIGPLOT_PARAMS"] = str(ligplus_root / "lib" / "params")
        components_cif = self._resolve_components_cif(ligplot_bin)
        if components_cif:
            env["HET_GROUP_DICTIONARY"] = str(components_cif)

        if hbadd_bin and components_cif:
            subprocess.run([hbadd_bin, str(complex_pdb.resolve()), str(components_cif), "-wkdir", str(out_dir)],
                           check=False, cwd=str(out_dir), env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if hbplus_bin:
            work_pdb = out_dir / complex_pdb.name
            if work_pdb.resolve() != complex_pdb.resolve():
                shutil.copy2(complex_pdb, work_pdb)
            subprocess.run([hbplus_bin, work_pdb.name],
                           check=False, cwd=str(out_dir), env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        cmd = [
            ligplot_bin,
            str(complex_pdb.resolve()),
            ligand_resname,
            str(resnum),
            str(resnum),
            chain
        ]
        result = subprocess.run(
            cmd,
            check=False,
            cwd=str(out_dir),
            env=env,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            error = result.stderr.strip() or result.stdout.strip() or "unknown error"
            logger.warning("LigPlot failed for {}: {}", complex_pdb.name, error)
            return None
        ps_file = out_dir / "ligplot.ps"
        if not ps_file.exists():
            return None
        png_file = out_dir / "ligplot.png"
        if not self._ligplot_ps_to_png(ps_file, png_file):
            return None
        return png_file

    def _pml_quote(self, path: Path) -> str:
        normalized = str(path).replace("\\", "/")
        return f"\"{normalized}\""

    def _prepare_pymol_overlay(
        self,
        case: Dict[str, object]
    ) -> Optional[Path]:
        output_file = Path(str(case["output_file"]))
        case_dir = Path(str(case["case_dir"]))
        viewer_dir = Path(str(case["viewer_dir"])) / "pymol"
        viewer_dir.mkdir(parents=True, exist_ok=True)

        receptor_pdb = self._ensure_receptor_pdb(case_dir)
        if receptor_pdb is None or not receptor_pdb.exists():
            return None

        poses, scores = self._load_poses_and_scores(output_file)
        if not poses:
            return None

        crystal_ligand = case.get("crystal_ligand_pdb")
        native_ligand = Path(str(crystal_ligand)) if crystal_ligand else None
        if native_ligand is not None and not native_ligand.exists():
            native_ligand = None

        ref_mol = self._load_reference_mol(native_ligand) if native_ligand is not None else None
        if ref_mol is None:
            best_score_idx = None
            if scores:
                best_score_idx = min(
                    (i for i, s in enumerate(scores) if s is not None),
                    key=lambda i: scores[i],
                    default=None
                )
            if best_score_idx is None:
                best_score_idx = 0

            pose_file = viewer_dir / "best_score.pdb"
            self._write_ligand_pdb(poses[best_score_idx], str(case["display_name"]), pose_file)
            pml_lines = [
                "reinitialize",
                f"load {self._pml_quote(receptor_pdb)}, receptor",
                f"load {self._pml_quote(pose_file)}, docked",
                "hide everything",
                "show cartoon, receptor",
                "color grey80, receptor",
                "show sticks, docked",
                "color blue, docked",
                "set stick_radius, 0.22",
                "select binding_site, byres (receptor within 6 of docked)",
                "show sticks, binding_site",
                "color white, binding_site",
                "set bg_rgb, white",
                "zoom docked, 12",
                "set label_color, black",
                "set label_outline_color, white",
                "set label_size, 18",
                "pseudoatom legend_docked, pos=[0,0,0], label='Docked pose (blue); receptor cartoon grey; nearby residues white'",
                "hide everything, legend_docked",
                "show labels, legend_docked",
                "set label_screen_point, 1, legend_docked",
                "set label_position, [2,2,0], legend_docked"
            ]
            pml_path = viewer_dir / "overlay.pml"
            pml_path.write_text("\n".join(pml_lines) + "\n")
            return pml_path

        rmsd_values: List[Optional[float]] = [
            self._pose_rmsd(ref_mol, pose) for pose in poses
        ]
        best_rmsd_idx = None
        best_rmsd = None
        for idx, rmsd in enumerate(rmsd_values):
            if rmsd is None:
                continue
            if best_rmsd is None or rmsd < best_rmsd:
                best_rmsd = rmsd
                best_rmsd_idx = idx
        if best_rmsd_idx is None:
            best_rmsd_idx = 0

        best_score_idx = None
        if scores:
            best_score_idx = min(
                (i for i, s in enumerate(scores) if s is not None),
                key=lambda i: scores[i],
                default=None
            )
        if best_score_idx is None:
            best_score_idx = best_rmsd_idx

        pose_indices = []
        pose_entries: List[Tuple[int, str]] = []
        if best_rmsd_idx is not None:
            pose_entries.append((best_rmsd_idx, "best_rmsd"))
        if best_score_idx is not None:
            pose_entries.append((best_score_idx, "best_score"))

        pose_files = []
        for idx, label in pose_entries:
            pose_file = viewer_dir / f"{label}.pdb"
            aligned_pose = self._align_pose_to_ref(ref_mol, poses[idx])
            self._write_ligand_pdb(aligned_pose, str(case["display_name"]), pose_file)
            pose_files.append((label, pose_file))

        pml_lines = [
            "reinitialize",
            f"load {self._pml_quote(receptor_pdb)}, receptor",
            f"load {self._pml_quote(native_ligand)}, native",
            "hide everything",
            "show cartoon, receptor",
            "color grey80, receptor",
            "show sticks, native",
            "color green, native",
            "set stick_radius, 0.2",
            "select binding_site, byres (receptor within 6 of native)",
            "show sticks, binding_site",
            "color white, binding_site"
        ]

        color_map = {
            "best_rmsd": "red",
            "best_score": "blue"
        }

        for label, pose_file in pose_files:
            obj = label
            color = color_map.get(label, "orange")
            pml_lines.append(f"load {self._pml_quote(pose_file)}, {obj}")
            pml_lines.append(f"show sticks, {obj}")
            pml_lines.append(f"color {color}, {obj}")

        pml_lines.extend([
            "set bg_rgb, white",
            "zoom native, 12",
            "set label_color, black",
            "set label_outline_color, white",
            "set label_size, 18",
            "pseudoatom legend_native, pos=[0,0,0], label='Native (green)'",
            "pseudoatom legend_best_rmsd, pos=[0,0,0], label='Best RMSD (red)'",
            "pseudoatom legend_best_score, pos=[0,0,0], label='Best Score (blue)'",
            "group legend, legend_native legend_best_rmsd legend_best_score",
            "set label_position, [0,2.0,0], legend_native",
            "set label_position, [0,1.0,0], legend_best_rmsd",
            "set label_position, [0,0.0,0], legend_best_score",
            "hide everything, legend",
            "show labels, legend",
            "set label_screen_point, 1, legend",
            "set label_position, [2,2,0], legend",
            "set label_connector, 0"
        ])

        pml_path = viewer_dir / "overlay.pml"
        pml_path.write_text("\n".join(pml_lines) + "\n")
        return pml_path

    def _select_best_pose(
        self,
        crystal_ligand_pdb: Optional[Path],
        docked_file: Path,
        selection_mode: str = "best_rmsd"
    ) -> Optional[Tuple[Optional[Chem.Mol], Chem.Mol, Optional[float], Optional[float], int, int]]:
        poses, scores = self._load_poses_and_scores(docked_file)
        if not poses:
            return None

        ref_mol = self._load_reference_mol(crystal_ligand_pdb) if crystal_ligand_pdb else None
        if selection_mode == "best_score" or ref_mol is None:
            best_idx = None
            if scores:
                best_idx = min(
                    (i for i, s in enumerate(scores) if s is not None),
                    key=lambda i: scores[i],
                    default=None
                )
            if best_idx is None:
                best_idx = 0
            best_score = scores[best_idx] if scores and best_idx < len(scores) else None
            return ref_mol, poses[best_idx], None, best_score, best_idx, len(poses)

        if selection_mode != "best_rmsd":
            raise ValueError(f"Unknown pose selection mode: {selection_mode}")

        best_idx = None
        best_rmsd = None
        for idx, pose in enumerate(poses):
            rmsd = self._pose_rmsd(ref_mol, pose)
            if rmsd is None:
                continue
            if best_rmsd is None or rmsd < best_rmsd:
                best_rmsd = rmsd
                best_idx = idx

        if best_idx is None:
            best_idx = 0

        best_score = None
        if scores and best_idx < len(scores):
            best_score = scores[best_idx]

        return ref_mol, poses[best_idx], best_rmsd, best_score, best_idx, len(poses)

    @staticmethod
    def _pose_results_csv(results_path: Path) -> Path:
        results_path = Path(results_path)
        return (
            results_path if results_path.suffix.lower() == ".csv"
            else results_path.with_name("redock_results.csv")
        )

    def _show_pose_viewer(self, results_path: Path) -> None:
        csv_path = self._pose_results_csv(results_path)
        if not csv_path.exists():
            messagebox.showwarning("Pose viewer", "Results CSV not found.")
            return

        df = pd.read_csv(csv_path)
        if df.empty or "output_file" not in df.columns:
            messagebox.showwarning("Pose viewer", "No docked poses found in results.")
            return

        cases = []
        seen = set()
        viewer_root = Path(results_path).parent / "pose_viewer"
        for _, row in df.iterrows():
            output_file_val = row.get("output_file")
            if not isinstance(output_file_val, str):
                continue
            output_file = Path(output_file_val)
            if not output_file.exists():
                continue
            pdb_id = str(row.get("pdb_id", "")).upper()
            ligand = str(row.get("ligand_resname", "")).upper()
            dock_name = str(row.get("dock_name", "")).strip()
            display_name = dock_name or ligand
            key = (pdb_id, display_name, str(output_file))
            if key in seen:
                continue
            seen.add(key)
            case_dir = self._case_dir_from_output_file(output_file)
            crystal_ligand_pdb = case_dir / "crystal_ligand.pdb"
            if not crystal_ligand_pdb.exists():
                crystal_ligand_pdb = None
            cases.append(
                {
                    "pdb_id": pdb_id,
                    "ligand": ligand,
                    "display_name": display_name,
                    "output_file": output_file,
                    "crystal_ligand_pdb": crystal_ligand_pdb,
                    "case_dir": case_dir,
                    "viewer_dir": viewer_root / self._safe_case_id(f"{pdb_id}_{display_name}"),
                }
            )

        if not cases:
            messagebox.showwarning("Pose viewer", "No valid pose files found.")
            return

        case_labels = [f"{case['pdb_id']}_{case['display_name']}" for case in cases]
        case_label_map = {label.upper(): idx for idx, label in enumerate(case_labels)}

        dialog = tk.Toplevel(self)
        dialog.title("Docking Pose Viewer")
        dialog.geometry("980x640")
        dialog.transient(self)

        state = {"index": 0, "request_id": 0, "closed": False}

        def _on_close():
            state["closed"] = True
            dialog.destroy()

        dialog.protocol("WM_DELETE_WINDOW", _on_close)

        header = tk.Frame(dialog)
        header.pack(fill="x", padx=10, pady=10)

        info_var = tk.StringVar(value="Loading pose...")
        info_label = tk.Label(header, textvariable=info_var, anchor="w")
        info_label.pack(side="left", fill="x", expand=True)

        nav_frame = tk.Frame(header)
        nav_frame.pack(side="right")

        prev_btn = tk.Button(nav_frame, text="Prev", width=8)
        next_btn = tk.Button(nav_frame, text="Next", width=8)
        prev_btn.pack(side="left", padx=4)
        next_btn.pack(side="left", padx=4)
        pymol_btn = tk.Button(nav_frame, text="Open PyMOL", width=12)
        pymol_btn.pack(side="left", padx=4)
        ligplot_btn = tk.Button(nav_frame, text="LigPlot", width=10)
        ligplot_btn.pack(side="left", padx=4)

        jump_frame = tk.Frame(dialog)
        jump_frame.pack(fill="x", padx=10, pady=(0, 6))
        tk.Label(jump_frame, text="Jump to").pack(side="left")
        case_select_var = tk.StringVar(value=case_labels[0])
        case_select = ttk.Combobox(
            jump_frame,
            textvariable=case_select_var,
            values=case_labels,
            width=30
        )
        case_select.pack(side="left", padx=6)
        jump_btn = tk.Button(jump_frame, text="Go", width=6)
        jump_btn.pack(side="left")

        content = tk.Frame(dialog)
        content.pack(fill="both", expand=True, padx=10, pady=10)
        content.grid_columnconfigure(0, weight=1)
        content.grid_columnconfigure(1, weight=1)
        content.grid_rowconfigure(1, weight=1)

        tk.Label(content, text="Reference pose / status", font=("Helvetica", 12, "bold")).grid(
            row=0, column=0, sticky="w", padx=5, pady=(0, 6)
        )
        tk.Label(content, text="Best docked ligand (2D)", font=("Helvetica", 12, "bold")).grid(
            row=0, column=1, sticky="w", padx=5, pady=(0, 6)
        )

        native_label = tk.Label(content, bd=1, relief="solid")
        docked_label = tk.Label(content, bd=1, relief="solid")
        native_label.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        docked_label.grid(row=1, column=1, sticky="nsew", padx=5, pady=5)

        native_contact_frame = tk.Frame(content)
        docked_contact_frame = tk.Frame(content)
        native_contact_frame.grid(row=2, column=0, sticky="nsew", padx=5, pady=5)
        docked_contact_frame.grid(row=2, column=1, sticky="nsew", padx=5, pady=5)
        native_contact_frame.grid_columnconfigure(0, weight=1)
        docked_contact_frame.grid_columnconfigure(0, weight=1)

        tk.Label(native_contact_frame, text="Contacts <= 4.0Å").grid(row=0, column=0, sticky="w")
        tk.Label(docked_contact_frame, text="Contacts <= 4.0Å").grid(row=0, column=0, sticky="w")

        native_contact_text = tk.Text(native_contact_frame, height=8, wrap="none")
        docked_contact_text = tk.Text(docked_contact_frame, height=8, wrap="none")
        native_contact_text.grid(row=1, column=0, sticky="nsew")
        docked_contact_text.grid(row=1, column=0, sticky="nsew")
        native_scroll = ttk.Scrollbar(native_contact_frame, orient="vertical", command=native_contact_text.yview)
        docked_scroll = ttk.Scrollbar(docked_contact_frame, orient="vertical", command=docked_contact_text.yview)
        native_contact_text.configure(yscrollcommand=native_scroll.set)
        docked_contact_text.configure(yscrollcommand=docked_scroll.set)
        native_scroll.grid(row=1, column=1, sticky="ns")
        docked_scroll.grid(row=1, column=1, sticky="ns")

        detail_var = tk.StringVar(value="")
        detail_label = tk.Label(dialog, textvariable=detail_var, anchor="w")
        detail_label.pack(fill="x", padx=10, pady=(0, 10))

        def _render_case(idx: int) -> None:
            if idx < 0 or idx >= len(cases):
                return
            state["index"] = idx
            state["request_id"] += 1
            request_id = state["request_id"]
            case = cases[idx]
            case_select_var.set(case_labels[idx])
            info_var.set(
                f"Case {idx + 1}/{len(cases)}: {case['pdb_id']}_{case['display_name']} (loading)"
            )
            detail_var.set("")
            native_label.config(image="", text="Loading...")
            docked_label.config(image="", text="Loading...")
            native_contact_text.config(state="normal")
            docked_contact_text.config(state="normal")
            native_contact_text.delete("1.0", "end")
            docked_contact_text.delete("1.0", "end")
            native_contact_text.insert("1.0", "Loading...")
            docked_contact_text.insert("1.0", "Loading...")
            native_contact_text.config(state="disabled")
            docked_contact_text.config(state="disabled")

            def _worker():
                selected = self._select_best_pose(
                    case["crystal_ligand_pdb"],
                    case["output_file"]
                )
                if not selected:
                    return None
                receptor_pdb = self._ensure_receptor_pdb(case["case_dir"])
                if receptor_pdb is None:
                    return selected, None, None
                viewer_dir = case["viewer_dir"]
                viewer_dir.mkdir(parents=True, exist_ok=True)
                native_lig = viewer_dir / "native_ligand.pdb"
                docked_lig = viewer_dir / "docked_best_ligand.pdb"
                ref_mol, best_pose, best_rmsd, best_score, best_idx, pose_count = selected
                if ref_mol is not None:
                    self._write_ligand_pdb(ref_mol, case["display_name"], native_lig)
                self._write_ligand_pdb(best_pose, case["display_name"], docked_lig)
                native_contacts = self._compute_contact_summary(receptor_pdb, ref_mol) if ref_mol is not None else None
                docked_contacts = self._compute_contact_summary(receptor_pdb, best_pose)
                return selected, native_contacts, docked_contacts

            def _apply(payload):
                if state["closed"] or request_id != state["request_id"]:
                    return
                if not payload:
                    info_var.set(
                        f"Case {idx + 1}/{len(cases)}: {case['pdb_id']}_{case['display_name']} (no poses)"
                    )
                    native_label.config(text="No native pose", image="")
                    docked_label.config(text="No docked pose", image="")
                    native_contact_text.config(state="normal")
                    docked_contact_text.config(state="normal")
                    native_contact_text.delete("1.0", "end")
                    docked_contact_text.delete("1.0", "end")
                    native_contact_text.insert("1.0", "No contacts")
                    docked_contact_text.insert("1.0", "No contacts")
                    native_contact_text.config(state="disabled")
                    docked_contact_text.config(state="disabled")
                    return

                selected, native_contacts, docked_contacts = payload
                ref_mol, best_pose, best_rmsd, best_score, best_idx, pose_count = selected
                native_img = (
                    self._mol_to_interaction_photoimage(ref_mol, native_contacts, (420, 420))
                    if ref_mol is not None else None
                )
                docked_img = self._mol_to_interaction_photoimage(best_pose, docked_contacts, (420, 420))
                if native_img is None and ref_mol is not None:
                    native_img = self._mol_to_photoimage(ref_mol, (420, 420))
                if docked_img is None:
                    docked_img = self._mol_to_photoimage(best_pose, (420, 420))

                if native_img:
                    native_label.config(image=native_img, text="")
                    native_label.image = native_img
                elif ref_mol is None:
                    native_label.config(text="No crystal ligand\nscreening pose only", image="")
                else:
                    native_label.config(text="Native image unavailable", image="")

                if docked_img:
                    docked_label.config(image=docked_img, text="")
                    docked_label.image = docked_img
                else:
                    docked_label.config(text="Docked image unavailable", image="")

                rmsd_text = f"{best_rmsd:.2f}" if best_rmsd is not None else "N/A"
                score_text = f"{best_score:.2f}" if isinstance(best_score, (int, float)) else "N/A"
                info_var.set(
                    f"Case {idx + 1}/{len(cases)}: {case['pdb_id']}_{case['display_name']}"
                )
                detail_var.set(
                    f"Best pose {best_idx + 1}/{pose_count} | RMSD {rmsd_text} | Score {score_text}"
                )
                native_contact_text.config(state="normal")
                docked_contact_text.config(state="normal")
                native_contact_text.delete("1.0", "end")
                docked_contact_text.delete("1.0", "end")
                if ref_mol is None:
                    native_contact_text.insert("1.0", "No native/reference ligand for screening rows")
                elif native_contacts:
                    native_contact_text.insert(
                        "1.0",
                        "\n".join(
                            f"{item['chain']} {item['resname']} {int(item['resnum']):>4}  "
                            f"{float(item['dist']):5.2f}Å"
                            for item in native_contacts[:40]
                        )
                    )
                else:
                    native_contact_text.insert("1.0", "No contacts <= 4.0Å")
                if docked_contacts:
                    docked_contact_text.insert(
                        "1.0",
                        "\n".join(
                            f"{item['chain']} {item['resname']} {int(item['resnum']):>4}  "
                            f"{float(item['dist']):5.2f}Å"
                            for item in docked_contacts[:40]
                        )
                    )
                else:
                    docked_contact_text.insert("1.0", "No contacts <= 4.0Å")
                native_contact_text.config(state="disabled")
                docked_contact_text.config(state="disabled")

            def _run():
                try:
                    payload = _worker()
                except Exception as exc:
                    logger.warning("Pose viewer failed: {}", exc)
                    payload = None
                self._run_on_ui(lambda: _apply(payload))

            threading.Thread(target=_run, daemon=True).start()

        def _prev():
            new_idx = max(0, state["index"] - 1)
            _render_case(new_idx)

        def _next():
            new_idx = min(len(cases) - 1, state["index"] + 1)
            _render_case(new_idx)

        def _jump_to_label(label: str) -> None:
            raw = label.strip()
            if not raw:
                return
            key = raw.upper().replace(" ", "_")
            idx = case_label_map.get(key)
            if idx is None:
                for i, case_label in enumerate(case_labels):
                    if case_label.upper().startswith(key):
                        idx = i
                        break
            if idx is None:
                info_var.set(f"No match for '{raw}'")
                return
            _render_case(idx)

        def _on_jump(_event=None):
            _jump_to_label(case_select_var.get())

        prev_btn.config(command=_prev)
        next_btn.config(command=_next)
        jump_btn.config(command=_on_jump)
        case_select.bind("<<ComboboxSelected>>", _on_jump)
        case_select.bind("<Return>", _on_jump)

        def _open_pymol():
            case = cases[state["index"]]
            pymol_bin = self._resolve_pymol_bin()
            if not pymol_bin:
                messagebox.showwarning(
                    "PyMOL not found",
                    "PyMOL executable not found in PATH. Install PyMOL or add it to PATH."
                )
                return

            info_var.set(f"Launching PyMOL for {case['pdb_id']}_{case['display_name']}...")

            def _worker():
                try:
                    pml_path = self._prepare_pymol_overlay(case)
                    if not pml_path:
                        return None
                    subprocess.Popen([pymol_bin, str(pml_path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    return pml_path
                except Exception as exc:
                    logger.warning("Failed to launch PyMOL: {}", exc)
                    return None

            def _apply(pml_path):
                if state["closed"]:
                    return
                if pml_path:
                    info_var.set(f"Opened PyMOL overlay for {case['pdb_id']}_{case['display_name']}")
                else:
                    info_var.set(f"Failed to open PyMOL for {case['pdb_id']}_{case['display_name']}")

            def _run():
                result = _worker()
                self._run_on_ui(lambda: _apply(result))

            threading.Thread(target=_run, daemon=True).start()

        pymol_btn.config(command=_open_pymol)

        def _open_ligplot():
            case = cases[state["index"]]
            ligplot_bin = self._resolve_ligplot_bin()
            if not ligplot_bin:
                messagebox.showwarning(
                    "LigPlot not found",
                    "LigPlot executable not found. Please check installation."
                )
                return

            info_var.set(f"Generating LigPlot for {case['pdb_id']}_{case['display_name']}...")

            def _worker():
                selected = self._select_best_pose(
                    case["crystal_ligand_pdb"],
                    case["output_file"],
                    selection_mode="best_score"
                )
                if not selected:
                    return None

                receptor_pdb = self._ensure_receptor_pdb(case["case_dir"])
                if receptor_pdb is None or not receptor_pdb.exists():
                    return None

                ref_mol, best_pose, _, _, _, _ = selected
                viewer_dir = Path(case["viewer_dir"])
                viewer_dir.mkdir(parents=True, exist_ok=True)

                # A stable generic name avoids collisions with unrelated CCD
                # entries when display names happen to share a three-letter code.
                ligand_resname = "LIG"
                ligand_chain = "L"
                ligand_resnum = 1

                native_lig = viewer_dir / "native_ligand.pdb"
                docked_lig = viewer_dir / "docked_best_ligand.pdb"
                self._write_ligand_pdb(best_pose, ligand_resname, docked_lig, chain=ligand_chain, resnum=ligand_resnum)

                native_dir = viewer_dir / "ligplot_native"
                docked_dir = viewer_dir / "ligplot_docked"
                docked_dir.mkdir(parents=True, exist_ok=True)
                docked_complex = docked_dir / "docked_complex.pdb"

                native_png = None
                if ref_mol is not None:
                    native_dir.mkdir(parents=True, exist_ok=True)
                    native_complex = native_dir / "native_complex.pdb"
                    self._write_ligand_pdb(ref_mol, ligand_resname, native_lig, chain=ligand_chain, resnum=ligand_resnum)
                    self._combine_complex(receptor_pdb, native_lig, native_complex)
                    native_png = self._run_ligplot(
                        ligplot_bin,
                        native_complex,
                        ligand_resname,
                        ligand_resnum,
                        ligand_chain,
                        native_dir
                    )

                self._combine_complex(receptor_pdb, docked_lig, docked_complex)

                docked_png = self._run_ligplot(
                    ligplot_bin,
                    docked_complex,
                    ligand_resname,
                    ligand_resnum,
                    ligand_chain,
                    docked_dir
                )
                return native_png, docked_png

            def _apply(result):
                if not result or state["closed"]:
                    info_var.set("LigPlot generation failed")
                    return
                native_png, docked_png = result
                if not docked_png:
                    info_var.set("LigPlot generation failed")
                    return

                dialog_lp = tk.Toplevel(dialog)
                dialog_lp.title("LigPlot Interaction Diagrams")
                dialog_lp.geometry("900x540")
                dialog_lp.transient(dialog)

                frame = tk.Frame(dialog_lp)
                frame.pack(fill="both", expand=True, padx=10, pady=10)
                frame.grid_columnconfigure(0, weight=1)
                frame.grid_columnconfigure(1, weight=1)
                frame.grid_rowconfigure(1, weight=1)

                header = tk.Frame(frame)
                header.grid(row=0, column=0, columnspan=2, sticky="ew", padx=5, pady=(0, 6))
                header.grid_columnconfigure(1, weight=1)
                tk.Label(header, text="Native", font=("Helvetica", 12, "bold")).grid(row=0, column=0, sticky="w")
                tk.Label(header, text="Docked", font=("Helvetica", 12, "bold")).grid(row=0, column=2, sticky="w")

                def _fit_image(img: tk.PhotoImage, max_w: int, max_h: int) -> tk.PhotoImage:
                    width = img.width()
                    height = img.height()
                    if width <= 0 or height <= 0:
                        return img
                    scale_w = math.ceil(width / max_w) if width > max_w else 1
                    scale_h = math.ceil(height / max_h) if height > max_h else 1
                    scale = max(scale_w, scale_h)
                    return img.subsample(scale) if scale > 1 else img

                native_img_raw = tk.PhotoImage(file=str(native_png)) if native_png else None
                docked_img_raw = tk.PhotoImage(file=str(docked_png))

                scale_steps = [50, 75, 100, 125, 150, 200]
                scale_map = {
                    50: (1, 2),
                    75: (3, 4),
                    100: (1, 1),
                    125: (5, 4),
                    150: (3, 2),
                    200: (2, 1),
                }
                scale_idx = {"value": scale_steps.index(100)}

                def _scale_image(img: tk.PhotoImage, percent: int) -> tk.PhotoImage:
                    zoom, subsample = scale_map.get(percent, (1, 1))
                    if zoom > 1:
                        img = img.zoom(zoom, zoom)
                    if subsample > 1:
                        img = img.subsample(subsample, subsample)
                    return img

                def _apply_scale(percent: int) -> None:
                    native_img = _scale_image(native_img_raw, percent) if native_img_raw else None
                    docked_img = _scale_image(docked_img_raw, percent)
                    if native_img:
                        native_lbl.config(image=native_img, text="")
                    else:
                        native_lbl.config(image="", text="No native/reference\nligand for screening")
                    docked_lbl.config(image=docked_img)
                    native_lbl.image = native_img
                    docked_lbl.image = docked_img
                    scale_var.set(f"{percent}%")

                native_lbl = tk.Label(frame, bd=1, relief="solid")
                docked_lbl = tk.Label(frame, bd=1, relief="solid")
                native_lbl.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
                docked_lbl.grid(row=1, column=1, sticky="nsew", padx=5, pady=5)

                controls = tk.Frame(frame)
                controls.grid(row=2, column=0, columnspan=2, sticky="ew", padx=5, pady=(6, 0))
                scale_var = tk.StringVar(value="100%")
                tk.Label(controls, text="Zoom").pack(side="left")
                tk.Button(
                    controls,
                    text="-",
                    width=3,
                    command=lambda: _zoom(-1)
                ).pack(side="left", padx=4)
                tk.Button(
                    controls,
                    text="+",
                    width=3,
                    command=lambda: _zoom(1)
                ).pack(side="left", padx=4)
                tk.Label(controls, textvariable=scale_var, width=6).pack(side="left", padx=(4, 12))
                tk.Button(
                    controls,
                    text="Fit",
                    command=lambda: _fit()
                ).pack(side="left")

                def _fit() -> None:
                    native_img = _fit_image(native_img_raw, 420, 420) if native_img_raw else None
                    docked_img = _fit_image(docked_img_raw, 420, 420)
                    if native_img:
                        native_lbl.config(image=native_img, text="")
                    else:
                        native_lbl.config(image="", text="No native/reference\nligand for screening")
                    docked_lbl.config(image=docked_img)
                    native_lbl.image = native_img
                    docked_lbl.image = docked_img
                    scale_var.set("Fit")

                def _zoom(delta: int) -> None:
                    idx = max(0, min(len(scale_steps) - 1, scale_idx["value"] + delta))
                    scale_idx["value"] = idx
                    _apply_scale(scale_steps[idx])

                _apply_scale(100)
                native_lbl.image_raw = native_img_raw
                docked_lbl.image_raw = docked_img_raw

                info_var.set(f"LigPlot generated for {case['pdb_id']}_{case['display_name']}")

            def _run():
                try:
                    result = _worker()
                except Exception as exc:
                    logger.warning("LigPlot failed: {}", exc)
                    result = None
                self._run_on_ui(lambda: _apply(result))

            threading.Thread(target=_run, daemon=True).start()

        ligplot_btn.config(command=_open_ligplot)

        _render_case(0)

    def _pose_rmsd(self, ref_mol: Chem.Mol, pose_mol: Chem.Mol) -> Optional[float]:
        from docking_platform_gui.utils.rmsd import coordinate_rmsd
        try:
            return coordinate_rmsd(ref_mol, pose_mol, use_symmetry=True)
        except Exception:
            return None

    def _align_pose_to_ref(self, ref_mol: Chem.Mol, pose_mol: Chem.Mol) -> Chem.Mol:
        aligned = Chem.Mol(pose_mol)
        try:
            rdMolAlign.AlignMol(aligned, ref_mol)
            return aligned
        except Exception:
            pass

        try:
            ref_noh = Chem.RemoveHs(ref_mol, sanitize=False)
            pose_noh = Chem.RemoveHs(aligned, sanitize=False)
            min_atoms = min(ref_noh.GetNumAtoms(), pose_noh.GetNumAtoms())
            if min_atoms < 3:
                return aligned
            mcs = rdFMCS.FindMCS(
                [ref_noh, pose_noh],
                ringMatchesRingOnly=False,
                completeRingsOnly=False,
                matchValences=False,
                bondCompare=rdFMCS.BondCompare.CompareAny,
                atomCompare=rdFMCS.AtomCompare.CompareElements
            )
            min_required = max(3, int(min_atoms * 0.7))
            if mcs.numAtoms < min_required:
                return aligned
            pattern = Chem.MolFromSmarts(mcs.smartsString)
            if pattern is None:
                return aligned
            ref_match = ref_noh.GetSubstructMatch(pattern)
            pose_match = pose_noh.GetSubstructMatch(pattern)
            if not ref_match or not pose_match:
                return aligned
            ref_map = []
            pose_map = []
            for ref_idx, pose_idx in zip(ref_match, pose_match):
                ref_atom = ref_noh.GetAtomWithIdx(ref_idx)
                pose_atom = pose_noh.GetAtomWithIdx(pose_idx)
                if ref_atom.HasProp("_origAtomIdx"):
                    ref_map.append(ref_atom.GetIntProp("_origAtomIdx"))
                else:
                    ref_map.append(ref_idx)
                if pose_atom.HasProp("_origAtomIdx"):
                    pose_map.append(pose_atom.GetIntProp("_origAtomIdx"))
                else:
                    pose_map.append(pose_idx)
            atom_map = list(zip(pose_map, ref_map))
            rdMolAlign.AlignMol(aligned, ref_mol, atomMap=atom_map)
        except Exception:
            return aligned
        return aligned

    def _score_rmsd_corr(
        self,
        scores: List[Optional[float]],
        rmsds: List[Optional[float]]
    ) -> Tuple[Optional[float], Optional[float]]:
        pairs = [
            (score, rmsd)
            for score, rmsd in zip(scores, rmsds)
            if score is not None and rmsd is not None
        ]
        if len(pairs) < 2:
            return None, None

        s_vals = np.array([p[0] for p in pairs], dtype=float)
        r_vals = np.array([p[1] for p in pairs], dtype=float)

        pearson = float(np.corrcoef(s_vals, r_vals)[0, 1])

        s_rank = s_vals.argsort().argsort()
        r_rank = r_vals.argsort().argsort()
        spearman = float(np.corrcoef(s_rank, r_rank)[0, 1])

        return pearson, spearman

    def _compute_enrichment_metrics(
        self,
        score_labels: List[Tuple[float, int]],
        min_fpr: float = 0.001
    ) -> Optional[dict]:
        if not score_labels:
            return None

        pos = sum(1 for _, label in score_labels if label == 1)
        neg = sum(1 for _, label in score_labels if label == 0)
        if pos == 0 or neg == 0:
            return None

        ranked = sorted(score_labels, key=lambda item: item[0], reverse=True)
        tpr = [0.0]
        fpr = [0.0]
        tp = 0
        fp = 0
        average_precision = 0.0
        index = 0
        while index < len(ranked):
            score = ranked[index][0]
            tied_labels = []
            while index < len(ranked) and ranked[index][0] == score:
                tied_labels.append(ranked[index][1])
                index += 1
            # Advance tied observations together. Ordering equal scores one by
            # one makes AUC depend on spreadsheet order rather than ranking.
            previous_tp = tp
            tp += sum(label == 1 for label in tied_labels)
            fp += sum(label == 0 for label in tied_labels)
            tpr.append(tp / pos)
            fpr.append(fp / neg)
            average_precision += ((tp - previous_tp) / pos) * (tp / (tp + fp))

        # np.trapz was removed in NumPy 2.0 (renamed to np.trapezoid).
        # This environment runs NumPy 2.5, where the old name raises
        # AttributeError and takes the whole summary down.
        _trap = getattr(np, "trapezoid", None) or np.trapz
        roc_auc = float(_trap(tpr, fpr))

        min_fpr = max(min_fpr, 1e-6)
        fpr_grid = np.logspace(np.log10(min_fpr), 0.0, num=200)
        tpr_interp = np.interp(fpr_grid, fpr, tpr)
        log_fpr = np.log10(fpr_grid)
        auc_log = float(_trap(tpr_interp, log_fpr))
        auc_log_random = float(_trap(fpr_grid, log_fpr))
        log_range = 0.0 - np.log10(min_fpr)
        log_auc = (auc_log - auc_log_random) / log_range * 100.0

        return {
            "roc_auc": roc_auc,
            "log_auc": log_auc,
            "average_precision": float(average_precision),
            "active_prevalence": float(pos / (pos + neg)),
            "actives": pos,
            "decoys": neg
        }

    @staticmethod
    def _assay_benchmark_chart_data(
        score_labels: List[Tuple[float, int]], bins: int = 20
    ) -> Optional[dict]:
        """Build tie-aware curves and class-normalized score distributions."""
        valid = [
            (float(score), int(label)) for score, label in score_labels
            if label in (0, 1) and np.isfinite(score)
        ]
        n_actives = sum(label == 1 for _, label in valid)
        n_inactives = sum(label == 0 for _, label in valid)
        if not valid or not n_actives or not n_inactives:
            return None

        ranked = sorted(valid, key=lambda item: item[0], reverse=True)
        roc_curve = [[0.0, 0.0]]
        precision_recall_curve = [[0.0, 1.0]]
        cumulative_recovery_curve = [[0.0, 0.0]]
        tp = fp = screened = 0
        index = 0
        while index < len(ranked):
            score = ranked[index][0]
            tied_labels = []
            while index < len(ranked) and ranked[index][0] == score:
                tied_labels.append(ranked[index][1])
                index += 1
            tp += sum(label == 1 for label in tied_labels)
            fp += sum(label == 0 for label in tied_labels)
            screened += len(tied_labels)
            recall = tp / n_actives
            roc_curve.append([fp / n_inactives, recall])
            precision_recall_curve.append([recall, tp / (tp + fp)])
            cumulative_recovery_curve.append([screened / len(ranked), recall])

        score_values = np.asarray([score for score, _ in valid], dtype=float)
        score_min = float(np.min(score_values))
        score_max = float(np.max(score_values))
        if math.isclose(score_min, score_max):
            score_min -= 0.5
            score_max += 0.5
        bin_count = max(1, min(int(bins), len(valid)))
        edges = np.linspace(score_min, score_max, bin_count + 1)
        active_scores = [score for score, label in valid if label == 1]
        inactive_scores = [score for score, label in valid if label == 0]
        active_counts, _ = np.histogram(active_scores, bins=edges)
        inactive_counts, _ = np.histogram(inactive_scores, bins=edges)

        return {
            "roc_curve": roc_curve,
            "precision_recall_curve": precision_recall_curve,
            "cumulative_recovery_curve": cumulative_recovery_curve,
            "score_histogram": {
                "edges": [float(value) for value in edges],
                "active_percent": [
                    float(100.0 * value / n_actives) for value in active_counts
                ],
                "inactive_percent": [
                    float(100.0 * value / n_inactives) for value in inactive_counts
                ],
            },
            "actives": n_actives,
            "inactives": n_inactives,
            "prevalence": float(n_actives / len(valid)),
        }

    @staticmethod
    def _tie_aware_enrichment_factor(
        score_labels: List[Tuple[float, int]],
        percent: float,
    ) -> Optional[dict]:
        """Return expected EF and its range when the cutoff splits a score tie."""
        if not score_labels:
            return None
        n_total = len(score_labels)
        n_actives = sum(label == 1 for _, label in score_labels)
        if n_actives == 0:
            return None

        ranked = sorted(score_labels, key=lambda item: item[0], reverse=True)
        n_select = max(1, int(math.ceil(n_total * percent / 100.0)))
        cutoff_score = ranked[n_select - 1][0]
        strict_labels = [label for score, label in ranked if score > cutoff_score]
        tied_labels = [label for score, label in ranked if score == cutoff_score]
        slots = n_select - len(strict_labels)
        strict_actives = sum(label == 1 for label in strict_labels)
        tied_actives = sum(label == 1 for label in tied_labels)
        tied_inactives = len(tied_labels) - tied_actives

        expected_actives = strict_actives + slots * tied_actives / len(tied_labels)
        minimum_actives = strict_actives + max(0, slots - tied_inactives)
        maximum_actives = strict_actives + min(slots, tied_actives)
        baseline_rate = n_actives / n_total

        def _ef(active_count: float) -> float:
            return float((active_count / n_select) / baseline_rate)

        return {
            "percent": float(percent),
            "selected": n_select,
            "cutoff_rank_score": float(cutoff_score),
            "cutoff_tie_size": len(tied_labels),
            "cutoff_slots": slots,
            "expected_actives": float(expected_actives),
            "enrichment_factor": _ef(expected_actives),
            "enrichment_factor_min": _ef(minimum_actives),
            "enrichment_factor_max": _ef(maximum_actives),
        }

    def _score_charge_corr(
        self,
        scores: List[float],
        charges: List[int]
    ) -> Tuple[Optional[float], Optional[float]]:
        if len(scores) < 2:
            return None, None
        s_vals = np.array(scores, dtype=float)
        c_vals = np.array(charges, dtype=float)
        if np.std(s_vals) == 0 or np.std(c_vals) == 0:
            return None, None
        pearson = float(np.corrcoef(s_vals, c_vals)[0, 1])
        s_rank = s_vals.argsort().argsort()
        c_rank = c_vals.argsort().argsort()
        spearman = float(np.corrcoef(s_rank, c_rank)[0, 1])
        return pearson, spearman

    def _rank_score_value(self, result: "RedockResult") -> Optional[float]:
        if result.rescore_cnn_affinity is not None:
            # GNINA CNNaffinity is a predicted pK; larger values rank better.
            return float(result.rescore_cnn_affinity)
        if result.rescore_cnn_score is not None:
            return float(result.rescore_cnn_score)
        if result.rescore_score is not None:
            return -float(result.rescore_score)
        if result.best_score is not None:
            return -float(result.best_score)
        return None

    @staticmethod
    def _selected_score_details(
        result: "RedockResult",
    ) -> Optional[Tuple[float, float, str, str]]:
        """Return normalized rank, displayed score, source, and direction."""
        if result.rescore_cnn_affinity is not None:
            value = float(result.rescore_cnn_affinity)
            return value, value, "GNINA CNN affinity", "higher"
        if result.rescore_cnn_score is not None:
            value = float(result.rescore_cnn_score)
            return value, value, "GNINA CNN score", "higher"
        if result.rescore_score is not None:
            value = float(result.rescore_score)
            method = str(result.rescore_method or "Smina rescore")
            if method.startswith("smina_score_only:"):
                method = f"{method.split(':', 1)[1]} (Smina score-only)"
            return -value, value, method, "lower"
        if result.best_score is not None:
            value = float(result.best_score)
            return -value, value, f"{result.engine or 'Docking'} score", "lower"
        return None

    @staticmethod
    def _property_matched(active: RedockResult, decoy: RedockResult) -> Optional[bool]:
        required = (
            active.molecular_weight, decoy.molecular_weight,
            active.logp, decoy.logp, active.tpsa, decoy.tpsa,
            active.rotatable_bonds, decoy.rotatable_bonds,
            active.ligand_charge, decoy.ligand_charge,
        )
        if any(value is None for value in required):
            return None
        mw_ratio = decoy.molecular_weight / active.molecular_weight
        return bool(
            0.75 <= mw_ratio <= 1.25
            and abs(decoy.logp - active.logp) <= 1.5
            and abs(decoy.tpsa - active.tpsa) <= 40.0
            and abs(decoy.rotatable_bonds - active.rotatable_bonds) <= 3
            and decoy.ligand_charge == active.ligand_charge
        )

    def _build_summary(self, results: List[RedockResult], threshold: float) -> dict:
        """
        Build summary statistics using the enhanced RedockAnalyzer.
        
        This properly separates actives (self-docking) from decoys (enrichment).
        
        Args:
            results: List of RedockResult objects
            threshold: RMSD threshold for success (default 2.0Å)
            
        Returns:
            Dictionary with comprehensive statistics separated by active/decoy
        """
        if not results:
            return {
                "total_cases": 0,
                "threshold": threshold,
                "n_actives": 0,
                "n_decoys": 0,
                "success_rate_best": 0.0,
                "mean_best_rmsd": None,
                "by_protocol": {},
                "by_engine": {},
            }
        
        # Convert RedockResult to DataFrame for analyzer
        results_data = []
        for r in results:
            results_data.append({
                'pdb_id': r.pdb_id,
                'ligand_resname': r.ligand_resname,
                'ligand_chain': r.ligand_chain,
                'mode': r.mode,
                'engine': r.engine,
                'protocol': r.protocol,
                'best_rmsd': r.best_rmsd,
                'success': r.success,
                'runtime_sec': r.runtime_sec,
                'output_file': r.output_file,
                'pose_count': r.pose_count,
                'best_score': r.best_score,
                'dock_name': r.dock_name,
                'top1_rmsd': r.top1_rmsd,
                'top5_rmsd': r.top5_rmsd,
                'top10_rmsd': r.top10_rmsd,
                'best_rmsd_rank': r.best_rmsd_rank,
                'rmsd_best_score': r.rmsd_best_score,
                'rmsd_mean': r.rmsd_mean,
                'rmsd_median': r.rmsd_median,
                'rmsd_std': r.rmsd_std,
                'near_native_fraction': r.near_native_fraction,
                'score_rmsd_pearson': r.score_rmsd_pearson,
                'score_rmsd_spearman': r.score_rmsd_spearman,
                'control_label': r.control_label,
                'ligand_charge': r.ligand_charge,
                'molecular_weight': r.molecular_weight,
                'logp': r.logp,
                'tpsa': r.tpsa,
                'rotatable_bonds': r.rotatable_bonds,
                'site_method': r.site_method,
            })
        
        df = pd.DataFrame(results_data)
        
        try:
            # Use the enhanced analyzer
            analyzer = RedockAnalyzer(rmsd_threshold=threshold)
            report = analyzer.generate_report(df)
            
            # Extract statistics
            active_stats = report['active_analysis']
            enrichment_stats = report['enrichment_analysis']
            
            # Build summary dict in format GUI expects
            summary = {
                # Basic info
                "total_cases": len(results),
                "threshold": threshold,
                
                # NEW: Classification
                "n_actives": report['summary']['n_actives'],
                "n_decoys": report['summary']['n_decoys'],
                
                # Active (self-docking) metrics - NOW CORRECT!
                "success_rate_best": (
                    active_stats['success_rate'] 
                    if active_stats['n_valid_rmsd'] > 0 
                    else 0.0
                ),
                "mean_best_rmsd": active_stats['mean_rmsd'],
                "median_best_rmsd": active_stats['median_rmsd'],
                "min_rmsd": active_stats['min_rmsd'],
                "max_rmsd": active_stats['max_rmsd'],
                "mean_runtime_sec": active_stats['mean_runtime'],
                
                # Pose-level metrics. The enhanced analyzer does not produce
                # these, but _compute_pose_metrics already stores them on each
                # RedockResult, so aggregate them from the per-case results
                # rather than reporting N/A.
                **self._aggregate_pose_metrics(results, threshold),
                
                # NEW: Enrichment metrics
                "control_actives": enrichment_stats['n_actives'],
                "control_decoys": enrichment_stats['n_decoys'],
                "roc_auc": enrichment_stats['roc_auc'],
                "ef_1_percent": enrichment_stats.get('ef_1_percent'),
                "ef_5_percent": enrichment_stats.get('ef_5_percent'),
                "ef_10_percent": enrichment_stats.get('ef_10_percent'),
                "mean_active_score": enrichment_stats.get('mean_active_score'),
                "mean_decoy_score": enrichment_stats.get('mean_decoy_score'),
                "score_separation": enrichment_stats.get('score_separation'),
                
                # Legacy fields
                "log_auc": None,
                "rescore_count": None,
                "rescore_failures": None,
                "rescore_methods": None,
                "charge_count": None,
                "mean_charge": None,
                "median_charge": None,
                "charge_frac_positive": None,
                "charge_frac_negative": None,
                "charge_frac_neutral": None,
                "charge_top10_frac_positive": None,
                "charge_top10_frac_negative": None,
                "charge_top10_frac_neutral": None,
                "score_charge_pearson": None,
                "score_charge_spearman": None,
                
                # NEW: Interpretation
                "interpretation": report['interpretation'],
                
                # Protocol/engine breakdown
                "by_protocol": {},
                "by_engine": {},
                "attempts_by_protocol": {}
            }
            
        except Exception as e:
            logger.error(f"Enhanced analysis failed: {e}, falling back to basic stats")
            # Fallback to basic stats if enhanced analysis fails
            summary = {
                "total_cases": len(results),
                "threshold": threshold,
                "success_rate_best": 0.0,
                "mean_best_rmsd": None,
                "by_protocol": {},
                "by_engine": {},
            }

        # Override/augment enrichment from explicit control labels.
        # This supports apo validation where RMSD is intentionally unavailable.
        score_labels: List[Tuple[float, int]] = []
        summary["n_actives"] = sum(r.control_label == 1 for r in results)
        summary["n_decoys"] = sum(r.control_label == 0 for r in results)
        summary["n_samples"] = sum(
            r.mode == "screening" and r.control_label is None for r in results
        )
        for result in results:
            if result.control_label not in (0, 1):
                continue
            rank_score = self._rank_score_value(result)
            if rank_score is None:
                continue
            score_labels.append((rank_score, int(result.control_label)))

        pose_results_by_target: Dict[str, List[RedockResult]] = {}
        for result in results:
            if result.control_label == 0 or result.best_rmsd is None or result.best_rmsd >= 900:
                continue
            target_name = result.target_name or result.pdb_id
            pose_results_by_target.setdefault(target_name, []).append(result)

        def _target_pose_rate(target_results: List[RedockResult], field: str) -> Optional[float]:
            values = [
                getattr(result, field) for result in target_results
                if getattr(result, field) is not None and getattr(result, field) < 900
            ]
            if not values:
                return None
            return float(100.0 * sum(value < threshold for value in values) / len(values))

        target_pose_rows = []
        for target_name, target_results in sorted(pose_results_by_target.items()):
            target_pose_rows.append({
                "target_name": target_name,
                "cases": len(target_results),
                "success_rate_best": _target_pose_rate(target_results, "best_rmsd"),
                "success_rate_top1": _target_pose_rate(target_results, "top1_rmsd"),
                "success_rate_top5": _target_pose_rate(target_results, "top5_rmsd"),
                "success_rate_top10": _target_pose_rate(target_results, "top10_rmsd"),
                "mean_best_rmsd": float(np.mean([
                    result.best_rmsd for result in target_results
                ])),
                "mean_top1_rmsd": (
                    float(np.mean([
                        result.top1_rmsd for result in target_results
                        if result.top1_rmsd is not None and result.top1_rmsd < 900
                    ]))
                    if any(
                        result.top1_rmsd is not None and result.top1_rmsd < 900
                        for result in target_results
                    ) else None
                ),
            })
        summary["per_target_pose_recovery"] = target_pose_rows

        # Per-structure enrichment.
        #
        # Docking scores are not comparable across different proteins: each
        # receptor has its own score scale, so pooling controls from many
        # targets into ONE ROC curve measures between-protein score offsets,
        # not the ability to tell binders from decoys. With a multi-target
        # screening template (one crystal ligand + its decoys per structure)
        # the pooled AUC can read ~0.70 while EVERY structure separates its
        # own active from its own decoys perfectly.
        #
        # The meaningful unit is one structure. Compute enrichment per
        # (pdb_id, ligand_resname) and report the distribution; keep the
        # pooled number too, but label it as diagnostic only.
        per_structure: Dict[tuple, List[Tuple[float, int]]] = {}
        per_structure_details: Dict[tuple, List[Tuple[float, int, str]]] = {}
        structure_targets: Dict[tuple, str] = {}
        per_target: Dict[str, List[Tuple[float, int]]] = {}
        for result in results:
            if result.control_label not in (0, 1):
                continue
            rank_score = self._rank_score_value(result)
            if rank_score is None:
                continue
            key = (result.pdb_id, result.ligand_resname)
            target_name = result.target_name or result.pdb_id
            structure_targets[key] = target_name
            per_target.setdefault(target_name, []).append(
                (rank_score, int(result.control_label))
            )
            per_structure.setdefault(key, []).append(
                (rank_score, int(result.control_label))
            )
            per_structure_details.setdefault(key, []).append(
                (rank_score, int(result.control_label), result.dock_name or result.ligand_resname)
            )

        # A matched-control run has one reference ligand per receptor structure.
        # Assay benchmarks such as LIT-PCBA have many independently labelled
        # actives and inactives at the same structure and require different
        # terminology and diagnostics.
        assay_benchmark = any(
            len({name for _, label, name in details if label == 1}) > 1
            for details in per_structure_details.values()
        )
        summary["enrichment_dataset_type"] = (
            "assay_benchmark" if assay_benchmark else "matched_controls"
        )
        summary["negative_class_label"] = "inactives" if assay_benchmark else "decoys"

        struct_rows = []
        for (pdb_id, resname), sl in sorted(per_structure.items()):
            m = self._compute_enrichment_metrics(sl)
            if not m:
                continue
            # Did the crystal ligand outscore every decoy at this structure?
            actives = [s for s, lab in sl if lab == 1]
            decoys = [s for s, lab in sl if lab == 0]
            top_active = max(actives) if actives else None
            details = per_structure_details[(pdb_id, resname)]
            best_active = max(
                ((score, name) for score, label, name in details if label == 1),
                default=(None, None),
                key=lambda item: item[0] if item[0] is not None else float("-inf")
            )
            best_decoy = max(
                ((score, name) for score, label, name in details if label == 0),
                default=(None, None),
                key=lambda item: item[0] if item[0] is not None else float("-inf")
            )
            active_rank = (
                1 + sum(score > top_active for score in decoys)
                if top_active is not None else None
            )
            struct_rows.append({
                "pdb_id": pdb_id,
                "target_name": structure_targets[(pdb_id, resname)],
                "ligand": resname,
                "roc_auc": m["roc_auc"],
                "actives": m["actives"],
                "decoys": m["decoys"],
                "active_rank": active_rank,
                "best_active": best_active[1],
                "active_score": -top_active if top_active is not None else None,
                "best_decoy": best_decoy[1],
                "best_decoy_score": -best_decoy[0] if best_decoy[0] is not None else None,
                "score_margin": (
                    top_active - best_decoy[0]
                    if top_active is not None and best_decoy[0] is not None else None
                ),
                "active_beats_all_decoys": (
                    bool(top_active is not None and decoys
                         and top_active > max(decoys))
                ),
            })

        if struct_rows:
            aucs = [r["roc_auc"] for r in struct_rows]
            validation_auc_threshold = 0.70
            n_clean = sum(1 for r in struct_rows if r["active_beats_all_decoys"])
            n_auc_pass = sum(
                1 for r in struct_rows if r["roc_auc"] >= validation_auc_threshold
            )
            summary["per_structure_enrichment"] = struct_rows
            summary["n_structures_with_controls"] = len(struct_rows)
            summary["mean_structure_roc_auc"] = float(np.mean(aucs))
            summary["median_structure_roc_auc"] = float(np.median(aucs))
            summary["n_structures_active_top"] = n_clean
            summary["frac_structures_active_top"] = float(n_clean / len(struct_rows))
            summary["validation_auc_threshold"] = validation_auc_threshold
            summary["n_structures_auc_pass"] = n_auc_pass
            summary["failed_control_structures"] = [
                f"{r['pdb_id']}/{r['ligand']}"
                for r in struct_rows if not r["active_beats_all_decoys"]
            ]
            target_rows = []
            for target_name, target_scores in sorted(per_target.items()):
                target_metrics = self._compute_enrichment_metrics(target_scores)
                target_structures = [
                    row for row in struct_rows if row["target_name"] == target_name
                ]
                if not target_metrics or not target_structures:
                    continue
                target_rows.append({
                    "target_name": target_name,
                    "structures": len(target_structures),
                    "actives": target_metrics["actives"],
                    "decoys": target_metrics["decoys"],
                    "macro_roc_auc": float(np.mean([
                        row["roc_auc"] for row in target_structures
                    ])),
                    "pooled_roc_auc": target_metrics["roc_auc"],
                    "structures_active_top": sum(
                        bool(row["active_beats_all_decoys"])
                        for row in target_structures
                    ),
                })
            summary["per_target_enrichment"] = target_rows
            if assay_benchmark:
                # Assay inactives are experimental class labels, not synthetic
                # decoys selected to match one crystal ligand. Applying the
                # pairwise decoy-matching test here produces a false failure.
                summary["control_property_diagnostics"] = []
                summary["control_property_match_passed"] = None
                summary["screening_validation"] = "benchmark_result"
            else:
                property_checks = []
                for key in sorted(per_structure):
                    structure_results = [
                        r for r in results
                        if (r.pdb_id, r.ligand_resname) == key
                        and r.control_label in (0, 1)
                    ]
                    active_results = [
                        r for r in structure_results if r.control_label == 1
                    ]
                    if not active_results:
                        continue
                    active = active_results[0]
                    for decoy in (
                        r for r in structure_results if r.control_label == 0
                    ):
                        matched = self._property_matched(active, decoy)
                        property_checks.append({
                            "pdb_id": key[0],
                            "ligand": key[1],
                            "decoy": decoy.dock_name,
                            "property_matched": matched,
                        })
                known_checks = [
                    c for c in property_checks if c["property_matched"] is not None
                ]
                summary["control_property_diagnostics"] = property_checks
                summary["control_property_match_passed"] = (
                    all(c["property_matched"] for c in known_checks)
                    if known_checks else None
                )
                if summary["control_property_match_passed"] is False:
                    summary["screening_validation"] = "invalid_decoy_matching"
                elif n_clean == len(struct_rows) and len(struct_rows) > 0:
                    summary["screening_validation"] = "passed_strict"
                elif n_auc_pass == len(struct_rows) and len(struct_rows) > 0:
                    summary["screening_validation"] = "passed_enrichment"
                else:
                    summary["screening_validation"] = "needs_review"

        enrichment = self._compute_enrichment_metrics(score_labels)
        if enrichment:
            # Pooled across all targets: only interpretable for a single
            # receptor. See per_structure_enrichment above.
            summary["pooled_enrichment_is_cross_target"] = len(per_structure) > 1
            summary["control_actives"] = enrichment["actives"]
            summary["control_decoys"] = enrichment["decoys"]
            summary["roc_auc"] = enrichment["roc_auc"]
            summary["log_auc"] = enrichment["log_auc"]
            summary["average_precision"] = enrichment["average_precision"]
            summary["active_prevalence"] = enrichment["active_prevalence"]
            quality = (
                "Excellent" if enrichment["roc_auc"] >= 0.7 else
                "Good" if enrichment["roc_auc"] >= 0.6 else
                "Fair" if enrichment["roc_auc"] >= 0.5 else "Poor"
            )
            summary.setdefault("interpretation", {})["enrichment_quality"] = quality
            if assay_benchmark:
                summary["interpretation"]["enrichment_message"] = (
                    f"Assay benchmark: ROC AUC = {enrichment['roc_auc']:.3f} "
                    f"({enrichment['actives']} actives, "
                    f"{enrichment['decoys']} inactives)."
                )
            else:
                summary["interpretation"]["enrichment_message"] = (
                    f"Explicit controls only: ROC AUC = {enrichment['roc_auc']:.3f} "
                    f"({enrichment['actives']} actives, {enrichment['decoys']} decoys)."
                )

            for pct, key in ((1.0, "ef_1_percent"), (5.0, "ef_5_percent"), (10.0, "ef_10_percent")):
                tie_metrics = self._tie_aware_enrichment_factor(score_labels, pct)
                if tie_metrics is None:
                    continue
                summary[key] = tie_metrics["enrichment_factor"]
                summary[f"{key}_min"] = tie_metrics["enrichment_factor_min"]
                summary[f"{key}_max"] = tie_metrics["enrichment_factor_max"]
                summary[f"{key}_details"] = tie_metrics

            if assay_benchmark and len(per_structure) == 1:
                summary["assay_benchmark_charts"] = (
                    self._assay_benchmark_chart_data(score_labels)
                )
        
        # Calculate protocol/engine breakdown (only for actives with valid RMSD)
        # Older redock results predate the explicit completion field. A valid
        # saved output is sufficient to classify those historical cases.
        completed = sum(self._result_is_completed(r) for r in results)
        failed = sum(self._result_is_failed(r) for r in results)
        summary["docking_completed"] = completed
        summary["docking_failed"] = failed

        completed_results = [
            result for result in results
            if self._result_is_completed(result)
        ]
        runtimes = [
            float(result.runtime_sec) for result in completed_results
            if result.runtime_sec is not None
            and np.isfinite(result.runtime_sec)
            and result.runtime_sec > 0
        ]
        summary["mean_runtime_sec"] = float(np.mean(runtimes)) if runtimes else None
        summary["median_runtime_sec"] = float(np.median(runtimes)) if runtimes else None

        rescore_results = [
            result for result in results
            if result.rescore_method
            or result.rescore_score is not None
            or result.rescore_cnn_score is not None
            or result.rescore_cnn_affinity is not None
            or result.rescore_error
        ]
        summary["rescore_count"] = sum(
            result.rescore_score is not None
            or result.rescore_cnn_score is not None
            or result.rescore_cnn_affinity is not None
            for result in rescore_results
        )
        summary["rescore_failures"] = sum(
            bool(result.rescore_error) for result in rescore_results
        )
        summary["rescore_methods"] = sorted({
            str(result.rescore_method) for result in rescore_results
            if result.rescore_method
        })

        screening_results = [
            result for result in results
            if result.mode == "screening" and result.control_label is None
        ]
        screening_scored = []
        for result in screening_results:
            score_details = self._selected_score_details(result)
            if score_details is not None:
                screening_scored.append((result, score_details))
        summary["screening_score_count"] = len(screening_scored)
        summary["screening_unscored_count"] = len(screening_results) - len(screening_scored)
        summary["screening_score_methods"] = sorted({
            details[2] for _, details in screening_scored
        })
        score_directions = {details[3] for _, details in screening_scored}
        summary["screening_score_direction"] = (
            next(iter(score_directions)) if len(score_directions) == 1 else "mixed"
        ) if score_directions else None

        structure_groups: Dict[Tuple[str, str], List[RedockResult]] = {}
        for result in screening_results:
            structure_groups.setdefault(
                (result.pdb_id, result.ligand_resname), []
            ).append(result)

        screening_structure_rows = []
        screening_top_hits = []
        for (pdb_id, ligand), structure_results in sorted(structure_groups.items()):
            structure_scored = []
            for result in structure_results:
                details = self._selected_score_details(result)
                if details is not None:
                    structure_scored.append((result, details))
            structure_scored.sort(key=lambda item: item[1][0], reverse=True)
            structure_completed = sum(
                self._result_is_completed(result) for result in structure_results
            )
            methods = sorted({details[2] for _, details in structure_scored})
            display_values = [details[1] for _, details in structure_scored]
            best_result = structure_scored[0] if structure_scored else None
            target_name = next(
                (result.target_name for result in structure_results if result.target_name),
                pdb_id,
            )
            screening_structure_rows.append({
                "target_name": target_name,
                "pdb_id": pdb_id,
                "ligand": ligand,
                "cases": len(structure_results),
                "completed": structure_completed,
                "completion_rate": (
                    100.0 * structure_completed / len(structure_results)
                    if structure_results else 0.0
                ),
                "scored": len(structure_scored),
                "score_source": ", ".join(methods) if methods else None,
                "best_compound": (
                    best_result[0].dock_name or best_result[0].ligand_resname
                    if best_result else None
                ),
                "best_score": best_result[1][1] if best_result else None,
                "median_score": (
                    float(np.median(display_values)) if display_values else None
                ),
            })
            for rank, (result, details) in enumerate(structure_scored[:5], 1):
                screening_top_hits.append({
                    "target_name": target_name,
                    "pdb_id": pdb_id,
                    "ligand": ligand,
                    "rank": rank,
                    "compound": result.dock_name or result.ligand_resname,
                    "score": details[1],
                    "score_source": details[2],
                })
        summary["per_structure_screening"] = screening_structure_rows
        summary["screening_top_hits"] = screening_top_hits
        summary["screening_failures"] = [
            {
                "pdb_id": result.pdb_id,
                "ligand": result.ligand_resname,
                "compound": result.dock_name or result.ligand_resname,
                "error": result.error_message or "No docked output was produced",
            }
            for result in screening_results
            if self._result_is_failed(result)
        ]

        actives = [
            r for r in results 
            if r.best_rmsd is not None and r.best_rmsd < 900
        ]
        if not actives:
            # Unknown screening compounds have no experimental reference pose.
            # A zero here looks like failed pose recovery, but RMSD is undefined.
            summary["success_rate_best"] = None
        
        by_protocol = {}
        for result in actives:
            proto = result.protocol or "N/A"
            stats = by_protocol.setdefault(
                proto, 
                {"count": 0, "success": 0, "mean_rmsd": []}
            )
            stats["count"] += 1
            if result.best_rmsd < threshold:
                stats["success"] += 1
            stats["mean_rmsd"].append(result.best_rmsd)
        
        for proto, stats in by_protocol.items():
            rmsds = stats.pop("mean_rmsd")
            stats["mean_rmsd"] = float(np.mean(rmsds)) if rmsds else None
        summary["by_protocol"] = by_protocol
        
        by_engine = {}
        for result in actives:
            eng = result.engine or "N/A"
            stats = by_engine.setdefault(
                eng, 
                {"count": 0, "success": 0, "mean_rmsd": []}
            )
            stats["count"] += 1
            if result.best_rmsd < threshold:
                stats["success"] += 1
            stats["mean_rmsd"].append(result.best_rmsd)
        
        for eng, stats in by_engine.items():
            rmsds = stats.pop("mean_rmsd")
            stats["mean_rmsd"] = float(np.mean(rmsds)) if rmsds else None
        summary["by_engine"] = by_engine
        
        return summary
    
    def _summary_to_markdown(self, summary: dict) -> str:
        assay_benchmark = summary.get("enrichment_dataset_type") == "assay_benchmark"
        lines = [
            "# Docking Analysis Summary",
            "",
            f"- Total cases: {summary.get('total_cases')}",
            f"- Docking completed: {summary.get('docking_completed', 0)}",
            f"- Docking failed: {summary.get('docking_failed', 0)}",
            f"- Screening samples: {summary.get('n_samples', 0)}",
        ]
        if summary.get("n_actives") or summary.get("n_decoys"):
            if assay_benchmark:
                lines.extend([
                    f"- Assay actives: {summary.get('n_actives', 0)}",
                    f"- Assay inactives: {summary.get('n_decoys', 0)}",
                ])
            else:
                lines.extend([
                    f"- Control actives: {summary.get('n_actives', 0)}",
                    f"- Control decoys: {summary.get('n_decoys', 0)}",
                ])
        if summary.get("mean_runtime_sec") is not None:
            lines.append(
                f"- Mean recorded case runtime (s): "
                f"{self._fmt(summary.get('mean_runtime_sec'))}"
            )
        if summary.get("median_runtime_sec") is not None:
            lines.append(
                f"- Median recorded case runtime (s): "
                f"{self._fmt(summary.get('median_runtime_sec'))}"
            )

        target_pose_rows = summary.get("per_target_pose_recovery") or []
        has_pose_recovery = (
            summary.get("mean_best_rmsd") is not None or bool(target_pose_rows)
        )
        if has_pose_recovery:
            lines.extend([
                "",
                "## Pose Recovery",
                "",
                f"- RMSD threshold: {summary.get('threshold')}",
                f"- Success rate (best pose): "
                f"{self._fmt(summary.get('success_rate_best'))}%",
                f"- Success rate (Top-1): "
                f"{self._fmt(summary.get('success_rate_top1'))}%",
                f"- Success rate (Top-5): "
                f"{self._fmt(summary.get('success_rate_top5'))}%",
                f"- Success rate (Top-10): "
                f"{self._fmt(summary.get('success_rate_top10'))}%",
                f"- Mean best RMSD: {self._fmt(summary.get('mean_best_rmsd'))}",
                f"- Median best RMSD: {self._fmt(summary.get('median_best_rmsd'))}",
                f"- Mean Top-1 RMSD: {self._fmt(summary.get('mean_top1_rmsd'))}",
                f"- Median Top-1 RMSD: {self._fmt(summary.get('median_top1_rmsd'))}",
                f"- Mean RMSD (best score pose): "
                f"{self._fmt(summary.get('mean_rmsd_best_score'))}",
                f"- Mean near-native fraction: "
                f"{self._fmt(summary.get('mean_near_native_fraction'))}",
                f"- Mean pose count: {self._fmt(summary.get('mean_pose_count'))}",
                f"- Mean score-RMSD Pearson: "
                f"{self._fmt(summary.get('mean_score_rmsd_pearson'))}",
                f"- Mean score-RMSD Spearman: "
                f"{self._fmt(summary.get('mean_score_rmsd_spearman'))}",
            ])
        elif summary.get("n_samples"):
            lines.extend([
                "",
                "RMSD and pose-recovery metrics are omitted because unknown screening "
                "compounds do not have experimental reference poses.",
            ])

        if target_pose_rows:
            lines.extend([
                "",
                "## Pose Recovery (per target)",
                "",
                "| Target | Cases | Best-pose success | Top-1 success | Top-5 success | "
                "Top-10 success | Mean best RMSD | Mean Top-1 RMSD |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
            ])
            for row in target_pose_rows:
                lines.append(
                    f"| {row['target_name']} | {row['cases']} | "
                    f"{self._fmt(row['success_rate_best'])}% | "
                    f"{self._fmt(row['success_rate_top1'])}% | "
                    f"{self._fmt(row['success_rate_top5'])}% | "
                    f"{self._fmt(row['success_rate_top10'])}% | "
                    f"{self._fmt(row['mean_best_rmsd'])} | "
                    f"{self._fmt(row['mean_top1_rmsd'])} |"
                )

        screening_rows = summary.get("per_structure_screening") or []
        if summary.get("n_samples"):
            score_methods = summary.get("screening_score_methods") or []
            direction = summary.get("screening_score_direction")
            direction_text = {
                "lower": "Lower scores rank better.",
                "higher": "Higher scores rank better.",
                "mixed": "Score direction varies by method; ranks are normalized internally.",
            }.get(direction, "")
            lines.extend([
                "",
                "## Screening score summary",
                "",
                f"- Scored samples: {summary.get('screening_score_count', 0)}/"
                f"{summary.get('n_samples', 0)}",
                f"- Samples without a score: {summary.get('screening_unscored_count', 0)}",
            ])
            if score_methods:
                lines.append(f"- Ranking source: {', '.join(score_methods)}")
            if direction_text:
                lines.append(f"- {direction_text}")
            if summary.get("mean_pose_count") is not None:
                lines.append(
                    f"- Mean saved pose count: {self._fmt(summary.get('mean_pose_count'))}"
                )
            if screening_rows:
                lines.extend([
                    "",
                    "Raw scores should be compared within a receptor structure, not pooled "
                    "across different structures.",
                    "",
                    "| Target | PDB | Ligand | Cases | Completed | Completion | Scored | "
                    "Score source | Best compound | Best score | Median score |",
                    "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: |",
                ])
                for row in screening_rows:
                    lines.append(
                        f"| {row['target_name']} | {row['pdb_id']} | {row['ligand']} | "
                        f"{row['cases']} | {row['completed']} | "
                        f"{self._fmt(row['completion_rate'])}% | {row['scored']} | "
                        f"{row.get('score_source') or '-'} | "
                        f"{row.get('best_compound') or '-'} | "
                        f"{self._fmt_score(row.get('best_score'))} | "
                        f"{self._fmt_score(row.get('median_score'))} |"
                    )

        screening_top_hits = summary.get("screening_top_hits") or []
        if screening_top_hits:
            lines.extend([
                "",
                "## Top-ranked compounds per structure",
                "",
                "| Target | PDB | Ligand | Rank | Compound | Score | Score source |",
                "| --- | --- | --- | ---: | --- | ---: | --- |",
            ])
            for row in screening_top_hits:
                lines.append(
                    f"| {row['target_name']} | {row['pdb_id']} | {row['ligand']} | "
                    f"{row['rank']} | {row['compound']} | "
                    f"{self._fmt_score(row['score'])} | {row['score_source']} |"
                )

        screening_failures = summary.get("screening_failures") or []
        if screening_failures:
            lines.extend([
                "",
                "## Failed screening cases",
                "",
                "| PDB | Ligand | Compound | Error |",
                "| --- | --- | --- | --- |",
            ])
            for row in screening_failures:
                error = str(row["error"]).replace("|", "/")
                lines.append(
                    f"| {row['pdb_id']} | {row['ligand']} | "
                    f"{row['compound']} | {error} |"
                )

        if summary.get("rescore_count"):
            methods = summary.get("rescore_methods") or []
            method_text = ", ".join(methods) if methods else "unknown"
            lines.extend([
                "",
                "## Rescoring",
                "",
                f"- Rescored cases: {summary.get('rescore_count')}",
                f"- Rescore failures: {summary.get('rescore_failures')}",
                f"- Methods: {method_text}"
            ])

        target_rows = summary.get("per_target_enrichment") or []
        if target_rows:
            target_heading = (
                "Assay Benchmark Enrichment (per target)"
                if assay_benchmark else "Control Enrichment (per target)"
            )
            negative_label = "Inactives" if assay_benchmark else "Decoys"
            top_label = "Best active top" if assay_benchmark else "Active top"
            lines.extend([
                "",
                f"## {target_heading}",
                "",
                "Macro AUC is the mean of the target's per-structure AUC values. When a "
                "target contains multiple receptor structures, macro AUC is preferred; "
                "target-pooled AUC then combines receptor-specific score scales and is "
                "diagnostic only.",
                "",
                f"| Target | Structures | Macro AUC | Target-pooled AUC | Actives | "
                f"{negative_label} | {top_label} |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ])
            for row in sorted(target_rows, key=lambda item: item["target_name"]):
                lines.append(
                    f"| {row['target_name']} | {row['structures']} | "
                    f"{self._fmt(row['macro_roc_auc'])} | "
                    f"{self._fmt(row['pooled_roc_auc'])} | {row['actives']} | "
                    f"{row['decoys']} | {row['structures_active_top']}/{row['structures']} |"
                )

        if summary.get("n_structures_with_controls"):
            n_s = summary["n_structures_with_controls"]
            n_top = summary.get("n_structures_active_top")
            failed = summary.get("failed_control_structures") or []
            validation_labels = {
                "passed_strict": "PASSED (STRICT)",
                "passed_enrichment": "PASSED (ENRICHMENT)",
                "needs_review": "NEEDS REVIEW",
                "invalid_decoy_matching": "INVALID DECOY MATCHING",
                "benchmark_result": "BENCHMARK RESULT",
            }
            validation = summary.get("screening_validation", "unknown")
            validation_label = validation_labels.get(validation, str(validation).upper())
            structure_heading = (
                "Assay Benchmark Enrichment (per structure)"
                if assay_benchmark else "Control Enrichment (per structure)"
            )
            assessment_label = (
                "Benchmark assessment" if assay_benchmark else "Screening validation"
            )
            top_rank_label = (
                "Best assay active ranks first"
                if assay_benchmark else "Crystal ligand ranks first (strict diagnostic)"
            )
            lines.extend([
                "",
                f"## {structure_heading}",
                "",
                f"- {'Benchmark structures' if assay_benchmark else 'Structures with controls'}: {n_s}",
                f"- {assessment_label}: {validation_label}",
                f"- Structures meeting AUC >= "
                f"{self._fmt(summary.get('validation_auc_threshold'))}: "
                f"{summary.get('n_structures_auc_pass')}/{n_s}",
                f"- {top_rank_label}: {n_top}/{n_s} "
                f"({self._fmt(summary.get('frac_structures_active_top'))})",
                f"- Mean per-structure ROC AUC: "
                f"{self._fmt(summary.get('mean_structure_roc_auc'))}",
                f"- Median per-structure ROC AUC: "
                f"{self._fmt(summary.get('median_structure_roc_auc'))}",
            ])
            if failed:
                if assay_benchmark:
                    lines.extend([
                        "",
                        "Structures where at least one assay inactive outscored the best "
                        "assay active. This rank-one diagnostic complements AUC and early "
                        "enrichment; it is not a decoy-matching validation test:",
                        "",
                    ])
                else:
                    lines.extend([
                        "",
                        "Structures where one or more decoys marginally or substantially "
                        "outscored the crystal ligand. Review AUC, active rank, and score "
                        "margin; this does not by itself invalidate useful enrichment:",
                        "",
                    ])
                lines.extend(f"- {s}" for s in failed[:25])
                if len(failed) > 25:
                    lines.append(f"- ... and {len(failed) - 25} more")

            rows = summary.get("per_structure_enrichment") or []
            if rows:
                if assay_benchmark:
                    lines.extend([
                        "",
                        "| PDB | Site ligand | AUC | Best active | Best-active rank | "
                        "Best-active score | Best inactive | Best-inactive score | "
                        "Margin | Active top |",
                        "| --- | --- | ---: | --- | ---: | ---: | --- | ---: | ---: | :---: |",
                    ])
                else:
                    lines.extend([
                        "",
                        "| PDB | Ligand | AUC | Active rank | Active score | Best decoy | "
                        "Best decoy score | Margin | Active top |",
                        "| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | :---: |",
                    ])
                for r in sorted(rows, key=lambda x: x["roc_auc"]):
                    if assay_benchmark:
                        lines.append(
                            f"| {r['pdb_id']} | {r['ligand']} | "
                            f"{self._fmt(r['roc_auc'])} | {r.get('best_active')} | "
                            f"{r.get('active_rank')} | "
                            f"{self._fmt(r.get('active_score'))} | "
                            f"{r.get('best_decoy')} | "
                            f"{self._fmt(r.get('best_decoy_score'))} | "
                            f"{self._fmt(r.get('score_margin'))} | "
                            f"{'yes' if r['active_beats_all_decoys'] else 'no'} |"
                        )
                    else:
                        lines.append(
                            f"| {r['pdb_id']} | {r['ligand']} | "
                            f"{self._fmt(r['roc_auc'])} | {r.get('active_rank')} | "
                            f"{self._fmt(r.get('active_score'))} | "
                            f"{r.get('best_decoy')} | "
                            f"{self._fmt(r.get('best_decoy_score'))} | "
                            f"{self._fmt(r.get('score_margin'))} | "
                            f"{'yes' if r['active_beats_all_decoys'] else 'no'} |"
                        )

        if summary.get("roc_auc") is not None:
            cross = summary.get("pooled_enrichment_is_cross_target")
            pooled_heading = (
                "Assay Benchmark Enrichment (pooled)"
                if assay_benchmark else "Control Enrichment (pooled)"
            )
            negative_label = "Inactives" if assay_benchmark else "Decoys"
            lines.extend([
                "",
                f"## {pooled_heading}",
                "",
                f"- Actives: {summary.get('control_actives')}",
                f"- {negative_label}: {summary.get('control_decoys')}",
                f"- ROC AUC: {self._fmt(summary.get('roc_auc'))}",
                f"- Average precision: {self._fmt(summary.get('average_precision'))}",
                f"- Random average-precision baseline: "
                f"{self._fmt(summary.get('active_prevalence'))}",
                f"- LogAUC (FPR 0.001-1): {self._fmt(summary.get('log_auc'))}"
            ])
            for label, key in (
                ("EF1%", "ef_1_percent"),
                ("EF5%", "ef_5_percent"),
                ("EF10%", "ef_10_percent"),
            ):
                if summary.get(key) is None:
                    continue
                lines.append(
                    f"- {label} (tie-aware expected; min-max): "
                    f"{self._fmt(summary.get(key))} "
                    f"({self._fmt(summary.get(f'{key}_min'))}-"
                    f"{self._fmt(summary.get(f'{key}_max'))})"
                )
            if assay_benchmark:
                lines.extend([
                    "",
                    "Assay inactives are experimental class labels, not property-matched "
                    "decoys for one crystal ligand. Pairwise decoy-matching validation is "
                    "therefore not applied to this benchmark.",
                ])
            if cross:
                lines.extend([
                    "",
                    "DIAGNOSTIC ONLY: these controls span multiple receptor structures. "
                    "Raw docking scores are not necessarily comparable across receptors, "
                    "so the pooled AUC can include receptor-specific score offsets. "
                    "Use the per-structure table above to judge whether "
                    "screening scores are trustworthy.",
                ])

        if summary.get("charge_count"):
            lines.extend([
                "",
                "## Charge Bias Diagnostics",
                "",
                f"- Charge count: {summary.get('charge_count')}",
                f"- Mean charge: {self._fmt(summary.get('mean_charge'))}",
                f"- Median charge: {self._fmt(summary.get('median_charge'))}",
                f"- Charge fraction (pos/neu/neg): "
                f"{self._fmt(summary.get('charge_frac_positive'))}% / "
                f"{self._fmt(summary.get('charge_frac_neutral'))}% / "
                f"{self._fmt(summary.get('charge_frac_negative'))}%",
                f"- Charge fraction Top-10% scores (pos/neu/neg): "
                f"{self._fmt(summary.get('charge_top10_frac_positive'))}% / "
                f"{self._fmt(summary.get('charge_top10_frac_neutral'))}% / "
                f"{self._fmt(summary.get('charge_top10_frac_negative'))}%",
                f"- Score-charge Pearson: {self._fmt(summary.get('score_charge_pearson'))}",
                f"- Score-charge Spearman: {self._fmt(summary.get('score_charge_spearman'))}"
            ])

        if summary.get("by_protocol"):
            lines.extend([
                "",
                "## By Protocol",
                "",
                "| Protocol | Count | Success | Mean RMSD |",
                "| --- | ---: | ---: | ---: |"
            ])
            for proto, stats in summary.get("by_protocol", {}).items():
                lines.append(
                    f"| {proto} | {stats['count']} | {stats['success']} | "
                    f"{self._fmt(stats.get('mean_rmsd'))} |"
                )
        if summary.get("by_engine"):
            lines.extend([
                "",
                "## By Engine",
                "",
                "| Engine | Count | Success | Mean RMSD |",
                "| --- | ---: | ---: | ---: |"
            ])
            for eng, stats in summary.get("by_engine", {}).items():
                lines.append(
                    f"| {eng} | {stats['count']} | {stats['success']} | "
                    f"{self._fmt(stats.get('mean_rmsd'))} |"
                )

        if summary.get("attempts_by_protocol"):
            lines.extend([
                "",
                "## Protocol Attempts",
                "",
                "| Protocol | Attempts | Success | Mean RMSD |",
                "| --- | ---: | ---: | ---: |"
            ])
            for proto, stats in summary.get("attempts_by_protocol", {}).items():
                lines.append(
                    f"| {proto} | {stats['attempts']} | {stats['success']} | {self._fmt(stats.get('mean_rmsd'))} |"
                )

        return "\n".join(lines) + "\n"

    def _fmt(self, value: Optional[float]) -> str:
        if value is None:
            return "N/A"
        return f"{value:.2f}"

    @staticmethod
    def _fmt_score(value: Optional[float]) -> str:
        return "-" if value is None else f"{value:.3f}"

    def _safe_call(self, func):
        def _wrapper(*args, **kwargs):
            try:
                logger.debug("GUI action: {}", func.__name__)
                return func(*args, **kwargs)
            except Exception as exc:
                logger.error("GUI action failed: {}", exc)
                messagebox.showerror("Action failed", str(exc))
                self._set_status("Action failed")
        return _wrapper

    def _register_busy_widget(self, widget: tk.Widget) -> None:
        if widget in self._busy_widgets:
            return
        self._busy_widgets.append(widget)
        try:
            self._busy_widget_states[widget] = widget.cget("state")
        except Exception:
            self._busy_widget_states[widget] = None

    def _set_busy(self, busy: bool, message: Optional[str] = None) -> None:
        if message:
            self._set_status(message)
        if self._run_button:
            try:
                self._run_button.configure(state="disabled" if busy else "normal")
            except Exception:
                pass
        cursor = "watch" if busy else ""
        try:
            self.configure(cursor=cursor)
            for widget in self._busy_widgets:
                try:
                    widget.configure(cursor=cursor)
                except Exception:
                    pass
        except Exception:
            pass

    def _set_rmsd_variant_available(self, available: bool) -> None:
        self._rmsd_variant_available = available
        if self._variant_all_rmsd_btn is not None:
            try:
                self._variant_all_rmsd_btn.configure(state="normal" if available else "disabled")
            except Exception:
                pass
        if not available and self.variant_mode_var.get() == "all_rmsd":
            self.variant_mode_var.set("all_score")
            self._set_status(
                "RMSD variant selection disabled: no co-crystal Ligand entries in the input sheet."
            )

    def _set_status(self, message: str) -> None:
        self.status_var.set(message)

    def _bring_to_front(self) -> None:
        """Ensure window is visible on macOS when launched from terminal."""
        try:
            self.update_idletasks()
            self.lift()
            self.attributes("-topmost", True)
            self.after(250, lambda: self.attributes("-topmost", False))
            self.focus_force()
        except Exception:
            pass

    def _edit_filters(self) -> None:
        """Edit additive/cofactor filters."""
        dialog = tk.Toplevel(self)
        dialog.title("Edit Filter Lists")
        dialog.geometry("700x500")
        dialog.transient(self)
        dialog.grab_set()

        dialog.grid_columnconfigure(0, weight=1)
        dialog.grid_rowconfigure(1, weight=1)
        dialog.grid_rowconfigure(3, weight=1)

        tk.Label(dialog, text="Known additives (one ID per line):").grid(row=0, column=0, sticky="w", padx=10, pady=(10, 5))
        additives_frame = tk.Frame(dialog)
        additives_frame.grid(row=1, column=0, sticky="nsew", padx=10)
        additives_frame.grid_columnconfigure(0, weight=1)
        additives_frame.grid_rowconfigure(0, weight=1)

        additives_text = tk.Text(additives_frame, height=10, wrap="none")
        additives_text.grid(row=0, column=0, sticky="nsew")
        additives_scroll = ttk.Scrollbar(additives_frame, orient="vertical", command=additives_text.yview)
        additives_scroll.grid(row=0, column=1, sticky="ns")
        additives_text.configure(yscrollcommand=additives_scroll.set)

        tk.Label(dialog, text="Cofactors (one ID per line):").grid(row=2, column=0, sticky="w", padx=10, pady=(10, 5))
        cofactors_frame = tk.Frame(dialog)
        cofactors_frame.grid(row=3, column=0, sticky="nsew", padx=10)
        cofactors_frame.grid_columnconfigure(0, weight=1)
        cofactors_frame.grid_rowconfigure(0, weight=1)

        cofactors_text = tk.Text(cofactors_frame, height=8, wrap="none")
        cofactors_text.grid(row=0, column=0, sticky="nsew")
        cofactors_scroll = ttk.Scrollbar(cofactors_frame, orient="vertical", command=cofactors_text.yview)
        cofactors_scroll.grid(row=0, column=1, sticky="ns")
        cofactors_text.configure(yscrollcommand=cofactors_scroll.set)

        additives_text.insert("1.0", self._format_ids(KNOWN_ADDITIVES))
        cofactors_text.insert("1.0", self._format_ids(COFACTORS))

        button_frame = tk.Frame(dialog)
        button_frame.grid(row=4, column=0, sticky="e", padx=10, pady=10)

        def _save():
            additives = self._parse_ids(additives_text.get("1.0", "end"))
            cofactors = self._parse_ids(cofactors_text.get("1.0", "end"))
            if not additives and not cofactors:
                messagebox.showwarning("Empty lists", "Both lists are empty.")
                return

            global KNOWN_ADDITIVES, COFACTORS, ADDITIVES_ONLY
            KNOWN_ADDITIVES = additives
            COFACTORS = cofactors
            ADDITIVES_ONLY = KNOWN_ADDITIVES - COFACTORS
            self._save_filter_config()
            self._update_pair_count()
            self._set_status("Filter lists updated")
            dialog.destroy()

        def _cancel():
            dialog.destroy()

        tk.Button(button_frame, text="Cancel", command=_cancel, width=10).pack(side="right", padx=(10, 0))
        tk.Button(button_frame, text="Save", command=_save, width=10).pack(side="right")

    def _format_ids(self, values: set) -> str:
        return "\n".join(sorted(values)) + ("\n" if values else "")

    def _parse_ids(self, text: str) -> set:
        ids = set()
        for raw in text.split():
            token = raw.strip().upper()
            if token:
                ids.add(token)
        return ids

    def _load_filter_config(self) -> None:
        if not FILTERS_PATH.exists():
            return
        try:
            data = json.loads(FILTERS_PATH.read_text())
            additives = set(item.strip().upper() for item in data.get("additives", []) if item)
            cofactors = set(item.strip().upper() for item in data.get("cofactors", []) if item)
            if additives or cofactors:
                global KNOWN_ADDITIVES, COFACTORS, ADDITIVES_ONLY
                if additives:
                    KNOWN_ADDITIVES = additives
                if cofactors:
                    COFACTORS = cofactors
                ADDITIVES_ONLY = KNOWN_ADDITIVES - COFACTORS
                logger.info("Loaded filter lists from {}", FILTERS_PATH)
        except Exception as exc:
            logger.warning("Failed to load filter config: {}", exc)

    def _save_filter_config(self) -> None:
        try:
            FILTERS_PATH.parent.mkdir(parents=True, exist_ok=True)
            payload = {
                "additives": sorted(KNOWN_ADDITIVES),
                "cofactors": sorted(COFACTORS)
            }
            FILTERS_PATH.write_text(json.dumps(payload, indent=2))
            logger.info("Saved filter lists to {}", FILTERS_PATH)
        except Exception as exc:
            logger.warning("Failed to save filter config: {}", exc)
    def _load_pairs_from_excel(
        self,
        excel_path: Path,
        exclude_additives: bool = False,
        exclude_cofactors: bool = False,
        use_smiles: bool = False,
        include_controls: bool = False
    ) -> Tuple[List[Dict[str, str]], dict]:
        if not excel_path.exists():
            raise ValueError("Excel file not found")

        df = pd.read_excel(excel_path)
        if df.empty:
            raise ValueError("Excel file is empty")

        col_map = {self._norm_col(c): c for c in df.columns}
        pdb_col = self._find_col(col_map, ["pdb", "pdbid", "pdb_id", "pdbcode", "pdb_code"])
        ligand_col = self._find_col(col_map, ["ligand", "ligandresname", "ligand_resname", "resname", "lig"])
        target_ligand_col = self._find_col(
            col_map,
            ["target_ligand", "targetligand", "target_ligand_name", "targetligandname"]
        )
        target_col = self._find_col(
            col_map, ["protein", "target", "target_name", "targetname", "protein_name"]
        )
        chain_col = self._find_col(col_map, ["chain", "ligandchain", "ligand_chain"])
        label_col = self._find_col(col_map, ["label", "class", "is_active", "isactive", "active", "actives"])
        smiles_col = self._find_col(col_map, ["smiles", "smile", "smiles_string", "smilesstring"])
        decoy_smiles_col = self._find_col(col_map, ["decoy_smiles", "decoysmiles", "decoy_smile", "decoysmile"])
        decoy_compound_col = self._find_col(col_map, ["decoy_compound", "decoycompound", "decoy"])
        center_x_col = self._find_col(col_map, ["pocket_center_x", "site_center_x", "grid_center_x", "center_x"])
        center_y_col = self._find_col(col_map, ["pocket_center_y", "site_center_y", "grid_center_y", "center_y"])
        center_z_col = self._find_col(col_map, ["pocket_center_z", "site_center_z", "grid_center_z", "center_z"])
        site_residues_col = self._find_col(
            col_map,
            ["site_residues", "siteresidues", "pocket_residues", "binding_site_residues", "residues"]
        )

        if pdb_col is None:
            if len(df.columns) >= 1:
                pdb_col = df.columns[0]
            else:
                raise ValueError("Could not detect PDB column in Excel file")
        if ligand_col is None and not use_smiles:
            if len(df.columns) >= 2:
                ligand_col = df.columns[1]
            else:
                raise ValueError("Could not detect ligand column in Excel file")

        invalid_tokens = {"NAN", "NONE", "NA", "N/A", ""}
        pairs = []
        seen = set()
        def _split_list(value: object) -> List[str]:
            if value is None:
                return []
            if isinstance(value, float) and np.isnan(value):
                return []
            text = str(value).strip()
            if not text or text.lower() == "nan":
                return []
            return [item.strip() for item in text.split(",") if item.strip()]

        def _add_pair(
            pdb_id: str,
            ligand: str,
            target_name: Optional[str],
            chain: Optional[str],
            smiles: Optional[str],
            control_label: Optional[int],
            dock_name: Optional[str],
            site_ligand: Optional[str],
            case_id: Optional[str],
            pocket_center: Optional[Tuple[float, float, float]],
            site_residues: Optional[str]
        ) -> None:
            is_control = control_label is not None
            if site_ligand and exclude_additives and site_ligand in ADDITIVES_ONLY and not (include_controls and is_control):
                return
            if site_ligand and exclude_cofactors and site_ligand in COFACTORS and not (include_controls and is_control):
                return
            key = (pdb_id, ligand, smiles, control_label, dock_name or "", site_ligand, pocket_center, site_residues)
            if key in seen:
                return
            seen.add(key)
            pairs.append({
                "pdb_id": pdb_id,
                "ligand": ligand,
                "target_name": target_name,
                "chain": chain,
                "control_label": control_label,
                "smiles": smiles,
                "dock_name": dock_name,
                "site_ligand": site_ligand,
                "case_id": case_id,
                "site_mode": "residues" if site_residues else ("cocrystal" if site_ligand else "prediction"),
                "pocket_center": pocket_center,
                "site_residues": site_residues
            })

        def _parse_center_value(value: object) -> Optional[float]:
            if value is None:
                return None
            if isinstance(value, float) and np.isnan(value):
                return None
            text = str(value).strip()
            if not text or text.lower() == "nan":
                return None
            try:
                return float(text)
            except ValueError:
                return None

        for _, row in df.iterrows():
            pdb_id = str(row[pdb_col]).strip().upper()
            ligand_text = ""
            if ligand_col:
                ligand_value = row.get(ligand_col)
                if not pd.isna(ligand_value):
                    ligand_text = str(ligand_value).strip().upper()
            if not pdb_id or pdb_id in invalid_tokens or len(pdb_id) != 4:
                continue

            site_ligand = None if (not ligand_text or ligand_text in invalid_tokens) else ligand_text
            ligand = site_ligand or "APO"
            target_name = None
            if target_col:
                target_value = row.get(target_col)
                if not pd.isna(target_value):
                    target_text = str(target_value).strip()
                    if target_text and target_text.upper() not in invalid_tokens:
                        target_name = target_text

            chain = None
            if chain_col:
                chain_val = str(row[chain_col]).strip()
                if chain_val and chain_val != "NAN":
                    chain = chain_val

            pocket_center = None
            if center_x_col and center_y_col and center_z_col:
                cx = _parse_center_value(row.get(center_x_col))
                cy = _parse_center_value(row.get(center_y_col))
                cz = _parse_center_value(row.get(center_z_col))
                if cx is not None and cy is not None and cz is not None:
                    pocket_center = (cx, cy, cz)
            site_residues = None
            if site_residues_col:
                residue_value = row.get(site_residues_col)
                if not pd.isna(residue_value):
                    residue_text = str(residue_value).strip()
                    if residue_text and residue_text.lower() != "nan":
                        site_residues = residue_text

            if decoy_smiles_col and label_col is None:
                active_name = ligand
                if target_ligand_col:
                    active_text = str(row[target_ligand_col]).strip()
                    if active_text and active_text.upper() not in invalid_tokens:
                        active_name = active_text
                active_smiles = None
                if use_smiles and smiles_col:
                    smiles_val = row[smiles_col]
                    if not pd.isna(smiles_val):
                        smiles_text = str(smiles_val).strip()
                        if smiles_text and smiles_text.lower() != "nan":
                            active_smiles = smiles_text

                decoy_smiles_list = _split_list(row.get(decoy_smiles_col))
                decoy_name_list = _split_list(row.get(decoy_compound_col)) if decoy_compound_col else []
                has_decoy = bool(decoy_smiles_list)

                if has_decoy:
                    _add_pair(
                        pdb_id=pdb_id,
                        ligand=ligand,
                        target_name=target_name,
                        chain=chain,
                        smiles=active_smiles,
                        control_label=1,
                        dock_name=active_name,
                        site_ligand=site_ligand,
                        case_id=f"{pdb_id}_{ligand}_{active_name}",
                        pocket_center=pocket_center,
                        site_residues=site_residues
                    )
                    for j, decoy_smiles in enumerate(decoy_smiles_list, 1):
                        decoy_name = decoy_name_list[j - 1] if j - 1 < len(decoy_name_list) else f"decoy_{j}"
                        _add_pair(
                            pdb_id=pdb_id,
                            ligand=ligand,
                            target_name=target_name,
                            chain=chain,
                            smiles=decoy_smiles,
                            control_label=0,
                            dock_name=decoy_name,
                            site_ligand=site_ligand,
                            case_id=f"{pdb_id}_{ligand}_{decoy_name}",
                            pocket_center=pocket_center,
                            site_residues=site_residues
                        )
                else:
                    # No decoy provided: treat as regular (non-control) docking case
                    _add_pair(
                        pdb_id=pdb_id,
                        ligand=ligand,
                        target_name=target_name,
                        chain=chain,
                        smiles=active_smiles,
                        control_label=None,
                        dock_name=active_name,
                        site_ligand=site_ligand,
                        case_id=f"{pdb_id}_{ligand}_{active_name}",
                        pocket_center=pocket_center,
                        site_residues=site_residues
                    )
                continue

            smiles = None
            if use_smiles and smiles_col:
                smiles_val = row[smiles_col]
                if not pd.isna(smiles_val):
                    smiles_text = str(smiles_val).strip()
                    if smiles_text and smiles_text.lower() != "nan":
                        smiles = smiles_text

            control_label = None
            if label_col:
                control_label = self._parse_control_label(row[label_col])

            dock_name = ligand
            if target_ligand_col:
                target_value = row.get(target_ligand_col)
                if not pd.isna(target_value):
                    target_text = str(target_value).strip()
                    if target_text and target_text.upper() not in invalid_tokens:
                        dock_name = target_text

            _add_pair(
                pdb_id=pdb_id,
                ligand=ligand,
                target_name=target_name,
                chain=chain,
                smiles=smiles,
                control_label=control_label,
                dock_name=dock_name,
                site_ligand=site_ligand,
                case_id=f"{pdb_id}_{ligand}_{dock_name}" if dock_name != ligand else f"{pdb_id}_{ligand}",
                pocket_center=pocket_center,
                site_residues=site_residues
            )

        return pairs, {
            "pdb_col": pdb_col,
            "ligand_col": ligand_col,
            "chain_col": chain_col,
            "label_col": label_col,
            "smiles_col": smiles_col,
            "decoy_smiles_col": decoy_smiles_col,
            "decoy_compound_col": decoy_compound_col,
            "target_ligand_col": target_ligand_col,
            "target_col": target_col,
            "center_x_col": center_x_col,
            "center_y_col": center_y_col,
            "center_z_col": center_z_col,
            "site_residues_col": site_residues_col
        }

    def _norm_col(self, name: str) -> str:
        return "".join(ch for ch in name.lower() if ch.isalnum())

    def _find_col(self, mapping: Dict[str, str], options: List[str]) -> Optional[str]:
        for opt in options:
            key = self._norm_col(opt)
            if key in mapping:
                return mapping[key]
        return None

    def _variant_config(self) -> Tuple[str, str]:
        """
        Get variant selection mode and selection criterion.
        
        Returns:
            (variant_mode, variant_select_by) tuple
            
        Modes:
            - "adaptive": Smart selection based on molecular flexibility (1-10 variants)
            - "best": Only lowest energy variant (1 variant)
            - "thorough": Comprehensive sampling (10-15 variants)
            - "all": All variants
            - "first": First N variants (legacy)
            
        Selection criteria:
            - "rmsd": Select best by RMSD (for actives)
            - "score": Select best by docking score
        """
        mode = self.variant_mode_var.get()
        
        if mode == "adaptive":
            return "adaptive", "rmsd"
        elif mode == "best":
            return "best", "rmsd"
        elif mode == "thorough":
            return "thorough", "rmsd"
        elif mode == "all_rmsd":
            return "all", "rmsd"
        elif mode == "all_score":
            return "all", "score"
        else:
            # Fallback to adaptive (safer than "first")
            logger.warning(f"Unknown variant mode '{mode}', defaulting to adaptive")
            return "adaptive", "rmsd"

    @staticmethod
    def _variant_selection_for_mode(run_mode: str, requested: str) -> str:
        """Prevent native-pose information from influencing screening ranks."""
        return "score" if run_mode == "screening" else requested

    def _safe_case_id(self, text: str) -> str:
        cleaned = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
        return cleaned or "case"

    def _parse_control_label(self, value: object) -> Optional[int]:
        if value is None:
            return None
        if isinstance(value, float) and np.isnan(value):
            return None
        text = str(value).strip().lower()
        if not text or text == "nan":
            return None
        if text in {"1", "true", "yes", "y", "active", "actives", "positive"}:
            return 1
        if text in {"0", "false", "no", "n", "decoy", "inactive", "negative"}:
            return 0
        try:
            numeric = int(float(text))
        except ValueError:
            return None
        if numeric in (0, 1):
            return numeric
        return None

    def _parse_size_override(self) -> Tuple[Optional[object], bool]:
        raw_vals = [
            self.size_x_var.get().strip(),
            self.size_y_var.get().strip(),
            self.size_z_var.get().strip()
        ]
        vals = [v for v in raw_vals if v]
        if not vals:
            return None, True
        if len(vals) == 1:
            try:
                size = [float(vals[0])] * 3
                logger.info("Box size override using cubic size: {}", size[0])
            except ValueError:
                logger.error("Box size override invalid: {}", raw_vals)
                messagebox.showerror("Box size", "Box size values must be numeric.")
                return None, False
        elif len(vals) == 3:
            try:
                size = [float(v) for v in raw_vals]
            except ValueError:
                logger.error("Box size override invalid: {}", raw_vals)
                messagebox.showerror("Box size", "Box size values must be numeric.")
                return None, False
        else:
            logger.error("Box size override incomplete: {}", raw_vals)
            messagebox.showerror("Box size", "Provide one value for cubic size or all three values.")
            return None, False
        import numpy as np
        return np.array(size, dtype=float), True

    def _parse_int(self, value: str, label: str) -> Optional[int]:
        try:
            return int(value)
        except ValueError:
            logger.error("Invalid integer for {}: '{}'", label, value)
            messagebox.showerror("Input error", f"{label} must be an integer.")
            return None

    def _parse_float(self, value: str, label: str, allow_blank: bool = False) -> Optional[float]:
        if allow_blank and not value.strip():
            return None
        try:
            return float(value)
        except ValueError:
            logger.error("Invalid float for {}: '{}'", label, value)
            messagebox.showerror("Input error", f"{label} must be a number.")
            return None

    def _parse_sample_size(self, silent: bool = False) -> Optional[int]:
        raw = self.sample_size_var.get().strip()
        if not raw:
            if silent:
                return None
            messagebox.showerror("Sample size", "Provide a sample size or disable random sampling.")
            return None
        try:
            value = int(raw)
        except ValueError:
            if not silent:
                messagebox.showerror("Sample size", "Sample size must be an integer.")
            logger.error("Invalid sample size: '{}'", raw)
            return None
        if value <= 0:
            if not silent:
                messagebox.showerror("Sample size", "Sample size must be > 0.")
            logger.error("Invalid sample size (<=0): '{}'", raw)
            return None
        return value

    def _parse_sample_seed(self, silent: bool = False) -> Optional[int]:
        raw = self.sample_seed_var.get().strip()
        if not raw:
            return None
        try:
            return int(raw)
        except ValueError:
            if not silent:
                messagebox.showerror("Sample seed", "Seed must be an integer or blank.")
            logger.error("Invalid sample seed: '{}'", raw)
            return None

    def _apply_random_sample(
        self,
        pairs: List[Dict[str, str]],
        size: int,
        seed: Optional[int],
        include_controls: bool = False
    ) -> List[Dict[str, str]]:
        rng = random.Random(seed)
        if include_controls:
            controls = [p for p in pairs if p.get("control_label") is not None]
            non_controls = [p for p in pairs if p.get("control_label") is None]
            remaining = min(size, len(non_controls))
            if remaining == len(non_controls):
                return pairs

            # Sample receptor-compound cases in rounds across structures. This
            # avoids concentrating a small global sample in one receptor.
            grouped: Dict[tuple, List[Dict[str, str]]] = {}
            for pair in non_controls:
                key = (pair.get("pdb_id"), pair.get("ligand"))
                grouped.setdefault(key, []).append(pair)
            group_keys = list(grouped)
            rng.shuffle(group_keys)
            for group in grouped.values():
                rng.shuffle(group)

            selected = []
            while len(selected) < remaining:
                added = False
                for key in group_keys:
                    group = grouped[key]
                    if group and len(selected) < remaining:
                        selected.append(group.pop())
                        added = True
                if not added:
                    break

            selected_ids = {id(pair) for pair in selected}
            sample = [
                pair for pair in pairs
                if pair.get("control_label") is not None or id(pair) in selected_ids
            ]
            logger.info(
                "Stratified random sample selected: {} total "
                "({} controls, {} non-controls across {} structures, seed={})",
                len(sample),
                len(controls),
                remaining,
                len({(p.get('pdb_id'), p.get('ligand')) for p in selected}),
                seed
            )
            return sample

        if size >= len(pairs):
            return pairs
        sample = rng.sample(pairs, size)
        logger.info("Random sample selected: {} of {} (seed={})", len(sample), len(pairs), seed)
        return sample
