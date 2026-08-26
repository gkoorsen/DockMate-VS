"""
Utility functions for GUI.

Provides helper functions for runtime estimation, format conversion,
and common UI operations.
"""

import os
import urllib.request
import urllib.error
from pathlib import Path
from typing import Tuple, Dict, Optional


def estimate_runtime(
    num_ligands: int,
    flexibility_mode: str,
    num_residues: int,
    exhaustiveness: int,
    cpu_cores: int,
    base_time_per_ligand: float = 120.0
) -> Tuple[float, float, float]:
    """
    Estimate docking campaign runtime.

    Args:
        num_ligands: Number of ligands to dock
        flexibility_mode: 'rigid' or 'flexible'
        num_residues: Number of flexible residues (if flexible mode)
        exhaustiveness: Docking exhaustiveness (1-128)
        cpu_cores: Number of CPU cores
        base_time_per_ligand: Base time per ligand in seconds (default: 120s)

    Returns:
        Tuple of (time_per_ligand, total_sequential, total_parallel) in seconds
    """
    # Adjust for exhaustiveness
    time_factor = exhaustiveness / 8.0

    # Adjust for flexibility
    if flexibility_mode == "rigid":
        flex_multiplier = 1.0
    elif flexibility_mode == "flexible":
        # Based on benchmarks from flexible_sidechain.py
        if num_residues <= 2:
            flex_multiplier = 3.0
        elif num_residues <= 4:
            flex_multiplier = 6.0
        else:
            flex_multiplier = 10.0
    else:
        flex_multiplier = 1.0

    # Calculate times
    time_per_ligand = base_time_per_ligand * time_factor * flex_multiplier
    total_sequential = time_per_ligand * num_ligands
    total_parallel = total_sequential / max(1, cpu_cores)

    return time_per_ligand, total_sequential, total_parallel


def format_time(seconds: float) -> str:
    """
    Format seconds into human-readable time string.

    Args:
        seconds: Time in seconds

    Returns:
        Formatted string (e.g., "2h 15m", "45m", "30s")
    """
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        minutes = (seconds % 3600) / 60
        if minutes < 5:
            return f"{hours:.1f}h"
        else:
            return f"{int(hours)}h {int(minutes)}m"


def get_cpu_count() -> int:
    """
    Get number of available CPU cores.

    Returns:
        Number of CPU cores
    """
    try:
        return os.cpu_count() or 1
    except Exception:
        return 1


def format_configuration_summary(config: Dict) -> str:
    """
    Format configuration dictionary as human-readable summary.

    Args:
        config: Configuration dictionary

    Returns:
        Formatted summary string
    """
    lines = []

    # Project section
    if 'project' in config:
        proj = config['project']
        lines.append(f"Campaign ID: {proj.get('campaign_id', 'N/A')}")
        lines.append(f"Output: {proj.get('output_dir', 'N/A')}")
        lines.append("")

    # Protein section
    if 'protein' in config:
        prot = config['protein']
        lines.append(f"Protein: {os.path.basename(prot.get('pdb_file', 'N/A'))}")
        lines.append(f"  - pH: {prot.get('ph', 'N/A')}")
        water_mode = prot.get('water_handling', 'N/A')
        if water_mode == 'selective':
            lines.append(f"  - Water: Selective retention")
        elif water_mode == 'remove_all':
            lines.append(f"  - Water: Remove all")
        elif water_mode == 'retain_all':
            lines.append(f"  - Water: Retain all")
        lines.append("")

    # Binding site section
    if 'binding_site' in config:
        bs = config['binding_site']
        resname = bs.get('ligand_resname', 'N/A')
        chain = bs.get('ligand_chain', 'N/A')
        lines.append(f"Binding Site: Co-crystal ({resname}, chain {chain})")
        lines.append(f"  - Margin: {bs.get('margin', 'N/A')} Å")
        if 'box_dimensions' in bs:
            dims = bs['box_dimensions']
            lines.append(f"  - Box: {dims.get('size_x', '?')}×{dims.get('size_y', '?')}×{dims.get('size_z', '?')} Å")
        lines.append("")

    # Ligands section
    if 'ligands' in config:
        lig = config['ligands']
        num_ligands = lig.get('num_ligands', 0)
        csv_file = lig.get('csv_file', 'N/A')
        lines.append(f"Ligands: {num_ligands} compounds from {os.path.basename(csv_file)}")
        if 'ph_range' in lig:
            ph_range = lig['ph_range']
            lines.append(f"  - pH range: {ph_range[0]}-{ph_range[1]}")
        lines.append(f"  - Max tautomers: {lig.get('max_tautomers', 'N/A')}")
        lines.append("")

    # Docking section
    if 'docking' in config:
        dock = config['docking']
        flex_mode = dock.get('flexibility_mode', 'rigid')
        lines.append(f"Docking: {flex_mode.replace('_', ' ').title()}")
        lines.append(f"  - Exhaustiveness: {dock.get('exhaustiveness', 'N/A')}")

        if flex_mode == 'flexible' and 'flexible_residues' in dock:
            residues = dock['flexible_residues']
            residue_str = ', '.join(residues[:3])
            if len(residues) > 3:
                residue_str += f", ... ({len(residues)} total)"
            lines.append(f"  - Flexible residues: {len(residues)} ({residue_str})")

        lines.append(f"  - Scoring: {dock.get('scoring', 'vina')}")
        lines.append(f"  - CPU cores: {dock.get('cpu', 'N/A')}")
        lines.append("")

    # Validation section
    if 'validation' in config and config['validation'] is not None:
        val = config['validation']
        enabled_tiers = val.get('enabled_tiers', [])
        lines.append(f"Validation: {len(enabled_tiers)} tier(s) enabled")

        if 1 in enabled_tiers and 'tier1' in val:
            tier1 = val['tier1']
            filters = []
            if tier1.get('pains_filter', False):
                filters.append("PAINS")
            if tier1.get('lipinski_filter', False):
                filters.append("Lipinski")
            lines.append(f"  - Tier 1: {', '.join(filters) if filters else 'None'}")
            lines.append(f"    Max MW: {tier1.get('max_molecular_weight', 500)} Da")
            lines.append(f"    Min H-bonds: {tier1.get('min_hbonds', 1)}")

        if 2 in enabled_tiers:
            lines.append(f"  - Tier 2: ML Rescoring [Not Implemented]")

        if 3 in enabled_tiers:
            lines.append(f"  - Tier 3: Interaction Profiling [Not Implemented]")

        if 4 in enabled_tiers:
            lines.append(f"  - Tier 4: MM-GBSA [Not Implemented]")

        lines.append("")
    else:
        lines.append("Validation: Disabled")
        lines.append("")

    # Runtime estimate section
    if 'runtime_estimate' in config:
        rt = config['runtime_estimate']
        lines.append("Estimated Runtime:")
        lines.append(f"  - Per ligand: ~{format_time(rt.get('time_per_ligand', 0))}")
        lines.append(f"  - Total sequential: {format_time(rt.get('total_sequential', 0))}")
        lines.append(f"  - Parallel ({rt.get('cpu_cores', 1)} cores): ~{format_time(rt.get('total_parallel', 0))}")

    return '\n'.join(lines)


def create_tooltip(widget, text: str):
    """
    Create a tooltip for a widget.

    Args:
        widget: Tkinter widget
        text: Tooltip text
    """
    def on_enter(event):
        tooltip = Tooltip(widget, text)
        widget.tooltip = tooltip

    def on_leave(event):
        if hasattr(widget, 'tooltip'):
            widget.tooltip.destroy()
            del widget.tooltip

    widget.bind('<Enter>', on_enter)
    widget.bind('<Leave>', on_leave)


class Tooltip:
    """Simple tooltip widget."""

    def __init__(self, widget, text):
        """Create tooltip."""
        import tkinter as tk

        self.widget = widget
        self.text = text
        self.tooltip_window = None

        # Create tooltip window
        self.tooltip_window = tk.Toplevel(widget)
        self.tooltip_window.wm_overrideredirect(True)

        # Position near widget
        x = widget.winfo_rootx() + 20
        y = widget.winfo_rooty() + widget.winfo_height() + 5
        self.tooltip_window.wm_geometry(f"+{x}+{y}")

        # Create label
        label = tk.Label(
            self.tooltip_window,
            text=text,
            background="#ffffe0",
            relief='solid',
            borderwidth=1,
            padx=5,
            pady=3
        )
        label.pack()

    def destroy(self):
        """Destroy tooltip."""
        if self.tooltip_window:
            self.tooltip_window.destroy()
            self.tooltip_window = None


def download_pdb_structure(pdb_code: str, output_dir: Optional[str] = None) -> Tuple[bool, str, Optional[str]]:
    """
    Download PDB structure from RCSB PDB.

    Args:
        pdb_code: 4-character PDB code (e.g., '1ABC')
        output_dir: Optional output directory (defaults to temp dir)

    Returns:
        Tuple of (success, message, file_path)
    """
    from loguru import logger
    import tempfile

    # Validate PDB code format
    if not pdb_code or len(pdb_code) != 4:
        return False, "PDB code must be exactly 4 characters", None

    if not pdb_code.isalnum():
        return False, "PDB code must be alphanumeric", None

    # Normalize to uppercase
    pdb_code = pdb_code.upper()

    # Determine output path
    if output_dir:
        output_path = Path(output_dir) / f"{pdb_code}.pdb"
    else:
        temp_dir = tempfile.gettempdir()
        output_path = Path(temp_dir) / f"{pdb_code}.pdb"

    # Download URL
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"

    try:
        logger.info(f"Downloading PDB structure {pdb_code} from RCSB...")

        # Download file
        urllib.request.urlretrieve(url, str(output_path))

        # Verify file was downloaded and is not empty
        if not output_path.exists():
            return False, f"Failed to download {pdb_code}", None

        if output_path.stat().st_size == 0:
            output_path.unlink()  # Delete empty file
            return False, f"Downloaded file is empty (PDB code {pdb_code} may not exist)", None

        # Check if it's a valid PDB file (basic check)
        with open(output_path, 'r') as f:
            first_line = f.readline()
            if not (first_line.startswith('HEADER') or
                    first_line.startswith('ATOM') or
                    first_line.startswith('HETATM') or
                    first_line.startswith('REMARK')):
                output_path.unlink()
                return False, f"Downloaded file is not a valid PDB file", None

        logger.info(f"Successfully downloaded {pdb_code} to {output_path}")
        return True, f"Successfully downloaded {pdb_code}", str(output_path)

    except urllib.error.HTTPError as e:
        if e.code == 404:
            # Recent depositions and large structures are released as mmCIF only,
            # with no legacy PDB file. Convert CIF -> PDB so those entries still work.
            try:
                import gemmi
            except ImportError:
                return False, (
                    f"PDB code {pdb_code} has no legacy PDB file (mmCIF-only entry); "
                    "install gemmi to enable automatic CIF conversion"
                ), None
            try:
                cif_path = output_path.with_suffix(".cif")
                logger.info(f"No legacy PDB for {pdb_code}; retrying as mmCIF...")
                urllib.request.urlretrieve(
                    f"https://files.rcsb.org/download/{pdb_code}.cif", str(cif_path)
                )
                st = gemmi.read_structure(str(cif_path))
                st.setup_entities()
                st.write_pdb(str(output_path))
                logger.info(f"Converted {pdb_code} from mmCIF to {output_path}")
                return True, f"Successfully downloaded {pdb_code} (via mmCIF)", str(output_path)
            except Exception as conv_exc:
                return False, (
                    f"PDB code {pdb_code} not found as PDB and mmCIF conversion "
                    f"failed: {conv_exc}"
                ), None
        else:
            return False, f"HTTP error {e.code}: {e.reason}", None

    except urllib.error.URLError as e:
        return False, f"Network error: {e.reason}", None

    except Exception as e:
        return False, f"Error downloading PDB: {str(e)}", None
