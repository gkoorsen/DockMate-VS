"""
Configuration schema for docking platform using Pydantic models.

This module defines the complete configuration structure with validation.
"""

from enum import Enum
from pathlib import Path
from typing import List, Optional, Union, Literal, Dict, Any

try:
    # Pydantic v2
    from pydantic import BaseModel, Field, field_validator, model_validator
    PYDANTIC_V2 = True
except ImportError:
    # Pydantic v1
    from pydantic import BaseModel, Field, validator, root_validator
    PYDANTIC_V2 = False
else:
    # Provide v1-compatible decorator names when running on Pydantic v2.
    def validator(*fields, **kwargs):
        if "pre" in kwargs:
            kwargs["mode"] = "before" if kwargs.pop("pre") else "after"
        # v2 doesn't support always/each_item in the same way; ignore to keep compatibility.
        kwargs.pop("always", None)
        kwargs.pop("each_item", None)
        return field_validator(*fields, **kwargs)

    def root_validator(*, pre=False, **kwargs):
        kwargs.pop("skip_on_failure", None)
        return model_validator(mode="before" if pre else "after", **kwargs)


# Enums for configuration options

class DockingEngine(str, Enum):
    """Available docking engines."""
    SMINA = "smina"
    VINA = "vina"


class FlexibilityMode(str, Enum):
    """Receptor flexibility modes."""
    RIGID = "rigid"
    FLEXIBLE_SIDECHAIN = "flexible_sidechain"
    ENSEMBLE = "ensemble"
    HYBRID = "hybrid"


class BindingSiteMode(str, Enum):
    """Binding site definition methods."""
    COCRYSTAL = "cocrystal"
    PREDICTION = "prediction"


class WaterHandling(str, Enum):
    """Water molecule handling strategies."""
    REMOVE_ALL = "remove_all"
    RETAIN_ALL = "retain_all"
    SELECTIVE = "selective"


class SelectionStrategy(str, Enum):
    """Ensemble pose selection strategies."""
    BEST_SCORE = "best_score"
    CONSENSUS = "consensus"
    BAYESIAN = "bayesian"


# Configuration models

class WaterRetentionConfig(BaseModel):
    """Configuration for selective water retention."""

    distance_to_protein: float = Field(
        default=3.0,
        ge=0.0,
        le=10.0,
        description="Maximum distance from protein (Angstroms)"
    )
    distance_to_ligand_site: float = Field(
        default=3.0,
        ge=0.0,
        le=10.0,
        description="Maximum distance from ligand binding site (Angstroms)"
    )
    min_contacts: int = Field(
        default=2,
        ge=1,
        le=10,
        description="Minimum contacts with both protein and ligand site"
    )
    max_bfactor: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Maximum B-factor for water retention (optional)"
    )
    min_occupancy: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Minimum occupancy for water retention"
    )


class LigandPreparationConfig(BaseModel):
    """Configuration for ligand preparation."""

    ph_range: tuple[float, float] = Field(
        default=(7.0, 7.4),
        description="pH range for protonation enumeration"
    )
    max_tautomers: int = Field(
        default=8,
        ge=1,
        le=20,
        description="Maximum tautomers per compound"
    )
    max_conformers: int = Field(
        default=10,
        ge=1,
        le=10,
        description="Maximum conformers to generate (typically 1 for docking)"
    )
    energy_window: float = Field(
        default=10.0,
        ge=0.0,
        description="Energy window for conformer generation (kcal/mol)"
    )
    use_etkdg_v3: bool = Field(
        default=True,
        description="Use ETKDGv3 algorithm for conformer generation"
    )
    mmff_minimize: bool = Field(
        default=True,
        description="Perform MMFF94s force field minimization"
    )

    @validator("ph_range")
    def validate_ph_range(cls, v, values):
        """Validate pH range is reasonable."""
        if v[0] >= v[1]:
            raise ValueError("ph_range[0] must be less than ph_range[1]")
        if v[0] < 0 or v[1] > 14:
            raise ValueError("pH values must be between 0 and 14")
        return v


class ProteinPreparationConfig(BaseModel):
    """Configuration for protein preparation."""

    add_missing_atoms: bool = Field(
        default=True,
        description="Add missing heavy atoms using PDBFixer"
    )
    model_missing_loops: bool = Field(
        default=True,
        description="Model missing loops using PDBFixer"
    )
    ph: float = Field(
        default=7.4,
        ge=0.0,
        le=14.0,
        description="pH for protonation state assignment"
    )
    water_handling: WaterHandling = Field(
        default=WaterHandling.REMOVE_ALL,
        description="Strategy for handling crystallographic waters"
    )
    water_retention: Optional[WaterRetentionConfig] = Field(
        default=None,
        description="Parameters for selective water retention (only used if water_handling='selective')"
    )
    remove_hetero_residues: Optional[List[str]] = Field(
        default=None,
        description=(
            "Residue names to strip from the receptor before docking, e.g. the "
            "co-crystal ligand defining the site. Leaving the site ligand in "
            "place blocks its own pocket, so poses are pushed to the periphery "
            "and redocking cannot reproduce the reference geometry."
        )
    )
    optimize_hbonds: bool = Field(
        default=True,
        description="Optimize hydrogen bonds using Reduce"
    )
    validate_structure: bool = Field(
        default=True,
        description="Run MolProbity validation"
    )

    @root_validator
    def validate_water_config(cls, values):
        """Ensure water_retention is provided if selective mode is chosen."""
        if values.get('water_handling') == WaterHandling.SELECTIVE and values.get('water_retention') is None:
            # Provide default selective retention parameters
            values['water_retention'] = WaterRetentionConfig()
        return values


class BindingSiteConfig(BaseModel):
    """Configuration for binding site definition."""

    mode: BindingSiteMode = Field(
        description="Method for defining binding site"
    )

    # Co-crystal mode parameters
    ligand_resname: Optional[str] = Field(
        default=None,
        description="Residue name of co-crystal ligand (e.g., 'LIG', 'ATP')"
    )
    ligand_chain: Optional[str] = Field(
        default=None,
        description="Chain ID of co-crystal ligand (e.g., 'A')"
    )
    margin: float = Field(
        default=12.0,
        ge=5.0,
        le=25.0,
        description="Margin around ligand for box definition (Angstroms)"
    )

    # General parameters
    min_box_size: float = Field(
        default=20.0,
        ge=10.0,
        le=50.0,
        description="Minimum box size in each dimension (Angstroms)"
    )
    max_box_size: float = Field(
        default=50.0,
        ge=20.0,
        le=100.0,
        description="Maximum box size in each dimension (Angstroms)"
    )

    # Validation parameters
    validate_with_selfdocking: bool = Field(
        default=True,
        description="Perform self-docking validation before campaign"
    )
    selfdocking_rmsd_threshold: float = Field(
        default=2.0,
        ge=0.5,
        le=5.0,
        description="RMSD threshold for successful self-docking (Angstroms)"
    )

    @root_validator
    def validate_cocrystal_params(cls, values):
        """Validate co-crystal mode has required parameters."""
        if values.get('mode') == BindingSiteMode.COCRYSTAL:
            if not values.get('ligand_resname'):
                raise ValueError("ligand_resname required for cocrystal mode")
            if not values.get('ligand_chain'):
                raise ValueError("ligand_chain required for cocrystal mode")
        return values


class EnsembleConfig(BaseModel):
    """Configuration for ensemble docking."""

    conformations: List[str] = Field(
        description="Paths to protein conformations (PDB files)"
    )
    selection_strategy: SelectionStrategy = Field(
        default=SelectionStrategy.BEST_SCORE,
        description="Strategy for selecting best pose across ensemble"
    )
    cluster_rmsd_threshold: float = Field(
        default=2.0,
        ge=0.5,
        le=5.0,
        description="RMSD threshold for pose clustering (Angstroms)"
    )

    @validator("conformations")
    def validate_conformations(cls, v, values):
        """Validate conformation files exist."""
        if len(v) < 2:
            raise ValueError("Ensemble mode requires at least 2 conformations")
        if len(v) > 10:
            raise ValueError("Ensemble mode supports maximum 10 conformations")
        return v


class DockingConfig(BaseModel):
    """Configuration for docking execution."""

    engine: DockingEngine = Field(
        default=DockingEngine.SMINA,
        description="Docking engine to use"
    )
    flexibility_mode: FlexibilityMode = Field(
        default=FlexibilityMode.RIGID,
        description="Receptor flexibility mode"
    )

    # Docking parameters
    exhaustiveness: int = Field(
        default=16,
        ge=1,
        le=128,
        description="Search exhaustiveness (higher = more thorough but slower)"
    )
    num_modes: int = Field(
        default=9,
        ge=1,
        le=20,
        description="Number of binding modes to generate"
    )
    cpu: int = Field(
        default=8,
        ge=1,
        le=128,
        description="Number of CPU cores to use"
    )
    seed: int = Field(
        default=42,
        description="Random seed for reproducibility"
    )

    # Flexible side-chain parameters
    flexible_residues: Optional[List[str]] = Field(
        default=None,
        description="Residues to make flexible (e.g., ['A:315', 'A:318'])"
    )

    # Ensemble parameters
    ensemble: Optional[EnsembleConfig] = Field(
        default=None,
        description="Ensemble docking configuration"
    )

    # Scoring function
    scoring: str = Field(
        default="vina",
        description="Scoring function (e.g., 'vina', 'vinardo')"
    )

    @root_validator
    def validate_flexibility_config(cls, values):
        """Validate flexibility mode has required parameters."""
        if values.get('flexibility_mode') == FlexibilityMode.FLEXIBLE_SIDECHAIN:
            if not values.get('flexible_residues') or len(values.get('flexible_residues', [])) == 0:
                raise ValueError("flexible_residues required for flexible_sidechain mode")
            if len(values.get('flexible_residues', [])) > 10:
                raise ValueError("Maximum 10 flexible residues supported")

        if values.get('flexibility_mode') == FlexibilityMode.ENSEMBLE:
            if not values.get('ensemble'):
                raise ValueError("ensemble config required for ensemble mode")

        return values


class Tier1Config(BaseModel):
    """Configuration for Tier 1 validation (PAINS/Lipinski)."""

    pains_filter: bool = Field(default=True, description="Apply PAINS filter")
    lipinski_filter: bool = Field(default=True, description="Apply Lipinski Rule of Five")
    min_hbonds: int = Field(default=1, ge=0, description="Minimum H-bonds required")
    max_molecular_weight: float = Field(default=500.0, ge=0.0, description="Maximum molecular weight (Da)")


class Tier2Config(BaseModel):
    """Configuration for Tier 2 validation (ML rescoring)."""

    rescoring_method: str = Field(default="rfscore", description="Rescoring method")
    top_n: int = Field(default=1000, ge=1, description="Number of top compounds to keep")
    use_ensemble: bool = Field(default=False, description="Use ensemble rescoring")


class Tier3Config(BaseModel):
    """Configuration for Tier 3 validation (Interaction profiling)."""

    interaction_profiling: bool = Field(default=True, description="Generate interaction fingerprints")
    min_hbonds: int = Field(default=2, ge=0, description="Minimum H-bonds")
    min_hydrophobic: int = Field(default=3, ge=0, description="Minimum hydrophobic contacts")
    require_pi_stacking: bool = Field(default=False, description="Require π-stacking interaction")
    top_n: int = Field(default=200, ge=1, description="Number of top compounds to keep")


class Tier4Config(BaseModel):
    """Configuration for Tier 4 validation (MM-GBSA)."""

    method: str = Field(default="mmgbsa", description="Free energy method")
    run_md: bool = Field(default=False, description="Run short MD before MM-GBSA")
    md_time_ns: float = Field(default=2.0, ge=0.1, le=10.0, description="MD simulation time (ns)")
    epsilon_in: float = Field(default=2.0, ge=1.0, le=4.0, description="Interior dielectric constant")
    gb_model: str = Field(default="GBOBC", description="GB model (GBOBC or GBGBn1)")
    decompose_energy: bool = Field(default=False, description="Perform per-residue decomposition")


class ValidationConfig(BaseModel):
    """Configuration for validation pipeline."""

    enabled_tiers: List[int] = Field(
        default=[1, 2, 3],
        description="Validation tiers to enable (1-4)"
    )
    tier1: Tier1Config = Field(default_factory=Tier1Config)
    tier2: Tier2Config = Field(default_factory=Tier2Config)
    tier3: Tier3Config = Field(default_factory=Tier3Config)
    tier4: Tier4Config = Field(default_factory=Tier4Config)

    @validator("enabled_tiers")
    def validate_tiers(cls, v, values):
        """Validate tier numbers."""
        if not all(1 <= t <= 4 for t in v):
            raise ValueError("Tier numbers must be between 1 and 4")
        if len(v) != len(set(v)):
            raise ValueError("Duplicate tier numbers not allowed")
        return sorted(v)


class OutputConfig(BaseModel):
    """Configuration for output and export."""

    export_formats: List[str] = Field(
        default=["csv", "sdf", "json"],
        description="Output formats for results"
    )
    save_all_poses: bool = Field(
        default=False,
        description="Save all poses (not just top pose per ligand)"
    )
    generate_visualizations: bool = Field(
        default=True,
        description="Generate plots and visualizations"
    )
    compress_results: bool = Field(
        default=True,
        description="Compress large result files"
    )


class LoggingConfig(BaseModel):
    """Configuration for logging."""

    level: str = Field(default="INFO", description="Logging level")
    file: str = Field(default="campaign.log", description="Log file path")
    console: bool = Field(default=True, description="Enable console logging")
    structured: bool = Field(default=True, description="Enable structured JSON logging")
    rotation: str = Field(default="100 MB", description="Log rotation size")


class CheckpointConfig(BaseModel):
    """Configuration for checkpoint/resume."""

    enabled: bool = Field(default=True, description="Enable checkpointing")
    interval: int = Field(default=100, ge=1, description="Checkpoint every N ligands")
    compression: bool = Field(default=True, description="Compress checkpoint files")


class PhaseConfig(BaseModel):
    """Configuration for a single phase in hybrid mode."""

    name: str = Field(description="Phase name")
    flexibility_mode: FlexibilityMode = Field(description="Flexibility mode for this phase")
    exhaustiveness: int = Field(default=8, ge=1, le=128)
    select_top_n: int = Field(ge=1, description="Number of top compounds to pass to next phase")
    flexible_residues: Optional[List[str]] = None
    ensemble: Optional[EnsembleConfig] = None


class HybridConfig(BaseModel):
    """Configuration for hybrid multi-phase docking."""

    phases: List[PhaseConfig] = Field(
        description="List of phases to execute sequentially"
    )

    @validator("phases")
    def validate_phases(cls, v, values):
        """Validate phase configuration."""
        if len(v) < 2:
            raise ValueError("Hybrid mode requires at least 2 phases")
        if len(v) > 5:
            raise ValueError("Maximum 5 phases supported")
        return v


class ReceptorConfig(BaseModel):
    """Configuration for receptor (protein)."""

    pdb_file: Union[str, Path] = Field(description="Path to PDB file")
    preparation: ProteinPreparationConfig = Field(default_factory=ProteinPreparationConfig)

    @validator("pdb_file")
    def validate_pdb_file(cls, v, values):
        """Validate PDB file exists."""
        path = Path(v)
        if not path.exists():
            raise ValueError(f"PDB file not found: {v}")
        if path.suffix not in [".pdb", ".gz"]:
            raise ValueError(f"Invalid PDB file extension: {path.suffix}")
        return str(path)


class LigandsConfig(BaseModel):
    """Configuration for ligands."""

    input_file: Union[str, Path] = Field(description="Path to ligand file")
    format: str = Field(default="smiles", description="Input format (smiles, sdf, mol2)")
    preparation: LigandPreparationConfig = Field(default_factory=LigandPreparationConfig)

    @validator("input_file")
    def validate_input_file(cls, v, values):
        """Validate input file exists."""
        path = Path(v)
        if not path.exists():
            raise ValueError(f"Ligand file not found: {v}")
        return str(path)


class CampaignConfig(BaseModel):
    """Complete configuration for a docking campaign."""

    # Basic settings
    output_dir: str = Field(description="Output directory for campaign")
    name: Optional[str] = Field(default=None, description="Campaign name")
    description: Optional[str] = Field(default=None, description="Campaign description")

    # Main configuration sections
    receptor: ReceptorConfig
    ligands: LigandsConfig
    binding_site: BindingSiteConfig
    docking: DockingConfig
    validation: ValidationConfig = Field(default_factory=ValidationConfig)

    # Additional settings
    output: OutputConfig = Field(default_factory=OutputConfig)
    logging: LoggingConfig = Field(default_factory=LoggingConfig)
    checkpoint: CheckpointConfig = Field(default_factory=CheckpointConfig)

    # Hybrid mode (optional)
    hybrid: Optional[HybridConfig] = Field(default=None)

    @root_validator
    def validate_config(cls, values):
        """Validate overall configuration consistency."""
        # Create output directory if it doesn't exist
        Path(values['output_dir']).mkdir(parents=True, exist_ok=True)

        # Set default name if not provided
        if not values.get('name'):
            values['name'] = Path(values['output_dir']).name

        # Validate hybrid mode
        if values.get('docking', {}).get('flexibility_mode') == FlexibilityMode.HYBRID:
            if not values.get('hybrid'):
                raise ValueError("hybrid config required for hybrid flexibility mode")

        return values

    class Config:
        """Pydantic configuration."""
        json_schema_extra = {
            "example": {
                "output_dir": "./my_campaign",
                "receptor": {
                    "pdb_file": "protein.pdb",
                    "preparation": {
                        "ph": 7.4,
                        "water_handling": "selective",
                        "water_retention": {
                            "distance_to_protein": 3.0,
                            "distance_to_ligand_site": 3.0,
                            "min_contacts": 2
                        }
                    }
                },
                "ligands": {
                    "input_file": "ligands.smi",
                    "format": "smiles"
                },
                "binding_site": {
                    "mode": "cocrystal",
                    "ligand_resname": "LIG",
                    "ligand_chain": "A"
                },
                "docking": {
                    "engine": "smina",
                    "flexibility_mode": "rigid",
                    "exhaustiveness": 16
                },
                "validation": {
                    "enabled_tiers": [1, 2, 3]
                }
            }
        }
