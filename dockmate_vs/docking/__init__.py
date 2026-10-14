"""Built-in engine registrations; add explicit contributor registrations here."""

from .registry import EngineCapabilities, register_engine
from .smina import SminaDockingEngine
from .vina import VinaDockingEngine
from .rdock import RDockEngine

register_engine(
    "smina", SminaDockingEngine,
    EngineCapabilities(binary_names=("smina",), scoring_functions=("vina", "vinardo", "ad4_scoring", "dkoes_fast", "dkoes_scoring"),
                       score_units="kcal/mol", docking_score_label="Smina", color="#177E89"),
    binary_config_keys=("smina_binary",), required_by_doctor=True,
)
register_engine(
    "vina", VinaDockingEngine,
    EngineCapabilities(binary_names=("vina",), scoring_functions=("vina",),
                       score_units="kcal/mol", docking_score_label="Vina", color="#D97732"),
    binary_config_keys=("vina_binary",), required_by_doctor=True,
)
register_engine(
    "rdock", RDockEngine,
    EngineCapabilities(supports_exhaustiveness=False, requires_cavity_preparation=True,
                       site_definition_label="crystal_or_center_seed_v3", binary_names=("rbdock", "rbcavity"),
                       receptor_format="pdb", output_suffix=".sd", docking_score_label="Rdock",
                       color="#5B8E3E"),
    root_config_key="rdock_root", version_args=("-h",),
)
