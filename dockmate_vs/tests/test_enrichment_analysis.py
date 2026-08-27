import numpy as np
import pandas as pd
import pytest

import dockmate_vs.analysis.enrichment_analysis as enrichment_module
from dockmate_vs.analysis.enrichment_analysis import EnrichmentAnalysis


def _perfect_ranking():
    return EnrichmentAnalysis(
        pd.DataFrame(
            {
                "compound_id": ["active-1", "active-2", "decoy-1", "decoy-2"],
                "is_active": [True, True, False, False],
                "docking_score": [-10.0, -9.0, -5.0, -4.0],
            }
        )
    )


def _manual_trapezoid(y, x):
    y = np.asarray(y)
    x = np.asarray(x)
    return np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2.0)


@pytest.mark.parametrize("api_name", ["trapezoid", "trapz"])
def test_calculate_pauc_supports_current_and_legacy_numpy_apis(monkeypatch, api_name):
    calls = []

    def integration_api(y, x):
        calls.append(True)
        return _manual_trapezoid(y, x)

    unavailable_api = "trapz" if api_name == "trapezoid" else "trapezoid"
    monkeypatch.setattr(enrichment_module.np, api_name, integration_api, raising=False)
    monkeypatch.delattr(enrichment_module.np, unavailable_api, raising=False)

    assert _perfect_ranking().calculate_pauc(fraction=0.5) == pytest.approx(1.0)
    assert calls
