"""Pytest session fixtures for sympy medial pipeline."""

from __future__ import annotations

import pytest


@pytest.fixture(scope="session", autouse=True)
def _warm_energy_cache():
    from compose import compose_energy_abstract, get_graph

    compose_energy_abstract()
    get_graph("quadric")
