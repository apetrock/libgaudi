"""Shape operator W(G, H) — surface-agnostic primitive."""

from __future__ import annotations

from typing import Tuple

import sympy as sp
from sympy import Matrix


def abstract_gh_symbols() -> Tuple[Matrix, Matrix]:
    g = Matrix(sp.symbols("g0:3", real=True))
    H = Matrix(3, 3, lambda i, j: sp.symbols(f"H{i}{j}", real=True))
    return g, H


def shape_operator(g: Matrix, H: Matrix) -> Matrix:
    """Weingarten map: W = PHP / |g| (regular surface points only)."""
    s = sp.sqrt(g.dot(g))
    n = g / s
    P = sp.eye(3) - n * n.T
    return (P * H * P) / s
