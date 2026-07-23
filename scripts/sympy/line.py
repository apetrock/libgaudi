"""Line search along a fixed ray: f(t) = p0 + t * n0."""

from __future__ import annotations

from typing import Dict

import sympy as sp
from sympy import Matrix, symbols

t = sp.Symbol("t", real=True)
p0 = Matrix(symbols("p00:3", real=True))
n0 = Matrix(symbols("n00:3", real=True))


def f_expr() -> Matrix:
    return p0 + t * n0


def fp_expr() -> Matrix:
    return n0


def fpp_expr() -> Matrix:
    return Matrix(3, 1, [0, 0, 0])


def line_symbols() -> Dict[str, sp.Symbol]:
    out: Dict[str, sp.Symbol] = {str(s): s for s in list(p0) + list(n0)}
    out[str(t)] = t
    return out
