"""Quadric implicit surface — matches albers::quadric (vec10)."""

from __future__ import annotations

from sympy import Matrix, symbols

from implicit import gradient, hessian, x, y, z

Q = symbols("Q0:10", real=True)


def D_expr():
    return (
        Q[0] * x**2
        + Q[1] * y**2
        + Q[2] * z**2
        + Q[3] * x * y
        + Q[4] * x * z
        + Q[5] * y * z
        + 2 * Q[6] * x
        + 2 * Q[7] * y
        + 2 * Q[8] * z
        + Q[9]
    )


def G_expr() -> Matrix:
    return gradient(D_expr())


def H_expr() -> Matrix:
    return hessian(D_expr())


def sphere_coefficients(R: float) -> list:
    return [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -R**2]


def ellipsoid_coefficients(a: float, b: float, c: float) -> list:
    return [1 / a**2, 1 / b**2, 1 / c**2, 0, 0, 0, 0, 0, 0, -1.0]
