"""Darboux cyclide — matches albers::eval_darboux (vec14)."""

from __future__ import annotations

from sympy import Matrix, symbols

from implicit import gradient, hessian, x, y, z

C = symbols("C0:14", real=True)
A, B, Cc, Dc, E, F, G, H, I, J, lam, mu, nu, kap = C


def D_expr():
    X = x**2 + y**2 + z**2
    L = mu * x + nu * y + kap * z
    Q_quad = (
        A * x**2
        + B * y**2
        + Cc * z**2
        + 2 * Dc * x * y
        + 2 * E * x * z
        + 2 * F * y * z
        + 2 * G * x
        + 2 * H * y
        + 2 * I * z
        + J
    )
    return lam * X**2 + L * X + Q_quad


def G_expr() -> Matrix:
    return gradient(D_expr())


def H_expr() -> Matrix:
    return hessian(D_expr())


def sphere_coefficients(R: float) -> list:
    return [1, 1, 1, 0, 0, 0, 0, 0, 0, -R**2, 0, 0, 0, 0]


def darboux_row_at(pos) -> list:
    """D(x) linear coeffs w.r.t. vec14 Q at position pos."""
    D = D_expr().subs({x: pos[0], y: pos[1], z: pos[2]})
    return [float(D.coeff(C[i])) for i in range(14)]


def fit_darboux_homogeneous(samples) -> list:
    """Fit Q from (position, normal) samples — homogeneous value rows (SVD)."""
    import numpy as np

    rows = [darboux_row_at(p) for p, _ in samples]
    A = np.asarray(rows, dtype=float)
    _, _, vh = np.linalg.svd(A)
    Q = vh[-1]
    if Q[0] < 0:
        Q = -Q
    return Q.tolist()


def torus_point(u: float, v: float, major_R: float, minor_r: float):
    import numpy as np

    rho = major_R + minor_r * np.cos(v)
    return np.array([rho * np.cos(u), rho * np.sin(u), minor_r * np.sin(v)])


def torus_normal(p, major_R: float):
    import numpy as np

    radial = np.array([p[0], p[1], 0.0])
    if np.linalg.norm(radial) < 1e-12:
        radial = np.array([1.0, 0.0, 0.0])
    radial /= np.linalg.norm(radial)
    tube = major_R * radial
    n = p - tube
    return n / np.linalg.norm(n)


def fit_canonical_torus(major_R: float = 1.25, minor_r: float = 0.35) -> list:
    import numpy as np

    samples = []
    for i in range(48):
        u = 2.0 * np.pi * i / 48.0
        for j in range(24):
            v = 2.0 * np.pi * j / 24.0
            p = torus_point(u, v, major_R, minor_r)
            n = torus_normal(p, major_R)
            samples.append((p, n))
    return fit_darboux_homogeneous(samples)


def canonical_torus_coefficients() -> list:
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parent / "test_data" / "canonical_torus_Q.json"
    if path.exists():
        data = json.loads(path.read_text())
        return data["Q"]
    Q = fit_canonical_torus()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "major_R": 1.25,
                "minor_r": 0.35,
                "Q": Q,
            },
            indent=2,
        )
    )
    return Q
