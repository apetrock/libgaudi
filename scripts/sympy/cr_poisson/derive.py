"""Symbolic checks for MATH.md (one triangle / one diamond).

1. The CR parallel derivative of the coordinate functions is the unit edge
   vector, so `N_f . g_par = 0` and the A.2 `v_par` term vanishes.
2. `K_f` (4 x 6): the read-out `s_f = K_f n_f` on one face, symbolic.
3. The closed-form closest point `n* = n + K^+ (s - K n)` satisfies the
   constraint and is the least-norm correction.

Run: python3 derive.py
"""

from __future__ import annotations

import sympy as sp


def parallel_derivative_is_edge_direction():
    """Both faces of a diamond give the same parallel derivative: the edge."""
    a = sp.Matrix(sp.symbols("a0:3", real=True))
    b = sp.Matrix(sp.symbols("b0:3", real=True))
    c = sp.Matrix(sp.symbols("c0:3", real=True))
    d = sp.Matrix(sp.symbols("d0:3", real=True))
    li = sp.sqrt((b - a).dot(b - a))
    A1, A2 = sp.symbols("A1 A2", positive=True)
    # Face (a, b, c) holds a -> b (s = +1); face (b, a, d) holds b -> a
    # (s = -1, local tail b, head a). Stein A.2 parallel term per face.
    par = (A1 / (3 * li)) * (b - a) + (-1) * (A2 / (3 * li)) * (a - b)
    g = sp.simplify(par / ((A1 + A2) / 3))
    same = sp.simplify(g - (b - a) / li) == sp.zeros(3, 1)
    n1 = (b - a).cross(c - a)
    normal_zero = sp.simplify(n1.dot(g)) == 0
    return same, normal_zero


def read_out_matrix():
    """`K_f` for one face in its frame: `T` along edge 0, `B = N x T`.
    Edge directions at angles `theta_i` from `T`; `v_perp = N x v`."""
    th = sp.symbols("theta0:3", real=True)
    s = sp.symbols("s0:3")
    l = sp.symbols("l0:3", positive=True)
    A = sp.symbols("A", positive=True)
    cols = []
    for i in range(3):
        v = sp.Matrix([sp.cos(th[i]), sp.sin(th[i])])
        vp = sp.Matrix([-sp.sin(th[i]), sp.cos(th[i])])
        k = -s[i] * l[i] / A
        for M in (vp * v.T, vp * vp.T):
            cols.append(k * sp.Matrix([M[0, 0], M[0, 1], M[1, 0], M[1, 1]]))
    return sp.Matrix.hstack(*cols)


def closest_point_check():
    """Full-row-rank K (3 x 6): n* satisfies K n* = s and the correction is
    orthogonal to the null space of K (least norm)."""
    K = sp.Matrix(3, 6, lambda i, j: sp.Rational((i + 1) * (j + 2) % 7 + 1, 3))
    n = sp.Matrix(sp.symbols("n0:6"))
    s = sp.Matrix(sp.symbols("t0:3"))
    Kp = K.T * (K * K.T).inv()
    ns = n + Kp * (s - K * n)
    ok_constraint = sp.simplify(K * ns - s) == sp.zeros(3, 1)
    ok_least = all(sp.simplify((ns - n).dot(z)) == 0 for z in K.nullspace())
    return ok_constraint, ok_least


if __name__ == "__main__":
    same, normal_zero = parallel_derivative_is_edge_direction()
    print("g_par is the unit edge vector:", same, "| N_f . g_par = 0:",
          normal_zero)
    print("K_f (4 x 6), columns (par, perp) per edge:")
    sp.pprint(sp.simplify(read_out_matrix()))
    print("closest point: constraint, least norm:", closest_point_check())
