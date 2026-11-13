"""Gap of a ring point against the posed osculating Dupin cyclide, and its
gradient in the seven fit parameters, for `osculating_dupin_ring.hpp`.

The model mirrors `land_on_posed_dupin` exactly. In the foot's chart
(T, B, N) = (x, y, z), a ring point `d` is measured from the posed foot:

    n = normalize(N + tilt_t T + tilt_b B)
    t = normalize(T - (T.n) n),   b = n x t
    e_hi = hx t + hy b,   e_lo = lx t + ly b
    P = (p.e_hi, p.e_lo, p.n),   p = d - offset n

then pulled back by the special conformal map, dropped on the torus graph,
and pushed forward again:

    q = S_{-beta}(P),   beta = (b_hi, b_lo, 0)
    h(u, v) = sag(k_hi, u) + sag(k_lo / (1 - k_lo sag(k_hi, u)), v)
    land = S_beta(q_u, q_v, h),   normal = normalize(dS_beta (-h_u, -h_v, 1))
    gap = (P - land) . normal

with S_b(x) = (x + b |x|^2) / (1 + 2 b.x + |b|^2 |x|^2) and
sag(k, t) = k t^2 / (1 + sqrt(1 - k^2 t^2)).

Differentiating that composition in one piece does not finish, so it is
split into four stages, each emitted with its own Jacobian, and the caller
chains them (`osculating_dupin_ring.hpp`, `dupin_gap_and_gradient`):

    pose:      (offset, tilt_t, tilt_b)          -> P            3x3
    pullback:  (P, b_hi, b_lo)                   -> (q_u, q_v)   2x5
    torus:     (q_u, q_v, k_hi, k_lo)            -> (h, h_u, h_v) 3x4
    gap:       (P, b_hi, b_lo, q_u, q_v, h, h_u, h_v) -> gap     1x10

Parameters of the fit, in order: offset, tilt_t, tilt_b, k_hi, k_lo, b_hi, b_lo.

Run: python3 osculating_dupin_jacobian.py          # verify
     python3 osculating_dupin_jacobian.py --emit   # also write the C++ header
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.printing.cxx import CXX17CodePrinter

HEADER = Path(__file__).resolve().parent / "generated" / "osculating_dupin_gap_generated.hpp"


def sag(k, t):
    return k * t**2 / (1 + sp.sqrt(1 - k**2 * t**2))


def special_conformal(b, x):
    x2 = x.dot(x)
    den = 1 + 2 * b.dot(x) + b.dot(b) * x2
    return (x + b * x2) / den


def stages():
    """Each stage: (inputs, outputs). Jacobians are taken per stage."""
    # Pose: P from the chart offset of the ring point.
    d = sp.symbols("d0 d1 d2", real=True)
    hx, hy, lx, ly = sp.symbols("hx hy lx ly", real=True)
    offset, tilt_t, tilt_b = sp.symbols("offset tilt_t tilt_b", real=True)
    T = sp.Matrix([1, 0, 0])
    n_raw = sp.Matrix([tilt_t, tilt_b, 1])
    n = n_raw / sp.sqrt(n_raw.dot(n_raw))
    t_raw = T - T.dot(n) * n
    t = t_raw / sp.sqrt(t_raw.dot(t_raw))
    bt = n.cross(t)
    e_hi = hx * t + hy * bt
    e_lo = lx * t + ly * bt
    p = sp.Matrix(d) - offset * n
    pose = {
        "consts": list(d) + [hx, hy, lx, ly],
        "inputs": [offset, tilt_t, tilt_b],
        "outputs": [p.dot(e_hi), p.dot(e_lo), p.dot(n)],
    }

    # Pullback: q = S_{-beta}(P), tangential part only.
    P = sp.symbols("P0 P1 P2", real=True)
    b_hi, b_lo = sp.symbols("b_hi b_lo", real=True)
    beta = sp.Matrix([b_hi, b_lo, 0])
    q = special_conformal(-beta, sp.Matrix(P))
    pullback = {"consts": [], "inputs": list(P) + [b_hi, b_lo], "outputs": [q[0], q[1]]}

    # Torus graph: height and gradient.
    u, v, k_hi, k_lo = sp.symbols("u v k_hi k_lo", real=True)
    ht = sag(k_hi, u)
    h = ht + sag(k_lo / (1 - k_lo * ht), v)
    torus = {
        "consts": [],
        "inputs": [u, v, k_hi, k_lo],
        "outputs": [h, sp.diff(h, u), sp.diff(h, v)],
    }

    # Gap: push the torus point forward and measure along the cyclide normal.
    hh, hu, hv = sp.symbols("h h_u h_v", real=True)
    qs = sp.Matrix([u, v, hh])
    y = sp.Matrix(sp.symbols("y0 y1 y2", real=True))
    Sy = special_conformal(beta, y)
    at_qs = dict(zip(y, qs))
    land = Sy.xreplace(at_qs)
    n_c_raw = Sy.jacobian(y).xreplace(at_qs) * sp.Matrix([-hu, -hv, 1])
    n_c = n_c_raw / sp.sqrt(n_c_raw.dot(n_c_raw))
    gap = (sp.Matrix(P) - land).dot(n_c)
    gap_stage = {
        "consts": [],
        "inputs": list(P) + [b_hi, b_lo, u, v, hh, hu, hv],
        "outputs": [gap],
    }
    return {"pose": pose, "pullback": pullback, "torus": torus, "gap": gap_stage}


def jacobians(st):
    for s in st.values():
        s["jac"] = [[sp.diff(o, i) for i in s["inputs"]] for o in s["outputs"]]
    return st


def lambdas(st):
    out = {}
    for name, s in st.items():
        args = s["consts"] + s["inputs"]
        out[name] = (
            sp.lambdify(args, s["outputs"], "numpy"),
            sp.lambdify(args, s["jac"], "numpy"),
        )
    return out


def evaluate(L, d, axes, x, with_grad=True):
    """Chain the stages numerically, as the C++ does."""
    offset, tilt_t, tilt_b, k_hi, k_lo, b_hi, b_lo = x
    P = np.array(L["pose"][0](*d, *axes, offset, tilt_t, tilt_b), dtype=float)
    q = np.array(L["pullback"][0](*P, b_hi, b_lo), dtype=float)
    th = np.array(L["torus"][0](q[0], q[1], k_hi, k_lo), dtype=float)
    g = float(L["gap"][0](*P, b_hi, b_lo, q[0], q[1], *th)[0])
    if not with_grad:
        return g
    dP = np.zeros((3, 7))
    dP[:, 0:3] = np.array(L["pose"][1](*d, *axes, offset, tilt_t, tilt_b), dtype=float)
    Jq = np.array(L["pullback"][1](*P, b_hi, b_lo), dtype=float)
    dq = Jq[:, 0:3] @ dP
    dq[:, 5:7] += Jq[:, 3:5]
    Jt = np.array(L["torus"][1](q[0], q[1], k_hi, k_lo), dtype=float)
    dt = Jt[:, 0:2] @ dq
    dt[:, 3:5] += Jt[:, 2:4]
    Jg = np.array(L["gap"][1](*P, b_hi, b_lo, q[0], q[1], *th), dtype=float)[0]
    grad = Jg[0:3] @ dP + Jg[5:7] @ dq + Jg[7:10] @ dt
    grad[5:7] += Jg[3:5]
    return g, grad


def verify(L):
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(16):
        a = rng.uniform(-np.pi, np.pi)
        axes = [np.cos(a), np.sin(a), -np.sin(a), np.cos(a)]
        d = list(rng.uniform(-0.15, 0.15, 3))
        x = np.array([0.01, 0.04, -0.03, 2.5, -0.6, 0.12, -0.18]) + rng.uniform(-0.01, 0.01, 7)
        _, an = evaluate(L, d, axes, x)
        fd = np.zeros(7)
        for k in range(7):
            e = 1e-6
            xp, xm = x.copy(), x.copy()
            xp[k] += e
            xm[k] -= e
            fd[k] = (evaluate(L, d, axes, xp, False) - evaluate(L, d, axes, xm, False)) / (2 * e)
        worst = max(worst, float(np.max(np.abs(an - fd)) / max(1.0, np.max(np.abs(fd)))))
    print(f"chained gradient vs central differences, worst relative error: {worst:.2e}")
    assert worst < 1e-6


DOCS = {
    "pose": "Posed point `P = (p.e_hi, p.e_lo, p.n)` from the chart offset `d` and\n"
    "/// the principal coefficients, and `dP / d(offset, tilt_t, tilt_b)`.",
    "pullback": "Tangential part of `S_{-b}(P)`, and its Jacobian in `(P, b_hi, b_lo)`.",
    "torus": "Torus graph height and gradient `(h, h_u, h_v)` at `(u, v)`, and their\n"
    "/// Jacobian in `(u, v, k_hi, k_lo)`.",
    "gap": "Gap of `P` against the torus point `(u, v, h)` pushed forward by `S_b`,\n"
    "/// along the cyclide normal, and its gradient in\n"
    "/// `(P, b_hi, b_lo, u, v, h, h_u, h_v)`.",
}


class FastPowPrinter(CXX17CodePrinter):
    """C++ printer that writes small integer and half-integer powers as
    products and square roots. `std::pow` is a library call; after `cse` the
    bases are mostly temporaries, so products are a few multiplies."""

    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        b = self.parenthesize(base, 1000)
        if exp.is_Integer and 1 <= abs(int(exp)) <= 8:
            prod = " * ".join([b] * abs(int(exp)))
            return f"({prod})" if exp > 0 else f"(1.0 / ({prod}))"
        if exp.is_Rational and exp.q == 2 and abs(exp.p) <= 7:
            k = abs(int(exp.p))
            root = f"std::sqrt({self._print(base)})"
            parts = [b] * (k // 2) + [root]
            prod = " * ".join(parts)
            return f"({prod})" if exp > 0 else f"(1.0 / ({prod}))"
        return super()._print_Pow(expr)


PRINTER = FastPowPrinter()


def emit_stage(name, s):
    exprs = list(s["outputs"]) + [e for row in s["jac"] for e in row]
    repl, red = sp.cse(exprs, symbols=sp.numbered_symbols("t"), optimizations="basic")
    n_in = len(s["inputs"])
    n_out = len(s["outputs"])
    consts = s["consts"]
    sub = {c: sp.Symbol(f"c[{i}]") for i, c in enumerate(consts)}
    sub.update({x: sp.Symbol(f"x[{i}]") for i, x in enumerate(s["inputs"])})

    def cxx(e):
        return PRINTER.doprint(e.xreplace(sub))

    lines = [f"  const real {sym} = {cxx(e)};" for sym, e in repl]
    for o in range(n_out):
        lines.append(f"  y[{o}] = {cxx(red[o])};")
    for o in range(n_out):
        for i in range(n_in):
            lines.append(f"  J({o}, {i}) = {cxx(red[n_out + o * n_in + i])};")
    cparam = f"const Eigen::Matrix<real, {len(consts)}, 1> &c, " if consts else ""
    ops = sum(sp.count_ops(e) for _, e in repl) + sum(sp.count_ops(e) for e in red)
    print(f"  {name}: {len(repl)} temporaries, ~{ops} ops")
    body = "\n".join(lines)
    return f"""/// {DOCS[name]}
inline void dupin_{name}({cparam}const Eigen::Matrix<real, {n_in}, 1> &x,
                         Eigen::Matrix<real, {n_out}, 1> &y, Eigen::Matrix<real, {n_out}, {n_in}> &J) {{
{body}
}}
"""


def emit(st):
    t0 = time.time()
    funcs = "\n".join(emit_stage(name, s) for name, s in st.items())
    print(f"cse in {time.time() - t0:.1f}s")
    text = f"""// AUTO-GENERATED by scripts/sympy/osculating_dupin_jacobian.py --emit — do not edit
#pragma once

#include <cmath>
#include <Eigen/Dense>

#include "gaudi/common.h"

namespace gaudi {{
namespace hepworth {{
namespace block {{
namespace dupin_generated {{

// Stages of the ring point's gap against the posed osculating Dupin cyclide.
// `dupin_gap_and_gradient` in `osculating_dupin_ring.hpp` chains them. Every
// stage assumes the point lands, both sags inside their circles; the caller
// checks with `land_on_posed_dupin` first.

{funcs}
}} // namespace dupin_generated
}} // namespace block
}} // namespace hepworth
}} // namespace gaudi
"""
    HEADER.write_text(text)
    print(f"wrote {HEADER} ({len(text.splitlines())} lines)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    st = jacobians(stages())
    print(f"stage Jacobians in {time.time() - t0:.1f}s")
    verify(lambdas(st))
    if args.emit:
        emit(st)


if __name__ == "__main__":
    main()
