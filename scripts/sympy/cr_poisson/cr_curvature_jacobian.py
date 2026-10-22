"""Symbolic per-face CR shape operator and its Jacobian (Stein et al. 2020,
App. A.2, with the CR derivative recomputed from the stencil, so the
read-out follows the metric).

Face f, stencil P0..P5: the face (P0, P1, P2) in winding order and, for local
edge i = (P_i, P_{i+1}), the opposite vertex of the neighboring face in
P_{3+i}. `h_i` is 1 when that neighbor exists, 0 on a boundary edge (then the
caller puts a reflected dummy, P_{3+i} = P_i + P_{i+1} - P_{i+2}, in the slot
so the masked terms stay finite).

    s_f = (k11, k12, k22) of S_f in the face frame T = (P1 - P0) / |P1 - P0|,
    B = N x T, with
    S_f = K sum_i l_i / A (N . g_perp,i) v_perp,i v_perp,i^T,
    g_perp,i = (D_f + h_i D_g)_i P / M_e,  M_e = (A + h_i A_g) / 3.

The orientation sign of each edge appears twice (in `g_perp` and in the
read-out) and cancels, so the face's own winding is used throughout.

Run: python3 cr_curvature_jacobian.py          # verify
     python3 cr_curvature_jacobian.py --emit   # also write the C++ header
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import sympy as sp

K_CR_SHAPE_SIGN = -1

HEADER = (Path(__file__).resolve().parents[1] / "generated"
          / "cr_curvature_stencil_generated.hpp")


def symbolic():
    x, h, s, _ = symbolic_parts()
    return x, h, s


def symbolic_edge_frames():
    """Angle of each local edge in the face frame, as `(cos, sin)`:
    `c_i = T . v_i / |v_i|`, `s_i = B . v_i / |v_i|`, functions of the face's
    three corners only (9 coordinates). The edge reading of the face's
    operator is the operator rotated by that angle; differentiating only the
    angles and composing with the face read-out's Jacobian keeps the
    expressions small (differentiating the 9 readings symbolically does not
    finish)."""
    x = sp.symbols("x0:9", real=True)
    P = [sp.Matrix(x[3 * k:3 * k + 3]) for k in range(3)]
    n = (P[1] - P[0]).cross(P[2] - P[0])
    N = n / sp.sqrt(n.dot(n))
    T = (P[1] - P[0]) / sp.sqrt((P[1] - P[0]).dot(P[1] - P[0]))
    B = N.cross(T)
    out = []
    for i in range(3):
        v = P[(i + 1) % 3] - P[i]
        l = sp.sqrt(v.dot(v))
        out += [T.dot(v) / l, B.dot(v) / l]
    return x, sp.Matrix(out)  # (c0, s0, c1, s1, c2, s2)


def edge_readings_numeric(fs, fJ, fa, fJa):
    """Readings `R (9)` and `J (9 x 18)` composed from the face read-out and
    the edge angles (the same composition the C++ header emits)."""
    def rot(s, c, sn):
        k11, k12, k22 = s
        R = np.array([c * c * k11 + 2 * c * sn * k12 + sn * sn * k22,
                      -c * sn * k11 + (c * c - sn * sn) * k12 + c * sn * k22,
                      sn * sn * k11 - 2 * c * sn * k12 + c * c * k22])
        dRds = np.array([[c * c, 2 * c * sn, sn * sn],
                         [-c * sn, c * c - sn * sn, c * sn],
                         [sn * sn, -2 * c * sn, c * c]])
        dRdc = np.array([2 * c * k11 + 2 * sn * k12, -sn * k11 + 2 * c * k12 + sn * k22,
                         2 * c * k22 - 2 * sn * k12])
        dRdsn = np.array([2 * c * k12 + 2 * sn * k22, -c * k11 - 2 * sn * k12 + c * k22,
                          2 * sn * k11 - 2 * c * k12])
        return R, dRds, dRdc, dRdsn

    def f(P, hh):
        s, Js = fs(P, hh), fJ(P, hh)
        a, Ja = fa(P[:3]), fJa(P[:3])
        R = np.zeros(9)
        J = np.zeros((9, 18))
        for i in range(3):
            Ri, dRds, dRdc, dRdsn = rot(s, a[2 * i], a[2 * i + 1])
            R[3 * i:3 * i + 3] = Ri
            J[3 * i:3 * i + 3] = dRds @ Js
            J[3 * i:3 * i + 3, :9] += np.outer(dRdc, Ja[2 * i]) + np.outer(dRdsn, Ja[2 * i + 1])
        return R, J
    return f


def symbolic_parts():
    x = sp.symbols("x0:18", real=True)
    h = sp.symbols("h0:3", real=True)
    P = [sp.Matrix(x[3 * k:3 * k + 3]) for k in range(6)]

    def norm(v):
        return sp.sqrt(v.dot(v))

    n = (P[1] - P[0]).cross(P[2] - P[0])
    nn = norm(n)
    A = nn / 2
    N = n / nn
    T = (P[1] - P[0]) / norm(P[1] - P[0])
    B = N.cross(T)
    s = sp.zeros(3, 1)
    for i in range(3):
        t, hd, o, o2 = i, (i + 1) % 3, (i + 2) % 3, 3 + i
        v = P[hd] - P[t]
        l = norm(v)
        ejk = (P[o] - P[hd]).dot(P[o] - P[hd])
        eki = (P[t] - P[o]).dot(P[t] - P[o])
        k = -1 / (12 * l)
        row = (k * (l**2 - (eki - ejk)) * P[t] + k * (l**2 + (eki - ejk)) * P[hd]
               + l / 6 * P[o])
        Ag = norm((P[t] - P[hd]).cross(P[o2] - P[hd])) / 2
        ejk2 = (P[o2] - P[t]).dot(P[o2] - P[t])
        eki2 = (P[hd] - P[o2]).dot(P[hd] - P[o2])
        k2 = 1 / (12 * l)
        row += h[i] * (k2 * (l**2 - (eki2 - ejk2)) * P[hd]
                       + k2 * (l**2 + (eki2 - ejk2)) * P[t] - l / 6 * P[o2])
        Me = (A + h[i] * Ag) / 3
        gperp = row / Me
        vperp = N.cross(v / l)  # unit: v lies in the face plane
        kk = K_CR_SHAPE_SIGN * l / A * N.dot(gperp)
        tv, bv = T.dot(vperp), B.dot(vperp)
        s += kk * sp.Matrix([tv * tv, tv * bv, bv * bv])
    return x, h, s, {"P": P, "T": T, "B": B}


def build(verbose=True):
    t0 = time.time()
    x, h, s = symbolic()
    J = s.jacobian(x)
    if verbose:
        print(f"symbolic S and J: {time.time() - t0:.1f}s")
    return x, h, s, J


def numeric(x, h, s, J):
    args = list(x) + list(h)
    fs = sp.lambdify(args, s, "numpy", cse=True)
    fJ = sp.lambdify(args, J, "numpy", cse=True)
    return (lambda P, hh: np.asarray(fs(*P.ravel(), *hh), float).ravel(),
            lambda P, hh: np.asarray(fJ(*P.ravel(), *hh), float))


def stencils_with_dummies(st, X):
    """Stencil positions (m, 6, 3) and flags (m, 3); boundary slots get the
    reflected dummy."""
    P = X[st.ids].copy()
    hh = st.mask[:, 3:].copy()
    for i in range(3):
        b = hh[:, i] == 0
        P[b, 3 + i] = P[b, i] + P[b, (i + 1) % 3] - P[b, (i + 2) % 3]
    return P, hh


def verify(fs, fJ):
    import cr_ops as ops
    import face_loops as fl

    rng = np.random.default_rng(3)

    def meshes():
        X, F, _ = ops.torus(2.0, 0.7, 12, 8)
        yield "torus", X, F
        import closed_targets as ct
        X, F = ct.cube(3)
        yield "cube", X, F
        Pp, F = ops.strip_param(6, 4, 2.0, 1.0)
        yield "strip (boundary)", ops.cylinder_from_param(Pp, 1.5), F

    worst = {}
    for name, X0, F in meshes():
        m = ops.build_mesh(F, len(X0))
        st = fl.build_stencils(m, ops.cr_derivative(m, X0)[0])
        X = X0 + 0.03 * rng.normal(size=X0.shape) * np.ptp(X0, 0).max()
        P, hh = stencils_with_dummies(st, X)
        ref = fl.stencil_s3_live(st, X[st.ids])
        e_val = e_J = e_inv = 0.0
        for f in range(len(P)):
            s_sym = fs(P[f], hh[f])
            e_val = max(e_val, np.abs(s_sym - ref[f]).max() / (1 + np.abs(ref[f]).max()))
            Jf = fJ(P[f], hh[f])
            mask = np.repeat(np.r_[1, 1, 1, hh[f]], 3)
            Jfd = np.zeros((3, 18))
            eps = 1e-6
            for d in range(18):
                Pp, Pm = P[f].copy(), P[f].copy()
                Pp.ravel()[d] += eps
                Pm.ravel()[d] -= eps
                Jfd[:, d] = (fs(Pp, hh[f]) - fs(Pm, hh[f])) / (2 * eps)
            e_J = max(e_J, np.abs((Jf - Jfd) * mask).max() / (1 + np.abs(Jfd).max()))
            # Invariances: translation and rotation leave S unchanged; S is
            # homogeneous of degree -1 in the positions (J p = -s).
            Pf = P[f].ravel()
            tr = np.tile(rng.normal(size=3), 6)
            w = rng.normal(size=3)
            rot = np.cross(np.tile(w, (6, 1)), P[f]).ravel()
            e_inv = max(e_inv, np.abs(Jf @ tr).max(), np.abs(Jf @ rot).max(),
                        np.abs(Jf @ Pf + s_sym).max())
        worst[name] = (e_val, e_J, e_inv)
        print(f"{name:17s} value vs live read-out {e_val:.1e}  J vs FD {e_J:.1e}  "
              f"invariances {e_inv:.1e}  faces {len(P)}")
    return worst


def emit(x, h, s, J):
    exprs = list(s) + list(J)
    t0 = time.time()
    repl, red = sp.cse(exprs, symbols=sp.numbered_symbols("t"), optimizations="basic")
    ops_count = sum(sp.count_ops(e) for _, e in repl) + sum(sp.count_ops(e) for e in red)
    print(f"cse: {len(repl)} temporaries, ~{ops_count} ops, {time.time() - t0:.1f}s")
    sub = {x[k]: sp.Symbol(f"P[{k // 3}][{k % 3}]") for k in range(18)}
    sub.update({h[i]: sp.Symbol(f"h[{i}]") for i in range(3)})

    def cxx(e):
        return sp.cxxcode(e.xreplace(sub), standard="c++17")

    lines = []
    for sym, e in repl:
        lines.append(f"  const real {sym} = {cxx(e)};")
    for r in range(3):
        lines.append(f"  s[{r}] = {cxx(red[r])};")
    for r in range(3):
        for c in range(18):
            lines.append(f"  J({r}, {c}) = {cxx(red[3 + 18 * r + c])};")
    body = "\n".join(lines)
    text = f"""// AUTO-GENERATED by scripts/sympy/cr_poisson/cr_curvature_jacobian.py — do not edit
#pragma once

#include <array>
#include <cmath>
#include <Eigen/Dense>

#include "gaudi/common.h"

namespace gaudi {{
namespace hepworth {{
namespace block {{
namespace cr_generated {{

/// Per-face CR shape operator read-out and its Jacobian over the face's
/// stencil (Stein et al. 2020, App. A.2, with the CR derivative recomputed
/// from the stencil). `P[0..2]`: the face in winding order; `P[3 + i]`: the
/// vertex opposite local edge `(P[i], P[i + 1])` in the neighboring face.
/// `h[i]` is 1 when that neighbor exists, 0 on a boundary edge, where
/// `P[3 + i]` must be the reflected dummy `P[i] + P[i + 1] - P[i + 2]`.
/// `s = (k11, k12, k22)` in the frame `T = (P1 - P0) / |P1 - P0|`, `B = N x T`;
/// `J(r, 3 k + c) = d s_r / d P[k][c]`.
inline void cr_curvature_stencil(const std::array<vec3, 6> &P, const vec3 &h,
                                 vec3 &s, Eigen::Matrix<real, 3, 18> &J) {{
{body}
}}

}} // namespace cr_generated
}} // namespace block
}} // namespace hepworth
}} // namespace gaudi
"""
    HEADER.write_text(text)
    print(f"wrote {HEADER} ({len(text.splitlines())} lines)")


EDGE_HEADER = HEADER.parent / "cr_curvature_edge_generated.hpp"


def verify_edge(fRJ):
    """Edge readings against `shape_dirichlet.face_edge_readings` (the live
    Python read-out rotated into each edge frame) and J against central
    differences."""
    fR = lambda P, hh: fRJ(P, hh)[0]  # noqa: E731
    fJR = lambda P, hh: fRJ(P, hh)[1]  # noqa: E731
    import cr_ops as ops
    import face_loops as fl
    import shape_dirichlet as sd

    rng = np.random.default_rng(5)
    X, F, _ = ops.torus(2.0, 0.7, 12, 8)
    m = ops.build_mesh(F, len(X))
    st = fl.build_stencils(m, ops.cr_derivative(m, X)[0])
    X = X + 0.05 * rng.normal(size=X.shape)
    P, hh = stencils_with_dummies(st, X)
    ref = sd.face_edge_readings(st, X[st.ids]).reshape(len(P), 9)
    e_val = e_J = 0.0
    for f in range(len(P)):
        R = fR(P[f], hh[f])
        e_val = max(e_val, np.abs(R - ref[f]).max() / (1 + np.abs(ref[f]).max()))
        J = fJR(P[f], hh[f])
        Jfd = np.zeros((9, 18))
        for d in range(18):
            Pp, Pm = P[f].copy(), P[f].copy()
            Pp.ravel()[d] += 1e-6
            Pm.ravel()[d] -= 1e-6
            Jfd[:, d] = (fR(Pp, hh[f]) - fR(Pm, hh[f])) / 2e-6
        e_J = max(e_J, np.abs(J - Jfd).max() / (1 + np.abs(Jfd).max()))
    print(f"edge readings vs python {e_val:.1e}  J vs FD {e_J:.1e}")


def emit_edge(xa, a, Ja):
    exprs = list(a) + list(Ja)
    t0 = time.time()
    repl, red = sp.cse(exprs, symbols=sp.numbered_symbols("t"), optimizations="basic")
    print(f"edge-angle cse: {len(repl)} temporaries, {time.time() - t0:.1f}s")
    sub = {xa[k]: sp.Symbol(f"P[{k // 3}][{k % 3}]") for k in range(9)}

    def cxx(e):
        return sp.cxxcode(e.xreplace(sub), standard="c++17")

    lines = [f"  const real {sym} = {cxx(e)};" for sym, e in repl]
    lines += [f"  a[{r}] = {cxx(red[r])};" for r in range(6)]
    lines += [f"  Ja({r}, {c}) = {cxx(red[6 + 9 * r + c])};"
              for r in range(6) for c in range(9)]
    body = "\n".join(lines)
    text = f"""// AUTO-GENERATED by scripts/sympy/cr_poisson/cr_curvature_jacobian.py --emit-edge — do not edit
#pragma once

#include <array>
#include <cmath>
#include <Eigen/Dense>

#include "gaudi/common.h"
#include "cr_curvature_stencil_generated.hpp"

namespace gaudi {{
namespace hepworth {{
namespace block {{
namespace cr_generated {{

/// Angle of each local edge `(P[i], P[i + 1])` in the face frame
/// `T = (P1 - P0) / |P1 - P0|`, `B = N x T`, as `a = (c0, s0, c1, s1, c2, s2)`,
/// `c_i = T . v_i / |v_i|`, `s_i = B . v_i / |v_i|`; `Ja(r, 3 k + c) = d a_r / d P[k][c]`.
inline void cr_curvature_edge_angles(const std::array<vec3, 6> &P,
                                     Eigen::Matrix<real, 6, 1> &a,
                                     Eigen::Matrix<real, 6, 9> &Ja) {{
{body}
}}

/// The face's CR shape operator (`cr_curvature_stencil`) read in each of its
/// three edges' frames `(t_i, N x t_i)`, `t_i` along local edge `i`
/// (orientation free: `t` and `b` flip together). `R[3 i + c]`: edge `i`,
/// component `c` of `(k11, k12, k22)`; `J(r, 3 k + c) = d R_r / d P[k][c]`.
/// The reading is the face operator rotated by the edge angle, so `J` is
/// composed from the face read-out's Jacobian and the angles' Jacobian.
inline void cr_curvature_edge_readings(const std::array<vec3, 6> &P, const vec3 &h,
                                       Eigen::Matrix<real, 9, 1> &R,
                                       Eigen::Matrix<real, 9, 18> &J) {{
  vec3 s;
  Eigen::Matrix<real, 3, 18> Js;
  cr_curvature_stencil(P, h, s, Js);
  Eigen::Matrix<real, 6, 1> a;
  Eigen::Matrix<real, 6, 9> Ja;
  cr_curvature_edge_angles(P, a, Ja);
  const real k11 = s[0], k12 = s[1], k22 = s[2];
  for (int i = 0; i < 3; ++i) {{
    const real c = a[2 * i], sn = a[2 * i + 1];
    R[3 * i + 0] = c * c * k11 + 2 * c * sn * k12 + sn * sn * k22;
    R[3 * i + 1] = -c * sn * k11 + (c * c - sn * sn) * k12 + c * sn * k22;
    R[3 * i + 2] = sn * sn * k11 - 2 * c * sn * k12 + c * c * k22;
    Eigen::Matrix<real, 3, 3> dRds;
    dRds << c * c, 2 * c * sn, sn * sn,
        -c * sn, c * c - sn * sn, c * sn,
        sn * sn, -2 * c * sn, c * c;
    const vec3 dRdc(2 * c * k11 + 2 * sn * k12, -sn * k11 + 2 * c * k12 + sn * k22,
                    2 * c * k22 - 2 * sn * k12);
    const vec3 dRdsn(2 * c * k12 + 2 * sn * k22, -c * k11 - 2 * sn * k12 + c * k22,
                     2 * sn * k11 - 2 * c * k12);
    J.block<3, 18>(3 * i, 0) = dRds * Js;
    J.block<3, 9>(3 * i, 0) += dRdc * Ja.row(2 * i) + dRdsn * Ja.row(2 * i + 1);
  }}
}}

}} // namespace cr_generated
}} // namespace block
}} // namespace hepworth
}} // namespace gaudi
"""
    EDGE_HEADER.write_text(text)
    print(f"wrote {EDGE_HEADER} ({len(text.splitlines())} lines)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", action="store_true")
    ap.add_argument("--emit-edge", action="store_true")
    args = ap.parse_args()
    if args.emit_edge:
        x, h, s, J = build()
        fs, fJ = numeric(x, h, s, J)
        xa, a = symbolic_edge_frames()
        Ja = a.jacobian(xa)
        fa = sp.lambdify(list(xa), a, "numpy", cse=True)
        fJa = sp.lambdify(list(xa), Ja, "numpy", cse=True)
        verify_edge(edge_readings_numeric(
            fs, fJ, lambda P3: np.asarray(fa(*P3.ravel()), float).ravel(),
            lambda P3: np.asarray(fJa(*P3.ravel()), float)))
        emit_edge(xa, a, Ja)
    else:
        x, h, s, J = build()
        fs, fJ = numeric(x, h, s, J)
        verify(fs, fJ)
        if args.emit:
            emit(x, h, s, J)
