"""Stein's face read-out in the discrete fundamental forms (MATH.md section 16).

Face f's stencil (the face and the three vertices opposite its edges) is
rebuilt from 9 edge lengths and 3 dihedrals,

    ell = (l0, l1, l2, lt0, lh0, lt1, lh1, lt2, lh2),  theta = (th0, th1, th2),

l_i the face's local edge i = (P_i, P_{i+1}); lt_i, lh_i the wing's edges
from that edge's tail and head to the opposite vertex P_{3+i}; th_i the
dihedral across edge i in the sync convention (face_loops.dihedrals: the turn
about the edge from the + face's normal to the - face's, the same seen from
either face). The face lies in the plane z = 0 (P0 at the origin, P1 on +x),
each wing is placed opposite the face's third vertex and turned about the edge
by -(th + pi). The read-out is rigid-invariant, so

    J_ff = d s / d (log ell, theta) = J_P . dP / d(ell, theta) . diag(ell, 1)

with J_P the generated stencil Jacobian (cr_curvature_jacobian.py). Only the
rebuild and its Jacobian are new. Boundary edges (no wing) are the caller's:
put the reflected dummy in the slot, h_i = 0, and drop that wing's columns.

Run: python3 cr_ff_jacobian.py          # verify
     python3 cr_ff_jacobian.py --emit   # also write the C++ header
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import sympy as sp

HEADER = (Path(__file__).resolve().parents[1] / "generated" / "cr_ff_rebuild_generated.hpp")


def symbolic():
    """Rebuilt stencil P (18) and dP / d(ell, theta) (18 x 12)."""
    ell = sp.symbols("l0:9", positive=True)
    th = sp.symbols("th0:3", real=True)
    l0, l1, l2 = ell[:3]
    x = (l0 ** 2 + l2 ** 2 - l1 ** 2) / (2 * l0)
    y = sp.sqrt(l2 ** 2 - x ** 2)
    P = [sp.Matrix([0, 0, 0]), sp.Matrix([l0, 0, 0]), sp.Matrix([x, y, 0])]
    for i in range(3):
        Pt, Ph = P[i], P[(i + 1) % 3]
        li = ell[i]
        eh = (Ph - Pt) / li  # |Ph - Pt| = l_i by construction
        dt, dh = ell[3 + 2 * i], ell[4 + 2 * i]
        xw = (dt ** 2 + li ** 2 - dh ** 2) / (2 * li)
        yw = sp.sqrt(dt ** 2 - xw ** 2)
        perp = sp.Matrix([-eh[1], eh[0], 0])  # points to the face's side (counterclockwise face)
        v = xw * eh - yw * perp
        c, s = -sp.cos(th[i]), sp.sin(th[i])  # turn by -(th + pi)
        v = v * c + eh.cross(v) * s + eh * (eh.dot(v)) * (1 - c)
        P.append(Pt + v)
    Pf = sp.Matrix([P[k][c] for k in range(6) for c in range(3)])
    vars_ = list(ell) + list(th)
    return vars_, Pf, Pf.jacobian(vars_)


def numeric(vars_, Pf, dP):
    fP = sp.lambdify(vars_, Pf, "numpy", cse=True)
    fdP = sp.lambdify(vars_, dP, "numpy", cse=True)
    return (lambda v: np.asarray(fP(*v), float).reshape(6, 3),
            lambda v: np.asarray(fdP(*v), float))


def forms_from_stencil(P):
    """(ell (9), theta (3), h (3)) of one face's live stencil P (6, 3); missing
    wings (h_i = 0) get ell = 1, theta = 0 placeholders."""
    ell = np.ones(9)
    th = np.zeros(3)
    h = np.ones(3)
    for i in range(3):
        ell[i] = np.linalg.norm(P[(i + 1) % 3] - P[i])
    n = np.cross(P[1] - P[0], P[2] - P[0])
    Nf = n / np.linalg.norm(n)
    for i in range(3):
        Pt, Ph, Po = P[i], P[(i + 1) % 3], P[3 + i]
        if not np.all(np.isfinite(Po)):
            h[i] = 0
            continue
        ell[3 + 2 * i] = np.linalg.norm(Po - Pt)
        ell[4 + 2 * i] = np.linalg.norm(Po - Ph)
        ax = (Ph - Pt) / np.linalg.norm(Ph - Pt)
        Ng = np.cross(Po - Ph, Pt - Ph)
        Ng /= np.linalg.norm(Ng)
        phi = np.arctan2(np.cross(Nf, Ng) @ ax, Nf @ Ng)
        th[i] = (phi - np.pi + np.pi) % (2 * np.pi) - np.pi
    return ell, th, h


def verify(fP, fdP):
    import cr_curvature_jacobian as cj
    import cr_ops as ops
    import face_loops as fl

    x, hs, s, J = cj.build(verbose=False)
    fs, fJ = cj.numeric(x, hs, s, J)
    rng = np.random.default_rng(5)

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
        ref = fl.stencil_s3_live(st, X[st.ids])
        e_val = e_J = 0.0
        for f in range(len(m.F)):
            Pl = X[st.ids[f]].copy()
            for i in range(3):
                if st.mask[f, 3 + i] == 0:
                    Pl[3 + i] = np.nan
            ell, th, h = forms_from_stencil(Pl)
            v = np.r_[ell, th]
            P = fP(v)
            for i in range(3):  # boundary: the reflected dummy
                if h[i] == 0:
                    P[3 + i] = P[i] + P[(i + 1) % 3] - P[(i + 2) % 3]
            s_reb = fs(P, h)
            e_val = max(e_val, np.abs(s_reb - ref[f]).max() / (1 + np.abs(ref[f]).max()))
            Jp = fJ(P, h) * np.repeat(np.r_[1, 1, 1, h], 3)[None, :]
            Jff = Jp @ fdP(v)
            Jff[:, :9] *= ell[None, :]  # log lengths
            keep = np.r_[np.ones(3), np.repeat(h, 2), h] > 0
            # central differences of the rebuilt read-out in (log ell, theta)
            Jfd = np.zeros((3, 12))
            eps = 1e-6
            for d in range(12):
                if not keep[d]:
                    continue
                vp, vm = v.copy(), v.copy()
                if d < 9:
                    vp[d] *= np.exp(eps)
                    vm[d] *= np.exp(-eps)
                else:
                    vp[d] += eps
                    vm[d] -= eps

                def sval(w):
                    Q = fP(w)
                    for i in range(3):
                        if h[i] == 0:
                            Q[3 + i] = Q[i] + Q[(i + 1) % 3] - Q[(i + 2) % 3]
                    return fs(Q, h)
                Jfd[:, d] = (sval(vp) - sval(vm)) / (2 * eps)
            e_J = max(e_J, np.abs((Jff - Jfd)[:, keep]).max() / (1 + np.abs(Jfd).max()))
        worst[name] = (e_val, e_J)
        print(f"{name:17s} rebuilt read-out vs live {e_val:.1e}  J_ff vs FD {e_J:.1e}  faces {len(m.F)}")
    return worst


def emit(vars_, Pf, dP):
    exprs = list(Pf) + list(dP)
    t0 = time.time()
    repl, red = sp.cse(exprs, symbols=sp.numbered_symbols("t"), optimizations="basic")
    ops_count = sum(sp.count_ops(e) for _, e in repl) + sum(sp.count_ops(e) for e in red)
    print(f"cse: {len(repl)} temporaries, ~{ops_count} ops, {time.time() - t0:.1f}s")
    sub = {vars_[k]: sp.Symbol(f"l[{k}]") for k in range(9)}
    sub.update({vars_[9 + i]: sp.Symbol(f"th[{i}]") for i in range(3)})

    def cxx(e):
        return sp.cxxcode(e.xreplace(sub), standard="c++17")

    lines = [f"  const real {sym} = {cxx(e)};" for sym, e in repl]
    for k in range(6):
        for c in range(3):
            lines.append(f"  P[{k}][{c}] = {cxx(red[3 * k + c])};")
    for r in range(18):
        for c in range(12):
            lines.append(f"  dP({r}, {c}) = {cxx(red[18 + 12 * r + c])};")
    body = "\n".join(lines)
    text = f"""// AUTO-GENERATED by scripts/sympy/cr_poisson/cr_ff_jacobian.py --emit — do not edit
#pragma once

#include <array>
#include <cmath>
#include <Eigen/Dense>

#include "gaudi/common.h"

namespace gaudi {{
namespace hepworth {{
namespace block {{
namespace cr_generated {{

/// A face's 6-point stencil rebuilt from its discrete fundamental forms:
/// `l = (l0, l1, l2, lt0, lh0, lt1, lh1, lt2, lh2)` (the face's local edges
/// `(P[i], P[i + 1])`, then each wing's edges from that edge's tail and head
/// to `P[3 + i]`) and `th` the dihedral across each local edge in the sync
/// convention. The face lies in z = 0 with `P[0]` at the origin and `P[1]`
/// on +x; `dP(3 k + c, j) = d P[k][c] / d (l, th)_j`. Compose with
/// `cr_curvature_stencil` for `d s / d (l, th)`. Boundary slots are the
/// caller's (reflected dummy, `h_i = 0`, drop that wing's columns).
inline void cr_ff_rebuild(const std::array<real, 9> &l, const vec3 &th,
                          std::array<vec3, 6> &P, Eigen::Matrix<real, 18, 12> &dP) {{
{body}
}}

}} // namespace cr_generated
}} // namespace block
}} // namespace hepworth
}} // namespace gaudi
"""
    HEADER.write_text(text)
    print(f"wrote {HEADER} ({len(text.splitlines())} lines)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    vars_, Pf, dP = symbolic()
    print(f"symbolic rebuild and dP: {time.time() - t0:.1f}s")
    fP, fdP = numeric(vars_, Pf, dP)
    verify(fP, fdP)
    if args.emit:
        emit(vars_, Pf, dP)
