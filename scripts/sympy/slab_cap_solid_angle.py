#!/usr/bin/env python3
"""Cap wedge solid-angle + barrier energies for slab pierce fairness.

Hub-local spokes e0, e1 (from foot/pierce to ring verts). Smooth solid angle
at the hub for the spherical lune between e0 and e1:

  Ω(e0, e1) = 2 atan2(|e0×e1|, e0·e1 + |e0||e1|)

Optional trivector (existing wedge volume proxy):

  τ = (e0 × e1) · d,   d = c_pierce − c_foot

Barrier / fairness terms (per wedge, rest Ω_rest baked at init):

  E_rest   = (Ω − Ω_rest)²
  E_log    = −log(Ω + ε)
  E_inv    = 1 / (Ω + ε)
  E_inv_m  = (1/(Ω+ε) − 1/(Ω_rest+ε))²

Combined (try in GD / local GN):

  E = w (E_rest + λ_b E_barrier) + w_τ (τ − τ_rest)²

Chain rule to hub/spoke world points x0, x1, c with e0=x0−c, e1=x1−c.

Usage:
  python3 scripts/sympy/slab_cap_solid_angle.py
  python3 scripts/sympy/slab_cap_solid_angle.py --fd-check
  python3 scripts/sympy/slab_cap_solid_angle.py --emit
"""

from __future__ import annotations

import argparse
from typing import Iterable

import sympy as sp
from sympy import Matrix, symbols

# Hub-relative spokes
e0x, e0y, e0z, e1x, e1y, e1z = symbols("e0x e0y e0z e1x e1y e1z", real=True)
dx, dy, dz = symbols("dx dy dz", real=True)

e0 = Matrix([e0x, e0y, e0z])
e1 = Matrix([e1x, e1y, e1z])
d = Matrix([dx, dy, dz])

Omega_rest, tau_rest = symbols("Omega_rest tau_rest", real=True, positive=True)
Omega_min = symbols("Omega_min", real=True, positive=True)
w, w_tau, lam_barrier = symbols("w w_tau lam_barrier", real=True, nonnegative=True)
eps = symbols("eps", real=True, positive=True)
lam = symbols("lam", real=True, positive=True)


def solid_angle_expr(e0_m: Matrix, e1_m: Matrix) -> sp.Expr:
    """Smooth hub solid angle between directions e0, e1 (Van Oosterom lune form)."""
    cross = e0_m.cross(e1_m)
    num = sp.sqrt(sp.Max(cross.dot(cross), 0))  # for codegen use raw norm in C++
    # SymPy path: use sqrt(dot) — C++ uses cross.norm()
    num_sym = sp.sqrt(cross.dot(cross))
    n0 = sp.sqrt(e0_m.dot(e0_m))
    n1 = sp.sqrt(e1_m.dot(e1_m))
    den = e0_m.dot(e1_m) + n0 * n1
    return 2 * sp.atan2(num_sym, den)


def tau_expr(e0_m: Matrix, e1_m: Matrix, d_m: Matrix) -> sp.Expr:
    return e0_m.cross(e1_m).dot(d_m)


def barrier_exprs(omega: sp.Expr) -> dict[str, sp.Expr]:
    return {
        "rest": (omega - Omega_rest) ** 2,
        "log": -sp.log(omega + eps),
        "inv": 1 / (omega + eps),
        "inv_match": (1 / (omega + eps) - 1 / (Omega_rest + eps)) ** 2,
    }


def combined_energy(
    include_tau: bool = True,
    barrier_key: str = "log",
) -> sp.Expr:
    omega = solid_angle_expr(e0, e1)
    barriers = barrier_exprs(omega)
    e = w * (barriers["rest"] + lam_barrier * barriers[barrier_key])
    if include_tau:
        tau = tau_expr(e0, e1, d)
        e += w_tau * (tau - tau_rest) ** 2
    return sp.simplify(e)


def grad_wrt(vec: Matrix, energy: sp.Expr) -> Matrix:
    return Matrix([sp.diff(energy, v) for v in vec])


def hub_spoke_grads(energy: sp.Expr) -> dict[str, Matrix]:
    """e0=p0-c, e1=p1-c → grads w.r.t. c, p0, p1."""
    hc, hy, hz = symbols("hc hy hz", real=True)
    p0x, p0y, p0z, p1x, p1y, p1z = symbols("p0x p0y p0z p1x p1y p1z", real=True)
    hub = Matrix([hc, hy, hz])
    p0 = Matrix([p0x, p0y, p0z])
    p1 = Matrix([p1x, p1y, p1z])
    e0_w = p0 - hub
    e1_w = p1 - hub
    e_sub = energy.subs(
        {
            e0x: e0_w[0],
            e0y: e0_w[1],
            e0z: e0_w[2],
            e1x: e1_w[0],
            e1y: e1_w[1],
            e1z: e1_w[2],
        }
    )
    return {
        "c": grad_wrt(hub, e_sub),
        "p0": grad_wrt(p0, e_sub),
        "p1": grad_wrt(p1, e_sub),
        "d": grad_wrt(d, e_sub),
    }


def _symbol_map() -> dict:
    hc, hy, hz = symbols("hc hy hz", real=True)
    p0x, p0y, p0z, p1x, p1y, p1z = symbols("p0x p0y p0z p1x p1y p1z", real=True)
    return {
        e0x: "e0[0]",
        e0y: "e0[1]",
        e0z: "e0[2]",
        e1x: "e1[0]",
        e1y: "e1[1]",
        e1z: "e1[2]",
        dx: "d[0]",
        dy: "d[1]",
        dz: "d[2]",
        Omega_rest: "Omega_rest",
        tau_rest: "tau_rest",
        w: "w",
        w_tau: "w_tau",
        lam_barrier: "lam_barrier",
        lam: "lam",
        Omega_min: "omega_min",
        eps: "eps",
        hc: "c[0]",
        hy: "c[1]",
        hz: "c[2]",
        p0x: "x0[0]",
        p0y: "x0[1]",
        p0z: "x0[2]",
        p1x: "x1[0]",
        p1y: "x1[1]",
        p1z: "x1[2]",
    }


def _emit_vec3(name: str, args: str, comps: Iterable[sp.Expr], printer) -> str:
    body = ", ".join(printer.doprint(c) for c in comps)
    return f"inline vec3 {name}({args}) {{\n  return vec3({body});\n}}\n"


def _cse_body(exprs: list, printer) -> tuple[str, list[str]]:
    repls, reduced = sp.cse(
        exprs, optimizations="basic", symbols=sp.numbered_symbols("cse")
    )
    lines = []
    for sym, sub in repls:
        lines.append(f"  const real {sym} = {printer.doprint(sub)};")
    return "\n".join(lines), [printer.doprint(r) for r in reduced]


def emit_cap_solid_angle_generated() -> str:
    from albers_printer import AlbersPrinter

    printer = AlbersPrinter(symbol_map=_symbol_map())

    omega = solid_angle_expr(e0, e1)
    tau = tau_expr(e0, e1, d)
    barriers = barrier_exprs(omega)

    e_log = w * (barriers["rest"] + lam_barrier * barriers["log"]) + w_tau * (
        tau - tau_rest
    ) ** 2
    e_inv = w * (barriers["rest"] + lam_barrier * barriers["inv"]) + w_tau * (
        tau - tau_rest
    ) ** 2
    e_inv_m = w * (barriers["rest"] + lam_barrier * barriers["inv_match"]) + w_tau * (
        tau - tau_rest
    ) ** 2

    g_e0 = grad_wrt(e0, e_log)
    g_e1 = grad_wrt(e1, e_log)
    world = hub_spoke_grads(e_log)

    lam = symbols("lam", real=True, positive=True)
    Omega = solid_angle_expr(e0, e1)
    # Pure log barrier on Ω: −λ log(Ω+ε). Repels Ω→0 (tol is activation only in C++).
    e_barrier = -lam * sp.log(Omega + eps)
    g_bar_e0 = grad_wrt(e0, e_barrier)
    g_bar_e1 = grad_wrt(e1, e_barrier)

    def fn_energy(name: str, e_expr: sp.Expr) -> str:
        cse_lines, red = _cse_body([e_expr], printer)
        return f"""
inline real {name}(const vec3 &e0, const vec3 &e1, const vec3 &d,
                   real Omega_rest, real tau_rest, real w, real w_tau,
                   real lam_barrier, real eps) {{
{cse_lines}
  return {red[0]};
}}
"""

    def fn_grad_e(name: str) -> str:
        cse_lines, red = _cse_body(list(g_e0) + list(g_e1), printer)
        return f"""
inline void {name}(const vec3 &e0, const vec3 &e1, const vec3 &d,
                   real Omega_rest, real tau_rest, real w, real w_tau,
                   real lam_barrier, real eps, vec3 &grad_e0, vec3 &grad_e1) {{
{cse_lines}
  grad_e0 = vec3({red[0]}, {red[1]}, {red[2]});
  grad_e1 = vec3({red[3]}, {red[4]}, {red[5]});
}}
"""

    def fn_grad_world(name: str) -> str:
        cse_lines, red = _cse_body(
            list(world["c"]) + list(world["p0"]) + list(world["p1"]),
            printer,
        )
        return f"""
inline void {name}(const vec3 &c, const vec3 &x0, const vec3 &x1, const vec3 &d,
                   real Omega_rest, real tau_rest, real w, real w_tau,
                   real lam_barrier, real eps, vec3 &grad_c, vec3 &grad_x0,
                   vec3 &grad_x1) {{
{cse_lines}
  grad_c = vec3({red[0]}, {red[1]}, {red[2]});
  grad_x0 = vec3({red[3]}, {red[4]}, {red[5]});
  grad_x1 = vec3({red[6]}, {red[7]}, {red[8]});
}}
"""

    def fn_barrier_only(name_e: str, name_g: str) -> str:
        cse_e, red_e = _cse_body([e_barrier], printer)
        cse_g, red_g = _cse_body(list(g_bar_e0) + list(g_bar_e1), printer)
        return f"""
/// E = −λ log(Ω+ε) — log barrier repelling Ω→0.
inline real {name_e}(const vec3 &e0, const vec3 &e1, real lam, real eps) {{
{cse_e}
  return {red_e[0]};
}}

inline void {name_g}(const vec3 &e0, const vec3 &e1, real lam, real eps,
                     vec3 &grad_e0, vec3 &grad_e1) {{
{cse_g}
  grad_e0 = vec3({red_g[0]}, {red_g[1]}, {red_g[2]});
  grad_e1 = vec3({red_g[3]}, {red_g[4]}, {red_g[5]});
}}
"""

    header = """// AUTO-GENERATED by scripts/sympy/slab_cap_solid_angle.py — do not edit
#pragma once

#include <cmath>
#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace block {
namespace cap_generated {

/// Smooth hub solid angle Ω for cap wedge (spherical lune between e0, e1).
inline real cap_wedge_solid_angle(const vec3 &e0, const vec3 &e1) {
  const vec3 cr = e0.cross(e1);
  const real num = cr.norm();
  const real den = e0.dot(e1) + e0.norm() * e1.norm();
  return 2.0 * std::atan2(num, den);
}

/// Wedge trivector coefficient τ = (e0×e1)·d.
inline real cap_wedge_tau(const vec3 &e0, const vec3 &e1, const vec3 &d) {
  return e0.cross(e1).dot(d);
}

enum class cap_barrier_kind : int {
  none = 0,
  log = 1,
  inv = 2,
  inv_match = 3,
};

inline real cap_barrier_term(cap_barrier_kind kind, real omega, real Omega_rest,
                             real eps) {
  switch (kind) {
  case cap_barrier_kind::none:
    return 0.0;
  case cap_barrier_kind::log:
    return -std::log(omega + eps);
  case cap_barrier_kind::inv:
    return 1.0 / (omega + eps);
  case cap_barrier_kind::inv_match: {
    const real a = 1.0 / (omega + eps);
    const real b = 1.0 / (Omega_rest + eps);
    const real r = a - b;
    return r * r;
  }
  }
  return 0.0;
}

"""

    footer = """
} // namespace cap_generated
} // namespace block
} // namespace hepworth
} // namespace gaudi
"""

    return (
        header
        + fn_barrier_only("cap_log_barrier_energy", "cap_log_barrier_grad")
        + fn_energy("cap_wedge_fairness_energy_log", e_log)
        + fn_energy("cap_wedge_fairness_energy_inv", e_inv)
        + fn_energy("cap_wedge_fairness_energy_inv_match", e_inv_m)
        + fn_grad_e("cap_wedge_fairness_grad_log")
        + fn_grad_world("cap_wedge_fairness_grad_world_log")
        + footer
    )


def fd_check_energy(energy: sp.Expr, vars_: list, subs_base: dict, eps_fd: float = 1e-6) -> float:
    import numpy as np

    err = 0.0
    grad = Matrix([sp.diff(energy, v) for v in vars_])
    for v in vars_:
        s_plus = dict(subs_base)
        s_minus = dict(subs_base)
        s_plus[v] = float(subs_base[v]) + eps_fd
        s_minus[v] = float(subs_base[v]) - eps_fd
        fd = (float(energy.subs(s_plus)) - float(energy.subs(s_minus))) / (2 * eps_fd)
        sym = float(grad[vars_.index(v)].subs(subs_base))
        err = max(err, abs(fd - sym) / max(abs(sym), 1e-12))
    return err


def run_fd_checks() -> None:
    e = combined_energy(include_tau=True, barrier_key="log")
    base = {
        e0x: 0.1,
        e0y: 0.05,
        e0z: 0.12,
        e1x: -0.08,
        e1y: 0.11,
        e1z: 0.09,
        dx: 0.0,
        dy: 0.0,
        dz: 0.2,
        Omega_rest: 0.15,
        tau_rest: 1e-4,
        w: 2.0,
        w_tau: 5.0,
        lam_barrier: 0.1,
        eps: 1e-6,
    }
    err = fd_check_energy(e, [e0x, e0y, e0z, e1x, e1y, e1z, dx, dy, dz], base)
    print(f"FD vs sympy grad (log barrier) rel err max: {err:.3e}")

    lam = symbols("lam", real=True, positive=True)
    e_bar = -lam * sp.log(solid_angle_expr(e0, e1) + eps)
    base_bar = {e0x: 0.01, e0y: 0.005, e0z: 0.012, e1x: -0.008, e1y: 0.011, e1z: 0.009,
                lam: 1.0, eps: 1e-6}
    err_bar = fd_check_energy(e_bar, [e0x, e0y, e0z, e1x, e1y, e1z], base_bar)
    print(f"FD vs sympy grad (log barrier on Ω) rel err max: {err_bar:.3e}")
    if err > 1e-5 or err_bar > 1e-5:
        raise SystemExit(1)

    omega = solid_angle_expr(e0, e1)
    print("Omega =", sp.simplify(omega))
    print("E_log =", combined_energy(barrier_key="log"))
    print("E_inv =", combined_energy(barrier_key="inv"))
    print("E_inv_match =", combined_energy(barrier_key="inv_match"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fd-check", action="store_true")
    parser.add_argument("--emit", action="store_true")
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="output header (default scripts/sympy/generated/...)",
    )
    args = parser.parse_args()

    if args.fd_check or not args.emit:
        run_fd_checks()

    if args.emit:
        from pathlib import Path

        body = emit_cap_solid_angle_generated()
        out = Path(args.out) if args.out else Path(__file__).resolve().parent / "generated" / "slab_cap_solid_angle_generated.hpp"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(body)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
