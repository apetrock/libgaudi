"""Slab field energy between Darboux states — symbolic E and grad for codegen.

Variational objective (see docs/divergence_shape_field_bending_tdd.md):

  E_bulk = ∫ |∇(ΔD)|² dV + λ_gap ∫ (ΔD)² dV

Discrete: sum over sample points. SymPy differentiates w.r.t. sample coordinates.

Usage:
  python3 scripts/sympy/slab_field_energy.py
  python3 scripts/sympy/slab_field_energy.py --fd-check
"""

from __future__ import annotations

import argparse

import sympy as sp
from sympy import Matrix, symbols

from darboux import C as Csym
from darboux import D_expr
from implicit import gradient, x, y, z

# Rest / live coefficient vectors (fixed during local x-solve).
Cr = symbols("Cr0:14", real=True)
Ck = symbols("Ck0:14", real=True)


def delta_D_at(pos: Matrix | None = None) -> sp.Expr:
    D_r = D_expr().subs({Csym[i]: Cr[i] for i in range(14)})
    D_k = D_expr().subs({Csym[i]: Ck[i] for i in range(14)})
    if pos is not None:
        pt = {x: pos[0], y: pos[1], z: pos[2]}
        D_r = D_r.subs(pt)
        D_k = D_k.subs(pt)
    return sp.simplify(D_r - D_k)


def bulk_energy_density(delta: sp.Expr, lam_gap: sp.Symbol) -> sp.Expr:
    g = gradient(delta)
    return g.dot(g) + lam_gap * delta * delta


def energy_at_sample(lam_gap: sp.Symbol | None = None) -> tuple[sp.Expr, Matrix]:
    pos = Matrix([x, y, z])
    delta = delta_D_at(pos)
    lam = lam_gap if lam_gap is not None else sp.Symbol("lam_gap", real=True, positive=True)
    E = bulk_energy_density(delta, lam)
    grad = Matrix([sp.diff(E, v) for v in (x, y, z)])
    return sp.simplify(E), sp.simplify(grad)


def green_split_integrand(delta: sp.Expr) -> dict[str, sp.Expr]:
    """Boundary flux density ΔD (n·∇ΔD) and bulk correction −ΔD Δ(ΔD)."""
    g = gradient(delta)
    nx, ny, nz = symbols("nx ny nz", real=True)
    n = Matrix([nx, ny, nz])
    flux = delta * g.dot(n)
    lap = sp.diff(g[0], x) + sp.diff(g[1], y) + sp.diff(g[2], z)
    correction = -delta * lap
    return {
        "flux_density": sp.simplify(flux),
        "laplacian_correction": sp.simplify(correction),
        "grad_sq": sp.simplify(g.dot(g)),
    }


def fd_check(E: sp.Expr, grad: Matrix, eps: float = 1e-6) -> float:
    import numpy as np

    lam = sp.Symbol("lam_gap", real=True, positive=True)
    subs = {
        x: 0.1,
        y: -0.2,
        z: 0.05,
        lam: 1.0,
        **{Cr[i]: (0.1 * (i + 1)) for i in range(14)},
        **{Ck[i]: (0.1 * (i + 1) + 0.01 * i) for i in range(14)},
    }
    err = 0.0
    for j, var in enumerate((x, y, z)):
        s_plus = dict(subs)
        s_minus = dict(subs)
        s_plus[var] = float(subs[var]) + eps
        s_minus[var] = float(subs[var]) - eps
        fd = (float(E.subs(s_plus)) - float(E.subs(s_minus))) / (2 * eps)
        sym = float(grad[j].subs(subs))
        err = max(err, abs(fd - sym) / max(abs(sym), 1e-12))
    return err


def _slab_symbol_map() -> dict:
    lam = sp.Symbol("lam_gap", real=True)
    m = {x: "x[0]", y: "x[1]", z: "x[2]", lam: "lam_gap"}
    for i in range(14):
        m[Cr[i]] = f"Qr[{i}]"
        m[Ck[i]] = f"Qk[{i}]"
    return m


def emit_slab_field_energy() -> str:
    """Body only — write_header adds banner/footer."""
    from albers_printer import AlbersPrinter
    from eigen_printer import emit_vec3_function

    E, grad = energy_at_sample()
    printer = AlbersPrinter(symbol_map=_slab_symbol_map())
    args = "const vec14 &Qr, const vec14 &Qk, const vec3 &x, real lam_gap"
    energy_fn = f"""
inline real slab_bulk_energy({args}) {{
  return {printer.doprint(E)};
}}
"""
    grad_fn = emit_vec3_function(
        "slab_bulk_energy_grad", args, [grad[0], grad[1], grad[2]], printer
    )
    return (
        "// Slab bulk energy |∇(D_r-D_k)|² + lam_gap (D_r-D_k)² at sample x.\n"
        + energy_fn
        + grad_fn
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fd-check", action="store_true")
    args = parser.parse_args()

    E, grad = energy_at_sample()
    print("E =", E)
    print("grad E =", grad)

    parts = green_split_integrand(delta_D_at(Matrix([x, y, z])))
    print("Green flux density =", parts["flux_density"])
    print("Green bulk correction =", parts["laplacian_correction"])

    if args.fd_check:
        err = fd_check(E, grad)
        print(f"FD vs sympy grad rel err max: {err:.3e}")
        if err > 1e-5:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
