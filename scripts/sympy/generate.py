#!/usr/bin/env python3
"""Sympy derive + emit medial search kernels in albers C++ style."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import sympy as sp

from albers_printer import emit_all, write_header
from compose import compose_energy_abstract

INSPECT_KEYS = {
    "W": "W",
    "W_dot": "W_dot",
    "W_ddot": "W_ddot",
    "s2": "s2",
    "E": "E",
    "E_prime": "E_prime",
    "E_pprime": "E_pprime",
    "g_dot": "g_dot",
    "H_dot": "H_dot",
    "g_ddot": "g_ddot",
    "H_ddot": "H_ddot",
}


def inspect_symbols(name: Optional[str], use_cse: bool) -> None:
    energy = compose_energy_abstract()
    keys = [name] if name else list(INSPECT_KEYS.keys())
    for key in keys:
        if key not in energy and key not in INSPECT_KEYS:
            raise SystemExit(f"unknown symbol: {key} (choose from {list(INSPECT_KEYS)})")
        expr = energy.get(key, energy[INSPECT_KEYS.get(key, key)])
        print(f"\n=== {key} ===")
        if use_cse:
            if hasattr(expr, "shape"):
                flat = [expr[i, j] for i in range(expr.shape[0]) for j in range(expr.shape[1])]
                repls, reduced = sp.cse(flat, optimizations="basic")
                for sym, sub in repls:
                    print(f"{sym} = {sub}")
                idx = 0
                for i in range(expr.shape[0]):
                    for j in range(expr.shape[1]):
                        print(f"[{i},{j}] = {reduced[idx]}")
                        idx += 1
            else:
                repls, reduced = sp.cse([expr], optimizations="basic")
                for sym, sub in repls:
                    print(f"{sym} = {sub}")
                print(reduced[0])
        else:
            sp.pprint(expr, use_unicode=True)


def emit_targets(target: str, out_dir: Path) -> None:
    for filename, body in emit_all(target):
        path = out_dir / filename
        write_header(path, body)
        print(f"Wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Derive medial shape energy and emit C++")
    parser.add_argument("--inspect", nargs="?", const="__all__", metavar="SYMBOL")
    parser.add_argument("--cse", action="store_true", help="show CSE-reduced form with --inspect")
    parser.add_argument(
        "--emit",
        choices=["kernels", "quadric", "darboux", "cyclide_smooth", "all"],
        default=None,
        help="write generated headers (default: kernels if no --inspect)",
    )
    parser.add_argument("--out", type=Path, default=ROOT / "generated")
    parser.add_argument(
        "--profile-ray",
        choices=["sphere", "torus"],
        default=None,
        help="dump ray profile table (t, xyz, s2, 1/s2)",
    )
    parser.add_argument("--R", type=float, default=2.5, help="sphere radius for --profile-ray sphere")
    parser.add_argument("--t-max", type=float, default=None, help="max t for --profile-ray")
    parser.add_argument("--steps", type=int, default=25, help="sample count for --profile-ray")
    parser.add_argument("--eps", type=float, default=1e-12, help="denominator mollifier for 1/trace(W²)")
    parser.add_argument(
        "--inspect-march",
        choices=["quadric", "darboux"],
        default=None,
        help="pprint composed medial_energy_graph (E, E', E'')",
    )
    parser.add_argument(
        "--inspect-cyclide-smooth",
        nargs="?",
        const="E_neighbor",
        metavar="KEY",
        help="pprint cyclide_smooth jet energy (E_anchor, E_neighbor, grad_*, hess_*, M_G_i, ...)",
    )
    args = parser.parse_args()

    if args.inspect_cyclide_smooth is not None:
        from cyclide_smooth import INSPECT_KEYS, compose_cyclide_smooth, get_inspect_expr

        bundle = compose_cyclide_smooth()
        key = args.inspect_cyclide_smooth
        if key not in INSPECT_KEYS:
            raise SystemExit(f"unknown key: {key} (choose from {list(INSPECT_KEYS)})")
        expr = get_inspect_expr(bundle, key)
        print(f"\n=== {key} ===")
        if args.cse and hasattr(expr, "shape"):
            flat = [expr[i, j] for i in range(expr.shape[0]) for j in range(expr.shape[1])]
            repls, reduced = sp.cse(flat, optimizations="basic")
            for sym, sub in repls:
                print(f"{sym} = {sub}")
            idx = 0
            for i in range(expr.shape[0]):
                for j in range(expr.shape[1]):
                    print(f"[{i},{j}] = {reduced[idx]}")
                    idx += 1
        elif args.cse:
            repls, reduced = sp.cse([expr], optimizations="basic")
            for sym, sub in repls:
                print(f"{sym} = {sub}")
            print(reduced[0])
        else:
            sp.pprint(expr, use_unicode=True)
        return

    if args.inspect_march is not None:
        from compose import get_graph

        graph = get_graph(args.inspect_march)
        for label, expr in [
            ("s2", graph.s2),
            ("E", graph.E),
            ("E_prime", graph.E_prime),
            ("E_pprime", graph.E_pprime),
        ]:
            print(f"\n=== {label} ===")
            sp.pprint(expr, use_unicode=True)
        return

    if args.profile_ray is not None:
        if args.profile_ray == "sphere":
            from profile_ray import profile_sphere_ray

            profile_sphere_ray(R=args.R, t_max=args.t_max, steps=args.steps, eps=args.eps)
        elif args.profile_ray == "torus":
            from profile_ray import profile_torus_ray

            profile_torus_ray(t_max=args.t_max, steps=args.steps, eps=args.eps)
        return

    if args.inspect is not None:
        name = None if args.inspect == "__all__" else args.inspect
        inspect_symbols(name, args.cse)
        return

    target = args.emit or "kernels"
    emit_targets(target, args.out)


if __name__ == "__main__":
    main()
