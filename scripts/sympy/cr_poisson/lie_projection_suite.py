"""Suite for `lie_shape_projection` (MATH.md 19.16): reconstruction on the
closed cases by position route and projection setting, with cost."""

import sys
import time

import numpy as np

import cr_ops as ops
import lie_edge_net as lsn
import lie_edges as le
import lie_shape_projection as lp
from lie_compare import cases
from live_compare import shape_error

VARIANTS = [
    ("raw", dict(w_couple=0.0)),
    ("central", dict(w_couple=1.0, use_cyclide=False)),
    ("cyclide", dict(w_couple=1.0, use_cyclide=True)),
]

if __name__ == "__main__":
    routes = [a.split("=")[1] for a in sys.argv[1:] if a.startswith("--routes=")]
    routes = routes[0].split(",") if routes else ["spoke", "cnet", "mcp", "envelope"]
    only = [a for a in sys.argv[1:] if not a.startswith("-")]
    outer = 8
    for a in sys.argv[1:]:
        if a.startswith("--outer="):
            outer = int(a.split("=")[1])
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        tg = lsn.read_spheres(m, le.hinges(m), Xt)
        print(f"{name}  (start {shape_error(X0, Xt, dist, size):.1e})", flush=True)
        for route in routes:
            for vname, pc in VARIANTS:
                cfg = lp.SolverConfig(route=route, outer=outer, inner=20, projection=lp.ProjectionConfig(**pc))
                t0 = time.time()
                try:
                    X, info = lp.LieShapeSolver(m, X0, tg, cfg).solve(X0.copy())
                    err = shape_error(X, Xt, dist, size) if np.all(np.isfinite(X)) else float("nan")
                except (np.linalg.LinAlgError, ValueError) as ex:
                    print(f"    {route:8s} {vname:8s} failed: {ex}", flush=True)
                    continue
                print(f"    {route:8s} {vname:8s} shape {err:.1e}  factorizations {info['factorizations']:3d}  "
                      f"backsolves {info['backsolves']:5d}  projection {info['t_projection']:.1f}s  "
                      f"positions {info['t_positions']:.1f}s", flush=True)
