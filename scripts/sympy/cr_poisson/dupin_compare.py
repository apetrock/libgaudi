"""Dupin Poisson solver on the shared validation suite (`edge_compare.cases`,
the closed band of `exponent_solver.parity_table`), with the same metrics
(`live_compare.shape_error`, `closed_targets.edge_strain`)."""

import sys
import time

import numpy as np

import closed_targets as ct
import cr_ops as ops
import dupin_poisson as dp
from edge_compare import cases
from exponent_solver import closed_band
from live_compare import shape_error


def all_cases():
    rows = list(cases())
    X1, Fb = closed_band(24, 6, 1.0)
    X13, _ = closed_band(24, 6, 1.3)
    rows.append(("closed band r 1->1.3", X1, Fb, X13, None))
    return rows


def run_case(name, X0, F, Xt, dist, cfg, verbose=False):
    m = ops.build_mesh(F, len(X0))
    size = np.ptp(Xt, axis=0).max()
    t = time.time()
    body = dp.pose_free_targets(m, Xt, mode=cfg.readout, frames=cfg.frames)
    t_read = time.time() - t
    solver = dp.DupinSolver(m, X0, body, cfg)
    # Fixed point: start at the target.
    fx = solver.run(Xt.copy(), iterations=5)
    solver.theta.clear()

    def log(k, X, move, res):
        if verbose and (k % 10 == 0):
            print(f"    it {k:3d} move {move:.1e} land {res:.1e} shape {shape_error(X, Xt, dist, size):.1e}",
                  flush=True)

    r = solver.run(X0.copy(), log=log)
    best_it = int(np.argmin([h[1] for h in r.history]))
    finite = [h[1] for h in r.history if np.isfinite(h[1])]
    r.history.append((0.0, finite[-1] if finite else np.inf))
    print(f"{name:22s} its {len(r.history):3d} (best {best_it})  shape {shape_error(X0, Xt, dist, size):.1e} -> "
          f"{shape_error(r.X, Xt, dist, size):.1e}  (fixed point {shape_error(fx.X, Xt, dist, size):.1e})  "
          f"strain {ct.edge_strain(m, r.X, Xt).mean():.3f}  land {r.history[-1][1]:.1e}  "
          f"t {r.seconds:.1f}s ({r.seconds / len(r.history) * 1e3:.0f} ms/it, read-out {t_read:.1f}s)", flush=True)
    return r


if __name__ == "__main__":
    only = [a for a in sys.argv[1:] if not a.startswith("-")]
    verbose = "-v" in sys.argv
    cfg = dp.SolverConfig()
    for a in sys.argv[1:]:
        if a.startswith("--"):
            k, v = a[2:].split("=")
            cur = getattr(cfg, k)
            setattr(cfg, k, type(cur)(v) if not isinstance(cur, bool) else v in ("1", "true", "True"))
    for name, X0, F, Xt, dist in all_cases():
        if only and not any(o in name for o in only):
            continue
        run_case(name, X0, F, Xt, dist, cfg, verbose)
