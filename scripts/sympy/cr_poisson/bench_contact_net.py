"""Benchmarks for the contact element net (MATH.md 19.11, 19.12).

    python3 bench_contact_net.py cold [key=value ...]   # the shared suite from scratch
    python3 bench_contact_net.py warm [key=value ...]   # 20 warm frames per case

`cold`: iterates to merit / start 1e-3, 1e-6, 1e-9, factorizations,
backsolves and the final shape error, per case. `warm`: the target is reached
in 20 frames (each frame's target interpolates the start and final read-outs),
each frame warm-started from the last; per frame iterates and backsolves to
merit 1e-9 relative to the frame's start, for GN and for the projection solve
with an accelerator.
"""

import sys
import time

import numpy as np

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
from lie_compare import cases
from live_compare import shape_error


def parse(argv):
    cfg = cn.NetConfig()
    solver = "gn"
    only = []
    for a in argv:
        if "=" not in a:
            only.append(a)
            continue
        k, v = a.split("=")
        if k == "solver":
            solver = v
            continue
        cur = getattr(cfg, k)
        setattr(cfg, k, v in ("1", "true", "True") if isinstance(cur, bool) else type(cur)(v))
    return cfg, solver, only


def setup(F, X0):
    m = ops.build_mesh(F, len(X0))
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    return m, pairs


def first_below(h, t):
    i = np.nonzero(np.asarray(h) < t * h[0])[0]
    return int(i[0]) if len(i) else None


def make(solver, m, X0, tgt, cfg, pairs, Xt=None, w_plane=1.0):
    if solver == "gn":
        return cn.ContactNetGN(m, X0, tgt, cfg, pairs)
    if solver in ("sphere", "sphere_point"):
        sph = cn.read_sphere_targets(m, pairs, Xt)
        return cn.ContactNetSphereGN(m, X0, tgt, sph, cfg, pairs, w_plane=0.0 if solver == "sphere_point" else w_plane)
    return cn.ContactNetSolver(m, X0, tgt, cfg, pairs)


def cold(cfg, solver, only):
    print(f"cold, solver {solver}, accel {cfg.accel}, w_memory {cfg.w_memory}, fade {cfg.memory_fade}, "
          f"gn_refresh {cfg.gn_refresh}")
    print(f"  {'case':22s} {'1e-3':>5s} {'1e-6':>5s} {'1e-9':>5s} {'fact':>5s} {'solves':>7s}  shape")
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        m, pairs = setup(F, X0)
        tgt = cn.read_targets(m, pairs, Xt)
        s = make(solver, m, X0, tgt, cfg, pairs, Xt)
        r = s.run(X0.copy(), iterations=cfg.iterations)
        h = [q[1] for q in r.history]
        f = [first_below(h, t) for t in (1e-3, 1e-6, 1e-9)]
        size = np.ptp(Xt, axis=0).max()
        print(f"  {name:22s} " + " ".join(f"{('-' if v is None else v):>5}" for v in f) +
              f" {r.factorizations:5d} {r.backsolves:7d}  {shape_error(r.X, Xt, dist, size):.1e}", flush=True)


def warm(cfg, solver, only, frames=20):
    print(f"warm ({frames} frames), solver {solver}, accel {cfg.accel}, w_memory {cfg.w_memory}, "
          f"gn_refresh {cfg.gn_refresh}")
    print(f"  {'case':22s} {'its/frame':>10s} {'max':>5s} {'fact':>6s} {'solves':>7s}  shape   time")
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        m, pairs = setup(F, X0)
        t0 = cn.read_targets(m, pairs, X0)
        t1 = cn.read_targets(m, pairs, Xt)
        X, N = X0.copy(), dp.vertex_normals(m, X0)
        its, fact, solves = [], 0, 0
        tt = time.time()
        pd = None
        if solver == "pd":  # one factor for all frames: k_ref = the start's curvatures
            pd = cn.ContactNetSolver(m, X0, t0, cfg, pairs, k_ref=t0.k)
            fact += 1
        for f in range(1, frames + 1):
            a = f / frames
            tgt = cn.NetTargets((1 - a) * t0.k + a * t1.k, (1 - a) * t0.tau + a * t1.tau,
                                (1 - a) * t0.c + a * t1.c)
            if pd is not None:
                pd.set_targets(tgt)
                r = pd.run(X, N, iterations=cfg.iterations)
                r.factorizations = 0
            else:
                s = make(solver, m, X, tgt, cfg, pairs)
                r = s.run(X, N, iterations=cfg.iterations)
            h = [q[1] for q in r.history]
            k = first_below(h, 1e-9)
            its.append(k if k is not None else len(h) - 1)
            fact += r.factorizations
            solves += r.backsolves
            X, N = r.X, r.N
        size = np.ptp(Xt, axis=0).max()
        print(f"  {name:22s} {np.mean(its):10.1f} {max(its):5d} {fact:6d} {solves:7d}  "
              f"{shape_error(X, Xt, dist, size):.1e}  {time.time() - tt:.1f}s", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    cfg, solver, only = parse(sys.argv[2:])
    (cold if mode == "cold" else warm)(cfg, solver, only)
