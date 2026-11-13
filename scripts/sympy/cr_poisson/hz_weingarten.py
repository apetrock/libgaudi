"""The vertex (H, z) store under the Weingarten pathway (MATH.md 19.21): the
`shape_site_e2e.py` protocol (the CR edge read-out of the target as the
goal, `cr_goal.goal_gn_loop` with `site="edge_avg"`, 8 frames of 30), with the
goal taken from the edge store or projected onto the vertex store
(`hz_store.VertexStore.project`, frames at the target or at the start), clean
and with 10% noise (against the CR spin-2 fairing of `shape_site_e2e`)."""

import time

import numpy as np

import cr_ops as ops
import face_loops as fl
import hz_store as hz
import shape_dirichlet as sd
from edge_compare import align
from shape_site_e2e import LAMBDAS, cases, fair, recover

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    for name, X0, F, Xt, dist in cases():
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        ell = np.mean(np.linalg.norm(Xt[m.E[:, 1]] - Xt[m.E[:, 0]], axis=1))
        h = sd.edge_average_readings(m, fl.build_stencils(m, ops.cr_derivative(m, Xt)[0]), Xt)
        M = np.repeat(hz.edge_dual_area(m, Xt), 3)
        L = sd.edge_dirichlet(m, Xt, 2)
        vt, v0 = hz.VertexStore(m, Xt), hz.VertexStore(m, X0)
        noise = 0.1 * np.sqrt(np.mean(h ** 2)) * rng.normal(size=h.shape)
        print(f"{name}: unknowns {3 * m.n} (vertex) vs {3 * len(m.E)} (edge)", flush=True)

        def run(label, goal):
            t = time.time()
            X = recover(m, X0, goal, "edge_avg")
            shape = float(np.abs(dist(align(X, Xt))).mean() / size)
            print(f"  {label:34s} goal err {hz.relative(goal, h, m, Xt):.3f}  shape {shape:.4f}  "
                  f"t {time.time() - t:.0f}s", flush=True)
            return shape

        run("edge store, clean", h)
        run("vertex store (target frames), clean", vt.project(h))
        run("vertex store (start frames), clean", v0.project(h))
        run("edge store, noisy", h + noise)
        run("vertex store, noisy", vt.project(h + noise))
        best = None
        for lam in LAMBDAS:
            g = fair(h + noise, M, L, lam * ell ** 2)
            s = run(f"edge store, noisy, CR fair {lam}", g)
            best = s if best is None else min(best, s)
        print(f"  best CR-faired shape {best:.4f}", flush=True)
