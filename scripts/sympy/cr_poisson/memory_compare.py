"""Where the ARAP memory lives: vertex rotations pushed to edges (polar or
blend), the vertex ARAP energy itself, or per-face rotations. Isometric
cylinder (PD loop and goal GN), sphere cap, torus tube change, cube to sphere."""

import sys
import time

import numpy as np

import base_loop as bl
import closed_targets as ct
import cr_goal as cg
import cr_ops as ops
from run_goal import isometric_cylinder

MEMORIES = (("vertex_tangent", "polar"), ("vertex_tangent", "blend"),
            ("vertex_arap", "polar"), ("face_tangent", "polar"))


def report(label, r, m, Xt, size, t):
    E = np.array(r.energy)
    s = ct.edge_strain(m, r.X, Xt)
    print(f"  {label:32s} its {len(E):3d} E {E[0]:.1e} -> {E[-1]:.1e} "
          f"err/size {bl.rigid_error(r.X, Xt) / size:.2e} "
          f"strain vs target mean {s.mean():.4f} t {time.time() - t:.0f}s",
          flush=True)


def goal_runs(name, m, X0, Xt, rebuilds, w=1e-2, iterations=150):
    s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
    size = np.ptp(Xt, axis=0).max()
    print(f"{name}: rest vs target strain mean "
          f"{ct.edge_strain(m, X0, Xt).mean():.4f}", flush=True)
    for mem, push in MEMORIES:
        for rb in rebuilds:
            t = time.time()
            r = cg.goal_gn_loop(m, X0, X0, s4, cg.GoalGNConfig(
                iterations=iterations, w_arap=w, rebuild=rb, memory=mem,
                push=push))
            report(f"gn {mem}/{push} {rb}", r, m, Xt, size, t)


if __name__ == "__main__":
    which = sys.argv[1:] or ["cyl", "cap", "torus", "cube"]
    if "cyl" in which:
        P, F = ops.strip_param(17, 6, 2.0, 1.0)
        Xf = ops.flat_from_param(P)
        m = ops.build_mesh(F, len(Xf))
        Xi = isometric_cylinder(P, 1.5)
        s4 = bl.readout(m, ops.cr_derivative(m, Xf)[0], Xi)
        print("isometric cylinder, PD loop (frozen), flat start", flush=True)
        for base in ("arap_vertex", "arap_vertex_blend", "arap_face"):
            t = time.time()
            r = cg.tangent_loop(m, Xf, Xf, s4, cg.TangentConfig(
                iterations=150, nr_iterations=3, base=base, rebuild="never"))
            report(f"pd {base}", r, m, Xi, 3.0, t)
        goal_runs("isometric cylinder, goal GN, flat start", m, Xf, Xi, ("none",))
    if "cap" in which:
        P, F = ops.strip_param(13, 9, 1.6, 1.2)
        Xf = ops.flat_from_param(P)
        m = ops.build_mesh(F, len(Xf))
        rr = np.linalg.norm(P, axis=1)
        d = P / np.maximum(rr, 1e-12)[:, None]
        a = rr / 2.0
        Xs = np.stack([2 * np.sin(a) * d[:, 0], 2 * (1 - np.cos(a)),
                       2 * np.sin(a) * d[:, 1]], 1)
        goal_runs("sphere cap r 2", m, Xf, Xs, ("none", "D+G0"))
    if "torus" in which:
        X0, F, _ = ops.torus(2.0, 0.7, 24, 12)
        m = ops.build_mesh(F, len(X0))
        Xt, _, _ = ops.torus(2.0, 0.9, 24, 12)
        goal_runs("torus (2, 0.7) -> (2, 0.9)", m, X0, Xt, ("none", "D+G0"))
    if "cube" in which:
        X0, F = ct.cube(6)
        m = ops.build_mesh(F, len(X0))
        R = np.sqrt(24.0 / (4 * np.pi))
        goal_runs("cube -> sphere", m, X0, ct.radial_ellipsoid(X0, (R, R, R)),
                  ("none", "D+G0"))
