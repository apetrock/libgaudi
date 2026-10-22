"""Non-isometric targets on closed meshes with `cr_goal.goal_gn_loop`:
cube to sphere / ellipsoid (radial projection, same connectivity) and torus
radius changes. The target `s_g` is the target mesh's own read-out (its own
CR derivative), so it is exactly realizable, by a mesh with another metric."""

import sys
import time

import numpy as np

import base_loop as bl
import cr_goal as cg
import cr_ops as ops


def cube(k: int):
    """Cube [-1, 1]^3, k x k quads per side, two triangles each, outward."""
    index, X, F = {}, [], []

    def vid(p):
        key = tuple(np.round(p, 9))
        if key not in index:
            index[key] = len(X)
            X.append(np.array(p, dtype=float))
        return index[key]

    g = np.linspace(-1, 1, k + 1)
    for axis in range(3):
        for side in (-1.0, 1.0):
            u, v = [a for a in range(3) if a != axis]
            for i in range(k):
                for j in range(k):
                    q = []
                    for di, dj in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        p = np.zeros(3)
                        p[axis] = side
                        p[u], p[v] = g[i + di], g[j + dj]
                        q.append(vid(p))
                    for tri in ((q[0], q[1], q[2]), (q[0], q[2], q[3])):
                        F.append(tri)
    X, F = np.array(X), np.array(F)
    for f, (a, b, c) in enumerate(F):
        n = np.cross(X[b] - X[a], X[c] - X[a])
        if n @ (X[a] + X[b] + X[c]) < 0:
            F[f] = (a, c, b)
    return X, F


def radial_ellipsoid(X, abc):
    d = X / np.linalg.norm(X, axis=1)[:, None]
    t = 1.0 / np.sqrt(np.sum((d / np.asarray(abc)) ** 2, axis=1))
    return d * t[:, None]


def edge_strain(m, X, Y):
    """`|L_X / L_Y - 1|` per edge: metric mismatch of `X` against `Y`."""
    LX = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
    LY = np.linalg.norm(Y[m.E[:, 1]] - Y[m.E[:, 0]], axis=1)
    return np.abs(LX / LY - 1)


def run_case(name, m, X0, Xt, rebuilds=("none", "D", "D+G0"), w=1e-2,
             iterations=150, inner=30):
    s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
    size = np.ptp(Xt, axis=0).max()
    s0 = edge_strain(m, X0, Xt)
    print(f"{name}: n {m.n} faces {len(m.F)}; rest vs target strain "
          f"max {s0.max():.3f} mean {s0.mean():.3f}", flush=True)
    for rb in rebuilds:
        t = time.time()
        r = cg.goal_gn_loop(m, X0, X0, s4, cg.GoalGNConfig(
            iterations=iterations, inner=inner, w_arap=w, rebuild=rb))
        E = np.array(r.energy)
        s = edge_strain(m, r.X, Xt)
        print(f"  {rb:5s} its {len(E):3d} outer {r.outer} refactors {r.refactors} "
              f"E {E[0]:.1e} -> {E[-1]:.1e} "
              f"err/size {bl.rigid_error(r.X, Xt) / size:.2e} "
              f"strain vs target max {s.max():.3f} mean {s.mean():.4f} "
              f"t {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["torus", "cube"]
    if "torus" in which:
        X0, F, _ = ops.torus(2.0, 0.7, 24, 12)
        m = ops.build_mesh(F, len(X0))
        for R, r in ((2.0, 0.9), (2.5, 0.7), (1.6, 0.7)):
            Xt, _, _ = ops.torus(R, r, 24, 12)
            run_case(f"torus (2, 0.7) -> ({R}, {r})", m, X0, Xt)
    if "cube" in which:
        X0, F = cube(6)
        m = ops.build_mesh(F, len(X0))
        R = np.sqrt(24.0 / (4 * np.pi))
        run_case(f"cube -> sphere r {R:.3f}", m, X0, radial_ellipsoid(X0, (R, R, R)))
        run_case("cube -> ellipsoid (1.8, 1.2, 0.9)", m, X0,
                 radial_ellipsoid(X0, (1.8, 1.2, 0.9)))
