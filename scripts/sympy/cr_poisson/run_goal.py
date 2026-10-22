"""Strip to cylinder with `cr_goal.tangent_loop`: targets, ARAP advection,
PD schedule (frozen versus rebuilt every `inner` iterates), starts."""

import time

import numpy as np

import base_loop as bl
import cr_goal as cg
import cr_ops as ops


def isometric_cylinder(P, radius):
    """Columns on the circle with chord = flat spacing: an exact isometry of
    the strip (each column band stays planar)."""
    s = np.unique(P[:, 0])
    dx = s[1] - s[0]
    dth = 2 * np.arcsin(dx / (2 * radius))
    idx = np.round((P[:, 0] - s[0]) / dx).astype(int)
    t = (idx - (len(s) - 1) / 2) * dth - np.pi / 2
    return np.stack([radius * np.cos(t), radius * np.sin(t) + radius, P[:, 1]], 1)


if __name__ == "__main__":
    P, F = ops.strip_param(17, 6, 2.0, 1.0)
    Xf = ops.flat_from_param(P)
    m = ops.build_mesh(F, len(Xf))
    D0, _ = ops.cr_derivative(m, Xf)
    targets = (("isometric", isometric_cylinder(P, 1.5)),
               ("arc", ops.cylinder_from_param(P, 1.5)))
    starts = (("flat", Xf), ("inverted", ops.cylinder_from_param(P, -1.5)))
    for tname, Xt in targets:
        s4 = bl.readout(m, D0, Xt)
        for base in ("arap_face", "arap"):
            for rebuild in ("never", "L+G0"):
                for sname, X0 in starts:
                    t = time.time()
                    r = cg.tangent_loop(m, X0, Xf, s4, cg.TangentConfig(
                        iterations=150, nr_iterations=3, base=base,
                        rebuild=rebuild, inner=30))
                    E = np.array(r.energy)
                    print(f"{tname:9s} {base:9s} {rebuild:5s} {sname:8s} "
                          f"E30 {E[29]:.1e} E60 {E[59]:.1e} Eend {E[-1]:.1e} "
                          f"err {bl.rigid_error(r.X, Xt):.2e} "
                          f"rebuilds {r.rebuilds} t {time.time() - t:.0f}s",
                          flush=True)
