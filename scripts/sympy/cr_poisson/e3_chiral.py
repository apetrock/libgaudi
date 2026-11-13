"""Chiral curvature drive on the E3 net: Armon's incompatible strip (MATH.md 19.22).

Armon, Efrati, Kupferman, Sharon 2011 (`reference/armon2011-chiral-seed-pods.pdf`):
a strip with a flat reference metric and reference curvature
`b = k0 (e+ e+^T - e- e-^T)`, `e+-` at 45 degrees to its axis, is a twisted
helicoid (straight centerline, `H ~ 0`) when narrow and a cylindrical helix
(`K ~ 0`) when wide, sharply at `w sqrt(k0 / t) ~ 2.5`.

Here: `E3` with
  rest    edge lengths held at frame 0 (`w_rest`), the flat metric;
  chiral  per directed spoke the target normal curvature and twist of `b`
          (`h_ref` + saddle, in the vertex frame, the axis the carried frame
          tangent `T` seeded along the strip), in chord form
          `(-dn . t - k_t L) / L0`, `(-dn . b - tau_t L) / L0` as the C++
          spoke targets.
The bending modulus is the target weight, the stretching one `w_rest / h^2`,
so `t_eff = h sqrt(w_chiral / w_rest)`.

Run: `python e3_chiral.py [sweep] [k0=..] [frames=..]`.
"""

import sys
import time

import numpy as np

import cr_ops as ops
import e3_native as e3


def quat_from_frame(T, B, N):
    """`(w, x, y, z)` rows of the rotations with columns `T, B, N`."""
    R = np.stack([T, B, N], axis=2)
    n = len(T)
    q = np.zeros((n, 4))
    for i in range(n):
        m = R[i]
        tr = np.trace(m)
        if tr > 0:
            s = 2 * np.sqrt(tr + 1)
            q[i] = [0.25 * s, (m[2, 1] - m[1, 2]) / s, (m[0, 2] - m[2, 0]) / s, (m[1, 0] - m[0, 1]) / s]
        else:
            k = int(np.argmax(np.diag(m)))
            a, b = (k + 1) % 3, (k + 2) % 3
            s = 2 * np.sqrt(1 + m[k, k] - m[a, a] - m[b, b])
            v = np.zeros(3)
            v[k] = 0.25 * s
            v[a] = (m[a, k] + m[k, a]) / s
            v[b] = (m[b, k] + m[k, b]) / s
            q[i] = [(m[b, a] - m[a, b]) / s, *v]
    return q


class Chiral(e3.E3):
    def __init__(self, m, X, cfg, axis=np.array([1.0, 0.0, 0.0])):
        cfg = dict(cfg)
        self.k0 = cfg.pop("k0", 1.0)
        self.angle = cfg.pop("angle", np.pi / 4)
        self.w_chiral = cfg.pop("chiral", 1.0)
        self.w_rest = cfg.pop("rest", 1.0)
        self.mode = cfg.pop("mode", "held")
        super().__init__(m, X, dict(memory_len=0.0, **cfg))
        N = e3.unit(np.asarray(__import__("dupin_poisson").vertex_normals(m, X)))
        T = e3.unit(axis[None, :] - e3.dot(np.broadcast_to(axis, N.shape), N)[:, None] * N)
        self.Y = quat_from_frame(T, np.cross(N, T), N)
        self.l00 = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)

    def start(self):
        super().start()
        tp = self.tp
        t, _, _ = e3.spoke_curvature(tp, self.X0, self.N0)
        v = tp.sv
        a, b = e3.dot(t, self.T0[v]), e3.dot(t, self.B0[v])
        r = np.hypot(a, b)
        a, b = a / r, b / r
        # Saddle in the frame: k0 (e+ e+^T - e- e-^T), e+- at +-angle from T.
        S = np.zeros((len(v), 3))
        for sgn, ang in ((1, self.angle), (-1, -self.angle)):
            c, s = np.cos(ang), np.sin(ang)
            S += sgn * self.k0 * np.column_stack([c * c, c * s, s * s])
        if self.mode == "rate":
            S = S + self.s0[v]
        self.kt = a * a * S[:, 0] + 2 * a * b * S[:, 1] + b * b * S[:, 2]
        self.taut = -a * b * S[:, 0] + (a * a - b * b) * S[:, 1] + a * b * S[:, 2]
        dx = self.X0[tp.sj] - self.X0[v]
        nb = e3.unit(self.N0[v] + self.N0[tp.sj])
        self.Ls0 = np.linalg.norm(dx - e3.dot(dx, nb)[:, None] * nb, axis=1)

    def residual(self, z):
        base = super().residual(z)
        X, Y = self.unpack(z)
        T, B, N = e3.qframe(Y)
        m, tp = self.m, self.tp
        out = [base]
        L = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
        out.append(np.sqrt(self.w_rest) * (L - self.l00) / self.hl)
        v, j = tp.sv, tp.sj
        dx, dn = X[j] - X[v], N[j] - N[v]
        nb = e3.unit(N[v] + N[j])
        dxt = dx - e3.dot(dx, nb)[:, None] * nb
        Lt = np.linalg.norm(dxt, axis=1)
        t = dxt / Lt[:, None]
        bb = np.cross(nb, t)
        s = np.sqrt(self.w_chiral) * self.hl / self.Ls0
        out.append(s * (-e3.dot(dn, t) - self.kt * Lt))
        out.append(s * (-e3.dot(dn, bb) - self.taut * Lt))
        return np.concatenate(out)


# ------------------------------------------------------------- read-outs


def strip(length, width, h):
    ns = int(round(length / h)) + 1
    nz = max(3, int(round(width / h)) + 1) | 1  # odd: a centerline row
    P, F = ops.strip_param(ns, nz, length, width)
    X = ops.flat_from_param(P)
    return X, ops.build_mesh(F, len(X)), ns, nz


def readout(c, ns, nz):
    X = c.X
    mid = X.reshape(ns, nz, 3)[:, nz // 2]
    seg = np.diff(mid, axis=0)
    arc = np.linalg.norm(seg, axis=1).sum()
    chord = np.linalg.norm(mid[-1] - mid[0])
    tng = e3.unit(seg)
    turn = np.arccos(np.clip(e3.dot(tng[1:], tng[:-1]), -1, 1)).sum()
    T, B, N = e3.qframe(c.Y)
    s = e3.vertex_tensors(c.tp, X, T, B, N)
    z, H = e3.zH(s)
    K = H * H - np.abs(z) ** 2
    grid = np.arange(len(X)).reshape(ns, nz)
    inner = grid[2:-2, 1:-1].ravel() if nz > 3 else grid[2:-2, 1].ravel()
    # Twist of the normal about the centerline, per arc length.
    Nm = N.reshape(ns, nz, 3)[:, nz // 2]
    tw = np.einsum("ij,ij->i", np.cross(Nm[:-1], Nm[1:]), tng).sum() / arc
    L = np.linalg.norm(X[c.m.E[:, 1]] - X[c.m.E[:, 0]], axis=1)
    return dict(chord=chord / arc, kc=turn / arc, H=np.abs(H[inner]).mean(), K=np.abs(K[inner]).mean(),
                tw=tw, strain=np.abs(L / c.l00 - 1).mean())


def run(width, k0=2.0, frames=6, h=0.1, length=4.0, cfg=None, verbose=False, seed=0):
    X, m, ns, nz = strip(length, width, h)
    rng = np.random.default_rng(seed)
    X = X + 1e-3 * h * rng.standard_normal(X.shape)
    c = Chiral(m, X, dict(k0=k0, inertia=1e-3, **(cfg or {})))
    t0 = time.time()
    for f in range(frames):
        st = c.step(iterations=12)
        if verbose:
            r = readout(c, ns, nz)
            print(f"   [{f}] chord {r['chord']:.3f} kc {r['kc']:.3f} |H| {r['H']:.3f} |K| {r['K']:.3f} "
                  f"tw {r['tw']:+.3f} its {st['its']} ({time.time() - t0:.0f}s)", flush=True)
    return readout(c, ns, nz), m.n


def sweep(k0s=(0.5, 1.0), widths=(0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6), frames=6, cfg=None):
    h = 0.1
    teff = h * np.sqrt((cfg or {}).get("chiral", 1.0) / (cfg or {}).get("rest", 1.0))
    for k0 in k0s:
        print(f"k0 {k0}: t_eff {teff:.3f}, Armon w_c ~ {2.5 * np.sqrt(teff / k0):.2f}", flush=True)
        for w in widths:
            t = time.time()
            r, n = run(w, k0, frames, h, cfg=cfg)
            wt = w * np.sqrt(k0 / teff)
            print(f"  w {w:.2f} (w~ {wt:4.1f}, {n} v)  chord {r['chord']:.3f}  kc {r['kc']:.3f}  "
                  f"|H|/k0 {r['H'] / k0:.3f}  |K|/k0^2 {r['K'] / k0 ** 2:.3f}  tw/k0 {r['tw'] / k0:+.3f}  strain {r['strain']:.4f}  "
                  f"({time.time() - t:.0f}s)", flush=True)


if __name__ == "__main__":
    args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    frames = int(args.pop("frames", 6))
    k0 = float(args.pop("k0", 2.0))
    width = float(args.pop("width", 0.5))
    cfg = {k: (float(v) if v.replace(".", "").replace("-", "").replace("e", "").isdigit() else v)
           for k, v in args.items()}
    if "sweep" in sys.argv[1:]:
        sweep(frames=frames, cfg=cfg)
    else:
        r, n = run(width, k0, frames, cfg=cfg, verbose=True)
        print(r)
