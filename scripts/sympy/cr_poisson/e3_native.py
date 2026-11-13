"""E3: the contact net on vertex quaternion frames with native rows (MATH.md 19.23).

Unknowns per vertex: the position `x` and a quaternion `Y`. `Y` gives the
vertex frame `(T, B, N) = Y (e1, e2, e3) conj(Y) / |Y|^2`; `|Y|^2` is a free
scale (the conformal factor; the kite rows use it). No edge store.

Coupling (store-free contact net): per directed spoke `(v, j)` the frames'
normal curvature `k_vj = -(n_j - n_v) . t / |dx_t|` (Rodrigues, `t` the spoke
tangent) and the tangency `n_v . d / |d| = k_vj |d| / 2 + c0`, `c0` the frame
start's offset (increment form), so the positions' curvature follows the
frames'. The vertex tensor `S_v` (in `v`'s own frame) is the least-squares fit
of `k = t^T S t` and `tau = t^T S (N x t)` over `v`'s spokes; rows that compare
vertices read `S` along frozen world directions projected into each tangent
plane, or carry `z = (s11 - s22) / 2 + i s12` by the frames' relative turn,
so the turn of `Y` about `N` (the gauge) drops out.

Native rows (each optional, weight 0 off):
  memory     AMAP-style shape memory: `z_v` held at its frame-start value
             turned by the ring's best common turn; `H_v` held too
  tension    `H_v` to a per-frame constant (constant, frame_mean, volume)
  minimizer  across-line Dirichlet of `k1`, `k2` (frozen principal frames)
  dupin      along-line Dirichlet of `k1`, `k2`
  growth     normal curvature along a direction times `1 + x`, across `1 + y`
  det        each vertex's ring sum of face solid angles (Gauss image area)
  gradtl     Dirichlet on `z` with the frames' relative turn as transport
  twist      geodesic torsion along a winding direction: frame start plus `tau0`
             per frame (a uniform twist increment)
  kite       `|Y_v|^2 = e^{u_v}` and spin-scale edges `|d| = l0 |Y_i| |Y_j|`
  pgrowth    position growth: edge lengths `l0 (1 + g (d0 . a)^2)`

Run: `python e3_native.py [row=weight ...] [frames=N]`; no rows runs every
row alone for a few frames on the torus and prints the shape measures.
"""

import sys
import time

import numpy as np
import scipy.sparse as sp

import cr_ops as ops
import dupin_poisson as dp
import hz_edge_solvers as hes
from quat_store import quat_from_normal


# ------------------------------------------------------------- frames


def qframe(Y):
    """`T, B, N` rows of the normalized quaternions `(w, x, y, z)`."""
    q = Y / np.linalg.norm(Y, axis=1)[:, None]
    w, x, y, z = q.T
    T = np.column_stack([1 - 2 * (y * y + z * z), 2 * (x * y + w * z), 2 * (x * z - w * y)])
    B = np.column_stack([2 * (x * y - w * z), 1 - 2 * (x * x + z * z), 2 * (y * z + w * x)])
    N = np.column_stack([2 * (x * z + w * y), 2 * (y * z - w * x), 1 - 2 * (x * x + y * y)])
    return T, B, N


def dot(a, b):
    return np.einsum("ij,ij->i", a, b)


def unit(v):
    return v / np.linalg.norm(v, axis=1)[:, None]


def shortest_arc_apply(a, b, v):
    """Rotate `v` by the shortest arc taking unit `a` to unit `b` (rows)."""
    c = np.cross(a, b)
    d = dot(a, b)[:, None]
    return v * d + np.cross(c, v) + c * dot(c, v)[:, None] / (1 + d)


# ------------------------------------------------------------- mesh data


class Topo:
    def __init__(self, m):
        self.m = m
        ring = dp.rings(m)
        sv, sj = [], []
        for v in range(m.n):
            for j in ring[v]:
                sv.append(v)
                sj.append(j)
        self.sv, self.sj = np.array(sv), np.array(sj)
        self.E = m.E
        self.F = m.F


def vertex_areas(m, X):
    a, b, c = X[m.F[:, 0]], X[m.F[:, 1]], X[m.F[:, 2]]
    A = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    out = np.zeros(m.n)
    for k in range(3):
        np.add.at(out, m.F[:, k], A / 3)
    return out


def volume(m, X):
    a, b, c = X[m.F[:, 0]], X[m.F[:, 1]], X[m.F[:, 2]]
    return dot(a, np.cross(b, c)).sum() / 6


def spoke_curvature(tp, X, N):
    """Per directed spoke: unit tangent `t`, `k`, `tau` from the frames."""
    v, j = tp.sv, tp.sj
    dx, dn = X[j] - X[v], N[j] - N[v]
    nb = unit(N[v] + N[j])
    dxt = dx - dot(dx, nb)[:, None] * nb
    L = np.linalg.norm(dxt, axis=1)
    t = dxt / L[:, None]
    b = np.cross(nb, t)
    return t, -dot(dn, t) / L, -dot(dn, b) / L


def vertex_tensors(tp, X, T, B, N):
    """`(s11, s12, s22)` per vertex in its own `(T, B)` frame."""
    n = len(X)
    t, k, tau = spoke_curvature(tp, X, N)
    v = tp.sv
    a, b = dot(t, T[v]), dot(t, B[v])
    s = np.sqrt(a * a + b * b)
    a, b = a / s, b / s
    r1 = np.column_stack([a * a, 2 * a * b, b * b])
    r2 = np.column_stack([-a * b, a * a - b * b, a * b])
    M = np.zeros((n, 3, 3))
    rhs = np.zeros((n, 3))
    for r, y in ((r1, k), (r2, tau)):
        np.add.at(M, v, r[:, :, None] * r[:, None, :])
        np.add.at(rhs, v, r * y[:, None])
    M += 1e-12 * np.eye(3)
    return np.linalg.solve(M, rhs[:, :, None])[:, :, 0]


def along(s, T, B, w):
    """Normal curvature and geodesic torsion of `s` along world `w` projected."""
    a, b = dot(w, T), dot(w, B)
    r = np.sqrt(a * a + b * b) + 1e-300
    a, b = a / r, b / r
    k = a * a * s[:, 0] + 2 * a * b * s[:, 1] + b * b * s[:, 2]
    tau = -a * b * s[:, 0] + (a * a - b * b) * s[:, 1] + a * b * s[:, 2]
    return k, tau


def zH(s):
    return 0.5 * (s[:, 0] - s[:, 2]) + 1j * s[:, 1], 0.5 * (s[:, 0] + s[:, 2])


def relative_turn(Ti, Bi, Ni, Tj, Nj):
    """Angle of `T_j`, carried onto `i`'s plane, in `i`'s frame."""
    Tc = shortest_arc_apply(Nj, Ni, Tj)
    return np.arctan2(dot(Tc, Bi), dot(Tc, Ti))


def face_solid_angles(m, N):
    a, b, c = N[m.F[:, 0]], N[m.F[:, 1]], N[m.F[:, 2]]
    num = dot(a, np.cross(b, c))
    den = 1 + dot(a, b) + dot(b, c) + dot(c, a)
    return 2 * np.arctan2(num, den)


def ring_gauss(m, N):
    Om = face_solid_angles(m, N)
    out = np.zeros(m.n)
    for k in range(3):
        np.add.at(out, m.F[:, k], Om)
    return out


# ------------------------------------------------------------- the solver


DEFAULT = dict(
    coupling=10.0, inertia=0.02, memory_len=0.05, unit=1.0, gauge=1e-3,
    memory=0.0, memory_mean=True, tension=0.0, tension_mode="frame_mean", tension_h0=0.0, tension_kv=1.0,
    minimizer=0.0, dupin=0.0, growth=0.0, growth_x=0.3, growth_y=0.0, growth_angle=0.6,
    det=0.0, gradtl=0.0, twist=0.0, tau0=0.5, twist_angle=0.6,
    kite=0.0, kite_alpha=0.1, kite_max=0.1, pgrowth=0.0, pgrowth_g=0.1,
)


class E3:
    def __init__(self, m, X, cfg):
        self.m, self.tp = m, Topo(m)
        self.cfg = dict(DEFAULT, **cfg)
        self.X = X.copy()
        self.Y = quat_from_normal(dp.vertex_normals(m, X))
        self.n = m.n
        self.nvar = 7 * m.n
        self.V_ref = volume(m, X)
        self.pattern = None
        self.colors = None

    # -- frame start
    def start(self):
        c, m, tp = self.cfg, self.m, self.tp
        X, Y = self.X, self.Y
        self.X0, self.Y0 = X.copy(), Y.copy()
        T, B, N = qframe(Y)
        self.T0, self.B0, self.N0 = T, B, N
        self.A = vertex_areas(m, X)
        self.Am = self.A.mean()
        self.hl = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1).mean()
        self.l0 = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
        self.Ynorm0 = np.linalg.norm(Y, axis=1)
        s = vertex_tensors(tp, X, T, B, N)
        self.s0 = s
        self.z0, self.H0 = zH(s)
        # Principal directions in world, tight (larger |k|) first.
        th = 0.5 * np.arctan2(2 * s[:, 1], s[:, 0] - s[:, 2])
        e1 = np.cos(th)[:, None] * T + np.sin(th)[:, None] * B
        e2 = np.cross(N, e1)
        k1, _ = along(s, T, B, e1)
        k2, _ = along(s, T, B, e2)
        swap = np.abs(k2) > np.abs(k1)
        e1[swap], e2[swap] = e2[swap].copy(), -e1[swap].copy()
        self.e1, self.e2 = e1, e2
        # Coupling offsets.
        t, k, tau = spoke_curvature(tp, X, N)
        d = X[tp.sj] - X[tp.sv]
        L = np.linalg.norm(d, axis=1)
        self.c0 = dot(N[tp.sv], d) / L - 0.5 * k * L
        # Memory ring transport at the frame start.
        phi0 = relative_turn(T[tp.sv], B[tp.sv], N[tp.sv], T[tp.sj], N[tp.sj])
        self.zring0 = self.z0[tp.sj] * np.exp(2j * phi0)
        # Tension target.
        Hm = (self.A * self.H0).sum() / self.A.sum()
        mode = c["tension_mode"]
        if mode == "constant":
            self.Ht = c["tension_h0"]
        elif mode == "volume":
            V = volume(m, X)
            # Outward normals: convex `H < 0`; a volume deficit flattens the target.
            self.Ht = Hm + c["tension_kv"] * (self.V_ref - V) / abs(self.V_ref) * abs(Hm)
        else:
            self.Ht = Hm
        # Growth / twist directions.
        def winding(angle):
            return np.cos(angle) * e1 + np.sin(angle) * e2
        self.ga = winding(c["growth_angle"])
        self.gb = np.cross(N, self.ga)
        self.ka0, self.ta0 = along(s, T, B, self.ga)
        self.kb0, _ = along(s, T, B, self.gb)
        self.ta = winding(c["twist_angle"])
        _, self.tw0 = along(s, T, B, self.ta)
        self.Om0 = ring_gauss(m, N)
        # Kite: u from the frame-start mean curvature (convex positive).
        self.u = np.clip(c["kite_alpha"] * (-self.H0), -c["kite_max"], c["kite_max"])
        # Position growth: per edge, mean of the ends' growth axes.
        dhat0 = unit(X[m.E[:, 1]] - X[m.E[:, 0]])
        ae = unit(self.ga[m.E[:, 0]] + np.sign(dot(self.ga[m.E[:, 0]], self.ga[m.E[:, 1]]))[:, None] * self.ga[m.E[:, 1]])
        self.lg = self.l0 * (1 + c["pgrowth_g"] * dot(dhat0, ae) ** 2)

    # -- residual
    def unpack(self, z):
        n = self.n
        return z[:3 * n].reshape(n, 3), z[3 * n:].reshape(n, 4)

    def residual(self, z):
        c, m, tp = self.cfg, self.m, self.tp
        X, Y = self.unpack(z)
        T, B, N = qframe(Y)
        sa = np.sqrt(self.A / self.Am)
        out = []
        add = lambda w, r: out.append(np.sqrt(w) * np.ravel(r)) if w > 0 else None

        add(c["inertia"], (X - self.X0) / self.hl)
        i, j = m.E[:, 0], m.E[:, 1]
        L = np.linalg.norm(X[j] - X[i], axis=1)
        Yn = np.linalg.norm(Y, axis=1)
        if c["kite"] > 0:
            sc = (Yn[i] / self.Ynorm0[i]) * (Yn[j] / self.Ynorm0[j])
            add(c["kite"], (L - self.l0 * sc) / self.hl)
            add(c["kite"], Yn ** 2 / self.Ynorm0 ** 2 - np.exp(self.u))
        else:
            lt = self.lg if c["pgrowth"] > 0 else self.l0
            add(max(c["memory_len"], c["pgrowth"]), (L - lt) / self.hl)
            add(c["unit"], Yn - 1)
        # Coupling.
        t, k, tau = spoke_curvature(tp, X, N)
        d = X[tp.sj] - X[tp.sv]
        Ls = np.linalg.norm(d, axis=1)
        add(c["coupling"], dot(N[tp.sv], d) / Ls - 0.5 * k * Ls - self.c0)
        # Gauge: the turn of T about N against the frame start's, carried.
        Tc = shortest_arc_apply(self.N0, N, self.T0)
        add(c["gauge"], dot(np.cross(Tc, T), N))

        need_s = any(c[r] > 0 for r in ("memory", "tension", "minimizer", "dupin", "growth", "gradtl", "twist"))
        if need_s:
            s = vertex_tensors(tp, X, T, B, N)
            zc, H = zH(s)
            hk = self.hl
            if c["memory"] > 0:
                phi = relative_turn(T[tp.sv], B[tp.sv], N[tp.sv], T[tp.sj], N[tp.sj])
                zr = zc[tp.sj] * np.exp(2j * phi)
                acc = zc * np.conj(self.z0)
                np.add.at(acc, tp.sv, zr * np.conj(self.zring0))
                rot = np.exp(1j * np.angle(acc))
                g = rot * self.z0
                r = (zc - g) * sa
                add(c["memory"], np.column_stack([r.real, r.imag]) * hk)
                if c["memory_mean"]:
                    add(c["memory"], (H - self.H0) * sa * hk)
            if c["tension"] > 0:
                add(c["tension"], (H - self.Ht) * sa * hk)
            if c["minimizer"] > 0 or c["dupin"] > 0:
                v, w = tp.sv, tp.sj
                e1, e2 = self.e1[v], self.e2[v]
                k1v, _ = along(s[v], T[v], B[v], e1)
                k1w, _ = along(s[w], T[w], B[w], e1)
                k2v, _ = along(s[v], T[v], B[v], e2)
                k2w, _ = along(s[w], T[w], B[w], e2)
                dh = unit(X[w] - X[v])
                g1, g2 = dot(dh, e1), dot(dh, e2)
                ss = sa[v] * hk
                if c["minimizer"] > 0:
                    add(c["minimizer"], np.column_stack([g2 * (k1w - k1v), g1 * (k2w - k2v)]) * ss[:, None])
                if c["dupin"] > 0:
                    add(c["dupin"], np.column_stack([g1 * (k1w - k1v), g2 * (k2w - k2v)]) * ss[:, None])
            if c["growth"] > 0:
                ka, ta = along(s, T, B, self.ga)
                kb, _ = along(s, T, B, self.gb)
                r = np.column_stack([ka - (1 + c["growth_x"]) * self.ka0, kb - (1 + c["growth_y"]) * self.kb0, ta - self.ta0])
                add(c["growth"], r * (sa * hk)[:, None])
            if c["gradtl"] > 0:
                ii, jj = m.E[:, 0], m.E[:, 1]
                phi = relative_turn(T[ii], B[ii], N[ii], T[jj], N[jj])
                r = zc[ii] - zc[jj] * np.exp(2j * phi)
                se = np.sqrt(0.5 * (self.A[ii] + self.A[jj]) / self.Am) * hk
                add(c["gradtl"], np.column_stack([r.real, r.imag]) * se[:, None])
            if c["twist"] > 0:
                _, ta = along(s, T, B, self.ta)
                add(c["twist"], (ta - self.tw0 - c["tau0"]) * sa * hk)
        if c["det"] > 0:
            add(c["det"], ring_gauss(m, N) - self.Om0)
        return np.concatenate(out)

    # -- Jacobian by colored central differences
    def build_pattern(self, z):
        n = self.n
        r0 = self.residual(z)
        rows, cols = [], []
        for v in range(n):
            zp = z.copy()
            zp[3 * v:3 * v + 3] += 1e-3
            zp[3 * n + 4 * v:3 * n + 4 * v + 4] += np.array([1e-3, 2e-3, -1e-3, 1.5e-3])
            hit = np.nonzero(np.abs(self.residual(zp) - r0) > 0)[0]
            rows.append(hit)
            cols.append(np.full(len(hit), v))
        rows, cols = np.concatenate(rows), np.concatenate(cols)
        P = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(r0), n))
        C = (P.T @ P).tolil()
        color = -np.ones(n, dtype=int)
        for v in range(n):
            used = set(color[u] for u in C.rows[v] if color[u] >= 0)
            k = 0
            while k in used:
                k += 1
            color[v] = k
        self.pattern, self.colors = P.tocsc(), color
        self.nrows = len(r0)

    def jacobian(self, z):
        n, h = self.n, 1e-7
        P, color = self.pattern, self.colors
        rows, cols, vals = [], [], []
        for col in range(color.max() + 1):
            vs = np.nonzero(color == col)[0]
            sub = P[:, vs]
            for comp in range(7):
                idx = 3 * vs + comp if comp < 3 else 3 * n + 4 * vs + (comp - 3)
                zp, zm = z.copy(), z.copy()
                zp[idx] += h
                zm[idx] -= h
                D = (self.residual(zp) - self.residual(zm)) / (2 * h)
                coo = sub.tocoo()
                rows.append(coo.row)
                vcol = idx[coo.col]
                cols.append(vcol)
                vals.append(D[coo.row])
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.nrows, self.nvar))

    def step(self, iterations=8):
        self.start()
        z0 = np.concatenate([self.X.ravel(), self.Y.ravel()])
        if self.pattern is None:
            self.build_pattern(z0)
        z, st = hes.lm(z0, self.residual, self.jacobian, iterations, nx=3 * self.n)
        self.X, self.Y = (a.copy() for a in self.unpack(z))
        return st


# ------------------------------------------------------------- measures


def measure(m, X, Y, H0=None):
    rho = np.hypot(X[:, 0], X[:, 1])
    R = 0.5 * (rho.min() + rho.max())
    tube = np.hypot(rho - R, X[:, 2]).mean()
    a, b, c = X[m.F[:, 0]], X[m.F[:, 1]], X[m.F[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    T, B, N = qframe(Y)
    tp = Topo(m)
    s = vertex_tensors(tp, X, T, B, N)
    z, H = zH(s)
    rough = np.abs(H[m.E[:, 0]] - H[m.E[:, 1]]).mean()
    q = np.abs(z) / (np.abs(H) + np.abs(z) + 1e-12)
    Yn = np.linalg.norm(Y, axis=1)
    return dict(R=R, tube=tube, area=area, vol=volume(m, X), rough=rough, H=H.mean(), aniso=q.mean(),
                Ymin=Yn.min(), Ymax=Yn.max(), z=z, Hv=H)


def run(cfg, frames=3, res=(24, 12), label=""):
    X, F, _ = ops.torus(1.0, 0.35, *res)
    m = ops.build_mesh(F, len(X))
    e3 = E3(m, X, cfg)
    m0 = measure(m, e3.X, e3.Y)
    line = f"{label:22s} start  tube {m0['tube']:.4f} area {m0['area']:.3f} vol {m0['vol']:.4f} rough {m0['rough']:.4f} aniso {m0['aniso']:.3f}"
    print(line, flush=True)
    t0 = time.time()
    for f in range(frames):
        st = e3.step()
        mm = measure(m, e3.X, e3.Y)
        dz = np.abs(np.angle(mm["z"] * np.conj(m0["z"]))).mean() / 2
        print(f"{'':22s} [{f}] tube {mm['tube']:.4f} area {mm['area']:.3f} vol {mm['vol']:.4f} rough {mm['rough']:.4f} "
              f"aniso {mm['aniso']:.3f} H {mm['H']:+.3f} |Y| {mm['Ymin']:.3f}..{mm['Ymax']:.3f} turn {dz:.3f} "
              f"its {st['its']} ({time.time() - t0:.0f}s)", flush=True)
    return e3


SUITE = [
    ("idle", {}),
    ("memory", dict(memory=1.0, growth=1.0)),
    ("tension constant", dict(tension=1.0, tension_mode="constant", tension_h0=-1.5)),
    ("tension frame_mean", dict(tension=1.0)),
    ("tension volume", dict(tension=1.0, tension_mode="volume")),
    ("minimizer", dict(minimizer=1.0)),
    ("dupin", dict(dupin=1.0)),
    ("growth", dict(growth=1.0)),
    ("det", dict(det=1.0, growth=1.0)),
    ("gradtl", dict(gradtl=1.0)),
    ("twist", dict(twist=1.0)),
    ("twist + pgrowth", dict(twist=1.0, pgrowth=1.0)),
    ("kite", dict(kite=1.0)),
]


if __name__ == "__main__":
    args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    frames = int(args.pop("frames", 3))
    if args:
        cfg = {}
        for k, v in args.items():
            try:
                cfg[k] = float(v)
            except ValueError:
                cfg[k] = v
        run(cfg, frames, label="custom")
    else:
        for name, cfg in SUITE:
            run(cfg, frames, label=name)
