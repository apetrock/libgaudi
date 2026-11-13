"""Edge-native against chained curvature rows on the E3 net (MATH.md 19.21).

The E3 net carries the shape as vertex quaternion frames; curvature rows read
it either from each vertex's fitted tensor (chained: the fit over the ring, so
a row comparing two vertices couples both rings) or from per-spoke Rodrigues
readings (edge-native: a row touches the spoke's ends, plus the few spokes a
line reading averages). No new unknowns either way.

Edge-native line readings: per vertex, `k1(v)` the mean Rodrigues normal
curvature of the two spokes best aligned with `+-e1` (frame-start tight
direction), `k2(v)` with `+-e2`; the spoke choice is frozen at the frame start.
Rows: Dupin `k1(w) - k1(v)` for the `+e1` neighbour (and `k2` along `e2`);
minimizer `k1(w) - k1(v)` for the `+e2` neighbour (and `k2` along `e1`), or the
second difference `k1(w+) + k1(w-) - 2 k1(v)` (`order 2`).

Target rows for the shared suite (`lie_compare.cases()`):
  chained  each vertex's fitted tensor against the target mesh's fitted
           tensor, both as world tensors read along the live frame;
  edge     each directed spoke's `(k, tau)` against the target mesh's.

Run: `python e3_parity.py [torus] [suite] [cost]` (all by default).
"""

import sys
import time

import numpy as np

import cr_ops as ops
import dupin_poisson as dp
import e3_native as e3
from live_compare import shape_error
from quat_store import quat_from_normal


class Parity(e3.E3):
    """`E3` with edge-native Dupin / minimizer rows and target rows."""

    def __init__(self, m, X, cfg, Xt=None):
        cfg = dict(cfg)
        self.edge = dict(dupin=cfg.pop("edge_dupin", 0.0), minimizer=cfg.pop("edge_minimizer", 0.0),
                         order=cfg.pop("edge_order", 1), target=cfg.pop("edge_target", 0.0))
        self.w_target = cfg.pop("tensor_target", 0.0)
        super().__init__(m, X, cfg)
        self.Xt = Xt
        self.tp_list = None

    # -- frame start
    def start(self):
        super().start()
        tp, n = self.tp, self.n
        X, T, B, N = self.X0, self.T0, self.B0, self.N0
        d = e3.unit(X[tp.sj] - X[tp.sv])
        # Per vertex the spokes best aligned with +-e1 and +-e2.
        self.pick = np.full((n, 4), -1, dtype=int)  # +e1, -e1, +e2, -e2 (spoke index)
        best = np.full((n, 4), -2.0)
        for s in range(len(tp.sv)):
            v = tp.sv[s]
            c1, c2 = d[s] @ self.e1[v], d[s] @ self.e2[v]
            for slot, val in ((0, c1), (1, -c1), (2, c2), (3, -c2)):
                if val > best[v, slot]:
                    best[v, slot] = val
                    self.pick[v, slot] = s
        if self.Xt is not None:
            Nt = dp.vertex_normals(self.m, self.Xt)
            Yt = quat_from_normal(Nt)
            Tt, Bt, Nt = e3.qframe(Yt)
            st = e3.vertex_tensors(tp, self.Xt, Tt, Bt, Nt)
            self.Wt = st[:, 0, None, None] * np.einsum("ni,nj->nij", Tt, Tt) + \
                st[:, 1, None, None] * (np.einsum("ni,nj->nij", Tt, Bt) + np.einsum("ni,nj->nij", Bt, Tt)) + \
                st[:, 2, None, None] * np.einsum("ni,nj->nij", Bt, Bt)
            _, kt, taut = e3.spoke_curvature(tp, self.Xt, Nt)
            self.kt, self.taut = kt, taut

    def line_readings(self, k):
        """`k1, k2` per vertex from the picked spokes' normal curvature `k`."""
        p = self.pick
        k1 = 0.5 * (k[p[:, 0]] + k[p[:, 1]])
        k2 = 0.5 * (k[p[:, 2]] + k[p[:, 3]])
        return k1, k2

    def residual(self, z):
        base = super().residual(z)
        X, Y = self.unpack(z)
        T, B, N = e3.qframe(Y)
        sa = np.sqrt(self.A / self.Am)
        hk = self.hl
        out = [base]
        add = lambda w, r: out.append(np.sqrt(w) * np.ravel(r)) if w > 0 else None
        ed = self.edge
        if ed["dupin"] > 0 or ed["minimizer"] > 0 or ed["target"] > 0:
            _, k, tau = e3.spoke_curvature(self.tp, X, N)
            if ed["dupin"] > 0 or ed["minimizer"] > 0:
                k1, k2 = self.line_readings(k)
                p, tp = self.pick, self.tp
                nb = lambda slot: tp.sj[p[:, slot]]
                if ed["dupin"] > 0:
                    add(ed["dupin"], np.column_stack([k1[nb(0)] - k1, k2[nb(2)] - k2]) * (sa * hk)[:, None])
                if ed["minimizer"] > 0:
                    if ed["order"] == 2:
                        r = np.column_stack([k1[nb(2)] + k1[nb(3)] - 2 * k1, k2[nb(0)] + k2[nb(1)] - 2 * k2])
                    else:
                        r = np.column_stack([k1[nb(2)] - k1, k2[nb(0)] - k2])
                    add(ed["minimizer"], r * (sa * hk)[:, None])
            if ed["target"] > 0:
                sv = sa[self.tp.sv] * hk
                add(ed["target"], np.column_stack([(k - self.kt) * sv, (tau - self.taut) * sv]))
        if self.w_target > 0:
            s = e3.vertex_tensors(self.tp, X, T, B, N)
            W = s[:, 0, None, None] * np.einsum("ni,nj->nij", T, T) + \
                s[:, 1, None, None] * (np.einsum("ni,nj->nij", T, B) + np.einsum("ni,nj->nij", B, T)) + \
                s[:, 2, None, None] * np.einsum("ni,nj->nij", B, B)
            D = W - self.Wt
            dg = e3.unit(T + B)
            r = np.column_stack([np.einsum("ni,nij,nj->n", a, D, a) for a in (T, B, dg)])
            add(self.w_target, r * (sa * hk)[:, None])
        return np.concatenate(out)


def torus_mesh(res=(24, 12)):
    X, F, _ = ops.torus(1.0, 0.35, *res)
    return X, ops.build_mesh(F, len(X))


def run_torus():
    print("== torus: fixed point and smoothing (3 frames)")
    rng = np.random.default_rng(5)
    for noise in (0.0, 0.3):
        for label, cfg in (("chained dupin", dict(dupin=1.0)), ("edge dupin", dict(edge_dupin=1.0)),
                           ("chained minimizer", dict(minimizer=1.0)), ("edge minimizer", dict(edge_minimizer=1.0)),
                           ("edge minimizer o2", dict(edge_minimizer=1.0, edge_order=2))):
            X, m = torus_mesh()
            if noise > 0:
                hl = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1).mean()
                X = X + noise * hl * rng.normal(size=(len(X), 1)) * dp.vertex_normals(m, X)
            p = Parity(m, X, cfg)
            m0 = e3.measure(m, p.X, p.Y)
            t0 = time.time()
            for _ in range(3):
                p.step()
            mm = e3.measure(m, p.X, p.Y)
            nnz = p.pattern.nnz if p.pattern is not None else 0
            print(f"  noise {noise:.1f} {label:20s} tube {m0['tube']:.4f} -> {mm['tube']:.4f}  rough {m0['rough']:.4f} -> "
                  f"{mm['rough']:.4f}  pattern nnz {nnz:6d}  {time.time() - t0:.1f}s", flush=True)


def run_tiny():
    print("== a tiny weight leaves a drive alone (kite + twist, 3 frames)")
    for label, extra in (("none", {}), ("chained minimizer 1e-10", dict(minimizer=1e-10)),
                         ("edge minimizer 1e-10", dict(edge_minimizer=1e-10))):
        X, m = torus_mesh()
        p = Parity(m, X, dict(kite=1.0, twist=1.0, tau0=0.1, **extra))
        for _ in range(3):
            p.step()
        mm = e3.measure(m, p.X, p.Y)
        print(f"  {label:26s} tube {mm['tube']:.5f} area {mm['area']:.4f} rough {mm['rough']:.5f}", flush=True)


def run_suite(frames=5):
    from lie_compare import cases
    print("== suite: target rows, chained tensor against edge spoke readings")
    for name, X0, F, Xt, dist in cases():
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        line = f"  {name:22s} start {shape_error(X0, Xt, dist, size):.2e}"
        for label, cfg in (("chained", dict(tensor_target=1.0)), ("edge", dict(edge_target=1.0))):
            p = Parity(m, X0.copy(), dict(cfg, memory_len=0.01), Xt=Xt)
            t0 = time.time()
            for _ in range(frames):
                p.step()
            line += f"  {label} {shape_error(p.X, Xt, dist, size):.2e} ({time.time() - t0:.1f}s)"
        print(line, flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["torus", "tiny", "suite"]
    if "torus" in which:
        run_torus()
    if "tiny" in which:
        run_tiny()
    if "suite" in which:
        run_suite()
