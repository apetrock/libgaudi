"""Complex `(H, z)` edge store against the component store, in both joint
solvers (MATH.md 19.21).

The store stays on edges, 3 reals each, in either parametrization:

    components  s = (k11, k12, k22)
    complex     s = (H, Re z, Im z),   k = C s,  C = [[1, 1, 0], [0, 0, 1], [1, -1, 0]]

and the store's rows (Weingarten read-out coupling, store target) measure
the mismatch `dk` in either metric:

    component   |dk|^2 = dk11^2 + dk12^2 + dk22^2 = 2 dH^2 + 2 dRe(z)^2 + dIm(z)^2
    invariant   dk11^2 + 2 dk12^2 + dk22^2 = 2 dH^2 + 2 |dz|^2   (frame free)

Claim checked: the parametrization does not change the solve (same
residuals, a constant invertible change of unknowns per edge); the metric
does, where rows compete (noisy targets), not where they agree (clean).

1. Weingarten: unknowns `[x, s]`; rows `sqrt(w_c A_e) Wm (R_e(x) - C s_e)`
   (`R` the CR edge read-out, `shape_dirichlet.edge_average_readings`),
   `sqrt(w_t A_e) Wm (C s_e - g_e)`, fading length memory; Levenberg-Marquardt.
2. Contact net: `hz_contact_net.JointStore` with the store reparametrized
   (`k = (C s)_0`, `tau = (C s)_1`) and its target rows in either metric.
"""

import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import face_loops as fl
import hz_contact_net as hcn
import hz_store as hz
import lie_edge_net as lsn
import lie_edges as le
import shape_dirichlet as sd
from live_compare import shape_error

C_HZ = np.array([[1.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, -1.0, 0.0]])
C_K = np.eye(3)
W_COMP = np.eye(3)
W_INV = np.diag([1.0, np.sqrt(2.0), 1.0])


def lm(z0, residual, jacobian, iterations=40, lam0=1e-4, on_iter=None, tol=1e-12, nx=None):
    """Levenberg-Marquardt with `diag(J^T J)` damping (as `JointStore.run`)."""
    z = z0.copy()
    r = residual(z)
    cost = 0.5 * r @ r
    lam = lam0
    st = dict(its=0, factorizations=0, t_factor=0.0, fill=0)
    for it in range(iterations):
        if on_iter is not None:
            on_iter(z, cost)
            r = residual(z)
            cost = 0.5 * r @ r
        J = jacobian(z)
        H = (J.T @ J).tocsc()
        g = J.T @ r
        ok = False
        for _ in range(8):
            dg = H.diagonal()
            t = time.time()
            lu = spla.splu((H + sp.diags(lam * dg + 1e-12 * dg.mean())).tocsc())
            st["t_factor"] += time.time() - t
            st["factorizations"] += 1
            st["fill"] = lu.L.nnz + lu.U.nnz
            dz = -lu.solve(g)
            rn = residual(z + dz)
            if np.all(np.isfinite(rn)) and 0.5 * rn @ rn < cost:
                ok = True
                break
            lam *= 4
        if not ok:
            break
        z = z + dz
        r = rn
        cost = 0.5 * r @ r
        lam = max(lam / 3, 1e-12)
        st["its"] = it + 1
        if np.abs(dz[: (nx or len(dz))]).max() < tol:
            break
    return z, st


class JointWeingarten:
    def __init__(self, m, X0, g, C, Wm, w_c=1.0, w_t=10.0, w_m=0.01, fd=1e-7):
        self.m, self.g, self.C, self.Wm = m, g, C, Wm
        self.WC = Wm @ C
        self.w_c, self.w_t, self.w_m, self.fd = w_c, w_t, w_m, fd
        self.st = fl.build_stencils(m, ops.cr_derivative(m, X0)[0])
        A = hz.edge_dual_area(m, X0)
        self.sa = np.sqrt(A / A.mean())
        self.cnt = np.zeros(len(m.E))
        for i in range(3):
            np.add.at(self.cnt, m.fe[:, i], 1)
        self.hl = float(np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean())
        self.n, self.nE = m.n, len(m.E)
        self.memory_scale = 1.0
        self.L = None
        self.cost0 = None

    def unpack(self, z):
        return z[: 3 * self.n].reshape(-1, 3), z[3 * self.n:].reshape(-1, 3)

    def readout(self, X):
        Rf = sd.face_edge_readings(self.st, X[self.st.ids])
        acc = np.zeros((self.nE, 3))
        for i in range(3):
            np.add.at(acc, self.m.fe[:, i], Rf[:, i])
        return acc / self.cnt[:, None]

    def residual(self, z):
        X, S = self.unpack(z)
        K = S @ self.C.T
        E = self.m.E
        rc = np.sqrt(self.w_c) * self.sa[:, None] * ((self.readout(X) - K) @ self.Wm.T)
        rt = np.sqrt(self.w_t) * self.sa[:, None] * ((K - self.g) @ self.Wm.T)
        rm = np.sqrt(self.w_m * self.memory_scale) / self.hl * (np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1) - self.L)
        return np.concatenate([rc.ravel(), rt.ravel(), rm])

    def jacobian(self, z):
        X, S = self.unpack(z)
        m, st, n, nE = self.m, self.st, self.n, self.nE
        rows, cols, vals = [], [], []
        P = X[st.ids]
        sc = np.sqrt(self.w_c)
        for s in range(6):
            live = st.mask[:, s] > 0
            for c in range(3):
                Pp, Pm = P.copy(), P.copy()
                Pp[:, s, c] += self.fd
                Pm[:, s, c] -= self.fd
                dR = (sd.face_edge_readings(st, Pp) - sd.face_edge_readings(st, Pm)) / (2 * self.fd)
                dR = dR @ self.Wm.T
                for i in range(3):
                    e = m.fe[:, i]
                    coef = sc * self.sa[e] / self.cnt[e]
                    for r in range(3):
                        rows.append((3 * e + r)[live])
                        cols.append((3 * st.ids[:, s] + c)[live])
                        vals.append((coef * dR[:, i, r])[live])
        e = np.arange(nE)
        for blk_row0, w in ((0, -sc), (3 * nE, np.sqrt(self.w_t))):
            for r in range(3):
                for q in range(3):
                    if self.WC[r, q] != 0:
                        rows.append(blk_row0 + 3 * e + r)
                        cols.append(3 * n + 3 * e + q)
                        vals.append(w * self.sa * self.WC[r, q])
        E = m.E
        d = X[E[:, 1]] - X[E[:, 0]]
        u = d / np.linalg.norm(d, axis=1)[:, None]
        wm = np.sqrt(self.w_m * self.memory_scale) / self.hl
        r0 = 6 * nE
        for c in range(3):
            rows += [r0 + e, r0 + e]
            cols += [3 * E[:, 1] + c, 3 * E[:, 0] + c]
            vals += [wm * u[:, c], -wm * u[:, c]]
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(6 * nE + nE, 3 * n + 3 * nE))

    def run(self, X0, S0, iterations=40):
        E = self.m.E
        self.L = np.linalg.norm(X0[E[:, 1]] - X0[E[:, 0]], axis=1)
        z0 = np.concatenate([X0.ravel(), S0.ravel()])
        self.cost0 = None

        def on_iter(z, cost):
            X, _ = self.unpack(z)
            self.L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
            if self.cost0 is None:
                self.cost0 = max(cost, 1e-300)
            self.memory_scale = min(1.0, cost / self.cost0)

        z, st = lm(z0, self.residual, self.jacobian, iterations, on_iter=on_iter, nx=3 * self.n)
        return self.unpack(z), st


class JointContact(hcn.JointStore):
    """`JointStore` (edge store) with the store as `C s` and its target rows
    in the metric `Wm`."""

    def __init__(self, m, X0, h, c, C, Wm, **kw):
        super().__init__(m, X0, h, c, None, **kw)
        self.C, self.Wm = C, Wm

    def edge_kt(self, S, e_idx, V):
        K = V[:, :3] @ self.C.T
        return K[:, 0], K[:, 1]

    def blocks(self, X, N, S, L):
        b = super().blocks(X, N, S, L)
        w, _, Ve, starts = b[3]
        b[3] = (w, lambda V: self.hl * ((V @ self.C.T - self.h) @ self.Wm.T), Ve, starts)
        return b


VARIANTS = (("components, component metric", C_K, W_COMP),
            ("complex,    component metric", C_HZ, W_COMP),
            ("components, invariant metric", C_K, W_INV),
            ("complex,    invariant metric", C_HZ, W_INV))


def report(label, X, ref, Xt, dist, size, st, t):
    dx = "" if ref is None else f"  |X - X_ref| {np.abs(X - ref).max():.1e}"
    return (f"    {label}: shape {shape_error(X, Xt, dist, size):.2e}  its {st['its']:2d}  "
            f"facts {st['factorizations']:2d}  fill {st['fill']}  factor {st['t_factor']:.2f}s  "
            f"total {t:.1f}s{dx}")


def run_case(solver, name, X0, F, Xt, dist, rng):
    m = ops.build_mesh(F, len(X0))
    size = np.ptp(Xt, axis=0).max()
    if solver == "weingarten":
        g = sd.edge_average_readings(m, fl.build_stencils(m, ops.cr_derivative(m, Xt)[0]), Xt)
        h0 = sd.edge_average_readings(m, fl.build_stencils(m, ops.cr_derivative(m, X0)[0]), X0)
        extra = {}
    else:
        ring = dp.rings(m)
        pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        tg = cn.read_targets(m, pairs, Xt)
        g = lsn.read_spheres(m, le.hinges(m), Xt).H
        g[:, 0], g[:, 1] = tg.k, tg.tau
        h0 = lsn.read_spheres(m, le.hinges(m), X0).H
        extra = dict(c=tg.c, N0=dp.vertex_normals(m, X0))
    noise = 0.1 * np.sqrt(np.mean(g ** 2)) * rng.normal(size=g.shape)
    print(f"  {name}  (store unknowns {3 * len(m.E)} either way)", flush=True)
    for tlabel, tgt in (("clean", g), ("noisy 10%", g + noise)):
        print(f"   targets {tlabel}", flush=True)
        refs = {}
        for label, C, Wm in VARIANTS:
            S0 = h0 @ np.linalg.inv(C).T
            t = time.time()
            if solver == "weingarten":
                (X, S), st = JointWeingarten(m, X0, tgt, C, Wm).run(X0.copy(), S0)
            else:
                (X, N, S), st = JointContact(m, X0, tgt, extra["c"], C, Wm).run(X0.copy(), extra["N0"].copy(), S0)
            key = "inv" if Wm is W_INV else "comp"
            print(report(label, X, refs.get(key), Xt, dist, size, st, time.time() - t), flush=True)
            refs.setdefault(key, X)


if __name__ == "__main__":
    which = [a for a in sys.argv[1:] if a in ("weingarten", "contact")] or ["weingarten", "contact"]
    only = [a for a in sys.argv[1:] if a not in ("weingarten", "contact")]
    rng = np.random.default_rng(7)
    t0 = time.time()
    if "weingarten" in which:
        from shape_site_e2e import cases as w_cases
        print("1. Weingarten joint solve (positions + edge store)")
        for name, X0, F, Xt, dist in w_cases():
            if only and not any(o in name for o in only):
                continue
            run_case("weingarten", name, X0, F, Xt, dist, rng)
    if "contact" in which:
        from lie_compare import cases as c_cases
        print("2. Contact net joint solve (positions + normals + edge store)")
        for name, X0, F, Xt, dist in c_cases():
            if only and not any(o in name for o in only):
                continue
            run_case("contact", name, X0, F, Xt, dist, rng)
    print(f"({time.time() - t0:.0f}s)")
