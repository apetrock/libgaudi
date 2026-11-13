"""The exponent curvature solver (MATH.md section 16).

Each face carries an exponent xi_f = [omega_f]x + P_f L_f P_f^T: its
displacement gradient, a rotation omega (3) and an in-plane log-stretch L (3,
in the face frame T, B). Per restart, at the current X, one constant system:

  curvature rows    sqrt(W) ([J_theta C, J_ell E] xi - r),  r = s_g - s(X)
                    J_ff = [J_theta, J_ell]: Stein's read-out in the discrete
                    fundamental forms (cr_ff_jacobian: rebuild x generated
                    stencil Jacobian); C: dtheta_e = a_e . (omega_- - omega_+);
                    E: dlog l_e = mean over its faces of t^T L_f t
  compatibility     xi_f d_e = xi_g d_e per interior edge (the field is the
                    gradient of a displacement), stiff
  priors            kb |omega|^2 (also removes the rigid null space),
                    km |L_cur + dL|^2 (elastic: total stretch against a rest)

Each iterate: fresh residual, xi by one backsolve, F_f = expm(a xi_f) per
face, advected frames sum_f F_f g_f / M_e gathered into positions by Stein's
CR Dirichlet operator D^T M D (constant per restart); line search on
|r|_W^2 + km tau sum A_f |L_cur|^2. LU stands in for LDL^T.

Run: python3 exponent_solver.py           # parity table and plate buckling
"""

from __future__ import annotations

import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_curvature_jacobian as cj
import cr_ff_jacobian as ffj
import cr_loop as cl
import cr_ops as ops
import face_loops as fl

E_BASIS = (np.array([[1.0, 0.0], [0.0, 0.0]]), np.array([[0.0, 1.0], [1.0, 0.0]]),
           np.array([[0.0, 0.0], [0.0, 1.0]]))


def skew(w):
    return np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])


def sym2(v):
    return np.array([[v[0], v[1]], [v[1], v[2]]])


class Forms:
    """Topology of the fundamental-form Jacobian: per face, the global columns
    of its 9 lengths and 3 dihedrals, and the generated read-out pieces."""

    _fns = None

    def __init__(self, m):
        self.m = m
        D0, _ = ops.cr_derivative(m, np.random.default_rng(0).standard_normal((m.n, 3)))
        self.st = fl.build_stencils(m, D0)
        self.edges = fl._interior_edges(m)
        kidx = {e: k for k, (e, *_r) in enumerate(self.edges)}
        index = {(int(a), int(b)): e for e, (a, b) in enumerate(m.E)}
        nF = len(m.F)
        self.lcol = -np.ones((nF, 9), int)
        self.tcol = -np.ones((nF, 3), int)
        self.h = self.st.mask[:, 3:].copy()
        for f, t in enumerate(m.F):
            for i in range(3):
                e = m.fe[f, i]
                self.lcol[f, i] = e
                if self.h[f, i] == 0:
                    continue
                o = int(self.st.ids[f, 3 + i])
                tl, hd = int(t[i]), int(t[(i + 1) % 3])
                self.lcol[f, 3 + 2 * i] = index[tuple(sorted((tl, o)))]
                self.lcol[f, 4 + 2 * i] = index[tuple(sorted((hd, o)))]
                self.tcol[f, i] = kidx[e]
        if Forms._fns is None:
            x, hs, s, J = cj.build(verbose=False)
            fs, fJ = cj.numeric(x, hs, s, J)
            v, Pf, dP = ffj.symbolic()
            fP, fdP = ffj.numeric(v, Pf, dP)
            Forms._fns = (fs, fJ, fP, fdP)

    def readout(self, X):
        return fl.stencil_s3_live(self.st, X[self.st.ids])

    def jacobian(self, X):
        """J_ff (3 nF, nI + nE): d s / d(theta, log ell) at X."""
        fs, fJ, fP, fdP = Forms._fns
        m = self.m
        nF, nI, nE = len(m.F), len(self.edges), len(m.E)
        rows, cols, vals = [], [], []
        for f in range(nF):
            Pl = X[self.st.ids[f]].copy()
            for i in range(3):
                if self.h[f, i] == 0:
                    Pl[3 + i] = np.nan
            ell, th, h = ffj.forms_from_stencil(Pl)
            v = np.r_[ell, th]
            P = fP(v)
            for i in range(3):
                if h[i] == 0:
                    P[3 + i] = P[i] + P[(i + 1) % 3] - P[(i + 2) % 3]
            Jp = fJ(P, h) * np.repeat(np.r_[1, 1, 1, h], 3)[None, :]
            Jf = Jp @ fdP(v)
            Jf[:, :9] *= ell[None, :]
            for j in range(9):
                if self.lcol[f, j] >= 0:
                    for r in range(3):
                        rows.append(3 * f + r); cols.append(nI + self.lcol[f, j]); vals.append(Jf[r, j])
            for i in range(3):
                if self.tcol[f, i] >= 0:
                    for r in range(3):
                        rows.append(3 * f + r); cols.append(self.tcol[f, i]); vals.append(Jf[r, 9 + i])
        return sp.csr_matrix((vals, (rows, cols)), shape=(3 * nF, nI + nE))


def face_partials_frames(m, X, F, Me):
    """sum_f F_f (face f's part of each edge's CR frame at X) / Me: (2 nE, 3)."""
    nE = len(m.E)
    Fr = np.zeros((2 * nE, 3))
    for f, t in enumerate(m.F):
        Y = X[t]
        A = 0.5 * np.linalg.norm(np.cross(Y[1] - Y[0], Y[2] - Y[0]))
        L = [np.linalg.norm(Y[(i + 1) % 3] - Y[i]) for i in range(3)]
        for i in range(3):
            e = m.fe[f, i]
            tl, hd, op = i, (i + 1) % 3, (i + 2) % 3
            si = 1.0 if t[tl] < t[hd] else -1.0
            li = L[i]
            ejk, eki = L[(i + 1) % 3] ** 2, L[(i + 2) % 3] ** 2
            gp = si * A / (3 * li) * (Y[hd] - Y[tl])
            k = -si / (12 * li)
            gq = k * (li * li - (eki - ejk)) * Y[tl] + k * (li * li + (eki - ejk)) * Y[hd] + si * li / 6 * Y[op]
            Fr[e] += F[f] @ gp
            Fr[nE + e] += F[f] @ gq
    return Fr / np.concatenate([Me, Me])[:, None]


def current_log_stretch(m, X, X_rest):
    """Per face, log of the left stretch rest -> X in the current face frame."""
    T, B, _ = ops.face_frames(m, X)
    Tr, Br, _ = ops.face_frames(m, X_rest)
    out = np.zeros((len(m.F), 3))
    for f, t in enumerate(m.F):
        A = np.column_stack([Tr[f], Br[f]]).T @ np.column_stack([X_rest[t[1]] - X_rest[t[0]], X_rest[t[2]] - X_rest[t[0]]])
        Bm = np.column_stack([T[f], B[f]]).T @ np.column_stack([X[t[1]] - X[t[0]], X[t[2]] - X[t[0]]])
        Fd = Bm @ np.linalg.inv(A)
        w, v = np.linalg.eigh(Fd @ Fd.T)
        Lg = v @ np.diag(0.5 * np.log(np.maximum(w, 1e-30))) @ v.T
        out[f] = [Lg[0, 0], Lg[0, 1], Lg[1, 1]]
    return out


class ExponentSystem:
    """The constant per-restart system at X."""

    def __init__(self, fm: Forms, X, kb=1e-6, km=0.0, w_c=1e8, eps=1e-10, L_cur=None):
        m = fm.m
        self.m = m
        nF = len(m.F)
        edges = fm.edges
        nI, nE = len(edges), len(m.E)
        D0, Me = ops.cr_derivative(m, X)
        self.D0, self.Me = D0.tocsr(), Me
        self.M2 = sp.diags(np.concatenate([Me, Me]))
        self.gather = spla.factorized((self.D0.T @ self.M2 @ self.D0 + eps * sp.identity(m.n)).tocsc())
        T, B, _ = ops.face_frames(m, X)
        self.P = np.stack([T, B], 2)
        _, area = ops.face_normals(m, X)
        self.W = np.repeat(area, 3)
        J = fm.jacobian(X)
        Jt, Jl = J[:, :nI], J[:, nI:]
        rows, cols, vals = [], [], []
        for k, (e, fp, ip, fmn, im) in enumerate(edges):
            a = X[m.E[e, 1]] - X[m.E[e, 0]]
            a /= np.linalg.norm(a)
            for c in range(3):
                rows += [k, k]; cols += [3 * fmn + c, 3 * fp + c]; vals += [a[c], -a[c]]
        C = sp.csr_matrix((vals, (rows, cols)), shape=(nI, 3 * nF))
        rows, cols, vals = [], [], []
        cnt = np.bincount(m.fe.ravel(), minlength=nE).astype(float)
        for f, t in enumerate(m.F):
            for i in range(3):
                e = m.fe[f, i]
                d = X[t[(i + 1) % 3]] - X[t[i]]
                tv = self.P[f].T @ (d / np.linalg.norm(d))
                for j, Eb in enumerate(E_BASIS):
                    rows.append(e); cols.append(3 * f + j); vals.append(tv @ Eb @ tv / cnt[e])
        E = sp.csr_matrix((vals, (rows, cols)), shape=(nE, 3 * nF))
        Jc = sp.hstack([Jt @ C, Jl @ E]).tocsr()
        rows, cols, vals = [], [], []
        for k, (e, fp, ip, fmn, im) in enumerate(edges):
            d = X[m.E[e, 1]] - X[m.E[e, 0]]
            for f, sgn in ((fp, 1.0), (fmn, -1.0)):
                Dx = -skew(d)
                for r in range(3):
                    for c in range(3):
                        if Dx[r, c] != 0:
                            rows.append(3 * k + r); cols.append(3 * f + c); vals.append(sgn * Dx[r, c])
                for j, Eb in enumerate(E_BASIS):
                    v = self.P[f] @ Eb @ self.P[f].T @ d
                    for r in range(3):
                        rows.append(3 * k + r); cols.append(3 * nF + 3 * f + j); vals.append(sgn * v[r])
        Comp = sp.csr_matrix((vals, (rows, cols)), shape=(3 * nI, 6 * nF))
        JWJ = (Jc.T @ sp.diags(self.W) @ Jc).tocsr()
        self.tt = JWJ.diagonal()[:3 * nF].mean()
        CC = (Comp.T @ Comp).tocsr()
        sc = CC.diagonal().mean()
        prior = sp.diags(np.r_[np.full(3 * nF, kb * self.tt), np.full(3 * nF, (km + kb) * self.tt)])
        self.solve = spla.factorized((JWJ + (w_c * self.tt / max(sc, 1e-30)) * CC + prior).tocsc())
        self.Jc = Jc
        self.prior_rhs = np.zeros(6 * nF)
        if L_cur is not None and km > 0:
            self.prior_rhs[3 * nF:] = -km * self.tt * L_cur.ravel()

    def step(self, X, r, a):
        m = self.m
        nF = len(m.F)
        xi = self.solve(self.Jc.T @ (self.W * r.ravel()) + self.prior_rhs)
        om, L = xi[:3 * nF].reshape(nF, 3), xi[3 * nF:].reshape(nF, 3)
        F = np.array([sla.expm(a * (skew(om[f]) + self.P[f] @ sym2(L[f]) @ self.P[f].T)) for f in range(nF)])
        b = self.D0.T @ (self.M2 @ face_partials_frames(m, X, F, self.Me))
        return np.column_stack([self.gather(b[:, k]) for k in range(3)])


def run(m, X_init, s4_goal, iters=120, restart=30, kb=1e-6, km=0.0, w_c=1e8, tol=1e-13,
        elastic=False, X_rest=None, fm=None):
    """Returns (X, merit history)."""
    fm = fm or Forms(m)
    s3g = cl.sym_goal(s4_goal)
    X = X_init.copy()
    Xr = X_init.copy() if X_rest is None else X_rest
    sysm = None

    def Lc(X):
        return current_log_stretch(m, X, Xr) if elastic else None

    def merit(X):
        r = s3g - fm.readout(X)
        _, area = ops.face_normals(m, X)
        e = float(np.sum(area[:, None] * r * r))
        if elastic and km > 0 and sysm is not None:
            Lx = current_log_stretch(m, X, Xr)
            e += km * sysm.tt * float(np.sum(area[:, None] * Lx * Lx)) / area.mean()
        return r, float(np.sqrt(e))

    hist = []
    r, rn = merit(X)
    for it in range(iters):
        fresh = False
        if it % restart == 0 or sysm is None:
            sysm, fresh = ExponentSystem(fm, X, kb, km, w_c, L_cur=Lc(X)), True
            r, rn = merit(X)
        if not hist:
            hist.append(rn)
        ok = False
        for _attempt in range(2):
            a = 1.0
            for _ in range(12):
                Xn = sysm.step(X, r, a)
                rr, rnn = merit(Xn)
                if np.isfinite(rnn) and rnn < rn:
                    ok = True
                    break
                a *= 0.5
            if ok or fresh:
                break
            sysm, fresh = ExponentSystem(fm, X, kb, km, w_c, L_cur=Lc(X)), True
            r, rn = merit(X)
        if not ok:
            break
        X, r, rn = Xn, rr, rnn
        hist.append(rn)
        if rn < tol:
            break
    return X, np.array(hist)


def closed_band(ns, nz, radius, height=1.0):
    ang = 2 * np.pi * np.arange(ns) / ns
    z = np.linspace(-height / 2, height / 2, nz)
    X = np.array([[radius * np.cos(a), radius * np.sin(a), zz] for a in ang for zz in z])
    F = []
    for i in range(ns):
        for j in range(nz - 1):
            a, b = i * nz + j, ((i + 1) % ns) * nz + j
            F.append([a, b, b + 1])
            F.append([a, b + 1, a + 1])
    return X, np.array(F)


def parity_table():
    import base_loop as bl
    from edge_compare import cases
    from live_compare import shape_error
    rows = list(cases())
    X1, Fb = closed_band(24, 6, 1.0)
    X13, _ = closed_band(24, 6, 1.3)
    rows.append(("closed band r 1->1.3", X1, Fb, X13, None))
    for name, X0, F, Xt, dist in rows:
        m = ops.build_mesh(F, len(X0))
        s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
        size = np.ptp(Xt, axis=0).max()
        t = time.time()
        X, h = run(m, X0, s4, iters=90)
        print(f"{name:22s} its {len(h) - 1:3d}  residual {h[0]:.1e} -> {h[-1]:.1e}  "
              f"shape {shape_error(X0, Xt, dist, size):.1e} -> {shape_error(X, Xt, dist, size):.1e}  "
              f"t {time.time() - t:.0f}s", flush=True)


def plate_table():
    import base_loop as bl
    P, F = ops.strip_param(13, 13, 2.0, 2.0)
    Xf = ops.flat_from_param(P)
    m = ops.build_mesh(F, len(Xf))
    seed = Xf.copy()
    seed[:, 1] += 1e-3 * np.random.default_rng(3).standard_normal(len(Xf))
    fm = Forms(m)

    def cap(R):
        rr = np.linalg.norm(P, axis=1)
        d = P / np.maximum(rr, 1e-12)[:, None]
        a = rr / R
        return np.stack([R * np.sin(a) * d[:, 0], -R * (1 - np.cos(a)), R * np.sin(a) * d[:, 1]], 1)

    def ratio(X):
        c = X.mean(0)
        _, _, Vt = np.linalg.svd(X - c)
        u, v, h = (X - c) @ Vt[0], (X - c) @ Vt[1], (X - c) @ Vt[2]
        A = np.stack([0.5 * u * u, u * v, 0.5 * v * v, u, v, np.ones_like(u)], 1)
        a, b, cc = np.linalg.lstsq(A, h, rcond=None)[0][:3]
        q = np.linalg.eigvalsh(np.array([[a, b], [b, cc]]))
        return q[np.argmin(np.abs(q))] / q[np.argmax(np.abs(q))], np.abs(q).max(), np.abs(q).min()

    def strain(X):
        l0 = np.linalg.norm(Xf[m.E[:, 1]] - Xf[m.E[:, 0]], axis=1)
        return float(np.abs(np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1) / l0 - 1).max())

    for R in (1.5, 1.0):
        Xc = cap(R)
        s4 = bl.readout(m, ops.cr_derivative(m, Xc)[0], Xc)
        for km in (1.0, 1e2, 1e4):
            t = time.time()
            X, h = run(m, seed, s4, iters=90, km=km, elastic=True, X_rest=Xf, fm=fm)
            r, big, small = ratio(X)
            print(f"plate R {R} km {km:g}: ratio {r:+.3f}  |q| {big:.3f} {small:.3f}  max strain {strain(X):.1e}  "
                  f"its {len(h) - 1}  t {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    parity_table()
    plate_table()
