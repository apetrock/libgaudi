"""Curvature advection on a vertex (H, z) store (MATH.md 19.21).

Semi-Lagrangian transport of the curvature tensor along a tangent velocity
field on a triangle mesh:

- back trace: RK2, walking across faces (straightest path: the direction is
  turned about each crossed edge into the next face's plane);
- sample: each corner's tensor in ambient form `A = [T B] S [T B]^T`,
  interpolated barycentrically at the departure point (no frame angles);
- transport: `pullback` (Nabizadeh et al. 2022): `A <- J^T A J`, `J` the
  inverse flow map's Jacobian `exp(-dt G)`, `G` the one-ring velocity
  gradient with its normal column completed skew (exact for rigid motions),
  so the tensor turns and stretches with the flow;
  `parallel`: `J = I` (ambient components carried, not deformed); either
  way projected into the destination's tangent plane and read back as
  `(H, z)`; `frame`: `(H, z e^{-2 i psi})` interpolated as numbers against a
  reference tangent field (`set_reference`), so a field the flow only
  translates along the reference stays exact;
- schemes: `first` (Stam), `bfecc` (Kim et al. 2005: forward, backward, half
  the error corrected, forward again), `maccormack` (Selle et al. 2008:
  forward, backward for the error estimate, corrected, clamped to the range
  of the corner values the forward ray sampled).

Fields: `H` (n,), `z` (n,) complex in the vertex frames of `hz_store`.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

import cr_ops as ops
import dupin_poisson as dp
import lie_edge_net as lsn


class SurfaceAdvector:
    def __init__(self, m: ops.Mesh, X, N=None):
        self.m, self.X = m, X
        self.F = m.F
        self.Fn, _ = ops.face_normals(m, X)
        N = dp.vertex_normals(m, X) if N is None else N
        self.T, self.B, self.N = lsn.vertex_frames(N, lsn.reference_dirs(N))
        nF = len(m.F)
        # Face across each local edge k (between F[f, k] and F[f, k + 1]).
        self.nbr = -np.ones((nF, 3), dtype=int)
        for f in range(nF):
            for k in range(3):
                fs = m.edge_faces[m.fe[f, k]]
                for g in fs:
                    if g != f:
                        self.nbr[f, k] = g
        self.vfaces = [[] for _ in range(m.n)]
        for f, tri in enumerate(m.F):
            for v in tri:
                self.vfaces[v].append(f)

    # -- geometry on faces
    def bary(self, f, p):
        a, b, c = self.X[self.F[f]]
        M = np.column_stack([b - a, c - a])
        beta, gamma = np.linalg.lstsq(M, p - a, rcond=None)[0]
        return np.array([1 - beta - gamma, beta, gamma])

    def bary_dir(self, f, d):
        a, b, c = self.X[self.F[f]]
        M = np.column_stack([b - a, c - a])
        db, dg = np.linalg.lstsq(M, d, rcond=None)[0]
        return np.array([-db - dg, db, dg])

    def start_face(self, v, d):
        best, best_val = self.vfaces[v][0], -np.inf
        for f in self.vfaces[v]:
            n = self.Fn[f]
            dp_ = d - (d @ n) * n
            lam = self.bary(f, self.X[v] + 1e-6 * dp_ / max(np.linalg.norm(dp_), 1e-300))
            if lam.min() > best_val:
                best, best_val = f, lam.min()
        return best

    def walk(self, f, p, disp):
        """Move `disp` (3D) along the surface from `p` in face `f`."""
        n = self.Fn[f]
        disp = disp - (disp @ n) * n
        L = np.linalg.norm(disp)
        if not (L > 0):
            return f, p
        d = disp / L
        for _ in range(500):
            lam = self.bary(f, p)
            dl = self.bary_dir(f, d)
            t_exit, k_exit = np.inf, -1
            scale = np.abs(dl).max()
            for k in range(3):
                if dl[k] < -1e-15:
                    t = -lam[k] / dl[k]
                    # On a vertex or an edge of the face (`lam_k ~ 0`) a start
                    # heading into the face reads as an exit at `t ~ 0`; only a
                    # direction clearly leaving the face counts.
                    if t < 1e-12 and not (dl[k] < -1e-6 * scale and lam[k] < 1e-10):
                        continue
                    if t < 1e-12 and lam[k] >= 1e-10:
                        continue
                    if t < t_exit:
                        t_exit, k_exit = t, k
            if t_exit >= L:
                return f, p + L * d
            t_exit = max(t_exit, 0.0)
            p = p + t_exit * d
            L -= t_exit
            # Crossing the edge opposite corner k_exit: local edge (k_exit + 1) % 3.
            le_ = (k_exit + 1) % 3
            g = self.nbr[f, le_]
            if g < 0:
                return f, p
            q0, q1 = self.X[self.F[f, le_]], self.X[self.F[f, (le_ + 1) % 3]]
            e = (q1 - q0) / np.linalg.norm(q1 - q0)
            nf, ng = self.Fn[f], self.Fn[g]
            ang = np.arctan2(np.cross(nf, ng) @ e, nf @ ng)
            d = d * np.cos(ang) + np.cross(e, d) * np.sin(ang) + e * (e @ d) * (1 - np.cos(ang))
            d = d - (d @ ng) * ng
            d /= max(np.linalg.norm(d), 1e-300)
            f = g
            p = p + 1e-12 * d
        return f, p

    def velocity_at(self, f, p, V):
        lam = np.clip(self.bary(f, p), 0, 1)
        lam /= lam.sum()
        v = lam @ V[self.F[f]]
        n = self.Fn[f]
        return v - (v @ n) * n

    def trace(self, f, p, V, dt, sign):
        """RK2 from `(f, p)` for `sign * V * dt`."""
        v1 = self.velocity_at(f, p, V)
        fm, pm = self.walk(f, p, sign * 0.5 * dt * v1)
        vm = self.velocity_at(fm, pm, V)
        return self.walk(f, p, sign * dt * vm)

    # -- tensors
    def ambient(self, H, z):
        T, B, N = self.T, self.B, self.N
        Pt = np.eye(3)[None] - np.einsum("ni,nj->nij", N, N)
        return (H[:, None, None] * Pt + z.real[:, None, None] * (np.einsum("ni,nj->nij", T, T) - np.einsum("ni,nj->nij", B, B))
                + z.imag[:, None, None] * (np.einsum("ni,nj->nij", T, B) + np.einsum("ni,nj->nij", B, T)))

    def read(self, v, A):
        T, B = self.T[v], self.B[v]
        tt, bb, tb = T @ A @ T, B @ A @ B, T @ A @ B
        return 0.5 * (tt + bb), 0.5 * (tt - bb) + 1j * tb

    def velocity_gradient(self, V):
        """Per vertex the ambient gradient `G` of the velocity on the tangent
        plane (least squares over the one ring: `v_j - v_i ~ G (x_j - x_i)`)."""
        m = self.m
        ring = dp.rings(m)
        G = np.zeros((m.n, 3, 3))
        for i in range(m.n):
            d = self.X[ring[i]] - self.X[i]
            dv = V[ring[i]] - V[i]
            n = self.N[i]
            Pt = np.eye(3) - np.outer(n, n)
            # Solve G in the tangent basis: dv = G2 [d.T, d.B].
            uv = np.column_stack([d @ self.T[i], d @ self.B[i]])
            G2 = np.linalg.lstsq(uv, dv, rcond=None)[0].T  # 3 x 2
            G[i] = G2 @ np.vstack([self.T[i], self.B[i]])
            # The normal column: the normal turns with the flow, so complete
            # the normal-tangent block skew (`n . G e = -e . G n`); exact for a
            # rigid motion, whose gradient is skew.
            gn = -(n @ G2[:, 0]) * self.T[i] - (n @ G2[:, 1]) * self.B[i]
            G[i] += np.outer(gn, n)
        return G

    def set_reference(self, Rf):
        """A tangent reference field for `transport="frame"`: `z` is carried
        relative to it (its phase `psi` in each vertex frame), i.e. the
        reference field is treated as parallel."""
        self.psi = np.arctan2(np.einsum("ij,ij->i", Rf, self.B), np.einsum("ij,ij->i", Rf, self.T))

    def advect(self, H, z, V, dt, sign=-1, transport="pullback", bounds=False, G=None):
        """One semi-Lagrangian step: each vertex takes the tensor found by
        tracing `sign * V * dt` (back trace for `sign = -1`). The pullback's
        Jacobian is `exp(sign dt G)`, `G` the velocity gradient (smooth;
        finite differences of traces amplify the faceted path error)."""
        m = self.m
        Amb = self.ambient(H, z)
        Hn, zn = np.zeros(m.n), np.zeros(m.n, complex)
        lo = np.zeros((m.n, 3))
        hi = np.zeros((m.n, 3))
        if transport == "pullback" and G is None:
            G = self.velocity_gradient(V)
        for v in range(m.n):
            vel = V[v]
            f0 = self.start_face(v, sign * vel if np.linalg.norm(vel) > 0 else self.T[v])
            f, p = self.trace(f0, self.X[v], V, dt, sign)
            lam = np.clip(self.bary(f, p), 0, 1)
            lam /= lam.sum()
            corners = self.F[f]
            if transport == "frame":
                zr = z[corners] * np.exp(-2j * self.psi[corners])
                Hn[v] = lam @ H[corners]
                zn[v] = (lam @ zr) * np.exp(2j * self.psi[v])
                if bounds:
                    lo[v, :2] = H[corners].min(), np.abs(z[corners]).min()
                    hi[v, :2] = H[corners].max(), np.abs(z[corners]).max()
                continue
            A = np.einsum("k,kij->ij", lam, Amb[corners])
            J3 = None
            if transport == "pullback":
                J3 = sla.expm(sign * dt * G[v])
                Pt = np.eye(3) - np.outer(self.N[v], self.N[v])
                J3 = J3 @ Pt
                A = J3.T @ A @ J3
            Hn[v], zn[v] = self.read(v, A)
            if bounds:
                vals = []
                for c in corners:
                    Ac = Amb[c] if J3 is None else J3.T @ Amb[c] @ J3
                    hc, zc = self.read(v, Ac)
                    vals.append((hc, abs(zc)))
                vals = np.array(vals)
                lo[v, :2], hi[v, :2] = vals.min(0), vals.max(0)
        if bounds:
            return Hn, zn, lo, hi
        return Hn, zn

    def step(self, H, z, V, dt, scheme="maccormack", transport="pullback"):
        G = self.velocity_gradient(V) if transport == "pullback" else None
        if scheme == "first":
            return self.advect(H, z, V, dt, -1, transport, G=G)
        if scheme == "bfecc":
            H1, z1 = self.advect(H, z, V, dt, -1, transport, G=G)
            H2, z2 = self.advect(H1, z1, V, dt, +1, transport, G=G)
            Hb, zb = H + 0.5 * (H - H2), z + 0.5 * (z - z2)
            return self.advect(Hb, zb, V, dt, -1, transport, G=G)
        H1, z1, lo, hi = self.advect(H, z, V, dt, -1, transport, bounds=True, G=G)
        H2, z2 = self.advect(H1, z1, V, dt, +1, transport, G=G)
        Hc, zc = H1 + 0.5 * (H - H2), z1 + 0.5 * (z - z2)
        if scheme == "maccormack_free":
            return Hc, zc
        # Selle's clamp on `H`; on `z` only the modulus is bounded (by the
        # departure corners'), which keeps the limiter frame invariant.
        Hc = np.clip(Hc, lo[:, 0], hi[:, 0])
        r = np.abs(zc)
        zc = np.where(r > hi[:, 1], zc * hi[:, 1] / np.maximum(r, 1e-300), zc)
        return Hc, zc


def tangent(m, X, N, vec_fn):
    """Per-vertex tangent vectors from an ambient field, projected."""
    V = np.array([vec_fn(x) for x in X])
    return V - np.einsum("ij,ij->i", V, N)[:, None] * N
