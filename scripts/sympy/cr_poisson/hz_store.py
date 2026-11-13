"""Complex curvature store (MATH.md 19.21): mean `H` plus the trace-free part
as one complex number `z`, per vertex, related to the edge store by a fixed
map.

In any orthonormal tangent frame a symmetric 2-tensor is

    S = H I + [[a, b], [b, -a]],   H = (k11 + k22) / 2,   z = a + i b = (k11 - k22) / 2 + i k12,

and on a tangent vector written as a complex number `w`, `S w = H w + z conj(w)`.
Turning the frame by `phi` sends `z -> z e^{-2 i phi}` (components in the
turned frame) and leaves `H` alone.

The edge store keeps `(k11, k12, k22)` in each edge's frame `(t_e, b_e)`
(`t_e` the lower-to-higher vertex direction). A vertex store keeps `(H_i, z_i)`
in each vertex's tangent frame `(T_i, B_i)`. The map `B` (vertex to edge):
`H_e = (H_i + H_j) / 2`, `z_e = (z_i e^{-2 i phi_ie} + z_j e^{-2 i phi_je}) / 2`,
`phi_ve` the edge direction's angle in vertex `v`'s frame. `P` projects an edge
field onto the vertex space (least squares in the spin-2 metric times the
edge's dual area).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import dupin_poisson as dp
import lie_edge_net as lsn

# The spin-2 inner product on (k11, k12, k22): z . z' + tau tau' / 2.
Q = np.array([[0.75, 0.0, 0.25], [0.0, 1.0, 0.0], [0.25, 0.0, 0.75]])


def edge_to_hz(h):
    """(n, 3) `(k11, k12, k22)` -> `H` (n,), `z` (n,) complex."""
    return 0.5 * (h[:, 0] + h[:, 2]), 0.5 * (h[:, 0] - h[:, 2]) + 1j * h[:, 1]


def hz_to_edge(H, z):
    return np.column_stack([H + z.real, z.imag, H - z.real])


def apply(H, z, w):
    """`S w = H w + z conj(w)` for complex tangent vectors `w`."""
    return H * w + z * np.conj(w)


def edge_dual_area(m: ops.Mesh, X):
    _, area = ops.face_normals(m, X)
    A = np.zeros(len(m.E))
    for k in range(3):
        np.add.at(A, m.fe[:, k], area / 3)
    return A


class VertexStore:
    """Frames, edge angles, the map `B` and the projection `P` at `X`."""

    def __init__(self, m: ops.Mesh, X, N=None, R=None):
        self.m = m
        N = dp.vertex_normals(m, X) if N is None else N
        self.R = lsn.reference_dirs(N) if R is None else R
        self.T, self.Bv, self.N = lsn.vertex_frames(N, self.R)
        i, j = m.E[:, 0], m.E[:, 1]
        t = X[j] - X[i]
        self.phi = np.zeros((len(m.E), 2))
        for k, v in enumerate((i, j)):
            tt = t - np.einsum("ij,ij->i", t, self.N[v])[:, None] * self.N[v]
            self.phi[:, k] = np.arctan2(np.einsum("ij,ij->i", tt, self.Bv[v]), np.einsum("ij,ij->i", tt, self.T[v]))
        self.B = self._build_B()
        self.W = sp.kron(sp.diags(edge_dual_area(m, X)), sp.csr_matrix(Q)).tocsr()
        BtW = self.B.T @ self.W
        self._solve = spla.factorized((BtW @ self.B + 1e-12 * sp.eye(self.B.shape[1])).tocsc())
        self._BtW = BtW

    def _build_B(self):
        """Real (3 nE) x (3 nV): vertex `(H, Re z, Im z)` -> edge `(k11, k12, k22)`."""
        m = self.m
        rows, cols, vals = [], [], []
        for e, (i, j) in enumerate(m.E):
            for k, v in enumerate((i, j)):
                c = np.exp(-2j * self.phi[e, k])
                cr, ci = c.real, c.imag
                # k11 = H + Re z_e, k12 = Im z_e, k22 = H - Re z_e (each endpoint / 2).
                # Re(z c) = x cr - y ci, Im(z c) = x ci + y cr, z = x + i y.
                for r, (dH, dx, dy) in enumerate(((1, cr, -ci), (0, ci, cr), (1, -cr, ci))):
                    for col, val in ((0, dH), (1, dx), (2, dy)):
                        if val != 0:
                            rows.append(3 * e + r)
                            cols.append(3 * v + col)
                            vals.append(0.5 * val)
        return sp.csr_matrix((vals, (rows, cols)), shape=(3 * len(m.E), 3 * m.n))

    def to_edges(self, s):
        """Vertex store `(nV, 3)` `(H, Re z, Im z)` -> edge `(nE, 3)`."""
        return (self.B @ s.ravel()).reshape(-1, 3)

    def fit(self, h):
        """The vertex store whose edges best match `h` (nE, 3)."""
        return self._solve(self._BtW @ h.ravel()).reshape(-1, 3)

    def project(self, h):
        return self.to_edges(self.fit(h))


def relative(a, b, m=None, X=None):
    """Relative difference of two edge fields in the area-weighted spin-2 norm."""
    d = (a - b)
    if m is None:
        return float(np.linalg.norm(d) / max(np.linalg.norm(b), 1e-300))
    A = edge_dual_area(m, X)
    num = np.einsum("e,ei,ij,ej->", A, d, Q, d)
    den = np.einsum("e,ei,ij,ej->", A, b, Q, b)
    return float(np.sqrt(num / max(den, 1e-300)))


if __name__ == "__main__":
    import lie_edges as le
    import shape_dirichlet as sd
    import face_loops as fl
    from lie_compare import cases

    rng = np.random.default_rng(0)
    fails = 0

    def expect(ok, what, v=None):
        global fails
        fails += 0 if ok else 1
        print(f"  {'ok  ' if ok else 'FAIL'} {what}" + ("" if v is None else f"  ({v:.2e})"))

    print("Algebra")
    k = rng.normal(size=(6, 3))
    H, z = edge_to_hz(k)
    w = rng.normal(size=6) + 1j * rng.normal(size=6)
    Sw = np.array([np.array([[a, b], [b, c]]) @ np.array([ww.real, ww.imag]) for (a, b, c), ww in zip(k, w)])
    got = apply(H, z, w)
    expect(np.abs(Sw[:, 0] + 1j * Sw[:, 1] - got).max() < 1e-12, "S w = H w + z conj(w)",
           np.abs(Sw[:, 0] + 1j * Sw[:, 1] - got).max())
    expect(np.abs(hz_to_edge(H, z) - k).max() < 1e-14, "edge <-> (H, z) round trip")
    phi = 0.7
    c, s = np.cos(phi), np.sin(phi)
    Rm = np.array([[c, -s], [s, c]])
    worst = 0.0
    for (a, b, cc), zz in zip(k, z):
        Sp = Rm.T @ np.array([[a, b], [b, cc]]) @ Rm  # components in the frame turned by phi
        zp = 0.5 * (Sp[0, 0] - Sp[1, 1]) + 1j * Sp[0, 1]
        worst = max(worst, abs(zp - zz * np.exp(-2j * phi)))
    expect(worst < 1e-12, "frame turned by phi: z -> z e^{-2 i phi}", worst)

    print("Representational loss on the suite targets (edge store -> vertex store -> edges)")
    for name, X0, F, Xt, dist in cases():
        m = ops.build_mesh(F, len(X0))
        vs = VertexStore(m, Xt)
        st = fl.build_stencils(m, ops.cr_derivative(m, Xt)[0])
        h_cr = sd.edge_average_readings(m, st, Xt)
        h_cn = le_sph = lsn.read_spheres(m, le.hinges(m), Xt).H
        lcr = relative(vs.project(h_cr), h_cr, m, Xt)
        lcn = relative(vs.project(h_cn), h_cn, m, Xt)
        # Frame invariance: a different reference direction gives the same edges.
        vs2 = VertexStore(m, Xt, R=np.cos(0.6) * vs.T + np.sin(0.6) * vs.Bv)
        inv = relative(vs2.project(h_cr), vs.project(h_cr))
        print(f"  {name:22s} CR read-out loss {lcr:.3f}   contact read-out loss {lcn:.3f}   "
              f"frame invariance {inv:.1e}   unknowns {3 * m.n} vs {3 * len(m.E)}")
        if inv > 1e-8:
            fails += 1
    # A sphere: curvature I/(-r) everywhere, exact in the vertex store.
    name, X0, F, Xt, dist = [c for c in cases() if c[0].startswith("cube->sphere")][0]
    m = ops.build_mesh(F, len(Xt))
    vs = VertexStore(m, Xt)
    r = np.linalg.norm(Xt - Xt.mean(0), axis=1).mean()
    h = np.tile([-1 / r, 0.0, -1 / r], (len(m.E), 1))
    expect(relative(vs.project(h), h) < 1e-10, "umbilic sphere field exact under B", relative(vs.project(h), h))
    print("FAILED" if fails else "all passed")
