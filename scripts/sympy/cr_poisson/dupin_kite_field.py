"""Integrable Moebius kites from one direction field (MATH.md 19.20).

The demo's pathway, in Python: a tangent direction field (here the
principal direction `arg(z) / 2` of the vertex `(H, z)` store, sign-aligned
to the previous frame's field, standing in for the Dirichlet edge tangents),
turned by `angle` (0: along the tight family, trapezoids; 45 degrees: on the
cell diagonals, kites), scaled by `alpha |z|` (fades where the direction is
undefined), then one Poisson solve `L u = div g` for the log-scale `u`
(only the integrable part survives, so neighbouring cells agree). Each
cell's Moebius shift comes from the local `grad u` (`a = rho_c grad u / 2`
in the cell circle's unit-disk coordinates), matched against the frame-start
cell by an area-held Kabsch fit (variational).

Compared: `field` (the above), `direct` (each cell's shift from `g` itself,
no Poisson: independent kites), `checker` (`direct` with a checkerboard
sign). Torus with a bump (curvature varies both ways), 32 x 16 and 64 x 32.
"""

import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import advect_geometry as ag
import cr_ops as ops
import dupin_circle_study as st
import dupin_mobius_shift as ms
import hz_store as hz


def face_gradient(m, X):
    """(3 nF) x n gradient of a vertex scalar per face, and face areas."""
    rows, cols, vals = [], [], []
    A = np.zeros(len(m.F))
    for f, (i, j, k) in enumerate(m.F):
        a, b, c = X[i], X[j], X[k]
        n = np.cross(b - a, c - a)
        area2 = np.linalg.norm(n)
        A[f] = 0.5 * area2
        n /= area2
        for v, (p, q) in zip((i, j, k), ((b, c), (c, a), (a, b))):
            gv = np.cross(n, q - p) / area2  # grad of the hat function at v
            for d in range(3):
                rows.append(3 * f + d)
                cols.append(v)
                vals.append(gv[d])
    return sp.csr_matrix((vals, (rows, cols)), shape=(3 * len(m.F), m.n)), A


def poisson_u(m, X, g):
    """`argmin sum_f A_f |grad u - g_f|^2`, mean zero (`g` per vertex)."""
    Gr, A = face_gradient(m, X)
    gf = g[m.F].mean(1)
    N = np.cross(X[m.F[:, 1]] - X[m.F[:, 0]], X[m.F[:, 2]] - X[m.F[:, 0]])
    N /= np.linalg.norm(N, axis=1)[:, None]
    gf -= np.einsum("ij,ij->i", gf, N)[:, None] * N
    W = sp.diags(np.repeat(A, 3))
    K = (Gr.T @ W @ Gr).tocsc()
    b = Gr.T @ (W @ gf.ravel())
    K = K + 1e-9 * K.diagonal().mean() * sp.eye(m.n)
    u = spla.spsolve(K, b)
    return u - u.mean(), (Gr @ u).reshape(-1, 3)


def direction_field(m, X, prev=None, angle=0.0):
    """Principal direction of the vertex store, sign-aligned to `prev`,
    turned by `angle`; and `|z|`."""
    vs = hz.VertexStore(m, X)
    s = ag.cr_field(m, X, vs)
    z = s[:, 1] + 1j * s[:, 2]
    th = 0.5 * np.angle(z) + angle
    t = np.cos(th)[:, None] * vs.T + np.sin(th)[:, None] * vs.Bv
    if prev is not None:
        t *= np.sign(np.einsum("ij,ij->i", t, prev) + 1e-30)[:, None]
    return t, np.abs(z)


def cell_targets(Xf, C, a_vec, scale=None):
    """Per cell the shifted, centred target (a_vec: per-cell 3D shift
    vector, `a` = its in-plane part in unit-disk units), its area the
    frame-start area times `scale^2` (1: area held)."""
    T = []
    for g, ids in enumerate(C):
        P = Xf[ids]
        f = ms.circle_frame(P)
        if f is None:
            T.append(P - P.mean(0))
            continue
        c, e1, e2, r = f
        z = ((P - c) @ e1 + 1j * ((P - c) @ e2)) / r
        z /= np.abs(z)
        a = (a_vec[g] @ e1 + 1j * (a_vec[g] @ e2))
        zp = (z - a) / (1 - np.conj(a) * z) if abs(a) > 0 else z
        t = c + r * (np.outer(zp.real, e1) + np.outer(zp.imag, e2))
        t -= t.mean(0)
        s0, s1 = np.sqrt(ms.quad_area(P)), np.sqrt(ms.quad_area(t))
        k = 1.0 if scale is None else scale[g]
        T.append(t * (k * s0 / s1) if s1 > 0 else t)
    return T


def angles(X, C):
    P = X[C]
    out = np.zeros((len(C), 4))
    for k in range(4):
        u = P[:, (k - 1) % 4] - P[:, k]
        v = P[:, (k + 1) % 4] - P[:, k]
        out[:, k] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", u, v) /
                                                 (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1)))
    return out


def run(nu, nv, mode, alpha=0.5, angle=np.pi / 4, frames=6, iters=20, mu=0.05, base=None):
    """`base`: the angles of the same run at `alpha = 0` (the circle rows
    alone); the realized kite is measured against it."""
    m, X0 = ag.make(nu, nv)
    F = m.F
    C = ms.cells(nu, nv)
    n = m.n
    w = np.ones(len(C))
    solve = st.system(n, C, w, mu)
    X = X0.copy()
    prev = None
    a0 = angles(X0, C)
    pred = np.zeros_like(a0)
    checker = np.array([(-1) ** ((g // nv) + (g % nv)) for g in range(len(C))])
    V0, A0 = st.volume_area(X0, F)
    for fr in range(frames):
        Xf = X.copy()
        t, zabs = direction_field(m, Xf, prev, angle)
        prev = t
        g = alpha * zabs[:, None] * t
        rho_c = np.array([np.linalg.norm(Xf[ids] - Xf[ids].mean(0), axis=1).mean() for ids in C])
        scale = None
        if mode in ("field", "field_scaled"):
            u, gu = poisson_u(m, Xf, g)
            # Cell gradient: mean of its two triangles' gradients.
            fidx = {tuple(sorted(f)): k for k, f in enumerate(F)}
            gcell = np.zeros((len(C), 3))
            for k, ids in enumerate(C):
                tri = [fidx.get(tuple(sorted((ids[0], ids[1], ids[2])))), fidx.get(tuple(sorted((ids[0], ids[2], ids[3]))))]
                tri = [q for q in tri if q is not None]
                gcell[k] = gu[tri].mean(0) if tri else 0
            a_vec = 0.5 * rho_c[:, None] * gcell
            if mode == "field_scaled":
                # The conformal growth's own scale: cell sizes follow e^u.
                scale = np.exp(u[C].mean(1))
        else:
            gcell = g[C].mean(1)
            if mode == "checker":
                gcell = gcell * checker[:, None]
            a_vec = 0.5 * rho_c[:, None] * gcell
        T = cell_targets(Xf, C, a_vec, scale)
        T0 = cell_targets(Xf, C, np.zeros_like(a_vec))

        def tangles(TT):
            out = np.zeros((len(C), 4))
            for k, tk in enumerate(TT):
                for q in range(4):
                    u_ = tk[(q - 1) % 4] - tk[q]
                    v_ = tk[(q + 1) % 4] - tk[q]
                    out[k, q] = np.degrees(np.arccos(np.clip(u_ @ v_ / np.linalg.norm(u_) / np.linalg.norm(v_), -1, 1)))
            return out

        # Predicted kite this frame: the shifted targets' angles against the
        # unshifted (circle-projected) ones, so the circle projection itself
        # does not count.
        pred += tangles(T) - tangles(T0)
        for _ in range(iters):
            Tm = np.array([ms.kabsch(X[ids] - X[ids].mean(0), T[k]) for k, ids in enumerate(C)])
            b = st.rhs(n, C, w, Tm) + mu * Xf
            X = np.column_stack([solve(b[:, c]) for c in range(3)])
    if base is None:
        base = a0
    real = angles(X, C) - base
    corr = float(np.corrcoef(real.ravel(), pred.ravel())[0, 1]) if np.std(pred) > 0 else np.nan
    frac = float(np.sum(real * pred) / max(np.sum(pred * pred), 1e-30))
    # Roughness: a cell's change against its four neighbours' mean.
    vid = lambda i, j: (i % nu) * nv + (j % nv)
    nb = np.array([[vid(i + 1, j), vid(i - 1, j), vid(i, j + 1), vid(i, j - 1)] for i in range(nu) for j in range(nv)])
    rough = float(np.abs(real - real[nb].mean(1)).mean() / max(np.abs(real).mean(), 1e-30))
    V, A = st.volume_area(X, F)
    h = np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean()
    if alpha == 0:
        return angles(X, C)
    return dict(real=float(np.abs(real).mean()), pred=float(np.abs(pred).mean()), corr=corr, frac=frac,
                rough=rough, dV=100 * (V / V0 - 1), dA=100 * (A / A0 - 1), h=h)


if __name__ == "__main__":
    t0 = time.time()
    res = sys.argv[1:] or ["32", "64"]
    for nu in [int(r) for r in res]:
        nv = nu // 2
        base = run(nu, nv, "field", alpha=0.0)
        for alpha, modes in ((0.5, ("field", "direct", "checker")), (0.02, ("field", "field_scaled", "direct"))):
            print(f"torus with a bump, {nu} x {nv}, alpha {alpha}, angle 45 deg, 6 frames (against alpha 0)", flush=True)
            for mode in modes:
                r = run(nu, nv, mode, alpha=alpha, base=base)
                print(f"  {mode:12s} kite real {r['real']:.3f} deg (asked {r['pred']:.3f}, realized {r['frac']:.2f}, "
                      f"corr {r['corr']:.2f})  roughness {r['rough']:.2f}  volume {r['dV']:+.2f}%  area {r['dA']:+.2f}%",
                      flush=True)
    print(f"({time.time() - t0:.0f}s)")
