"""Why the position-space Dupin circle rows shrink, and ratio / shear as
settings (MATH.md 19.20).

A replica of `hepworth/block/dupin_circle.hpp`'s `circle4` on a torus whose
grid lines are exact curvature lines: groups are four consecutive vertices
along a ring (`u`) or a tube (`v`) line, rows the four positions minus their
centroid, targets the closest points on the group's fitted circle minus
their centroid. Projective dynamics: per frame `iters` local-global steps of
`(sum w A^T A + mu M) x = sum w A^T p + mu M x_frame` (one constant factor),
`mu` a weak mass toward the frame start.

1. Shrink. The fit (`kasa`: PCA plane + algebraic circle, as the C++;
   `geometric`: the same plane, then Gauss-Newton on the true distances) and
   the size fix (`spread`: targets rescaled to the live RMS spread, as the
   C++; `none`), on a clean torus (a Dupin cyclide: must be a fixed point),
   with normal noise, and with tangential noise (the zig-zag of lines that do
   not follow edges). Reported per frame: volume, area, ring radius `R`, tube
   radius `r`, against the noise-free torus.
2. Ratio: each family's groups scaled to `sqrt(rho)` (tube lines) and
   `1 / sqrt(rho)` (ring lines) times their live spread per frame: a metric
   aspect drive in position space that also sets size, so it cannot shrink.
3. Shear: each group of one family slid along the other family's direction
   by `sigma` times its arc coordinate (turning that family against the
   other; curvature lines are orthogonal, so the net leaves them).
"""

import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops

R0, r0 = 2.0, 0.7


# ------------------------------------------------------------------ fits


def plane_frame(P, plane="pca"):
    """`pca`: the least-squares plane of the four points (as the C++).
    `mid`: the plane through the three chord midpoints, in which an
    alternating sideways zig-zag cancels; PCA when they are collinear."""
    m = P.mean(0)
    C = (P - m).T @ (P - m)
    w, V = np.linalg.eigh(C)
    nu, e1 = V[:, 0], V[:, 2]
    if plane == "mid":
        M = 0.5 * (P[1:] + P[:-1])
        nn = np.cross(M[1] - M[0], M[2] - M[0])
        L = np.linalg.norm(nn)
        if L > 1e-9 * np.linalg.norm(M[2] - M[0]) ** 2:
            nu = nn / L
            e1 = P[-1] - P[0]
            e1 = e1 - (e1 @ nu) * nu
            e1 /= np.linalg.norm(e1)
    elif plane == "sag":
        # The chord and the sagitta vector (inner pair's mean off the ends'
        # mean, perpendicular to the chord): all four points, symmetric; an
        # alternating zig-zag cancels in the sagitta.
        d = P[3] - P[0]
        d = d / np.linalg.norm(d)
        s = 0.5 * (P[1] + P[2] - P[0] - P[3])
        s = s - (s @ d) * d
        if np.linalg.norm(s) > 1e-9 * np.linalg.norm(P[3] - P[0]):
            e1 = d
            nu = np.cross(d, s)
            nu /= np.linalg.norm(nu)
    e2 = np.cross(nu, e1)
    return m, nu, e1, e2


def fit_circle(P, kind="kasa", plane="pca"):
    """`(centre, normal, radius)` or None (a line)."""
    m, nu, e1, e2 = plane_frame(P, plane)
    xy = np.column_stack([(P - m) @ e1, (P - m) @ e2])
    A = np.column_stack([xy, np.ones(len(P))])
    b = -(xy ** 2).sum(1)
    s = np.linalg.lstsq(A, b, rcond=None)[0]
    cx, cy = -s[0] / 2, -s[1] / 2
    r2 = cx * cx + cy * cy - s[2]
    size = np.linalg.norm(P - m, axis=1).max()
    if not (r2 > 0) or np.sqrt(r2) > 1e3 * size:
        return None
    c, r = np.array([cx, cy]), np.sqrt(r2)
    if kind == "geometric":
        for _ in range(20):
            d = xy - c
            L = np.linalg.norm(d, axis=1)
            res = L - r
            J = np.column_stack([-d / L[:, None], -np.ones(len(P))])
            step = np.linalg.lstsq(J, -res, rcond=None)[0]
            c, r = c + step[:2], r + step[2]
            if np.abs(step).max() < 1e-14:
                break
    return m + c[0] * e1 + c[1] * e2, nu, r


def closest_on(circ, q):
    c, nu, r = circ
    d = q - np.outer((q - c) @ nu, nu) - c
    return c + r * d / np.linalg.norm(d, axis=1)[:, None]


# ------------------------------------------------------------------ setup


def torus_grid(nu=32, nv=16):
    X, F, _ = ops.torus(R0, r0, nu, nv)
    return X, F, nu, nv


def groups(nu, nv):
    """(G, 4) vertex ids and (G,) family: 0 ring lines (along u, the loose
    family on the outside), 1 tube lines (along v, the tight family)."""
    vid = lambda i, j: (i % nu) * nv + (j % nv)
    G, fam = [], []
    for j in range(nv):
        for i in range(nu):
            G.append([vid(i + k, j) for k in range(4)])
            fam.append(0)
    for i in range(nu):
        for j in range(nv):
            G.append([vid(i, j + k) for k in range(4)])
            fam.append(1)
    return np.array(G), np.array(fam)


def volume_area(X, F):
    a, b, c = X[F[:, 0]], X[F[:, 1]], X[F[:, 2]]
    V = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6
    A = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    return abs(V), A


def torus_radii(X):
    rho = np.linalg.norm(X[:, :2], axis=1)
    R = rho.mean()
    r = np.hypot(rho - R, X[:, 2]).mean()
    return R, r


def system(n, G, w, mu):
    """`sum_g w_g C^T C + mu I` on one coordinate (rows identical per axis)."""
    Cm = np.eye(4) - 0.25
    rows, cols, vals = [], [], []
    for g, ids in enumerate(G):
        for a in range(4):
            for b in range(4):
                rows.append(ids[a])
                cols.append(ids[b])
                vals.append(w[g] ** 2 * (Cm.T @ Cm)[a, b])
    K = sp.csr_matrix((vals, (rows, cols)), shape=(n, n)) + mu * sp.eye(n)
    return spla.factorized(K.tocsc())


def targets(X, G, fam, fit="kasa", size="spread", rho=1.0, sigma=0.0, shear_family=1, N=None, Xref=None,
            plane="pca"):
    """Per group the centred targets (G, 4, 3). With `rho != 1` the size is
    the frame-start spread (`Xref`) times the family's factor (a per-frame
    drive, not compounded over iterates)."""
    T = np.zeros((len(G), 4, 3))
    s_coord = np.array([-1.5, -0.5, 0.5, 1.5])
    for g, ids in enumerate(G):
        P = X[ids]
        circ = fit_circle(P, fit, plane)
        t = P.copy() if circ is None else closest_on(circ, P)
        tm, pm = t.mean(0), P.mean(0)
        t = t - tm
        st = np.sqrt((t ** 2).sum())
        if rho != 1.0:
            Pr = Xref[ids]
            sr = np.sqrt(((Pr - Pr.mean(0)) ** 2).sum())
            t *= (sr / st if st > 0 else 1.0) * (np.sqrt(rho) if fam[g] == 1 else 1 / np.sqrt(rho))
        elif size == "spread":
            sx = np.sqrt(((P - pm) ** 2).sum())
            t *= sx / st if st > 0 else 1.0
        if sigma != 0.0 and fam[g] == shear_family:
            # The other family's direction at the group: the in-surface
            # normal to the line (N x line direction), slid by sigma * s * h.
            d = P[-1] - P[0]
            h = np.linalg.norm(d) / 3
            n = N[ids].mean(0)
            b = np.cross(n, d / np.linalg.norm(d))
            b /= np.linalg.norm(b)
            t = t + sigma * h * np.outer(s_coord, b)
        T[g] = t
    return T


def rhs(n, G, w, T):
    b = np.zeros((n, 3))
    Cm = np.eye(4) - 0.25
    for g, ids in enumerate(G):
        np.add.at(b, ids, w[g] * (Cm.T @ (w[g] * T[g])))
    return b


def run(X0, F, G, fam, frames=10, iters=20, mu=0.05, fit="kasa", size="spread", rho=1.0, sigma=0.0,
        normals=None, plane="pca"):
    n = len(X0)
    w = np.ones(len(G))
    solve = system(n, G, w, mu)
    X = X0.copy()
    hist = [(*volume_area(X, F), *torus_radii(X))]
    for f in range(frames):
        Xf = X.copy()
        for _ in range(iters):
            N = normals(X) if normals is not None else None
            T = targets(X, G, fam, fit, size, rho, sigma, N=N, Xref=Xf, plane=plane)
            b = rhs(n, G, w, T) + mu * Xf
            X = np.column_stack([solve(b[:, c]) for c in range(3)])
        hist.append((*volume_area(X, F), *torus_radii(X)))
    return X, np.array(hist)


def vertex_normals(F, n):
    def f(X):
        a, b, c = X[F[:, 0]], X[F[:, 1]], X[F[:, 2]]
        fn = np.cross(b - a, c - a)
        N = np.zeros((n, 3))
        for k in range(3):
            np.add.at(N, F[:, k], fn)
        return N / np.linalg.norm(N, axis=1)[:, None]
    return f


def noisy(X, F, kind, amp, rng):
    N = vertex_normals(F, len(X))(X)
    h = 2 * np.pi * r0 / 16
    if kind == "normal":
        return X + amp * h * rng.normal(size=len(X))[:, None] * N
    # Tangential: a random in-plane offset (lines no longer pass through
    # the vertices: the zig-zag).
    g = rng.normal(size=X.shape)
    g -= np.einsum("ij,ij->i", g, N)[:, None] * N
    return X + amp * h * g


def zigzag(X, nu, nv, amp):
    """Ring lines zig-zag sideways: vertex `(i, j)` offset along the tube
    direction by `+-amp h` alternating in `i` (tube lines only get uneven
    spacing along themselves), as when groups follow edges, not lines."""
    h = 2 * np.pi * r0 / nv
    Y = X.copy()
    for i in range(nu):
        for j in range(nv):
            k = i * nv + j
            v = 2 * np.pi * j / nv
            u = 2 * np.pi * i / nu
            e_tube = np.array([-np.sin(v) * np.cos(u), -np.sin(v) * np.sin(u), np.cos(v)])
            Y[k] = X[k] + (amp * h if i % 2 == 0 else -amp * h) * e_tube
    return Y


def plane_stats(X, F, G, fam, nrm):
    """Fraction of each family's groups whose fitted plane is closer to the
    tangent plane than to the normal plane (`|nu . N| > 0.7`), and the mean
    sagitta / lateral excursion ratio."""
    N = nrm(X)
    out = []
    for f in (0, 1):
        bad = 0
        idx = np.nonzero(fam == f)[0]
        for g in idx:
            m, nu, e1, e2 = plane_frame(X[G[g]])
            bad += abs(nu @ N[G[g]].mean(0)) > 0.7
        out.append(bad / len(idx))
    return out


def line(label, H, ref):
    V, A, R, r = H[-1]
    V0, A0, R0_, r0_ = ref
    return (f"  {label:34s} volume {100 * (V / V0 - 1):+6.2f}%  area {100 * (A / A0 - 1):+6.2f}%  "
            f"R {R:.4f}  r {r:.4f}")


if __name__ == "__main__":
    which = sys.argv[1:] or ["shrink", "ratio", "shear"]
    rng = np.random.default_rng(3)
    X, F, nu, nv = torus_grid()
    G, fam = groups(nu, nv)
    ref = (*volume_area(X, F), *torus_radii(X))
    t0 = time.time()
    if "shrink" in which:
        print(f"1. Shrink (10 frames x 20 iterates, mu 0.05; clean torus: V {ref[0]:.4f} A {ref[1]:.4f} R {ref[2]:.4f} r {ref[3]:.4f})")
        starts = [("clean", X)]
        for kind in ("normal", "tangential"):
            starts.append((f"{kind} noise 0.1 h", noisy(X, F, kind, 0.1, rng)))
        for slabel, Xs in starts:
            print(f" start {slabel}: " + line("", np.array([(*volume_area(Xs, F), *torus_radii(Xs))]), ref).strip())
            for fit in ("kasa", "geometric"):
                for size in ("spread", "none"):
                    _, H = run(Xs, F, G, fam, fit=fit, size=size)
                    print(line(f"{fit:9s} {size:6s}", H, ref), flush=True)
    if "zigzag" in which:
        print("1b. Sideways zig-zag of the ring lines (amp in edge lengths), kasa + spread, per resolution")
        for nu_, nv_ in ((32, 16), (64, 32)):
            Xg, Fg, _ = ops.torus(R0, r0, nu_, nv_)
            Gg, famg = groups(nu_, nv_)
            refg = (*volume_area(Xg, Fg), *torus_radii(Xg))
            nrm = vertex_normals(Fg, len(Xg))
            for amp in (0.1, 0.25, 0.4):
                Xs = zigzag(Xg, nu_, nv_, amp)
                print(f"   {nu_}x{nv_} amp {amp}: start " +
                      line("", np.array([(*volume_area(Xs, Fg), *torus_radii(Xs))]), refg).strip(), flush=True)
                for plane in ("pca", "mid"):
                    _, H = run(Xs, Fg, Gg, famg, frames=6, iters=15, plane=plane)
                    print(line(f"plane {plane}", H, refg), flush=True)
    if "ratio" in which:
        print("2. Ratio rho per frame (tube-line spacing x sqrt(rho), ring-line x 1/sqrt(rho)), clean start, kasa + spread")
        for rho in (1.0, 1.02, 1.05):
            _, H = run(X, F, G, fam, rho=rho)
            print(line(f"rho {rho}", H, ref), flush=True)
    if "shear" in which:
        print("3. Shear sigma on the tube lines (slid along the ring direction by sigma * arc), clean start")
        nrm = vertex_normals(F, len(X))
        for sigma in (0.0, 0.05, 0.1):
            Xe, H = run(X, F, G, fam, sigma=sigma, normals=nrm)
            # How far the tube lines turned: mean angle of a tube line's chord
            # against the meridian plane.
            ids = G[fam == 1]
            d = Xe[ids[:, 3]] - Xe[ids[:, 0]]
            mid = Xe[ids].mean(1)
            u = np.arctan2(mid[:, 1], mid[:, 0])
            e_ring = np.column_stack([-np.sin(u), np.cos(u), np.zeros_like(u)])
            ang = np.degrees(np.arcsin(np.abs(np.einsum("ij,ij->i", d, e_ring)) / np.linalg.norm(d, axis=1)))
            print(line(f"sigma {sigma}", H, ref) + f"  tube-line turn {ang.mean():.2f} deg", flush=True)
    print(f"({time.time() - t0:.0f}s)")
