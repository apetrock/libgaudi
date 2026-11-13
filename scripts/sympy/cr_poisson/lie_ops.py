"""Lie sphere geometry and the Dupin cyclide coordinate (numpy port).

Ports `include/gaudi/albers/lie_sphere.hpp` and `lie_dupin.hpp` (MATH.md
section 18), checked against `include/gaudi/test/lie_sphere_tests.hpp` by
`check_lie_ops.py`. Additions for the Dupin Poisson prototype (section 19):

- `dupin_conic_basis` / `dupin_contact`: a cyclide's contact elements are the
  pairs (u, v) of null vectors, u on the conic of U, v on the conic of
  U^perp (U and U^perp are eta-orthogonal, so span{u, v} is always a null
  2-plane), so the surface is parametrized by two angles.
- `dupin_land`: closest point on a cyclide by Gauss-Newton on those angles.
- `group_polar` / `lie_procrustes`: the generalized polar decomposition onto
  O(4,2) (Higham, Mackey, Mackey, Tisseur 2004) and the motion between two
  sets of contact elements.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

ETA = np.diag([-1.0, 1.0, 1.0, 1.0, 1.0, -1.0])


def lie_dot(a, b):
    return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4] - a[5] * b[5]


def sphere(c, r):
    c = np.asarray(c, float)
    c2 = c @ c
    return np.array([(1 + c2 - r * r) / 2, (1 - c2 + r * r) / 2, c[0], c[1], c[2], r])


def point(x):
    return sphere(x, 0.0)


def plane(n, h):
    return np.array([h, -h, n[0], n[1], n[2], 1.0])


def infinity():
    return np.array([1.0, -1.0, 0, 0, 0, 0])


def curvature_sphere(p, n, k):
    return k * point(p) + plane(n, n @ p)


def decode(s, tol=1e-12):
    """('sphere', c, r) | ('plane', n, h) | ('infinity',)."""
    scale = max(np.linalg.norm(s), 1e-300)
    w = s[0] + s[1]
    if abs(w) > tol * scale:
        return ("sphere", s[2:5] / w, s[5] / w)
    if abs(s[5]) > tol * scale:
        return ("plane", s[2:5] / s[5], s[0] / s[5])
    return ("infinity",)


def contact(p, n):
    p = np.asarray(p, float)
    n = np.asarray(n, float)
    return point(p), plane(n, n @ p)


def contact_decode(a, b, tol=1e-12):
    """(p, n) of the contact element spanned by `a`, `b`, or None."""
    pt = b[5] * a - a[5] * b
    pl = (b[0] + b[1]) * a - (a[0] + a[1]) * b
    dp, dt = decode(pt, tol), decode(pl, tol)
    if dp[0] != "sphere" or dt[0] != "plane":
        return None
    n = dt[1]
    if abs(np.linalg.norm(n) - 1) > 1e-6:
        return None
    return dp[1], n


def inverse(A):
    return ETA @ A.T @ ETA


def group_defect(A):
    return np.abs(A.T @ ETA @ A - ETA).max()


def algebra_defect(X):
    return np.abs(X.T @ ETA + ETA @ X).max()


def rotation(R):
    A = np.eye(6)
    A[2:5, 2:5] = R
    return A


def translation(t):
    t = np.asarray(t, float)
    t2 = 0.5 * t @ t
    A = np.eye(6)
    A[0, 0] += t2
    A[0, 1] = t2
    A[1, 0] = -t2
    A[1, 1] -= t2
    A[0, 2:5] = t
    A[1, 2:5] = -t
    A[2:5, 0] = t
    A[2:5, 1] = t
    return A


def rigid(R, t):
    return translation(t) @ rotation(R)


def scale(lam):
    a, b = 0.5 * (1 / lam + lam), 0.5 * (1 / lam - lam)
    A = np.eye(6)
    A[0, 0] = a
    A[0, 1] = b
    A[1, 0] = b
    A[1, 1] = a
    return A


def normal_offset(d):
    d2 = 0.5 * d * d
    A = np.eye(6)
    A[0, 0] -= d2
    A[0, 1] = -d2
    A[1, 0] = d2
    A[1, 1] += d2
    A[0, 5] = d
    A[1, 5] = -d
    A[5, 0] = -d
    A[5, 1] = -d
    return A


def orientation_flip():
    return np.diag([1.0, 1, 1, 1, 1, -1])


def reflection(v):
    return np.eye(6) - (2 / lie_dot(v, v)) * np.outer(v, ETA @ v)


def mirror_plane(n, h):
    return np.array([h, -h, n[0], n[1], n[2], 0.0])


def inversion_sphere(c, rho):
    c = np.asarray(c, float)
    c2 = c @ c
    return np.array([(1 + c2 - rho * rho) / 2, (1 - c2 + rho * rho) / 2, c[0], c[1], c[2], 0.0]) / rho


def algebra_generator(i, j):
    X = np.zeros((6, 6))
    X[i, j] = ETA[j, j]
    X[j, i] = -ETA[i, i]
    return X


GENERATORS = [algebra_generator(i, j) for i in range(6) for j in range(i + 1, 6)]


def algebra(w):
    return sum(wk * G for wk, G in zip(w, GENERATORS))


def algebra_coords(X):
    return np.array([X[i, j] * ETA[j, j] for i in range(6) for j in range(i + 1, 6)])


def cayley(X):
    I = np.eye(6)
    return np.linalg.solve(I - 0.5 * X, I + 0.5 * X)


def lexp(X):
    return sla.expm(X)


def llog(A):
    return np.real(sla.logm(A))


def algebra_form(X):
    return -0.5 * np.trace(X @ X)


def axis_rotation(axis, angle):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def rotation_from_z(axis):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    z = np.array([0.0, 0, 1])
    v = np.cross(z, a)
    s, c = np.linalg.norm(v), z @ a
    if s < 1e-15:
        return np.eye(3) if c > 0 else np.diag([1.0, -1, -1])
    return axis_rotation(v / s, np.arctan2(s, c))


# ---------------------------------------------------------------- Dupin


def projector(B):
    """eta-orthogonal projector onto the columns of `B` (6 x k)."""
    G = B.T @ ETA @ B
    return B @ np.linalg.solve(G, B.T @ ETA)


def from_subspace(B):
    return 2 * projector(B) - np.eye(6)


def torus(R, r):
    """Canonical torus about z, outward normals; U holds the tube spheres."""
    B = np.zeros((6, 3))
    B[:, 0] = [(1 + R * R - r * r) / 2, (1 - R * R + r * r) / 2, 0, 0, 0, -r]
    B[2, 1] = 1
    B[3, 2] = 1
    return from_subspace(B)


def osculating_torus(p, n, k1, k2, e1, kmin_ratio=1e-3, kmin_abs=1e-3):
    """The torus osculating `(p, n)` to second order: tube spheres of
    curvature `k1` (centres on a circle through `p + n / k1` in the plane of
    `n` and the `k2` direction), axis along `e1` through `p + n / k2`. U is
    spanned by three tube spheres (their centres are linear in
    (1, cos t, sin t)), so the orientation follows the signs of `k1`, `k2`.
    `k2` is clamped away from 0 (a cylinder is the limit) by `kmin_ratio |k1|`,
    and `k1` by `kmin_abs` (a plane is the limit; flat points)."""
    p, n, e1 = (np.asarray(a, float) for a in (p, n, e1))
    if abs(k1) < kmin_abs:
        k1 = kmin_abs if k1 >= 0 else -kmin_abs
    if abs(k1 - k2) < kmin_ratio * abs(k1):
        k2 = k1 * (1 - kmin_ratio)  # umbilic: split the spheres
    kmin = kmin_ratio * abs(k1)
    if abs(k2) < kmin:
        k2 = kmin if k2 >= 0 else -kmin
    c1 = p + n / k1
    c2 = p + n / k2
    a = c1 - c2
    Rs = np.linalg.norm(a)
    u1 = a / Rs
    u2 = np.cross(e1, u1)
    B = np.stack([sphere(c2 + Rs * (np.cos(t) * u1 + np.sin(t) * u2), 1 / k1)
                  for t in (0, 2 * np.pi / 3, 4 * np.pi / 3)], 1)
    return from_subspace(B)


def transform(A, sigma):
    return A @ sigma @ inverse(A)


def torus_posed(c, axis, R, r):
    return transform(rigid(rotation_from_z(axis), c), torus(R, r))


def defect(sigma):
    es = ETA @ sigma
    return max(np.abs(sigma @ sigma - np.eye(6)).max(), np.abs(es - es.T).max(), abs(np.trace(sigma)))


def signature(sigma, sign=1, tol=1e-9):
    K = ETA @ (0.5 * (np.eye(6) + sign * sigma))
    ev = np.linalg.eigvalsh(0.5 * (K + K.T))
    s = max(np.abs(ev).max(), 1)
    return int((ev > tol * s).sum()), int((ev < -tol * s).sum())


def contact_frame(p, n):
    """Unit basis (a0, a1) of the contact element and its eta-complement (6 x 4)."""
    P, T = contact(p, n)
    a = np.stack([P / np.linalg.norm(P), T / np.linalg.norm(T)])
    A = a @ ETA
    _, _, Vt = np.linalg.svd(A)
    return a, Vt[2:].T


def contact_rows(sigma, frame):
    a, perp = frame
    return ((ETA @ (sigma @ a.T)).T @ perp).ravel()


def contact_residual(sigma, p, n):
    return np.linalg.norm(contact_rows(sigma, contact_frame(p, n)))


def curvatures(sigma, p, n):
    """(k_u, k_v, res_u, res_v): the contact element's members in U and U^perp."""
    P, T = contact(p, n)
    out = []
    for s in (1.0, -1.0):
        Q = np.eye(6) - s * sigma
        M = np.stack([Q @ P, Q @ T], 1)
        _, sv, Vt = np.linalg.svd(M)
        a, b = Vt[1]
        k = a / b if abs(b) > 1e-300 else np.inf
        out.append((k, sv[1] / max(np.linalg.norm(P) + np.linalg.norm(T), 1e-300)))
    return out[0][0], out[1][0], out[0][1], out[1][1]


def normalize(S, iterations=30, tol=1e-14):
    """Nearest cyclide coordinate: eta-symmetrize, then the matrix sign."""
    X = 0.5 * (S + ETA @ S.T @ ETA)
    for _ in range(iterations):
        nxt = 0.5 * (X + np.linalg.inv(X))
        step = np.abs(nxt - X).max()
        X = nxt
        if step < tol:
            break
    return X


_IU = np.triu_indices(6)
_TR = np.array([ETA[i, i] if i == j else 0.0 for i, j in zip(*_IU)])


def coefficient_rows(a, w):
    """Rows `a^T M w` in the 21 upper-triangle entries of symmetric M, for
    a (k x 6) and w (l x 6): returns (k * l) x 21."""
    i, j = _IU
    out = []
    for ak in a:
        for wl in w:
            r = ak[i] * wl[j] + ak[j] * wl[i]
            r[i == j] *= 0.5
            out.append(r)
    return np.array(out)


def unpack21(m):
    M = np.zeros((6, 6))
    M[_IU] = m
    M.T[_IU] = m
    return M


def pack21(M):
    return M[_IU]


def fit_normal_matrix(ps, ns, weights=None):
    N = np.zeros((21, 21))
    for k, (p, n) in enumerate(zip(ps, ns)):
        a, perp = contact_frame(p, n)
        rows = coefficient_rows(a, perp.T)
        w = 1.0 if weights is None else weights[k]
        N += w * rows.T @ rows
    return N


def fit_solve(N, count, gap_tol=1e-3):
    """(sigma, gap, ok) from accumulated normal equations (`lie_dupin_fit`)."""
    if count < 4:
        return np.eye(6), 1.0, False
    A = N + np.trace(N) * np.outer(_TR, _TR)
    ev, V = np.linalg.eigh(A)
    l0, l1 = max(ev[0], 0), ev[1]
    gap = l0 / l1 if l1 > 0 else 1.0
    S = ETA @ unpack21(V[:, 0])
    s2 = np.trace(S @ S) / 6
    if not s2 > 0:
        return np.eye(6), gap, False
    sig = normalize(S / np.sqrt(s2))
    return sig, gap, gap < gap_tol and defect(sig) < 1e-8


def fit(ps, ns, weights=None, gap_tol=1e-3):
    return fit_solve(fit_normal_matrix(ps, ns, weights), len(ps), gap_tol)


def transvection(a, b):
    return np.real(sla.sqrtm(b @ a))


def dupin_log(a, b):
    return llog(transvection(a, b))


def dupin_exp(a, xi, t=1.0):
    g = lexp(t * xi)
    return g @ a @ inverse(g)


def geodesic(a, b, t):
    return dupin_exp(a, dupin_log(a, b), t)


def holonomy(gij, gjk, gki):
    return gki @ gjk @ gij


def tangent_basis(sigma):
    C = np.stack([(0.5 * (X - sigma @ X @ sigma)).ravel() for X in GENERATORS], 1)
    U, s, _ = np.linalg.svd(C, full_matrices=False)
    k = int((s > 1e-9 * s[0]).sum())
    return [U[:, b].reshape(6, 6) for b in range(k)]


def correct(start, ps, ns, w, prior, iterations=4):
    """Damped Gauss-Newton in the 9-dim tangent space, Cayley update."""
    d = start
    frames = [contact_frame(p, n) for p, n in zip(ps, ns)]
    for _ in range(iterations):
        basis = tangent_basis(d)
        m = len(basis)
        H = prior * np.eye(m)
        g = np.zeros(m)
        dsig = [2 * B @ d for B in basis]
        for k, f in enumerate(frames):
            r0 = contact_rows(d, f)
            J = np.stack([contact_rows(ds, f) for ds in dsig], 1)
            H += w[k] * J.T @ J
            g += w[k] * J.T @ r0
        x = np.linalg.solve(H, -g)
        xi = sum(xk * B for xk, B in zip(x, basis))
        c = cayley(xi)
        d = normalize(c @ d @ inverse(c))
    return d


# ------------------------------------------------- parametrization, landing


def eta_orthonormal_basis(Pi):
    """Columns (e_-, e_1, e_2): eta-orthonormal basis of range(Pi) (signature
    (2, 1)), e_- timelike."""
    U, s, _ = np.linalg.svd(Pi)
    B = U[:, :3]
    G = B.T @ ETA @ B
    ev, Q = np.linalg.eigh(0.5 * (G + G.T))
    E = B @ Q / np.sqrt(np.abs(ev))
    neg = [k for k in range(3) if ev[k] < 0]
    pos = [k for k in range(3) if ev[k] > 0]
    if len(neg) != 1 or len(pos) != 2:
        raise ValueError(f"subspace signature {len(pos)},{len(neg)} is not (2,1)")
    return E[:, [neg[0], pos[0], pos[1]]]


def dupin_conic_basis(sigma):
    """eta-orthonormal bases of U (sigma = +1) and U^perp (sigma = -1)."""
    I = np.eye(6)
    return eta_orthonormal_basis(0.5 * (I + sigma)), eta_orthonormal_basis(0.5 * (I - sigma))


def conic_point(E, th):
    return E[:, 0] + np.cos(th) * E[:, 1] + np.sin(th) * E[:, 2]


def dupin_contact(EU, EV, th):
    """Contact element (p, n) at conic angles th = (th_u, th_v), or None."""
    return contact_decode(conic_point(EU, th[0]), conic_point(EV, th[1]))


def conic_angle(E, s):
    """Angle of the conic point nearest the direction of `s` (any vector)."""
    c = np.linalg.lstsq(E, s, rcond=None)[0]
    if c[0] < 0:
        c = -c
    return np.arctan2(c[2], c[1])


def dupin_angles_near(EU, EV, sigma, p, n):
    """Initial angles from the contact element (p, n): its member nearest U
    and nearest U^perp, projected onto each conic."""
    P, T = contact(p, n)
    th = []
    for s, E in ((1.0, EU), (-1.0, EV)):
        Q = np.eye(6) - s * sigma
        M = np.stack([Q @ P, Q @ T], 1)
        _, _, Vt = np.linalg.svd(M)
        a, b = Vt[1]
        mem = a * P + b * T
        th.append(conic_angle(E, 0.5 * (np.eye(6) + s * sigma) @ mem))
    return np.array(th)


def dupin_land(sigma, x, th0=None, n0=None, iterations=12, conic=None):
    """Closest point on the cyclide `sigma` to `x`: Gauss-Newton on the two
    conic angles. Returns (p, n, th, ok)."""
    EU, EV = conic if conic is not None else dupin_conic_basis(sigma)
    if th0 is None:
        th0 = dupin_angles_near(EU, EV, sigma, x, n0 if n0 is not None else np.array([0.0, 0, 1]))
    th = np.array(th0, float)
    h = 1e-6

    def P(t):
        pn = dupin_contact(EU, EV, t)
        return (None, None) if pn is None else pn

    p, n = P(th)
    if p is None:
        return None, None, th, False
    for _ in range(iterations):
        J = np.zeros((3, 2))
        for k in range(2):
            d = np.zeros(2)
            d[k] = h
            pp, _ = P(th + d)
            pm, _ = P(th - d)
            if pp is None or pm is None:
                return p, n, th, False
            J[:, k] = (pp - pm) / (2 * h)
        r = p - x
        step = np.linalg.lstsq(J, -r, rcond=None)[0]
        a = 1.0
        e0 = r @ r
        for _ in range(8):
            q, m = P(th + a * step)
            if q is not None and (q - x) @ (q - x) <= e0:
                break
            a *= 0.5
        else:
            break
        th = th + a * step
        p, n = q, m
        if np.abs(a * step).max() < 1e-12:
            break
    return p, n, th, True


# ------------------------------------------------------- group Procrustes


def group_polar(A, iterations=50, tol=1e-14):
    """W of the generalized polar decomposition A = W S onto O(4,2):
    Newton `W <- (W + W^{-*}) / 2`, `W^* = eta W^T eta`."""
    W = A.copy()
    for _ in range(iterations):
        nxt = 0.5 * (W + ETA @ np.linalg.inv(W).T @ ETA)
        step = np.abs(nxt - W).max()
        W = nxt
        if step < tol * max(1, np.abs(W).max()):
            break
    return W


def procrustes_rows(p0, n0, p1, n1):
    """Rows `<A a, w> = 0` in the 36 entries of A (row-major): A maps the
    contact element (p0, n0) into (p1, n1)."""
    a, _ = contact_frame(p0, n0)
    _, perp = contact_frame(p1, n1)
    rows = []
    for ak in a:
        for wl in perp.T:
            rows.append(np.outer(ETA @ wl, ak).ravel())
    return np.array(rows)


def lie_procrustes(p0, n0, p1, n1, weights=None, prior=None):
    """The O(4,2) motion taking contact elements (p0, n0) to (p1, n1),
    least squares then group projection. Returns (A, residual).

    Contact elements on one cyclide sigma_0 leave the pencil {A, A sigma_0}
    (sigma_0 fixes each of them), both of which transport sigma_0 the same
    way. The prior (default identity) projected onto the two smallest
    eigenvectors picks one; the polar step lands on the nearer of
    +-A, +-A sigma_0."""
    N = np.zeros((36, 36))
    for k in range(len(p0)):
        r = procrustes_rows(p0[k], n0[k], p1[k], n1[k])
        N += (1.0 if weights is None else weights[k]) * r.T @ r
    ev, V = np.linalg.eigh(N)
    P = np.eye(6) if prior is None else prior
    V2 = V[:, :2]
    A = (V2 @ (V2.T @ P.ravel())).reshape(6, 6)
    c2 = np.trace(inverse(A) @ A) / 6
    A = A / np.sqrt(abs(c2))
    W = group_polar(A)
    return W, ev[0] / max(ev[2], 1e-300)
