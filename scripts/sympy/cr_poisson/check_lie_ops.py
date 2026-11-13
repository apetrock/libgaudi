"""`lie_ops.py` against `include/gaudi/test/lie_sphere_tests.hpp`, plus the
prototype's additions (conic parametrization, landing, group Procrustes)."""

import numpy as np

import lie_ops as lo

FAILS = []


def expect(cond, what):
    if not cond:
        FAILS.append(what)
        print("FAIL", what)


def torus_contact(R, r, th, ph):
    n = np.array([np.cos(ph) * np.cos(th), np.cos(ph) * np.sin(th), np.sin(ph)])
    return np.array([R * np.cos(th), R * np.sin(th), 0]) + r * n, n


def push_contact(A, p, n):
    P, T = lo.contact(p, n)
    return lo.contact_decode(A @ P, A @ T)


def sdist(a, b):
    return np.abs(a - b).max()


def usdist(a, b):
    return min(sdist(a, b), sdist(a, -b))


def cyclide_map():
    R = lo.axis_rotation([1, 2, -0.5], 0.7)
    return lo.rigid(R, [0.3, -0.2, 0.5]) @ lo.reflection(lo.inversion_sphere([0.4, 2.2, 1.1], 1.5))


def vectors_and_contact():
    c = np.array([0.3, -1.2, 2.0])
    s = lo.sphere(c, -0.7)
    expect(abs(lo.lie_dot(s, s)) < 1e-12, "null sphere")
    d = lo.decode(3.5 * s)
    expect(d[0] == "sphere" and np.linalg.norm(d[1] - c) < 1e-12 and abs(d[2] + 0.7) < 1e-12, "decode sphere")
    n = np.array([1, -2, 0.5]) / np.linalg.norm([1, -2, 0.5])
    t = lo.plane(n, 0.8)
    d = lo.decode(-2 * t)
    expect(d[0] == "plane" and np.linalg.norm(d[1] - n) < 1e-12 and abs(d[2] - 0.8) < 1e-12, "decode plane")
    c1, c2, r1, r2 = np.array([0.1, 0.2, 0.3]), np.array([1.0, -0.5, 0.7]), 0.4, -0.9
    e = 0.5 * ((r1 - r2) ** 2 - (c1 - c2) @ (c1 - c2))
    expect(abs(lo.lie_dot(lo.sphere(c1, r1), lo.sphere(c2, r2)) - e) < 1e-12, "contact form")
    p = np.array([0.4, -0.1, 1.3])
    n = np.array([-0.3, 0.8, 0.5]) / np.linalg.norm([-0.3, 0.8, 0.5])
    s = lo.curvature_sphere(p, n, 2.5)
    d = lo.decode(s)
    expect(np.linalg.norm(d[1] - (p + n / 2.5)) < 1e-12 and abs(d[2] - 0.4) < 1e-12, "curvature sphere")
    a = lo.curvature_sphere(p, n, -1.7)
    q, m = lo.contact_decode(a, s)
    expect(np.linalg.norm(q - p) < 1e-10 and np.linalg.norm(m - n) < 1e-10, "contact decode")


def group():
    R = lo.axis_rotation([0.3, -1, 0.4], 1.1)
    for A in [lo.rotation(R), lo.translation([0.5, -1.5, 2]), lo.rigid(R, [-0.2, 0.1, 0.9]), lo.scale(1.7),
              lo.normal_offset(0.35), lo.orientation_flip(), lo.reflection(lo.mirror_plane([0, 0.6, 0.8], 0.25)),
              lo.reflection(lo.inversion_sphere([1, 2, 3], 0.75))]:
        tol = 1e-14 * max(1, (A * A).sum())
        expect(lo.group_defect(A) < tol, "group defect")
        expect(np.abs(lo.inverse(A) @ A - np.eye(6)).max() < tol, "group inverse")
    w = 0.1 * np.sin(1.0 + 2.3 * np.arange(15))
    X = lo.algebra(w)
    expect(lo.algebra_defect(X) < 1e-15 and np.linalg.norm(lo.algebra_coords(X) - w) < 1e-15, "algebra")
    expect(lo.group_defect(lo.lexp(X)) < 1e-12 and lo.group_defect(lo.cayley(X)) < 1e-12, "exp/cayley")
    expect(np.abs(lo.llog(lo.lexp(X)) - X).max() < 1e-10, "log")
    expect(lo.algebra_form(lo.algebra_generator(2, 3)) > 0 and lo.algebra_form(lo.algebra_generator(0, 2)) < 0,
           "invariant form signs")
    c, r, t = np.array([0.3, 0.4, -0.5]), 0.6, np.array([1, -2, 0.5])
    Ry = lo.axis_rotation([0, 1, 0], 0.9)
    d = lo.decode(lo.rigid(Ry, t) @ lo.sphere(c, r))
    expect(np.linalg.norm(d[1] - (Ry @ c + t)) < 1e-12 and abs(d[2] - r) < 1e-12, "rigid acts")
    d = lo.decode(lo.normal_offset(0.2) @ lo.sphere(c, r))
    expect(abs(d[2] - (r - 0.2)) < 1e-12, "offset acts")
    x, n = np.array([0.7, -0.4, 1.9]), np.array([0, 0.6, 0.8])
    pm = push_contact(lo.normal_offset(0.4), x, n)
    expect(np.linalg.norm(pm[0] - (x + 0.4 * n)) < 1e-12, "offset slides contact")
    ic, rho = np.array([0.1, 0.2, -0.3]), 0.8
    d = lo.decode(lo.reflection(lo.inversion_sphere(ic, rho)) @ lo.point(x))
    xi = ic + rho * rho * (x - ic) / ((x - ic) @ (x - ic))
    expect(np.linalg.norm(d[1] - xi) < 1e-12, "inversion acts")


def dupin():
    R, r = 1.3, 0.4
    d = lo.torus(R, r)
    expect(lo.defect(d) < 1e-12 and lo.group_defect(d) < 1e-12, "torus defect")
    expect(lo.signature(d, 1) == (2, 1) and lo.signature(d, -1) == (2, 1), "torus signature")
    on = 0
    for i in range(12):
        for j in range(9):
            p, n = torus_contact(R, r, 0.5 * i, 0.7 * j + 0.1)
            on = max(on, lo.contact_residual(d, p, n))
    expect(on < 1e-12, "torus membership")
    p, n = torus_contact(R, r, 0.3, 0.9)
    expect(lo.contact_residual(d, p, -n) > 1e-2, "orientation matters")
    for ph in (0.0, 0.6, 1.4, 2.5, np.pi):
        p, n = torus_contact(R, r, 0.8, ph)
        ku, kv, ru, rv = lo.curvatures(d, p, n)
        expect(abs(ku + 1 / r) < 1e-10 and abs(kv + np.cos(ph) / (R + r * np.cos(ph))) < 1e-10, "curvatures")
    grown = lo.transform(lo.normal_offset(0.15), lo.torus(R, r))
    expect(sdist(grown, lo.torus(R, r + 0.15)) < 1e-12, "Laguerre offset grows the tube")
    expect(sdist(lo.transform(lo.scale(2), lo.torus(R, r)), lo.torus(2 * R, 2 * r)) < 1e-12, "scale")
    A = cyclide_map()
    cyc = lo.transform(A, d)
    expect(lo.defect(cyc) < 1e-10, "Moebius image")
    for which, truth in enumerate((d, cyc)):
        ps, ns = [], []
        for i in range(12):
            p, n = torus_contact(R, r, 0.53 * i, 0.31 + 0.77 * i)
            if which:
                p, n = push_contact(A, p, n)
            ps.append(p)
            ns.append(n)
        sig, gap, ok = lo.fit(ps, ns)
        expect(ok and gap < 1e-8 and usdist(sig, truth) < 1e-6, f"fit {which}")
    ps, ns = zip(*[torus_contact(R, r, 0.4 + 0.08 * np.cos(0.9 * i), 1.0 + 0.08 * np.sin(0.9 * i))
                   for i in range(7)])
    sig, _, _ = lo.fit(ps, ns, gap_tol=1e-2)
    expect(usdist(sig, d) < 1e-4, "local patch fit")
    noise = 1e-3 * np.sin(1.0 + 7.0 * np.arange(6)[:, None] + 3.0 * np.arange(6)[None, :])
    nd = lo.normalize(d + noise)
    expect(lo.defect(nd) < 1e-12 and sdist(nd, d) < 1e-2 and sdist(lo.normalize(d), d) < 1e-14, "normalize")
    w = 0.15 * np.cos(0.5 + 1.7 * np.arange(15))
    b = lo.transform(lo.lexp(lo.algebra(w)), d)
    g = lo.transvection(d, b)
    expect(lo.group_defect(g) < 1e-10 and sdist(lo.transform(g, d), b) < 1e-10, "transvection")
    xi = lo.dupin_log(d, b)
    expect(np.abs(d @ xi + xi @ d).max() < 1e-10, "log anticommutes")
    mid = lo.geodesic(d, b, 0.5)
    expect(sdist(lo.transform(mid, d), b) < 1e-9, "geodesic midpoint")
    expect(len(lo.tangent_basis(cyc)) == 9, "tangent dimension")
    ps, ns = zip(*[torus_contact(R, r, 0.61 * i, 0.27 + 0.83 * i) for i in range(10)])
    v = 0.02 * np.sin(2.0 + 1.3 * np.arange(15))
    start = lo.transform(lo.lexp(lo.algebra(v)), d)
    fixed = lo.correct(start, ps, ns, [1.0] * 10, 1e-12, 6)
    expect(sdist(fixed, d) < 1e-7, "correct")


def additions():
    R, r = 1.3, 0.4
    for sig, name in ((lo.torus(R, r), "torus"), (lo.transform(cyclide_map(), lo.torus(R, r)), "cyclide")):
        EU, EV = lo.dupin_conic_basis(sig)
        worst = 0
        for tu in np.linspace(0, 2 * np.pi, 7):
            for tv in np.linspace(0, 2 * np.pi, 7):
                pn = lo.dupin_contact(EU, EV, (tu, tv))
                if pn is not None:
                    worst = max(worst, lo.contact_residual(sig, *pn))
        expect(worst < 1e-10, f"conic parametrization lies on the {name}")
    sig = lo.torus(R, r)
    worst = 0
    for i in range(10):
        p, n = torus_contact(R, r, 0.7 * i, 0.4 + 0.9 * i)
        x = p + 0.05 * n + 0.02 * np.cross(n, [0.3, 0.2, 1.0])
        q, m, _, ok = lo.dupin_land(sig, x, n0=n)
        # Analytic closest point on the torus.
        rho = np.hypot(x[0], x[1])
        c = np.array([R * x[0] / rho, R * x[1] / rho, 0])
        exact = c + r * (x - c) / np.linalg.norm(x - c)
        worst = max(worst, np.linalg.norm(q - exact) if ok else 1.0)
    expect(worst < 1e-8, f"landing on the torus ({worst:.1e})")
    A = lo.rigid(lo.axis_rotation([1, 0.2, 0.3], 0.4), [0.1, -0.3, 0.2]) @ lo.normal_offset(0.1) @ lo.scale(1.2)
    for spread in (0.1, 0.3, 1.0):
        p0, n0, p1, n1 = [], [], [], []
        for i in range(7):
            p, n = torus_contact(R, r, 0.4 + spread * np.cos(0.9 * i), 1.0 + spread * np.sin(1.7 * i))
            q, m = push_contact(A, p, n)
            p0.append(p), n0.append(n), p1.append(q), n1.append(m)
        W, res = lo.lie_procrustes(p0, n0, p1, n1)
        expect(lo.group_defect(W) < 1e-9, "procrustes in group")
        # On cyclide data the motion is defined up to sigma_0; the shape is not.
        expect(sdist(lo.transform(W, sig), lo.transform(A, sig)) < 1e-6,
               f"procrustes transports the cyclide, spread {spread} "
               f"({sdist(lo.transform(W, sig), lo.transform(A, sig)):.1e})")


if __name__ == "__main__":
    vectors_and_contact()
    group()
    dupin()
    additions()
    print("all passed" if not FAILS else f"{len(FAILS)} failed")
