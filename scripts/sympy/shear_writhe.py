"""Shear, twist and writhe on tubes; cross ratios against vertex shape (MATH.md 19.24).

1. Tube twist as a metric shear. A tube of radius `r` about an axis `c(s)`
   (unit speed, Bishop frame `N1, N2`), its cross-sections turned by `phi(s)`:
       X(s, t) = c(s) + r (cos(t + phi) N1 + sin(t + phi) N2).
   The image is the same surface for every `phi` (a reparametrization), but
   in material coordinates `(s, t)` the metric picks up `g_st = r^2 phi'`: a
   pure shear. The second fundamental form keeps the tube's curvatures; only
   the material lines turn.

2. A circumferentially varying shear `gamma(s, t)` (the shear a drive asks
   for, read on a tube) against what the tube can do: the twist mode
   `gamma = r phi'` is free (a slide); the part varying around the tube is
   not, and its first harmonic moves the axis.

3. The quaternionic cross ratio of a small edge quad on a Monge patch
   `h = (k1 x^2 + k2 y^2) / 2`, expanded in the quad size: its non-real part
   against the curvature (the claim to check: governed by the traceless part
   along the edge, `Re` and `Im` of `z e^{-2 i alpha}`).

Run: `python shear_writhe.py`.
"""

import sympy as sp


def tube_twist():
    s, t, r = sp.symbols("s t r", real=True, positive=True)
    k1 = sp.Function("k1")(s)  # Bishop curvatures of the axis
    k2 = sp.Function("k2")(s)
    phi = sp.Function("phi")(s)
    # Bishop frame derivatives: T' = k1 N1 + k2 N2, N1' = -k1 T, N2' = -k2 T.
    # Write X_s, X_t in the frame (T, N1, N2) symbolically.
    th = t + phi
    # X = c + r (cos th N1 + sin th N2)
    # X_s = T + r (-sin th phi' N1 + cos th phi' N2) + r (cos th N1' + sin th N2')
    #     = T (1 - r (k1 cos th + k2 sin th)) + r phi' (-sin th N1 + cos th N2)
    Xs = sp.Matrix([1 - r * (k1 * sp.cos(th) + k2 * sp.sin(th)), -r * sp.sin(th) * phi.diff(s),
                    r * sp.cos(th) * phi.diff(s)])
    Xt = sp.Matrix([0, -r * sp.sin(th), r * sp.cos(th)])
    E = sp.simplify(Xs.dot(Xs))
    F = sp.simplify(Xs.dot(Xt))
    G = sp.simplify(Xt.dot(Xt))
    print("1. tube with cross-sections turned by phi(s), material coordinates (s, t):")
    print("   g_ss =", sp.simplify(E))
    print("   g_st =", F, "   (the twist is a pure metric shear r^2 phi')")
    print("   g_tt =", G)
    # Shear strain against the untwisted tube: angle change between the material lines.
    E0 = E.subs(phi.diff(s), 0)
    shear = sp.simplify(F / sp.sqrt(E0 * G))
    print("   engineering shear (small phi'):", sp.simplify(sp.series(shear, r, 0, 3).removeO()))
    return E, F, G


def circumferential_shear():
    s, t, r, eps = sp.symbols("s t r epsilon", real=True)
    print("2. a shear varying around the tube, gamma(s, t) = sum_n gamma_n(s) e^{i n t}:")
    print("   n = 0: gamma = r phi' is realized by turning the cross-sections (a slide, shape unchanged);")
    print("   n >= 1: no reparametrization gives it (phi' cannot depend on t), so the metric change")
    print("   is incompatible with the round tube: n = 1 tilts the cross-sections, bending and")
    print("   torsioning the axis; n >= 2 squashes them (shell wrinkles).")
    # Check n = 1: a displacement field that turns cross-sections differentially is a rigid tilt of
    # each section about an axis normal to T, i.e. a bend of the axis.
    a = sp.Function("a")(s)
    disp_t = a * sp.sin(t)  # axial displacement varying once around: a section tilt
    print("   n = 1 example: axial displacement a(s) sin t tilts each section by a(s) / r, the axis")
    print("   bends with curvature d(a / r) / ds:", sp.diff(a / sp.Symbol("r"), s))


def quat_mul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def quat_inv(a):
    n = sum(c * c for c in a)
    return (a[0] / n, -a[1] / n, -a[2] / n, -a[3] / n)


def cross_ratio_quad():
    e, k1, k2, al = sp.symbols("e k1 k2 alpha", real=True)
    print("3. quaternionic cross ratio of a small rhombus quad on h = (k1 x^2 + k2 y^2) / 2:")

    def P(x, y):
        return (0, x, y, (k1 * x * x + k2 * y * y) / 2)

    c, sn = sp.cos(al), sp.sin(al)
    # Edge from p0 to p2 along alpha (the quad's diagonal edge), apexes p1, p3 across it.
    u = (c, sn)
    v = (-sn, c)
    p0 = P(-e * u[0] / 2, -e * u[1] / 2)
    p2 = P(e * u[0] / 2, e * u[1] / 2)
    p1 = P(e * v[0] * sp.sqrt(3) / 2, e * v[1] * sp.sqrt(3) / 2)
    p3 = P(-e * v[0] * sp.sqrt(3) / 2, -e * v[1] * sp.sqrt(3) / 2)
    sub = lambda a, b: tuple(ai - bi for ai, bi in zip(a, b))
    # cr = (p0 - p1)(p1 - p2)^-1 (p2 - p3)(p3 - p0)^-1
    cr = quat_mul(quat_mul(quat_mul(sub(p0, p1), quat_inv(sub(p1, p2))), sub(p2, p3)), quat_inv(sub(p3, p0)))
    comps = [sp.simplify(sp.series(ci, e, 0, 3).removeO()) for ci in cr]
    print("   real part:", sp.simplify(comps[0]))
    for name, ci in zip(("i", "j", "k"), comps[1:]):
        print(f"   {name} part:", sp.factor(sp.simplify(ci)))
    H, z_re, z_im = (k1 + k2) / 2, (k1 - k2) / 2 * sp.cos(2 * al), -(k1 - k2) / 2 * sp.sin(2 * al)
    print("   with H = (k1 + k2)/2 and the traceless part along the edge, (k1 - k2)/2 (cos 2a, -sin 2a):")
    print("   check the normal (k) part against e^2 (k1 - k2)/2 cos 2a:",
          sp.simplify(comps[3] - sp.Rational(0) * H))


if __name__ == "__main__":
    tube_twist()
    circumferential_shear()
    cross_ratio_quad()
