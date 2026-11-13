"""Checks for `lie_cyclide_net.py` (MATH.md 19.15)."""

import numpy as np

import cr_ops as ops
import dupin_poisson as dp
import lie_cyclide_net as lcn
import lie_edge_net as lsn
import lie_edges as le
import lie_ops as lo
from lie_compare import cases

fails = 0


def expect(ok, what, value=None):
    global fails
    if not ok:
        fails += 1
    print(f"  {'ok  ' if ok else 'FAIL'} {what}" + ("" if value is None else f"  ({value:.2e})"))


def torus_element(R, r, u, v):
    c = np.array([np.cos(u), np.sin(u), 0.0])
    p = (R + r * np.cos(v)) * c + r * np.sin(v) * np.array([0, 0, 1.0])
    n = np.cos(v) * c + np.sin(v) * np.array([0, 0, 1.0])
    return p, n


R, r = 2.0, 0.7
print("Cyclide fit on the analytic torus")
truth = lo.torus(R, r)
worst = 0.0
for (u, v) in [(0.3, 0.4), (1.1, 2.0), (2.5, -1.2)]:
    p, n = torus_element(R, r, u, v)
    hstep = 0.15
    ring = [torus_element(R, r, u + hstep * np.cos(a), v + hstep * np.sin(a)) for a in np.linspace(0, 2 * np.pi, 7)[:-1]]
    rx, rn = np.array([q[0] for q in ring]), np.array([q[1] for q in ring])
    k_tube, k_ring = -1.0 / r, -np.cos(v) / (R + r * np.cos(v))
    e_tube = -np.sin(v) * np.array([np.cos(u), np.sin(u), 0]) + np.cos(v) * np.array([0, 0, 1.0])
    k1, k2 = (k_tube, k_ring) if abs(k_tube) >= abs(k_ring) else (k_ring, k_tube)
    start = lo.osculating_torus(p, n, k1, k2, e_tube)
    s = lcn.fit_vertex_cyclide(start, p, n, rx, rn, k1, k2, hstep, iterations=8)
    worst = max(worst, dp.field_error(s[None], truth)[0])
expect(worst < 1e-6, "fit recovers lo.torus(2, 0.7) from its contact elements and curvatures", worst)

print("Vectorized contact residual")
pts = [torus_element(R, r, u, v) for u, v in [(0.1, 0.2), (1.0, 2.5), (3.0, -0.7), (4.4, 1.3)]]
d = np.array([q[0] for q in pts])
nn = np.array([q[1] for q in pts])
x0 = np.array([0.3, -0.2, 0.1])
loc = lo.translation(-x0) @ truth @ lo.translation(x0)
res_on = np.abs(lcn.contact_residual_rows(np.repeat(loc[None], 4, 0), d - x0, nn)).max()
expect(res_on < 1e-12, "zero for contact elements on the cyclide", res_on)
res_off = np.abs(lcn.contact_residual_rows(np.repeat(loc[None], 4, 0), d - x0 + 0.05 * nn, nn)).max()
expect(res_off > 1e-3, "nonzero off it (offset 0.05 along the normal)", res_off)

print("Coupling: equivariance under a global Lie transformation")
name, X0, F, Xt, dist = [c for c in cases() if c[0].startswith("torus tube")][0]
m = ops.build_mesh(F, len(Xt))
hin = le.hinges(m)
net = lcn.LieCyclideNet(m, Xt, lsn.read_spheres(m, hin, Xt))
rng = np.random.default_rng(1)
sig = np.stack([lo.normalize(truth + 1e-2 * lo.algebra(rng.normal(size=15)) @ truth) for _ in range(m.n)])
ok = np.ones(m.n, bool)
A = lo.rigid(lo.axis_rotation([0.2, 1, -0.3], 0.7), [0.5, -0.3, 0.8]) @ lo.scale(1.3)
Ai = lo.inverse(A)
a = net.couple(np.einsum("ab,nbc,cd->nad", A, sig, Ai), ok)
b = np.einsum("ab,nbc,cd->nad", A, net.couple(sig, ok), Ai)
eq = np.abs(a - b).max() / np.abs(b).max()
expect(eq < 1e-8, "couple(A sigma A^-1) = A couple(sigma) A^-1", eq)
const = net.couple(np.repeat(truth[None], m.n, 0), ok)
expect(np.abs(const - truth).max() < 1e-10, "a constant field is a fixed point", np.abs(const - truth).max())

print("Joint rows")
z = net.initial(Xt)
nc = net.shape_block(z)
expect(nc > 0.9 * m.n, f"cyclides at the torus target: {nc} of {m.n}")
X, N, H, S = net.unpack(z)
L = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
blocks = net._blocks(X, N, H, S, L)
tens = np.sqrt(np.mean((blocks[-2][0] * blocks[-2][1](blocks[-2][2])) ** 2))
ring = np.sqrt(np.mean((blocks[-1][0] * blocks[-1][1](blocks[-1][2])) ** 2))
print(f"       at the target (rms): tensor rows {tens:.2e}  ring contact rows {ring:.2e}")
dz = rng.normal(size=net.nvar) * 1e-6
J = net.jacobian(z, L)
fd = net.residual(z + dz, L) - net.residual(z, L)
rel = np.linalg.norm(J @ dz - fd) / np.linalg.norm(fd)
expect(rel < 1e-3, "Jacobian against a directional difference", rel)

print("FAILED" if fails else "all passed")
