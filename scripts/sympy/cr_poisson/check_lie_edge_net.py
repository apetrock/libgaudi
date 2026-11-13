"""Checks for `lie_edge_net.py` (MATH.md 19.14)."""

import numpy as np

import cr_ops as ops
import dupin_poisson as dp
import lie_edge_net as len_
import lie_edges as le
import lie_ops as lo
from lie_compare import cases

fails = 0


def expect(ok, what, value=None):
    global fails
    tag = "ok  " if ok else "FAIL"
    if not ok:
        fails += 1
    print(f"  {tag} {what}" + ("" if value is None else f"  ({value:.2e})"))


def rot(axis, ang):
    a = np.asarray(axis, float)
    a /= np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


rng = np.random.default_rng(0)
print("Lie rows")
p, q = rng.normal(size=(5, 3)), rng.normal(size=(5, 3))
n = len_._nrm(rng.normal(size=(5, 3)))
k = rng.normal(size=5)
s = len_.curvature_sphere_rows(p, n, k)
ref = np.array([lo.curvature_sphere(p[r], n[r], k[r]) for r in range(5)])
expect(np.abs(s - ref).max() < 1e-12, "curvature_sphere_rows = lie_ops.curvature_sphere", np.abs(s - ref).max())
dots = len_.lie_dot_rows(s, len_.point_rows(q))
ref = np.array([lo.lie_dot(s[r], lo.point(q[r])) for r in range(5)])
expect(np.abs(dots - ref).max() < 1e-12, "lie_dot_rows = lie_ops.lie_dot", np.abs(dots - ref).max())
d = q - p
tang = np.einsum("ij,ij->i", n, d) - 0.5 * k * np.einsum("ij,ij->i", d, d)
expect(np.abs(dots - tang).max() < 1e-12, "incidence <kP+T, P(x_j)> = n . d - k |d|^2 / 2", np.abs(dots - tang).max())

print("Euler on the contact line")
S = rng.normal(size=(1, 3))
M = np.array([[S[0, 0], S[0, 1]], [S[0, 1], S[0, 2]]])
kk, V = np.linalg.eigh(M)
worst = 0.0
for phi in np.linspace(0, np.pi, 7):
    u = np.cos(phi) * V[:, 0] + np.sin(phi) * V[:, 1]
    kap = len_.tensor_in_dir(S, u[None], u[None])[0]
    s1 = lo.curvature_sphere(p[0], n[0], kk[0])
    s2 = lo.curvature_sphere(p[0], n[0], kk[1])
    s_phi = np.cos(phi) ** 2 * s1 + np.sin(phi) ** 2 * s2
    worst = max(worst, np.abs(s_phi - lo.curvature_sphere(p[0], n[0], kap)).max())
expect(worst < 1e-12, "s(phi) = cos^2 s1 + sin^2 s2 = (u^T S u) P + T", worst)

print("Read-outs and rows on the suite's targets")
for name, X0, F, Xt, dist in cases():
    if not any(t in name for t in ("tube", "sphere", "saddle")):
        continue
    m = ops.build_mesh(F, len(X0))
    hin = le.hinges(m)
    tg = len_.read_spheres(m, hin, Xt)
    net = len_.LieEdgeNet(m, Xt, tg)
    z = net.initial(Xt)
    X, N, H, Sv = net.unpack(z)
    L = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
    blocks = net._blocks(X, N, H, Sv, L)
    names = ["rodrigues", "incidence", "vertex", "couple", "unit", "target", "memory", "hinge"]
    norms = {nm: float(np.sqrt(np.mean((w * f(V)) ** 2))) for nm, (w, f, V, _) in zip(names, blocks)}
    expect(norms["rodrigues"] < 1e-10, f"{name}: Rodrigues exact at the target", norms["rodrigues"])
    expect(norms["hinge"] < 1e-10, f"{name}: hinge exact at the target", norms["hinge"])
    expect(norms["target"] < 1e-12 and norms["unit"] < 1e-12, f"{name}: targets and unit exact")
    print(f"       bias at the target (rms): incidence {norms['incidence']:.2e}  vertex {norms['vertex']:.2e}  "
          f"couple {norms['couple']:.2e}")
    # Jacobian against a directional difference.
    dz = rng.normal(size=net.nvar) * 1e-6
    J = net.jacobian(z, L)
    lin = J @ dz
    fd = net.residual(z + dz, L) - net.residual(z, L)
    rel = np.linalg.norm(lin - fd) / max(np.linalg.norm(fd), 1e-300)
    expect(rel < 1e-3, f"{name}: Jacobian against a directional difference", rel)
    # Rigid invariance: rotate and translate the target, rebuild, same merit.
    Rm = rot([0.3, -0.5, 0.8], 0.9)
    Xr = Xt @ Rm.T + np.array([0.4, -1.1, 2.0])
    netr = len_.LieEdgeNet(m, Xr, len_.read_spheres(m, hin, Xr))
    zr = netr.initial(Xr)
    a, b = net.merit(z, L), netr.merit(zr, L)
    expect(abs(a - b) <= 1e-8 * max(abs(a), 1e-12) + 1e-14, f"{name}: rigid invariance of the merit", abs(a - b))

print("FAILED" if fails else "all passed")
