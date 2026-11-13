"""The edge-sphere solver (`lie_edges.py`) on the shared suite: `edge_compare`'s
cases, the closed band, and a flat sheet to a saddle; same metrics as
section 16 (`live_compare.shape_error`, `closed_targets.edge_strain`)."""

import sys
import time

import numpy as np

import closed_targets as ct
import cr_ops as ops
import lie_edges as le
from dupin_compare import all_cases
from live_compare import shape_error


def saddle_case():
    P, F = ops.strip_param(13, 13, 1.6, 1.6)
    Xf = ops.flat_from_param(P)
    # The flat strip lies in y = 0 with -y normals; the target is
    # y = -a (x^2 - z^2) / 2 over the same parameters.
    a = 0.8
    Xs = Xf.copy()
    Xs[:, 1] = -0.5 * a * (P[:, 0] ** 2 - P[:, 1] ** 2)
    return ("saddle a 0.8", Xf, F, Xs, None)


def cases():
    return all_cases() + [saddle_case()]


def run_case(name, X0, F, Xt, dist, cfg, verbose=False, solver_kind="edges", net_cfg=None):
    m = ops.build_mesh(F, len(X0))
    size = np.ptp(Xt, axis=0).max()
    h = le.hinges(m)
    tgt = le.read_targets(m, h, Xt)
    if solver_kind == "lsn":
        import lie_edge_net as lsn
        solver = lsn.LieEdgeNet(m, X0, lsn.read_spheres(m, h, Xt), net_cfg)
        fx = solver.run(Xt.copy(), iterations=10, keep_best=False)
    elif solver_kind == "lcn":
        import lie_cyclide_net as lcn
        import lie_edge_net as lsn
        tg = lsn.read_spheres(m, h, Xt)
        fx = lcn.LieCyclideNet(m, Xt, tg, net_cfg).solve(Xt.copy())
        solver = lcn.LieCyclideNet(m, X0, tg, net_cfg)
    elif solver_kind in ("cnet", "cnet_gn"):
        import contact_net as cn
        import dupin_poisson as dp
        ring = dp.rings(m)
        pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        if solver_kind == "cnet_gn":
            solver = cn.ContactNetGN(m, X0, cn.read_targets(m, pairs, Xt), net_cfg, pairs,
                                     refresh=getattr(net_cfg, "gn_refresh", 1))
        else:
            solver = cn.ContactNetSolver(m, X0, cn.read_targets(m, pairs, Xt), net_cfg, pairs)
        fx = solver.run(Xt.copy(), iterations=10, keep_best=False)
    elif solver_kind == "net":
        import dupin_poisson as dp
        import lie_net as ln
        body = dp.pose_free_targets(m, Xt, mode="corrected", frames="principal")
        net_cfg.edge = cfg
        solver = ln.NetSolver(m, X0, body, tgt, net_cfg)
        fx = solver.run(Xt.copy(), iterations=5, keep_best=False)
    else:
        solver = le.EdgeSphereSolver(m, X0, tgt, cfg)
        fx = solver.run(Xt.copy(), iterations=10, keep_best=False)

    def log(k, X, move, mer):
        if verbose and k % 20 == 0:
            print(f"    it {k:4d} move {move:.1e} merit {mer:.2e} shape {shape_error(X, Xt, dist, size):.2e}",
                  flush=True)

    if solver_kind in ("cnet_gn", "lsn"):
        r = solver.run(X0.copy(), iterations=net_cfg.iterations, log=log)
    elif solver_kind == "lcn":
        r = solver.solve(X0.copy(), log=log)
        print(f"    cyclide vertices per sweep {r.cyclides} of {m.n}; merit per sweep "
              f"{' '.join(f'{q[1]:.2e}' for q in r.history)}")
    else:
        r = solver.run(X0.copy(), log=log)
    best_it = int(np.argmin([q[1] for q in r.history]))
    if hasattr(r, "backsolves"):
        h = np.array([q[1] for q in r.history])
        firsts = []
        for t in (1e-3, 1e-6, 1e-9):
            i = np.nonzero(h < t * h[0])[0]
            firsts.append(str(i[0]) if len(i) else "-")
        print(f"    iterates to merit/start 1e-3 / 1e-6 / 1e-9: {' / '.join(firsts)}   "
              f"factorizations {getattr(r, 'factorizations', 1)}  backsolves {r.backsolves}")
    print(f"{name:22s} its {len(r.history) - 1:4d} (best {best_it:4d})  shape {shape_error(X0, Xt, dist, size):.1e} -> "
          f"{shape_error(r.X, Xt, dist, size):.1e}  merit {r.history[0][1]:.1e} -> {r.best_merit:.1e}  "
          f"(fixed point shape {shape_error(fx.X, Xt, dist, size):.1e})  strain {ct.edge_strain(m, r.X, Xt).mean():.3f}  "
          f"t {r.seconds:.1f}s ({r.seconds / max(len(r.history) - 1, 1) * 1e3:.0f} ms/it)", flush=True)
    return r


if __name__ == "__main__":
    only = [a for a in sys.argv[1:] if not a.startswith("-")]
    verbose = "-v" in sys.argv
    cfg = le.EdgeConfig()
    import contact_net as cn
    import lie_net as ln
    solver_kind = "edges"
    for a in sys.argv[1:]:
        if a.startswith("--solver="):
            solver_kind = a.split("=")[1]
    if solver_kind == "lsn":
        import lie_edge_net as lsn
        net_cfg = lsn.EdgeNetConfig()
    elif solver_kind == "lcn":
        import lie_cyclide_net as lcn
        net_cfg = lcn.CyclideConfig()
    else:
        net_cfg = cn.NetConfig() if solver_kind in ("cnet", "cnet_gn") else ln.NetConfig()
    for a in sys.argv[1:]:
        if a.startswith("--"):
            k, v = a[2:].split("=")
            if k == "solver":
                continue
            else:
                target, name = (net_cfg, k[4:]) if k.startswith("net.") else (cfg, k)
                cur = getattr(target, name)
                setattr(target, name, v in ("1", "true", "True") if isinstance(cur, bool) else type(cur)(v))
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        run_case(name, X0, F, Xt, dist, cfg, verbose, solver_kind, net_cfg)
