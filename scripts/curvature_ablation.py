"""Curvature-row ablation over `principal_curvature_probe` (docs/TODO.md item 2).

Each curvature row runs alone on top of a base with every curvature row
off (the demo's position rows and coupling as configured), at a sane
weight, for a few frames on the torus target. Per row: the last frame's
R, tube, area, roughness, anisotropy (log r ellipsoidal / saddle), turn,
asymmetry, circle residual and frame time, their change against the base,
and the variational audit of the row's constraint at the first step.

    python3 scripts/curvature_ablation.py [frames] [row ...]
"""

import re
import subprocess
import sys
import time

PROBE = "build-test/tests/principal_curvature_probe"

# Every curvature row off (store, position-side curvature rows, projectors).
BASE_OFF = ("wc=0 wsshear=0 wwm=0 sarap=0 wgrad=0 wgradtl=0 wcmean=0 wcratio=0 wmin=0 wcmc=0 dupin=0 skew=0 "
            "wgrowth=0 wdet=0 wgauss=0 fshear=0 wotorus=0 wodupin=0 wodupinshape=0 wframefair=0 walbers=0 "
            "wkite=0 wdcircle=0").split()

# name -> (flags turning it on, the constraint name(s) in the audit table)
ROWS = {
    "anisotropy": ("wc=0.5", ["edge_shape_fair"]),
    "curvature_shear": ("wsshear=0.5", ["edge_shape_fair"]),
    "willmore": ("wwm=0.5", ["edge_shape_fair"]),
    "shape_arap": ("sarap=0.5", ["edge_shape_fair"]),
    "grad": ("wgrad=0.5", ["edge_shape_fair"]),
    "grad_traceless": ("wgradtl=0.5", ["edge_shape_fair"]),
    "curvature_mean": ("wcmean=0.5", ["edge_shape_fair"]),
    "curvature_ratio": ("wcratio=0.5", ["edge_shape_fair"]),
    "minimizer": ("wmin=0.5", ["edge_shape_fair"]),
    "tension": ("wcmc=0.15", ["edge_shape_fair"]),
    "dupin": ("dupin=0.5", ["edge_shape_fair"]),
    "skew": ("skew=0.5", ["edge_shape_fair"]),
    "growth": ("wgrowth=0.5", ["edge_shape_fair"]),
    "det": ("wdet=0.5", ["ring_det_conserve"]),
    "gauss": ("wgauss=0.5", ["vertex_gauss_conserve"]),
    "face_shear": ("fshear=0.5", ["face_shear", "face_stretch"]),
    "osculating_torus": ("wotorus=0.8", ["osculating_torus_ring"]),
    "osculating_dupin": ("wodupin=0.5", ["osculating_dupin_ring"]),
    "frame_fair": ("wframefair=0.5", ["frame_fair"]),
    "albers": ("walbers=0.5", ["albers_ring_fair"]),
    "dupin_circle": ("wdcircle=0.2", ["dupin_circle"]),
    "kite": ("wkite=0.05", ["edge_strain"]),
}

FRAME = re.compile(
    r"^frame (\d+).* R ([\d.e+-]+) tube ([\d.e+-]+)  area ([\d.e+-]+).*rough ([\d.e+-]+)\s+log r ell ([\d.e+-]+|nan) "
    r"sad ([\d.e+-]+|nan).*turn ([\d.e+-]+|nan).*asym ([\d.e+-]+|nan)\s+dcres ([\d.e+-]+|nan).*wall ([\d.e+-]+) s")
AUDIT = re.compile(r"^\s+(\S+)( \[GN\])? n=(\d+) flagged=(\d+) median=([\d.e+-]+) max=([\d.e+-]+)(?: instances: (.*))?")
KEYS = ["R", "tube", "area", "rough", "logr_ell", "logr_sad", "turn", "asym", "dcres", "wall"]


def run(flags, frames):
    cmd = [PROBE, str(frames), "torus", "1", "20", "audit"] + BASE_OFF + flags
    out = subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()
    last, audit, in_audit, audited = None, {}, False, False
    for line in out:
        m = FRAME.match(line)
        if m:
            last = dict(zip(KEYS, [float(v) for v in m.groups()[1:]]))
        if "variational audit" in line:
            in_audit = not audited
            audited = True
            continue
        if in_audit:
            a = AUDIT.match(line)
            if a:
                inst = [float(v) for v in a.group(7).split()] if a.group(7) else []
                audit[a.group(1)] = (int(a.group(3)), int(a.group(4)), float(a.group(5)), float(a.group(6)), inst)
            elif not line.startswith("  "):
                in_audit = False
    return last, audit


def fmt(v):
    return f"{v:.4g}" if v == v else "nan"


if __name__ == "__main__":
    args = sys.argv[1:]
    frames = int(args[0]) if args and args[0].isdigit() else 6
    names = [a for a in args if not a.isdigit()] or list(ROWS)
    t0 = time.time()
    base, base_audit = run([], frames)
    print(f"base (curvature rows off), {frames} frames: " + "  ".join(f"{k} {fmt(base[k])}" for k in KEYS), flush=True)
    print(f"{'row':18s} " + " ".join(f"{('d' + k):>10s}" for k in KEYS[:-1]) + f" {'wall':>7s}  audit (n, flagged, median)")
    for name in names:
        flag, cons = ROWS[name]
        res, aud = run(flag.split(), frames)
        if res is None:
            print(f"{name:18s} (no frame output)", flush=True)
            continue
        d = [res[k] - base[k] for k in KEYS[:-1]]
        audit_txt = []
        for c in cons:
            if c not in aud:
                continue
            n, fl, med, mx, inst = aud[c]
            b = base_audit.get(c)
            if b and b[0] == n:
                continue  # unchanged from the base: not this row
            if b and inst and b[4]:
                # Drop the base's own instances (matched by value) and report the rest.
                rest = list(inst)
                for v in b[4]:
                    j = min(range(len(rest)), key=lambda k: abs(rest[k] - v)) if rest else None
                    if j is not None and abs(rest[j] - v) <= 1e-3 * max(abs(v), 1e-12):
                        rest.pop(j)
                audit_txt.append(f"{c} row audit {' '.join(f'{v:.2g}' for v in rest) or '?'}")
            else:
                audit_txt.append(f"{c} (n {n}, flagged {fl}, median {med:.2g})")
        print(f"{name:18s} " + " ".join(f"{fmt(v):>10s}" for v in d) + f" {res['wall']:7.2f}  " +
              ("; ".join(audit_txt) or "-"), flush=True)
    print(f"({time.time() - t0:.0f}s)")
