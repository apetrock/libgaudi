#!/usr/bin/env python3
"""Replay the Albers minor-tube landing on the demo torus.

The mesh is load_torus(1, 0.33, 48, 24). Curvatures are the analytic ones
for the outward normal (positive curvature sits on -n). The seating and the
circle-versus-line choice match osculating_torus_from_shape_operator and
choose_minor_tube_landing.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = 1.0
MINOR = 0.33
NU = 48
NV = 24
DU = 2.0 * np.pi / NU
DV = 2.0 * np.pi / NV
K_LO_EPS = 1e-12
K_HI_EPS = 1e-8

OUT = Path(__file__).with_suffix(".png")


def surface(u: float, v: float) -> np.ndarray:
    rho = R + MINOR * np.cos(v)
    return np.array([rho * np.cos(u), rho * np.sin(u), MINOR * np.sin(v)])


def canal(u: float) -> np.ndarray:
    return np.array([R * np.cos(u), R * np.sin(u), 0.0])


def frame(u: float, v: float) -> dict:
    """Outward normal. Tube curvature is +1/r. Parallel curvature changes sign."""
    p = surface(u, v)
    n = (p - canal(u)) / MINOR
    e_u = np.array([-np.sin(u), np.cos(u), 0.0])
    e_v = np.array(
        [-np.sin(v) * np.cos(u), -np.sin(v) * np.sin(u), np.cos(v)]
    )
    denom = R + MINOR * np.cos(v)
    k_tube = 1.0 / MINOR
    k_par = np.cos(v) / denom if abs(denom) > 1e-14 else 0.0
    return {
        "p": p,
        "n": n,
        "e_u": e_u,
        "e_v": e_v,
        "k_tube": k_tube,
        "k_par": k_par,
    }


def seat(u: float, v: float) -> dict:
    """Same fields as osculating_torus_from_shape_operator."""
    fr = frame(u, v)
    k_tube = fr["k_tube"]
    k_par = fr["k_par"]
    if abs(k_tube) >= abs(k_par):
        k_hi, k_lo = k_tube, k_par
        e_lo = fr["e_u"]
    else:
        k_hi, k_lo = k_par, k_tube
        e_lo = fr["e_v"]
    out = {
        "valid": False,
        "tube": False,
        "fit": "none",
        "k_hi": k_hi,
        "k_lo": k_lo,
        "p": fr["p"],
        "n": fr["n"],
    }
    if abs(k_hi) <= K_HI_EPS:
        return out
    radius = 1.0 / abs(k_hi)
    sgn = 1.0 if k_hi >= 0.0 else -1.0
    s_prod = 1.0 if k_tube * k_par >= 0.0 else -1.0
    n = fr["n"]
    ruling = e_lo - n * np.dot(e_lo, n)
    ruling_n = np.linalg.norm(ruling)
    if ruling_n <= 1e-8:
        return out
    ruling = ruling / ruling_n
    spine = fr["p"] + sgn * radius * n
    out.update(
        {
            "tube": True,
            "r": radius,
            "s_prod": s_prod,
            "radial": -sgn * n,
            "ruling": ruling,
            "spine": spine,
            "R": 0.0,
            "center": np.zeros(3),
        }
    )
    if abs(k_lo) <= K_LO_EPS:
        return out
    inv_lo = 1.0 / abs(k_lo)
    major = max(0.0, inv_lo - s_prod * radius)
    axis = np.cross(n, e_lo)
    axis_n = np.linalg.norm(axis)
    if axis_n <= 1e-8 or not np.isfinite(major):
        return out
    center = fr["p"] + sgn * s_prod * inv_lo * n
    radial = spine - center
    radial_n = np.linalg.norm(radial)
    out["valid"] = True
    out["R"] = major
    out["axis"] = axis / axis_n
    out["center"] = center
    if radial_n > 1e-8:
        out["radial"] = radial / radial_n
    return out


def nearest_circle(tor: dict, p: np.ndarray) -> np.ndarray | None:
    if not tor["valid"] or not tor["tube"]:
        return None
    v = p - tor["spine"]
    tang = float(np.dot(v, tor["ruling"]))
    radial_off = float(np.dot(v, tor["radial"]))
    along = tor["R"] + radial_off
    rho = float(np.hypot(along, tang))
    if not np.isfinite(rho):
        return None
    if rho <= 1e-18:
        return tor["spine"].copy()
    scale = tor["R"] / rho
    denom = along + rho
    sag = (scale * tang * tang) / denom if abs(denom) > 1e-18 else 0.0
    return tor["spine"] - sag * tor["radial"] + (scale * tang) * tor["ruling"]


def nearest_line(tor: dict, p: np.ndarray) -> np.ndarray | None:
    if not tor["tube"]:
        return None
    v = p - tor["spine"]
    return tor["spine"] + float(np.dot(v, tor["ruling"])) * tor["ruling"]


def land(ps: np.ndarray | None, radius: float, p: np.ndarray) -> np.ndarray | None:
    if ps is None or not (radius > 1e-12):
        return None
    ray = p - ps
    q = float(np.linalg.norm(ray))
    if q <= 1e-12 or not np.isfinite(q):
        return None
    return ps + radius * (ray / q)


def choose(tor: dict, spokes: list[np.ndarray]) -> tuple[str, list[np.ndarray | None]]:
    """Smaller sum of squared moves wins. A tie keeps the circle."""
    n = len(spokes)
    tor_hit: list[np.ndarray | None] = [None] * n
    cyl_hit: list[np.ndarray | None] = [None] * n
    tor_l2 = 0.0
    cyl_l2 = 0.0
    reach = 0.0
    sep2 = 0.0
    both = 0
    for i, p in enumerate(spokes):
        reach = max(reach, float(np.linalg.norm(p)))
        if tor["valid"]:
            hit = land(nearest_circle(tor, p), tor.get("r", 0.0), p)
            if hit is not None:
                tor_hit[i] = hit
        if tor["tube"]:
            hit = land(nearest_line(tor, p), tor.get("r", 0.0), p)
            if hit is not None:
                cyl_hit[i] = hit
        if tor_hit[i] is not None and cyl_hit[i] is not None:
            tor_l2 += float(np.sum((tor_hit[i] - p) ** 2))
            cyl_l2 += float(np.sum((cyl_hit[i] - p) ** 2))
            sep2 = max(sep2, float(np.sum((tor_hit[i] - cyl_hit[i]) ** 2)))
            both += 1
    major_is_flat = (
        tor["valid"]
        and both > 0
        and reach > 1e-8
        and tor["R"] > reach
        and sep2 <= reach * reach * 1e-4
    )
    want_cyl = tor["tube"] and (
        not tor["valid"] or major_is_flat or (both > 0 and cyl_l2 < tor_l2)
    )
    if want_cyl:
        return "cylinder", cyl_hit
    if tor["valid"]:
        return "torus", tor_hit
    return "none", [None] * n


def torus_abs(p: np.ndarray) -> float:
    radial = float(np.hypot(p[0], p[1]))
    return abs(float(np.hypot(radial - R, p[2])) - MINOR)


def section_rho(p: np.ndarray) -> float:
    return float(np.copysign(np.hypot(p[0], p[1]), p[0]))


def main() -> None:
    rows = []
    for j in range(NV):
        v = DV * j
        tor = seat(0.0, v)
        spokes = [
            surface(DU, v),
            surface(-DU, v),
            surface(0.0, v + DV),
            surface(0.0, v - DV),
        ]
        fit, hits = choose(tor, spokes)
        hole = hits[0]
        moved = float(np.linalg.norm(hole - spokes[0])) if hole is not None else np.nan
        off = torus_abs(hole) if hole is not None else np.nan
        spine = tor["spine"] if tor["tube"] else None
        canal_p = canal(0.0)
        spine_err = (
            float(np.linalg.norm(spine - canal_p)) if spine is not None else np.nan
        )
        foot_to_spine = (
            float(np.linalg.norm(spine - tor["p"])) if spine is not None else np.nan
        )
        rows.append(
            {
                "v": v,
                "fit": fit,
                "k_lo": tor["k_lo"],
                "foot_to_spine": foot_to_spine,
                "spine_err": spine_err,
                "moved": moved,
                "off": off,
                "hit": hole,
                "spine": spine,
                "foot": tor["p"],
            }
        )

    print(f"demo torus R={R} r={MINOR}  nu={NU} nv={NV}")
    print(
        f"{'v/pi':>8} {'fit':>10} {'k_par':>10} {'spine':>8} "
        f"{'spine_err':>10} {'|hit-p|':>10} {'off torus':>10}"
    )
    for row in rows:
        print(
            f"{row['v'] / np.pi:8.3f} {row['fit']:>10} {row['k_lo']:10.4f} "
            f"{row['foot_to_spine']:8.4f} {row['spine_err']:10.4f} "
            f"{row['moved']:10.6f} {row['off']:10.6f}"
        )
    worst = max(rows, key=lambda row: row["off"] if np.isfinite(row["off"]) else -1.0)
    n_cyl = sum(row["fit"] == "cylinder" for row in rows)
    print(
        f"cylinder wins {n_cyl}/{NV} tube samples; "
        f"worst off-surface {worst['off']:.6f} at v/pi={worst['v'] / np.pi:.3f} "
        f"({worst['fit']})"
    )

    theta = np.linspace(0.0, 2.0 * np.pi, 256)
    true_rho = R + MINOR * np.cos(theta)
    true_z = MINOR * np.sin(theta)
    hit_rho = [section_rho(row["hit"]) for row in rows]
    hit_z = [row["hit"][2] for row in rows]
    spine_rho = [section_rho(row["spine"]) for row in rows]
    spine_z = [row["spine"][2] for row in rows]
    colors = ["C1" if row["fit"] == "cylinder" else "C0" for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))

    ax = axes[0]
    ax.plot(true_rho, true_z, color="0.35", lw=1.5, label="true section")
    ax.scatter(hit_rho, hit_z, c=colors, s=22, zorder=3, label="landed +u neighbor")
    ax.plot(spine_rho, spine_z, color="C3", lw=1.2, label="recovered spine")
    ax.scatter([R], [0.0], c="k", s=28, zorder=4, label="true canal")
    ax.set_aspect("equal")
    ax.set_xlabel("radial")
    ax.set_ylabel("z")
    ax.set_title("one landing, u = 0")
    ax.legend(fontsize=8, loc="best")

    ax = axes[1]
    vv = [row["v"] / np.pi for row in rows]
    off = [row["off"] for row in rows]
    ax.scatter(vv, off, c=colors, s=22)
    ax.plot(vv, off, color="0.7", lw=0.8, zorder=0)
    ax.set_xlabel("v / pi")
    ax.set_ylabel("distance of landed neighbor off the torus")
    ax.set_title("around-the-hole spoke")
    cyl = matplotlib.lines.Line2D(
        [0], [0], marker="o", color="none", markerfacecolor="C1", label="cylinder"
    )
    tor_mark = matplotlib.lines.Line2D(
        [0], [0], marker="o", color="none", markerfacecolor="C0", label="torus"
    )
    ax.legend(handles=[tor_mark, cyl], fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT, dpi=140)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
