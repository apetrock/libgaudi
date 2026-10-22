#!/usr/bin/env python3
"""Regenerate assets/knots/braids.json from Stoimenov braid.out, torus knots, and ht_braids.json."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "assets/knots/stoimenov_braid.out"
HT = ROOT / "assets/knots/ht_braids.json"
DST = ROOT / "assets/knots/braids.json"

PAT = re.compile(r"^\s*(\d+)\s+(\d+)\s+\{([^}]*)\}\s*$", re.M)


def parse_stoimenov(text: str) -> list[dict]:
    knots: list[dict] = []
    for m in PAT.finditer(text):
        n, k, body = int(m.group(1)), int(m.group(2)), m.group(3)
        word = [int(x.strip()) for x in body.split(",") if x.strip()]
        strands = max(abs(g) for g in word) + 1
        knots.append({"name": f"{n}_{k}", "strands": strands, "word": word})
    return knots


def torus_knots() -> list[dict]:
    out: list[dict] = []
    # (p, q) with gcd 1. Closure of (σ1…σ_{q-1})^p is T(p,q) on q strands.
    # The last four are denser than the 12-crossing table (8–14 strands).
    for p, q in [
        (5, 2), (7, 2), (8, 3), (5, 4), (8, 7), (9, 7),
        (9, 8), (11, 10), (13, 12), (15, 14),
    ]:
        if math.gcd(p, q) != 1:
            continue
        word: list[int] = []
        for _ in range(p):
            word.extend(range(1, q))
        out.append({"name": f"T{p}_{q}", "strands": q, "word": word})
    return out


def hoste_thistlethwaite_extras() -> list[dict]:
    """Named HT knots from Knot Atlas (not Stoimenov serial indices)."""
    # https://katlas.org/wiki/K11a359 — BR[8, {...}]
    k11a359_word = [
        1, -2, -3, 4, 5, 4, 3, 2, -1, -3, 2, 4, -5, 6, 5, 4, 4, 3, 4, 7, -6, -5,
        4, -3, -2, 4, 5, 4, 3, 4, 6, -5, 4, -7, -6,
    ]
    # Conway (K11n34) and Kinoshita–Terasaka (K11n42). Knot Atlas draws these
    # as BraidPart tiles, not a BR[...] list. The tile map reproduces
    # 3_1 = [-1,-1,-1] and 4_1 = [-1,2,-1,2]; both closures are one component
    # and match the Jones polynomial printed on their atlas pages.
    k11n34_word = [1, 1, 2, -3, 2, 1, -3, -2, -2, -3, -3]
    k11n42_word = [1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2]
    return [
        {
            "name": "K11a359",
            "strands": 8,
            "word": k11a359_word,
            "source": "https://katlas.org/wiki/K11a359",
        },
        {
            "name": "K11n34",
            "strands": 4,
            "word": k11n34_word,
            "source": "https://katlas.org/wiki/K11n34",
        },
        {
            "name": "K11n42",
            "strands": 4,
            "word": k11n42_word,
            "source": "https://katlas.org/wiki/K11n42",
        },
    ]


def load_ht_braids() -> list[dict]:
    """Spherogram braid words for Knot Atlas HT names. Empty if not generated."""
    if not HT.exists():
        return []
    data = json.loads(HT.read_text())
    return list(data["knots"])


def sort_key(entry: dict):
    name = entry["name"]
    ht = re.fullmatch(r"K(\d+)([an])(\d+)", name)
    if ht:
        return (1, int(ht.group(1)), ht.group(2), int(ht.group(3)))
    if name.startswith("T"):
        return (2, name)
    if name.startswith("K"):
        return (1, 0, name, 0)
    a, b = name.split("_", 1)
    return (0, int(a), int(b), 0)


def main() -> None:
    ht = load_ht_braids()
    ht_names = {k["name"] for k in ht}
    extras = [k for k in hoste_thistlethwaite_extras() if k["name"] not in ht_names]
    knots = parse_stoimenov(SRC.read_text()) + torus_knots() + extras + ht
    knots.sort(key=sort_key)
    payload = {
        "source": (
            "Stoimenov braid.out (prime knots through 12 crossings, "
            "minimal braid-index width) + torus knots T{p}_{q} "
            "(including dense T9_8, T11_10, T13_12, T15_14) + "
            "Knot Atlas Hoste–Thistlethwaite Knots11–15 via spherogram "
            f"({len(ht)} braid words)"
        ),
        "knots": knots,
    }
    # Pretty-print the small table. The HT dump is hundreds of thousands of
    # words; indent would roughly triple the file the demo parses at startup.
    indent = None if len(knots) > 10000 else 2
    DST.write_text(json.dumps(payload, indent=indent) + "\n")
    print(f"wrote {DST} ({len(knots)} braids)")


if __name__ == "__main__":
    main()
