#!/usr/bin/env python3
"""Regenerate assets/knots/braids.json from Stoimenov braid.out + torus knots."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "assets/knots/stoimenov_braid.out"
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
    for p, q in [(5, 2), (7, 2), (8, 3), (5, 4), (8, 7), (9, 7)]:
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
    return [{
        "name": "K11a359",
        "strands": 8,
        "word": k11a359_word,
        "source": "https://katlas.org/wiki/K11a359",
    }]


def sort_key(entry: dict):
    name = entry["name"]
    if name.startswith("T"):
        return (2, name)
    if name.startswith("K"):
        return (1, name)
    a, b = name.split("_", 1)
    return (0, int(a), int(b))


def main() -> None:
    knots = (
        parse_stoimenov(SRC.read_text())
        + torus_knots()
        + hoste_thistlethwaite_extras()
    )
    knots.sort(key=sort_key)
    payload = {
        "source": (
            "Stoimenov braid.out (prime knots through 12 crossings, "
            "minimal braid-index width) + torus knots T{p}_{q} + "
            "Knot Atlas Hoste-Thistlethwaite extras (K11a359, …)"
        ),
        "knots": knots,
    }
    DST.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {DST} ({len(knots)} braids)")


if __name__ == "__main__":
    main()
