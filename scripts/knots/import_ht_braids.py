#!/usr/bin/env python3
"""Braid words for Knot Atlas Hoste–Thistlethwaite knots.

Reads DT codes from assets/knots/katlas/ht_index.json (Knots11–15) and
converts each one with spherogram.Link.braid_word. Writes
assets/knots/ht_braids.json, then regenerates assets/knots/braids.json.

    python3 scripts/knots/fetch_katlas_rdf.py --index
    python3 scripts/knots/import_ht_braids.py

Requires: pip install spherogram
"""

from __future__ import annotations

import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INDEX = ROOT / "assets/knots/katlas/ht_index.json"
DST = ROOT / "assets/knots/ht_braids.json"

_Link = None


def _init_worker() -> None:
    global _Link
    from spherogram import Link

    _Link = Link


def _components(word: list[int]) -> int:
    n = max(abs(g) for g in word) + 1
    pos = list(range(n))
    for g in word:
        i = abs(g) - 1
        pos[i], pos[i + 1] = pos[i + 1], pos[i]
    seen = [False] * n
    cycles = 0
    for s in range(n):
        if seen[s]:
            continue
        cycles += 1
        j = s
        while not seen[j]:
            seen[j] = True
            j = pos[j]
    return cycles


def _convert(item: tuple[str, str]) -> tuple[str, int, list[int], str]:
    name, dt = item
    code = "DT: [(%s)]" % ", ".join(dt.split())
    try:
        word = [int(g) for g in _Link(code).braid_word()]
    except Exception as exc:  # noqa: BLE001 — library errors vary by diagram
        return name, 0, [], f"{type(exc).__name__}: {exc}"
    if not word:
        return name, 0, [], "empty braid word"
    if _components(word) != 1:
        return name, 0, [], "closure is not one component"
    strands = max(abs(g) for g in word) + 1
    return name, strands, word, ""


def main() -> None:
    if not INDEX.exists():
        raise SystemExit(
            f"missing {INDEX}; run scripts/knots/fetch_katlas_rdf.py --index"
        )
    index = json.loads(INDEX.read_text())
    items = [
        (name, rec["DT_Code"])
        for name, rec in index["knots"].items()
        if rec.get("DT_Code")
    ]
    items.sort()
    workers = max(1, min(8, mp.cpu_count() or 1))
    print(f"converting {len(items)} DT codes with spherogram ({workers} workers)")
    t0 = time.perf_counter()
    knots: list[dict] = []
    failed: list[tuple[str, str]] = []
    done = 0
    with mp.Pool(workers, initializer=_init_worker) as pool:
        for name, strands, word, err in pool.imap_unordered(
            _convert, items, chunksize=64
        ):
            done += 1
            if err:
                failed.append((name, err))
            else:
                knots.append(
                    {"name": name, "strands": strands, "word": word}
                )
            if done % 20000 == 0 or done == len(items):
                elapsed = time.perf_counter() - t0
                print(
                    f"  {done}/{len(items)}  ok {len(knots)}  "
                    f"fail {len(failed)}  {elapsed:.0f}s"
                )
    knots.sort(key=lambda k: k["name"])
    payload = {
        "source": (
            "Knot Atlas Hoste-Thistlethwaite DT codes "
            "(assets/knots/katlas/ht_index.json, Knots11–15) converted by "
            "spherogram.Link.braid_word"
        ),
        "knots": knots,
    }
    DST.write_text(json.dumps(payload) + "\n")
    print(f"wrote {DST} ({len(knots)} braids, {len(failed)} failed)")
    for name, err in failed[:20]:
        print(f"  {name}: {err}")
    if len(failed) > 20:
        print(f"  … {len(failed) - 20} more")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from import_stoimenov_braids import main as rebuild_catalog

    rebuild_catalog()


if __name__ == "__main__":
    main()
