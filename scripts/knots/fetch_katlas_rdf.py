#!/usr/bin/env python3
"""Fetch Knot Atlas take-home RDF dumps (Hoste–Thistlethwaite tables).

Source: https://katlas.org/wiki/The_Take_Home_Database
Data:   https://katlas.org/Data/{Rolfsen,Knots11..Knots15,Links,TorusKnots}.rdf.gz

Note: HT Knots11–15 RDF include Gauss/DT/PD etc., but generally not BraidWord
(Rolfsen does). Braid words for HT knots still come from wiki BR[…] / KnotTheory`.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "assets/knots/katlas"
BASE_URL = "https://katlas.org/Data"

DEFAULT_TABLES = [
    "Knots11",
    "Knots12",
    "Knots13",
    "Knots14",
    "Knots15",
]

# Lightweight invariants for HT lookup (skip PD_Presentation — huge HTML).
INDEX_KEYS = (
    "Gauss_Code",
    "DT_Code",
    "Determinant",
    "HyperbolicVolume",
    "Symmetry_Type",
    "BraidWord",
    "BraidIndex",
)

TRIPLE_RE = re.compile(
    r'<knot:([^>]+)>\s*<invariant:([^>]+)>\s*"((?:\\.|[^"\\])*)"'
)


def fetch_table(name: str, dest_dir: Path, force: bool = False) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    gz_path = dest_dir / f"{name}.rdf.gz"
    if gz_path.exists() and not force:
        print(f"skip {gz_path.name} (exists)")
        return gz_path
    url = f"{BASE_URL}/{name}.rdf.gz"
    print(f"fetch {url}")
    urllib.request.urlretrieve(url, gz_path)
    print(f"  -> {gz_path} ({gz_path.stat().st_size} bytes)")
    return gz_path


def parse_rdf_gz(gz_path: Path) -> dict[str, dict[str, str]]:
    """Map knot name -> selected invariants."""
    knots: dict[str, dict[str, str]] = {}
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = TRIPLE_RE.search(line)
            if not m:
                continue
            knot, key, val = m.group(1), m.group(2), m.group(3)
            if key not in INDEX_KEYS:
                continue
            # Unescape minimal RDF/N-Triples style backslashes.
            val = val.replace("\\n", "\n").replace('\\"', '"').replace("\\\\", "\\")
            knots.setdefault(knot, {})[key] = val
    return knots


def build_index(tables: list[str], dest_dir: Path) -> Path:
    index: dict[str, dict] = {
        "source": "https://katlas.org/Data/",
        "tables": {},
        "knots": {},
    }
    for name in tables:
        gz_path = dest_dir / f"{name}.rdf.gz"
        if not gz_path.exists():
            raise FileNotFoundError(gz_path)
        print(f"index {gz_path.name}")
        table_knots = parse_rdf_gz(gz_path)
        index["tables"][name] = {"count": len(table_knots)}
        # Later tables overwrite on name clash (shouldn't happen across 11–15).
        index["knots"].update(table_knots)
        print(f"  {len(table_knots)} knots")
    out = dest_dir / "ht_index.json"
    out.write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {out} ({len(index['knots'])} knots)")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tables",
        nargs="+",
        default=DEFAULT_TABLES,
        help="RDF table basenames (default: Knots11..Knots15)",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if .rdf.gz already present",
    )
    ap.add_argument(
        "--index",
        action="store_true",
        help="Build ht_index.json (name → Gauss/DT/PD/…)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=OUT_DIR,
        help=f"Output directory (default: {OUT_DIR})",
    )
    args = ap.parse_args()

    for name in args.tables:
        fetch_table(name, args.out, force=args.force)
    if args.index:
        build_index(args.tables, args.out)


if __name__ == "__main__":
    main()
