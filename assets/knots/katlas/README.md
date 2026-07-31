# Knot Atlas take-home RDF

Hoste–Thistlethwaite knot tables from the [Knot Atlas Take Home Database](https://katlas.org/wiki/The_Take_Home_Database).

Dumps are **gitignored** (~45MB compressed); fetch locally:

```bash
# Download Knots11..Knots15
python3 scripts/knots/fetch_katlas_rdf.py

# Also build ht_index.json (name → Gauss / DT / …)
python3 scripts/knots/fetch_katlas_rdf.py --index

# Optional extras
python3 scripts/knots/fetch_katlas_rdf.py --tables Rolfsen Knots11 --index
```

| File | Size (gz) | Knots (approx.) |
|------|-----------|-----------------|
| `Knots11.rdf.gz` | ~0.5MB | 552 |
| `Knots12.rdf.gz` | ~0.6MB | 2.2k |
| `Knots13.rdf.gz` | ~2.8MB | 10k |
| `Knots14.rdf.gz` | ~10MB | 47k |
| `Knots15.rdf.gz` | ~29MB | 253k |
| `ht_index.json` | regenerable | union of tables |

**Note:** HT dumps generally do **not** include `BraidWord` (Rolfsen does). For `K11a359` etc., use wiki `BR[…]` / `KnotTheory\``. Knotilus renders Gauss codes; it is not this RDF set.
