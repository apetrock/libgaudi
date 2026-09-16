#!/usr/bin/env bash
# Fetch known reference PDFs into reference/ (skips existing files).
# Extend the table when adding papers; keep MANIFEST.md in sync.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$ROOT/reference"
mkdir -p "$DEST"

fetch() {
  local url="$1"
  local out="$2"
  local path="$DEST/$out"
  if [[ -f "$path" ]]; then
    echo "skip (exists): $out"
    return 0
  fi
  echo "fetch: $out"
  curl -fsSL -L -o "$path" "$url"
  ls -la "$path"
}

# filename <- url
fetch "https://arxiv.org/pdf/1910.02926.pdf" "liu2019-cubic-stylization.pdf"
fetch "https://www.dgp.toronto.edu/projects/cubic-stylization/cubicStyle_high.pdf" \
  "liu2019-cubic-stylization-hires.pdf"
fetch "https://ddg.math.uni-goettingen.de/pub/vector-dirichlet.pdf" "stein2020-vector-dirichlet-energy.pdf"
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/DiscreteDevelopable/paper.pdf" \
  "stein2018-developability-triangle-meshes.pdf"
fetch "https://multires.caltech.edu/pubs/GeomElastic.pdf" \
  "chao2010-geometric-elastic-deformations.pdf"
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/RepulsiveCurves.pdf" \
  "crane2021-repulsive-curves.pdf"
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveSurfaces/RepulsiveSurfaces.pdf" \
  "crane2021-repulsive-surfaces.pdf"
fetch "https://markjgillespie.com/Research/harnack-tracing/HarnackTracing.pdf" \
  "gillespie-harnack-tracing.pdf"
fetch "https://markjgillespie.com/Research/integer-coordinates/integer_coordinates.pdf" \
  "gillespie-integer-coordinates.pdf"
fetch "https://nzfeng.github.io/research/PointsAsTori/PointsAsTori.pdf" \
  "feng-points-as-tori.pdf"
fetch "https://minarcik.com/SurfaceUntangling.pdf" \
  "minarcik-surface-untangling.pdf"
fetch "https://repository.gatech.edu/server/api/core/bitstreams/ebf3251d-156b-4cbc-8b49-0317de363ed7/content" \
  "turk1999-variational-implicit-surfaces.pdf"
fetch "https://arxiv.org/pdf/math-ph/0503012" \
  "dennis2005-geometry-calugareanu.pdf"

echo "done. Update reference/MANIFEST.md if you added new files."
