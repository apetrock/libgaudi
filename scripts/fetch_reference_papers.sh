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
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/RepulsiveCurves_Supplemental.pdf" \
  "crane2021-repulsive-curves-supplemental.pdf"
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
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/ConformalWillmoreFlow/paper.pdf" \
  "crane2013-conformal-willmore-flow.pdf"
fetch "https://www.yousufsoliman.com/projects/download/ConstrainedWillmoreSurfaces.pdf" \
  "soliman2021-constrained-willmore-surfaces.pdf"
fetch "https://www3.math.tu-berlin.de/geometrie/Lehre/SS12/MathVisII/Pinkall%20-%201985%20-%20Hopf%20tori%20in%20S3.pdf" \
  "pinkall1985-hopf-tori.pdf"
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/StripePatterns/paper.pdf" \
  "knoeppel2015-stripe-patterns.pdf"
fetch "https://diglib.eg.org/bitstream/handle/10.2312/8240/eigensatz.pdf?sequence=1" \
  "eigensatz2010-constraint-based-surface-processing.pdf"

fetch "https://www.cse.msu.edu/~ytong/DiscreteFundamentalForms.pdf" \
  "wang2012-discrete-fundamental-forms.pdf"

fetch "https://arxiv.org/pdf/1805.05715" \
  "peng2018-anderson-acceleration.pdf"

fetch "https://dgp.toronto.edu/~silex/Publications/ACWK04/ACWK04b.pdf" \
  "angelidis2006-swirling-sweepers.pdf"
fetch "https://graphics.stanford.edu/courses/cs468-07-winter/Papers/fts-vfbsd-06.pdf" \
  "vonfunck2006-vector-field-shape-deformations.pdf"
fetch "http://evasion.imag.fr/Publications/2002/CN02/pg.pdf" \
  "combaz2002-expansion-textures.pdf"
fetch "https://biomechanics.stanford.edu/me338/me338_project04.pdf" \
  "rodriguez1994-finite-growth.pdf"
fetch "https://pure.uva.nl/ws/files/3391622/26478_Merks.pdf" \
  "merks2003-branching-coral-thesis.pdf"

fetch "https://www.cs.cmu.edu/~kmcrane/Projects/VectorHeatMethod/paper.pdf" \
  "sharp2019-vector-heat-method.pdf"
fetch "https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/paper.pdf" \
  "crane2010-trivial-connections.pdf"

fetch "https://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/paper.pdf" \
  "knoeppel2013-globally-optimal-direction-fields.pdf"

echo "done. Update reference/MANIFEST.md if you added new files."
