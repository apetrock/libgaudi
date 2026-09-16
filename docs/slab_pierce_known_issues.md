# Slab pierce — known issues

Status: **contact reconstruction fixed** (null-frame `(s,u,w)` + axial contact point). Remaining items are tuning / visualization only.

## Fixed (2026-03)

Ray+slop triangle hits now return disk contact `(u*, w*, s*)` in the null frame and rebuild

`hit.point = ray.origin + (s* + h)·dir + u*·e1 + w*·e2`

instead of projecting `ray.at(s)` onto the triangle. Bunny audit: `slab_ok=34834/34834`, `pierce_fail=0`, `avg_travel≈2.8e-4`.

## Symptom (historical)

On closed meshes (bunny), some feet get a pierce **hit** (`slab_ok`, cyan segment in medial demo) that **visually skips the nearest inward surface**. The cast cylinder registers contacts on nearby caps (good `\|n_face · n_in\|`), but the chosen `hit.point` can lie on a distant portion of a tied triangle, so the foot→pierce chord tunnels through local geometry.

This is **not** a `pierce_miss` / BVH miss — the heap contains valid local entries.

## Suspected mechanism

1. Fat-cylinder entry assigns the **same** `t_ray` to several one-ring faces (t-tie at min-heap pop).
2. `ray_slop_intersects_triangle` sets `hit.point` via `va::distance_from_triangle(tri, ray.at(t), …)`, not the disk–triangle contact in 3D. On large shallow triangles, that projection can land far from the foot while `t` stays minimal.
3. Min-heap pop on `t_ray` alone picks one tied face; heap order can select the geometrically wrong sibling.

Example (centered bunny, `fat_ray_radius_scale=0.05`): foot **31619** — pop[0] travel **0.505** vs pop[1] at the same `t` travel **1.7e-5**.

## Outlier audit (recommended debug)

Run:

```bash
./tests/pierce_audit outliers          # centered bunny (tests)
./tests/pierce_audit outliers demo     # demo normalization (medial viz)
```

For each pierced foot, report:

| Metric | Meaning |
|--------|---------|
| `foot` | Vertex index |
| `travel` | `\|hit.point − foot\|` |
| `ray_gap` | Distance from `hit.point` to inward cast ray axis (`arp::ray_distance_to_point`) |
| `surface_gap` | BVH nearest triangle distance (`arp::test_point_tri`) — sanity that point is on mesh |
| `claim_gap` | Distance from `hit.point` to the **claimed** hit triangle |
| `t_tie_ratio` | First valid heap pop travel ÷ min valid travel among hits at min `t` |

Sort by `t_tie_ratio` and `ray_gap` to find the major outlier index without visual inspection.

## Related code

- [`slab_pierce_stencil.hpp`](../include/gaudi/hepworth/constraints/slab/slab_pierce_stencil.hpp) — heap pop + self filter
- [`ray_simplex.hpp`](../include/gaudi/arp/ray_simplex.hpp) — slop `hit.point` projection
- [`pairwise_tests.hpp`](../include/gaudi/arp/pairwise_tests.hpp) — `test_point_tri` / BVH nearest
- [`pierce_audit_main.cpp`](../tests/pierce_audit_main.cpp) — `outliers` mode
