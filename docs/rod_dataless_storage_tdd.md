# Rod Data-less Storage TDD

## Status

**Complete** — [`asawa::rod`](../include/gaudi/asawa/rod/rod.hpp) uses
shell-parity topology + [`asawa::datum`](../include/gaudi/asawa/datums.hpp) pool.
Remesh via `rod_split` / `rod_collapse`; accessors `x()`…`t()`; see
[`rod_datums.hpp`](../include/gaudi/asawa/rod/rod_datums.hpp).

**Upstream:** [`rod_residual_forest_tdd.md`](rod_residual_forest_tdd.md) Phase A.

**Model:** [`shell.hpp`](../include/gaudi/asawa/shell/shell.hpp) —
`insert_datum` / `get_datum` / `__data`; remesh via virtual hooks; pack via
`permute` + `resize`.

---

## Concept

Rod keeps **only** corner topology (`__corners_prev/next`). Attributes live in
typed `asawa::datum` instances. **Interpolation policy = datum subclass.**

| Datum type | Slot use | `rod_split` / `rod_collapse` |
|------------|----------|------------------------------|
| `vec3_catmull_datum` | `x` (0) | Catmull–Rom / midpoint fallback |
| `vec3_datum` | `v`, `M` | midpoint |
| `quat_slerp_datum` | `u`, `o` | slerp |
| `vec4_datum` | `J` | midpoint |
| `real_interp_l_datum` | `l0` | geometric rest-length split |
| `scalar_datum` | `t` | midpoint |

`dynamic` loops:

```cpp
for (auto &d : R.get_data())
  d->rod_split(R, c01, c00, c10, c11, cnew);
```

---

## Conventional slot ids

```cpp
namespace rod_datum_id {
  constexpr index_t x = 0, v = 1, u = 2, o = 3, M = 4, J = 5, l0 = 6, t = 7;
}
```

Accessors `x()`, `v()`, … return references into the corresponding datum.

---

## Tests

See [`rod_datum_storage_tests.hpp`](../include/gaudi/test/rod_datum_storage_tests.hpp):

- pool insert / size match
- split midpoint / slerp / interp_l / catmull
- collapse + pack permute all channels

---

## Out of scope

- Residual forest
- Changing interp math vs current `dynamic` behavior
