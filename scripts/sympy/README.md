# Medial shape search — sympy derive + C++ emit

Flat sympy modules; bundling lives in C++ (`include/gaudi/albers/`).

## Setup

```bash
cd scripts/sympy
python3 -m venv .venv && source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Usage

```bash
# Inspect sympy (before C++)
python3 generate.py --inspect
python3 generate.py --inspect E
python3 generate.py --inspect --cse E_prime

# Composed march graph
python3 generate.py --inspect-march quadric

# Ray profile dump (bidirectional march along ∇D)
python3 generate.py --profile-ray sphere --R 2.5 --t-max 2.5 --steps 25
python3 generate.py --profile-ray torus

# Emit C++ headers → generated/
python3 generate.py --emit kernels
python3 generate.py --emit quadric
python3 generate.py --emit darboux
python3 generate.py --emit slab_cap_solid_angle
python3 slab_cap_solid_angle.py --fd-check
python3 generate.py --emit all

python3 -m pytest test.py -v
```

## Energy

Default medial objective (minimize along ray `x(t) = p0 + t·∇D/|∇D|`, `t ∈ ℝ`):

```text
W  = (I - nnᵀ) H (I - nnᵀ) / |g|     pure Weingarten map (no σ)
s₂ = trace(W²) = κ₁² + κ₂²
E  = 1 / (s₂ + ε)
```

Ray is along `∇D/|∇D|`; medial may lie at `t* < 0` or `t* > 0` (bidirectional). Tests validate `|t*|`, medial position, and κ signature at `x(t*)`.

`ε` stabilizes `E` in flat stiff regions. Near negative singularities (`|g|→0`), search uses **hybrid damped Newton + bracket fallback + backoff** — not mollification of `W`.

Human-readable math spec: [`ila/medial_shape_energy.ila`](ila/medial_shape_energy.ila).

Symbolic bundles are disk-cached under `.cache/` (gitignored). Delete `.cache/` after changing `compose.py` or `shape_fcn.py`.

## Files

| File | Role |
|------|------|
| `darboux.py`, `quadric.py` | Layer 1: D, G, H sympy |
| `slab_field_energy.py` | Dual-field bulk energy `∫|∇(Q_r−Q_k)|²` + SymPy grad |
| `slab_cap_solid_angle.py` | Cap wedge Ω + barrier fairness energy + SymPy grad → C++ |
| `torus_analytic.py` | Parametric torus ground truth |
| `line.py` | f, fp, fpp |
| `shape_fcn.py` | W(g, H) — pure Weingarten |
| `compose.py` | Single DAG + jet: W → 1/trace(W²) |
| `eigen_printer.py` | Matrix-aware sympy → Eigen C++ |
| `eval_ray.py` | Numeric march + lambdify (`ray_from_foot`) |
| `search.py` | Hybrid Newton (damped + bracket + backoff) |
| `medial_expect.py` | κ inspection at medial point |
| `cache_utils.py` | Disk pickle cache |
| `profile_ray.py` | Ray profile table dump |
| `albers_printer.py` | C++ emit in albers style |
| `generate.py` | CLI |
| `ila/medial_shape_energy.ila` | Math documentation (optional I❤️LA parity) |

Layer 2: `shape_operator(g,H)` → `E = 1/(trace(W²)+ε)`. Composed graph in `compose.py` fans out to lambdify, structured `EigenPrinter` emit, and tests.

## Incompressible φ on sphere caps (Python iteration)

Discrete wedge-basis prototype matching C++ (`phi_div_free_basis` + quadric gauge):

```bash
cd scripts/sympy
python3 slab_phi_sphere_caps.py --separation 0.08 --nu-wedge 200 --gd-iters 40
python3 slab_phi_sphere_caps.py --sweep-nu --separation 0.08
python3 slab_phi_sphere_caps.py --sweep-mu --nu-wedge 200 --separation 0.08
python3 slab_phi_sphere_caps.py --two-phase --nu-wedge 200 --mu-gauge 1
```

Energy: `E = mu * E_gauge + nu * E_wedge` with `E_gauge = sum |f_r(phi(x))|^2` and `E_wedge = sum w_k (tau_k - tau_k^rest)^2`.

**Findings (2026-08):** `E_gauge` alone reduces field error but **increases** wedge volume error under cap separation. Adding `nu * E_wedge` (~200 on unit sphere) restores `tau` and `|V|` but `max|f_r|` stays ~0.6–0.7 — iso-surface gauge and discrete div-free volume memory conflict in the current basis. If joint tuning fails, revisit Hirota et al. volume-preserving FFD (`docs/slab_incompressible_phi_material.md`, DOI 10.1145/304012.304036).
