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
| `quadric.py`, `darboux.py` | Layer 1: D, G, H sympy |
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
