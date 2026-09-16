#!/usr/bin/env python3
"""Convert triangle meshes to OBJ for asawa loaders.

GLB/GLTF: trimesh (gpytoolbox does not read glTF).
OBJ/STL/PLY: gpytoolbox.read_mesh / write_mesh.

  pip install gpytoolbox trimesh
  python scripts/assets/convert_mesh.py INPUT.glb [-o assets/models/obj/NAME.obj]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_DIR = ROOT / "assets" / "models" / "obj"

GLTF_EXTS = {".glb", ".gltf"}
GPY_EXTS = {".obj", ".stl", ".ply"}


def _require(mod: str):
    try:
        return __import__(mod)
    except ImportError as e:
        sys.exit(f"Missing dependency '{mod}'. Try: pip install gpytoolbox trimesh\n{e}")


def load_vf(path: Path):
    ext = path.suffix.lower()
    if ext in GLTF_EXTS:
        trimesh = _require("trimesh")
        loaded = trimesh.load(path, force="scene")
        if isinstance(loaded, trimesh.Scene):
            mesh = trimesh.util.concatenate(
                tuple(g for g in loaded.geometry.values() if hasattr(g, "faces"))
            )
        else:
            mesh = loaded
        if mesh is None or len(getattr(mesh, "faces", [])) == 0:
            sys.exit(f"No triangle geometry in {path}")
        return mesh.vertices, mesh.faces
    if ext in GPY_EXTS:
        gpytoolbox = _require("gpytoolbox")
        return gpytoolbox.read_mesh(str(path))
    sys.exit(f"Unsupported input extension '{ext}' (want glb/gltf/obj/stl/ply)")


def write_obj(path: Path, V, F) -> None:
    gpytoolbox = _require("gpytoolbox")
    path.parent.mkdir(parents=True, exist_ok=True)
    gpytoolbox.write_mesh(str(path), V, F)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path, help="Input mesh path")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .obj (default: assets/models/obj/<stem>.obj)",
    )
    args = ap.parse_args()
    inp = args.input.expanduser().resolve()
    if not inp.is_file():
        sys.exit(f"Input not found: {inp}")

    out = args.output
    if out is None:
        out = DEFAULT_OUT_DIR / f"{inp.stem}.obj"
    else:
        out = out.expanduser()
        if not out.is_absolute():
            out = ROOT / out
        out = out.resolve()

    V, F = load_vf(inp)
    write_obj(out, V, F)
    print(f"wrote {out}  V={len(V)} F={len(F)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
