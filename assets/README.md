# Assets

## Layout

- `models/obj/` — triangle meshes for asawa loaders (`bunny.obj`, `vhm_skull.obj`, …)
- `knots/` — braid / knot catalog text & JSON (stays in main repo)

Loaders use repo-relative paths, e.g. `assets/models/obj/bunny.obj`
(`include/gaudi/asawa/shell/asset_loader.hpp`).

## Convert meshes

```bash
pip install gpytoolbox trimesh
python scripts/assets/convert_mesh.py INPUT.glb -o assets/models/obj/NAME.obj
```

GLB/GLTF via trimesh; OBJ write via gpytoolbox. See `scripts/assets/convert_mesh.py`.

## Future: submodule for bulky models

`assets/models/` is intended to become a **git submodule** (`libgaudi-assets`) so large OBJs are not diffed in libgaudi history. Until that remote exists, models live in-tree under `models/obj/` (gitignored `*.obj` with `!assets/models/obj/**` exceptions).
