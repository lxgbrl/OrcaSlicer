# Custom Scriptable Infill — pattern definitions

This folder holds the data that drives the **Custom (scripted)** sparse infill
pattern (`FillCustomScripted`). New infills are added here as files; no rebuild
of OrcaSlicer is required.

Lookup order for a `custom_infill_pattern_id`:
1. `<user data dir>/custom_infill/`  (your overrides win)
2. `<resources>/custom_infill/`      (this bundled folder)

## 2D tiles (`*.tile`) — used when `custom_infill_mode = 2D Tile`

Plain text, line oriented. Coordinates are tile-local millimetres.

```
TILE <width_mm> <height_mm>     # tile period
MODE serpentine                 # optional: reverse alternate rows for continuity
PATH x1,y1 x2,y2 x3,y3 ...      # one polyline; multiple PATH lines allowed
# lines starting with # are comments
```

The tile is repeated across the region bounding box; lines are clipped to the
region, density culls rows, and the global infill angle rotates the result.

## 3D volumetric patterns (`*.json` / `*.ini`) — used when `custom_infill_mode = 3D Volume`

```json
{ "type": "gyroid", "cell_size": 6.0, "level": 0.0, "thickness": 0.6, "density_scale": 1.0 }
```

`type` is one of `gyroid`, `schwarzp`. The implicit field is sampled on each
layer Z and the iso-contour `f(x,y,z)=level` becomes the infill curves. Print
density scales the effective cell size (denser = finer lattice).

The INI-like form is also accepted:

```
PATTERN gyroid
CELL_SIZE 6.0
LEVEL 0.0
THICKNESS 0.6
DENSITY_SCALE 1.0
```

## Authoring tiles from CAD

See `tools/blender_export_tile.py` and `tools/grasshopper_export_tile.py`.
