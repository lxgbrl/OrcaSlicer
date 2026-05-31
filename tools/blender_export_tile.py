"""Export selected Blender curves/edges as an OrcaSlicer custom-infill .tile file.

Usage (inside Blender's Scripting workspace, or `blender --python this.py`):
    1. Model a flat pattern in the XY plane inside a rectangle of the desired
       tile size (e.g. 10 x 10 mm). Work in millimetres (Blender units = mm).
    2. Select the curve/mesh objects that make up the pattern.
    3. Set TILE_W / TILE_H below (or leave 0 to auto-fit the selection bbox).
    4. Run. The .tile file is written to OUT_PATH.

Each curve spline / mesh edge-loop becomes one PATH line. Coordinates are
emitted tile-local (origin at the selection's min-XY corner) in millimetres.
"""

import bpy
import os

# ---- config ----------------------------------------------------------------
OUT_PATH = os.path.expanduser("~/orca_custom_infill/my_pattern.tile")
TILE_W = 0.0          # mm; 0 = auto from selection bounding box
TILE_H = 0.0          # mm
MODE_SERPENTINE = True
# ----------------------------------------------------------------------------


def _polylines_from_object(obj):
    """Return list of point-lists [(x,y), ...] in world XY (mm)."""
    polys = []
    mw = obj.matrix_world
    if obj.type == "CURVE":
        for spline in obj.data.splines:
            pts = []
            if spline.type == "BEZIER":
                for bp in spline.bezier_points:
                    co = mw @ bp.co
                    pts.append((co.x, co.y))
            else:
                for p in spline.points:
                    co = mw @ p.co.to_3d()
                    pts.append((co.x, co.y))
            if spline.use_cyclic_u and pts:
                pts.append(pts[0])
            if len(pts) >= 2:
                polys.append(pts)
    elif obj.type == "MESH":
        mesh = obj.data
        for edge in mesh.edges:
            a = mw @ mesh.vertices[edge.vertices[0]].co
            b = mw @ mesh.vertices[edge.vertices[1]].co
            polys.append([(a.x, a.y), (b.x, b.y)])
    return polys


def main():
    selected = [o for o in bpy.context.selected_objects if o.type in {"CURVE", "MESH"}]
    if not selected:
        raise RuntimeError("Select at least one curve or mesh object.")

    polys = []
    for obj in selected:
        polys.extend(_polylines_from_object(obj))
    if not polys:
        raise RuntimeError("No polylines extracted from selection.")

    xs = [x for p in polys for (x, _) in p]
    ys = [y for p in polys for (_, y) in p]
    min_x, min_y = min(xs), min(ys)
    w = TILE_W if TILE_W > 0 else (max(xs) - min_x)
    h = TILE_H if TILE_H > 0 else (max(ys) - min_y)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        f.write("# Exported from Blender by blender_export_tile.py\n")
        f.write("TILE %g %g\n" % (w, h))
        if MODE_SERPENTINE:
            f.write("MODE serpentine\n")
        for pts in polys:
            coords = " ".join("%g,%g" % (x - min_x, y - min_y) for (x, y) in pts)
            f.write("PATH " + coords + "\n")
    print("Wrote %d paths to %s (tile %g x %g mm)" % (len(polys), OUT_PATH, w, h))


if __name__ == "__main__":
    main()
