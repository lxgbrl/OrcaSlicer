"""Grasshopper GhPython component: export curves as an OrcaSlicer .tile file.

Place a GhPython component with these inputs:
    crvs   : List Access, Curve   -- the tile polylines (in a rectangular domain)
    tile_w : Item Access, float   -- tile width in mm  (e.g. 10)
    tile_h : Item Access, float   -- tile height in mm (e.g. 10)
    path   : Item Access, str     -- output .tile file path
    serp   : Item Access, bool    -- write MODE serpentine
    run    : Item Access, bool    -- set True to write the file

Curves are discretised into polylines; points are emitted tile-local (origin at
the min-XY corner of the combined bounding box) in millimetres.
"""

import Rhino.Geometry as rg
import os

TOL = 0.01  # mm; polyline discretisation tolerance


def _curve_to_points(crv):
    pl = crv.ToPolyline(TOL, 0.5, 0.0, 1e9) if hasattr(crv, "ToPolyline") else None
    if pl is None:
        # fallback: try direct polyline conversion
        ok, poly = crv.TryGetPolyline()
        if ok:
            return [(pt.X, pt.Y) for pt in poly]
        # sample by division
        params = crv.DivideByCount(64, True)
        return [(crv.PointAt(t).X, crv.PointAt(t).Y) for t in params]
    return [(pt.X, pt.Y) for pt in pl]


def run_export(crvs, tile_w, tile_h, path, serp):
    polys = [_curve_to_points(c) for c in crvs if c is not None]
    polys = [p for p in polys if len(p) >= 2]
    if not polys:
        return "No curves."

    xs = [x for p in polys for (x, _) in p]
    ys = [y for p in polys for (_, y) in p]
    min_x, min_y = min(xs), min(ys)
    w = tile_w if tile_w and tile_w > 0 else (max(xs) - min_x)
    h = tile_h if tile_h and tile_h > 0 else (max(ys) - min_y)

    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(path, "w") as f:
        f.write("# Exported from Grasshopper by grasshopper_export_tile.py\n")
        f.write("TILE %g %g\n" % (w, h))
        if serp:
            f.write("MODE serpentine\n")
        for pts in polys:
            coords = " ".join("%g,%g" % (x - min_x, y - min_y) for (x, y) in pts)
            f.write("PATH " + coords + "\n")
    return "Wrote %d paths to %s" % (len(polys), path)


# GhPython entry point
if "run" in dir() and run:
    info = run_export(crvs, tile_w, tile_h, path, serp)
else:
    info = "Set run = True to export."
