#include "../ClipperUtils.hpp"
#include "../MarchingSquares.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include "../TriangleMesh.hpp"
#include "../TriangleMeshSlicer.hpp"

#include <cmath>
#include <algorithm>

#include "FillBase.hpp"
#include "FillCustomScripted.hpp"
#include "CustomInfillPattern.hpp"

// ---------------------------------------------------------------------------
// 3D implicit-field marching-squares scaffolding.
//
// Generalizes the gyroid field used in FillGyroid.cpp to a small catalog of
// TPMS families. For a given layer Z we evaluate f(x,y,z) on a raster and
// extract the iso-contour f == level, which becomes the layer's infill curves.
//
//   Gyroid:    sin(fx)cos(fy) + sin(fy)cos(fz) + sin(fz)cos(fx)
//   Schwarz P: cos(fx) + cos(fy) + cos(fz)
// ---------------------------------------------------------------------------
namespace marchsq {
using namespace Slic3r;

using coordr_t = long;
using Pointf   = Vec2d;

struct TpmsField
{
    coord_t                rsize;       // raster step in scaled coords
    coordr_t               gsize;       // marching-squares window (in raster cells)
    Point                  size;
    Point                  offs;
    coordf_t               z;
    float                  freq;        // angular frequency = 2*pi / cell (periodic families)
    float                  isoval;
    VolumeType             type;
    double                 cell;        // cell size (mm), for Expr tiling
    double                 dmin, dmax;  // Expr per-cell coordinate domain
    double                 tval;        // Expr animation parameter
    const ImplicitExpr*    expr = nullptr;

    TpmsField(const BoundingBox bb, const coordf_t z, const double cell, const double level,
              const VolumePattern& vp)
        : size{bb.size()}, offs{bb.min}, z{z}, type{vp.type}, cell{cell},
          dmin{vp.domain_min}, dmax{vp.domain_max}, tval{vp.t}, expr{vp.expr.get()}
    {
        freq   = float(2.0 * PI) / float(std::max(cell, 1e-3));
        isoval = float(level);
        // Raster resolution: fine fixed step for the analytic families; for Expr
        // tie it to the cell so heavy formulas stay affordable (~24 samples/cell).
        const float rsizef = (type == VolumeType::Expr)
            ? float(std::max(cell / 24.0, 0.05)) : 0.004f;
        rsize = scaled(rsizef);
        gsize = std::max<coordr_t>(2, coordr_t(std::round(0.40f / rsizef)));
    }

    // Map a world coordinate (mm) into the formula's per-cell domain for Expr.
    inline double map_expr(double w) const
    {
        const double c = cell > 1e-6 ? cell : 1.0;
        const double loc = w - c * std::floor(w / c);
        return dmin + (loc / c) * (dmax - dmin);
    }

    float get_scalar(coordf_t x, coordf_t y, coordf_t z_arg) const
    {
        if (type == VolumeType::Expr && expr) {
            double v = expr->eval(map_expr(x), map_expr(y), map_expr(z_arg), tval);
            if (! std::isfinite(v)) v = isoval + 1.0; // push non-finite outside the iso band
            return float(v);
        }
        const float a = freq * float(x);
        const float b = freq * float(y);
        const float c = freq * float(z_arg);
        switch (type) {
        case VolumeType::SchwarzP:
            return std::cos(a) + std::cos(b) + std::cos(c);
        case VolumeType::Gyroid:
        default:
            return std::sin(a) * std::cos(b) + std::sin(b) * std::cos(c) + std::sin(c) * std::cos(a);
        }
    }

    float get_scalar(Coord p) const
    {
        Pointf pf = to_Pointf(p);
        return get_scalar(pf.x(), pf.y(), z);
    }

    inline coord_t  to_coord (const coordr_t& x) const { return x * rsize; }
    inline coordr_t to_coordr(const coord_t& x)  const { return x / rsize; }
    inline Point  to_Point (const Coord& p) const { return Point(to_coord(p.c) + offs.x(), to_coord(p.r) + offs.y()); }
    inline Pointf to_Pointf(const Point& p) const { return Pointf(unscaled(p.x()), unscaled(p.y())); }
    inline Pointf to_Pointf(const Coord& p) const { return to_Pointf(to_Point(p)); }
};

template<> struct _RasterTraits<TpmsField>
{
    using ValueType = float;
    static float  get (const TpmsField& sf, size_t row, size_t col) { return sf.get_scalar(Coord(row, col)); }
    static size_t rows(const TpmsField& sf) { return sf.to_coordr(sf.size.y()); }
    static size_t cols(const TpmsField& sf) { return sf.to_coordr(sf.size.x()); }
};

inline Polylines get_tpms_polylines(const TpmsField& sf, const double tolerance = SCALED_EPSILON)
{
    std::vector<Ring> rings = execute_with_policy(ex_tbb, sf, sf.isoval, {sf.gsize, sf.gsize});
    Polylines polys;
    polys.reserve(rings.size());
    for (const Ring& ring : rings) {
        Polyline poly;
        Points&  pts = poly.points;
        pts.reserve(ring.size() + 1);
        for (const Coord& crd : ring)
            pts.emplace_back(sf.to_Point(crd));
        if (pts.size() < 2)
            continue;
        pts.push_back(pts.front());
        if (tolerance >= 0.0)
            poly.simplify(tolerance);
        polys.emplace_back(poly);
    }
    return polys;
}

} // namespace marchsq

namespace Slic3r {

Polylines FillCustomScripted::fill_tile_2d(const FillParams& params, const ExPolygon& expolygon, const BoundingBox& bb)
{
    const PrintRegionConfig* cfg = params.config;
    const std::string id   = cfg ? cfg->custom_infill_pattern_id.value : std::string("rectilinear");
    const double      tw_mm = cfg ? cfg->custom_infill_tile_width.value  : 10.0;
    const double      th_mm = cfg ? cfg->custom_infill_tile_height.value : 10.0;

    const TilePattern& tile = PatternManager::instance().get_tile_pattern(id, tw_mm, th_mm);
    if (! tile.valid || tile.paths.empty())
        return {};

    const coord_t tw = std::max<coord_t>(scale_(tile.width),  scale_(0.1));
    const coord_t th = std::max<coord_t>(scale_(tile.height), scale_(0.1));

    // Density -> row culling: keep one row every `step` rows.
    const double density = std::clamp(double(params.density), 0.01, 1.0);
    const int    step    = std::max(1, int(std::lround(1.0 / density)));

    // Number of tiles needed to cover the (expanded) bounding box.
    const coord_t x0 = bb.min.x();
    const coord_t y0 = bb.min.y();
    const int cols = int((bb.size().x() + tw - 1) / tw) + 1;
    const int rows = int((bb.size().y() + th - 1) / th) + 1;

    Polylines out;
    out.reserve(size_t(cols) * rows * tile.paths.size());
    for (int r = 0; r < rows; ++r) {
        if (r % step != 0)
            continue; // density culling
        const bool reverse_row = tile.serpentine && (r & 1);
        for (int c = 0; c < cols; ++c) {
            const coord_t dx = x0 + coord_t(c) * tw;
            const coord_t dy = y0 + coord_t(r) * th;
            for (const Polyline& src : tile.paths) {
                Polyline pl = src;
                pl.translate(dx, dy);
                if (reverse_row)
                    pl.reverse();
                out.emplace_back(std::move(pl));
            }
        }
    }
    return out;
}

Polylines FillCustomScripted::fill_volume_3d(const FillParams& params, const ExPolygon& /*expolygon*/, const BoundingBox& bb)
{
    const PrintRegionConfig* cfg = params.config;
    const std::string id = cfg ? cfg->custom_infill_pattern_id.value : std::string("gyroid");

    const VolumePattern& vp = PatternManager::instance().get_volume_pattern(id);

    if (vp.type == VolumeType::Mesh)
        return fill_mesh_3d(params, bb);

    // Base cell: GUI "Custom cell size" overrides the pattern file when set (>0).
    double base_cell = vp.cell_size;
    if (cfg && cfg->custom_infill_volume_cell.value > 0.0)
        base_cell = cfg->custom_infill_volume_cell.value;

    // Density -> lattice fineness. Denser prints use a smaller effective cell so
    // more iso-contours cross the region. Bounded to keep raster cost sane.
    const double density = std::clamp(double(params.density), 0.01, 1.0);
    const double scale_f = std::clamp(std::sqrt(density * std::max(vp.density_scale, 1e-3) * 2.0), 0.25, 4.0);
    const double eff_cell = std::max(base_cell / scale_f, 0.2);

    const double level = (cfg ? double(cfg->custom_infill_level.value) : vp.level);

    marchsq::TpmsField sf(bb, this->z, eff_cell, level, vp);
    return marchsq::get_tpms_polylines(sf, SCALED_SPARSE_INFILL_RESOLUTION);
}

Polylines FillCustomScripted::fill_mesh_3d(const FillParams& params, const BoundingBox& bb)
{
    const PrintRegionConfig* cfg = params.config;
    const std::string id = cfg ? cfg->custom_infill_pattern_id.value : std::string();
    const VolumePattern& vp = PatternManager::instance().get_volume_pattern(id);
    if (vp.type != VolumeType::Mesh || ! vp.mesh)
        return {};

    double cell = (cfg && cfg->custom_infill_volume_cell.value > 0.0)
        ? cfg->custom_infill_volume_cell.value : vp.cell_size;
    if (cell <= 1e-3) cell = 10.0;

    const double mx0 = vp.mesh_min.x(), my0 = vp.mesh_min.y(), mz0 = vp.mesh_min.z();
    const double mw = std::max(vp.mesh_max.x() - mx0, 1e-6);
    const double mh = std::max(vp.mesh_max.y() - my0, 1e-6);
    const double md = std::max(vp.mesh_max.z() - mz0, 1e-6);

    // World layer Z -> mesh-local Z, tiling vertically with period = cell.
    const double localz = this->z - cell * std::floor(this->z / cell);
    double mesh_z = mz0 + (localz / cell) * md;
    mesh_z = std::clamp(mesh_z, mz0 + md * 1e-4, mz0 + md * (1.0 - 1e-4));

    std::vector<ExPolygons> slices = slice_mesh_ex(vp.mesh->its, std::vector<float>{ float(mesh_z) }, [](){});
    if (slices.empty() || slices[0].empty())
        return {};

    // Scale mesh XY (mm) into a cell of side `cell` (mm). Build the cell-local
    // polylines once, then tile them across the region bounding box.
    const double sx = cell / mw;
    const double sy = cell / mh;
    const coord_t cw = scale_(cell);
    const coord_t ch = scale_(cell);

    Polylines cell_polys;
    auto add_contour = [&](const Polygon& poly) {
        Polyline pl;
        pl.points.reserve(poly.points.size() + 1);
        for (const Point& p : poly.points)
            pl.points.emplace_back(Point(scale_((unscaled(p.x()) - mx0) * sx),
                                         scale_((unscaled(p.y()) - my0) * sy)));
        if (! pl.points.empty()) pl.points.push_back(pl.points.front());
        if (pl.points.size() >= 2) cell_polys.emplace_back(std::move(pl));
    };
    for (const ExPolygon& ex : slices[0]) {
        add_contour(ex.contour);
        for (const Polygon& h : ex.holes) add_contour(h);
    }
    if (cell_polys.empty())
        return {};

    const coord_t x0 = bb.min.x(), y0 = bb.min.y();
    const int cols = int((bb.size().x() + cw - 1) / cw) + 1;
    const int rows = int((bb.size().y() + ch - 1) / ch) + 1;

    Polylines out;
    out.reserve(size_t(cols) * rows * cell_polys.size());
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) {
            const coord_t dx = x0 + coord_t(c) * cw;
            const coord_t dy = y0 + coord_t(r) * ch;
            for (const Polyline& src : cell_polys) {
                Polyline pl = src;
                pl.translate(dx, dy);
                out.emplace_back(std::move(pl));
            }
        }
    return out;
}

void FillCustomScripted::_fill_surface_single(
    const FillParams                &params,
    unsigned int                     /*thickness_layers*/,
    const std::pair<float, Point>   &/*direction*/,
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    const PrintRegionConfig* cfg = params.config;
    const CustomInfillMode mode = cfg ? cfg->custom_infill_mode.value : cimTile2D;

    // Apply the infill rotation by rotating the sampling domain (spec 9.1), then
    // rotate the resulting polylines back. This keeps both the tile and implicit
    // pipelines axis-aligned internally.
    const float infill_angle = this->angle;
    if (std::abs(infill_angle) >= EPSILON)
        expolygon.rotate(-infill_angle);

    BoundingBox bb = expolygon.contour.bounding_box();
    // Expand to avoid edge artifacts (matches FillGyroid).
    const coord_t expand = 10 * coord_t(scale_(this->spacing));
    bb.offset(expand);

    Polylines polylines = (mode == cimVolume3D)
        ? fill_volume_3d(params, expolygon, bb)
        : fill_tile_2d(params, expolygon, bb);

    // Clip to the actual infill region.
    polylines = intersection_pl(std::move(polylines), expolygon);

    if (! polylines.empty()) {
        // Drop tiny fragments, but keep lines that bridge thin walls.
        const double minlength = scale_(0.8 * this->spacing);
        polylines.erase(
            std::remove_if(polylines.begin(), polylines.end(),
                [minlength](const Polyline& pl) { return pl.length() < minlength; }),
            polylines.end());
    }

    if (! polylines.empty()) {
        const size_t first_idx = polylines_out.size();
        chain_or_connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);
        // Rotate generated paths back into model space.
        if (std::abs(infill_angle) >= EPSILON)
            for (auto it = polylines_out.begin() + first_idx; it != polylines_out.end(); ++it)
                it->rotate(infill_angle);
    }
}

} // namespace Slic3r
