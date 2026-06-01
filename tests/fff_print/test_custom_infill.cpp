#include <catch2/catch_all.hpp>

#include <fstream>
#include <cstdio>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Fill/Fill.hpp"
#include "libslic3r/Fill/CustomInfillPattern.hpp"
#include "libslic3r/Fill/ImplicitExpr.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/libslic3r.h"

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

// A 50x50 mm square region to fill.
static ExPolygon square_region()
{
    ExPolygon ex;
    ex.contour = Polygon({ Point(scale_(0.), scale_(0.)),
                           Point(scale_(50.), scale_(0.)),
                           Point(scale_(50.), scale_(50.)),
                           Point(scale_(0.), scale_(50.)) });
    return ex;
}

// total_length(const Polylines&) is provided by libslic3r (Polyline.hpp).

TEST_CASE("CustomInfill: factory creates FillCustomScripted", "[CustomInfill]")
{
    std::unique_ptr<Fill> filler(Fill::new_from_type(ipCustomScripted));
    REQUIRE(filler != nullptr);
    REQUIRE(filler->clone() != nullptr);
}

TEST_CASE("CustomInfill: .tile parser", "[CustomInfill]")
{
    const std::string path = "/tmp/_orca_test.tile";
    {
        std::ofstream o(path);
        o << "# comment\nTILE 8 4\nMODE serpentine\nPATH 0,2 8,2\nPATH 0,0 0,4\n";
    }
    TilePattern tp;
    REQUIRE(PatternManager::parse_tile_file(path, tp));
    REQUIRE(tp.valid);
    REQUIRE_THAT(tp.width,  WithinAbs(8.0, 1e-6));
    REQUIRE_THAT(tp.height, WithinAbs(4.0, 1e-6));
    REQUIRE(tp.serpentine);
    REQUIRE(tp.paths.size() == 2);
    REQUIRE(tp.paths[0].points.size() == 2);
    std::remove(path.c_str());
}

TEST_CASE("CustomInfill: JSON volume config parser", "[CustomInfill]")
{
    const std::string path = "/tmp/_orca_test_vol.json";
    {
        std::ofstream o(path);
        o << R"({"type":"schwarzp","cell_size":4.5,"level":0.1,"thickness":0.7,"density_scale":1.3})";
    }
    VolumePattern vp;
    REQUIRE(PatternManager::parse_volume_config(path, vp));
    REQUIRE(vp.type == VolumeType::SchwarzP);
    REQUIRE_THAT(vp.cell_size, WithinAbs(4.5, 1e-6));
    REQUIRE_THAT(vp.level,     WithinAbs(0.1, 1e-6));
    std::remove(path.c_str());
}

TEST_CASE("CustomInfill: 2D tile mode produces polylines inside region", "[CustomInfill]")
{
    std::unique_ptr<Fill> filler(Fill::new_from_type(ipCustomScripted));
    filler->spacing = 0.45;
    filler->angle   = 0.f;

    PrintRegionConfig cfg;
    cfg.custom_infill_mode.value       = cimTile2D;
    cfg.custom_infill_pattern_id.value = "rectilinear"; // falls back to built-in tile if file absent
    cfg.custom_infill_tile_width.value  = 10.0;
    cfg.custom_infill_tile_height.value = 2.0;

    FillParams fp;
    fp.density = 0.5f;
    fp.config  = &cfg;

    ExPolygon region = square_region();
    Surface surface(stInternal, region);
    Polylines out = filler->fill_surface(&surface, fp);

    REQUIRE_FALSE(out.empty());
    REQUIRE(total_length(out) > scale_(10.0)); // some meaningful path length
    // all points within the 50x50 region bbox (allow small epsilon)
    BoundingBox bb = region.contour.bounding_box();
    for (const Polyline& p : out)
        for (const Point& pt : p.points) {
            REQUIRE(pt.x() >= bb.min.x() - SCALED_EPSILON);
            REQUIRE(pt.x() <= bb.max.x() + SCALED_EPSILON);
            REQUIRE(pt.y() >= bb.min.y() - SCALED_EPSILON);
            REQUIRE(pt.y() <= bb.max.y() + SCALED_EPSILON);
        }
}

TEST_CASE("CustomInfill: ImplicitExpr basic formula", "[CustomInfill]")
{
    ImplicitExpr e;
    std::string err;
    REQUIRE(e.compile_formula("x*x + y*y + z*z - 1", {}, err));
    REQUIRE(e.valid());
    REQUIRE_THAT(e.eval(0, 0, 0), WithinAbs(-1.0, 1e-9));
    REQUIRE_THAT(e.eval(1, 0, 0), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(e.eval(1, 1, 1), WithinAbs(2.0, 1e-9));
}

TEST_CASE("CustomInfill: ImplicitExpr consts, builtins, if/min", "[CustomInfill]")
{
    ImplicitExpr e;
    std::string err;
    REQUIRE(e.compile_formula("if(x<0, min(k, y), max(k, y))", {"k=3/2"}, err));
    REQUIRE_THAT(e.eval(-1, 5, 0), WithinAbs(1.5, 1e-9)); // x<0 -> min(1.5,5)
    REQUIRE_THAT(e.eval( 1, 0, 0), WithinAbs(1.5, 1e-9)); // x>=0 -> max(1.5,0)
}

TEST_CASE("CustomInfill: ImplicitExpr MathMod function composition + redefinition", "[CustomInfill]")
{
    // Schwarz(R(x),R(y),R(z)); R rebinds first arg to x. Then redefine to add 1.
    ImplicitExpr e;
    std::string err;
    std::vector<std::string> funct = {
        "R = 2*x",
        "Schwarz = cos(x)+cos(y)+cos(z)",
        "F = Schwarz(R(x,y,z,t), R(y,x,z,t), R(z,x,y,t), t)",
        "F = F(x,y,z,t) + 1"
    };
    REQUIRE(e.compile({}, funct, "F(x,y,z,t)", err));
    // F = cos(2x)+cos(2y)+cos(2z) + 1 ; at 0 -> 3 + 1 = 4
    REQUIRE_THAT(e.eval(0, 0, 0), WithinAbs(4.0, 1e-9));
}

TEST_CASE("CustomInfill: 3D expr volume mode produces polylines", "[CustomInfill]")
{
    std::unique_ptr<Fill> filler(Fill::new_from_type(ipCustomScripted));
    filler->spacing = 0.45;
    filler->z       = 1.0;

    PrintRegionConfig cfg;
    cfg.custom_infill_mode.value        = cimVolume3D;
    cfg.custom_infill_pattern_id.value  = "expr_gyroid"; // bundled type:expr preset (or fallback)
    cfg.custom_infill_volume_cell.value = 6.0;

    FillParams fp;
    fp.density = 0.3f;
    fp.config  = &cfg;

    ExPolygon region = square_region();
    Surface surface(stInternal, region);
    Polylines out = filler->fill_surface(&surface, fp);
    // expr_gyroid resolves only if resources are present; engine falls back to a
    // default gyroid otherwise. Either way we expect a non-empty result.
    REQUIRE_FALSE(out.empty());
}

TEST_CASE("CustomInfill: mesh cell loads from JSON", "[CustomInfill]")
{
    // Write a tiny 10mm ASCII-STL cube and a JSON referencing it.
    const std::string stl = "/tmp/_orca_cell.stl";
    const std::string js  = "/tmp/_orca_cell_mesh.json";
    {
        std::ofstream o(stl);
        auto tri = [&](double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz){
            o << "facet normal 0 0 0\nouter loop\n";
            o << "vertex " << ax << " " << ay << " " << az << "\n";
            o << "vertex " << bx << " " << by << " " << bz << "\n";
            o << "vertex " << cx << " " << cy << " " << cz << "\n";
            o << "endloop\nendfacet\n";
        };
        o << "solid c\n";
        // 12 triangles of a 0..10 cube
        double v[8][3]={{0,0,0},{10,0,0},{10,10,0},{0,10,0},{0,0,10},{10,0,10},{10,10,10},{0,10,10}};
        int f[12][3]={{0,3,2},{0,2,1},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
        for (auto& t : f) tri(v[t[0]][0],v[t[0]][1],v[t[0]][2], v[t[1]][0],v[t[1]][1],v[t[1]][2], v[t[2]][0],v[t[2]][1],v[t[2]][2]);
        o << "endsolid c\n";
    }
    { std::ofstream o(js); o << R"({"type":"mesh","file":"_orca_cell.stl","cell_size":10})"; }

    VolumePattern vp;
    REQUIRE(PatternManager::parse_volume_config(js, vp));
    REQUIRE(vp.type == VolumeType::Mesh);
    REQUIRE(vp.mesh != nullptr);
    REQUIRE_THAT(vp.mesh_max.x() - vp.mesh_min.x(), WithinAbs(10.0, 1e-3));
    REQUIRE_THAT(vp.mesh_max.z() - vp.mesh_min.z(), WithinAbs(10.0, 1e-3));
    std::remove(stl.c_str());
    std::remove(js.c_str());
}

TEST_CASE("CustomInfill: 3D gyroid volume mode produces polylines", "[CustomInfill]")
{
    std::unique_ptr<Fill> filler(Fill::new_from_type(ipCustomScripted));
    filler->spacing = 0.45;
    filler->angle   = 0.f;
    filler->z       = 1.0; // layer Z in unscaled mm

    PrintRegionConfig cfg;
    cfg.custom_infill_mode.value       = cimVolume3D;
    cfg.custom_infill_pattern_id.value = "gyroid"; // falls back to default gyroid params
    cfg.custom_infill_volume_cell.value = 6.0;

    FillParams fp;
    fp.density = 0.3f;
    fp.config  = &cfg;

    ExPolygon region = square_region();
    Surface surface(stInternal, region);
    Polylines out = filler->fill_surface(&surface, fp);

    REQUIRE_FALSE(out.empty());
    REQUIRE(total_length(out) > scale_(10.0));
}
