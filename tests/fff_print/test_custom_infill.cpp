#include <catch2/catch_all.hpp>

#include <fstream>
#include <cstdio>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Fill/Fill.hpp"
#include "libslic3r/Fill/CustomInfillPattern.hpp"
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
