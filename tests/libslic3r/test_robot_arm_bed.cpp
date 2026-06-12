#include <catch2/catch_all.hpp>

#include "libslic3r/Point.hpp"
#include "libslic3r/RobotArmBed.hpp"

using namespace Slic3r;

static double polygon_area(const Pointfs &pts)
{
    double a = 0.;
    for (size_t i = 0; i < pts.size(); ++i) {
        const Vec2d &p = pts[i];
        const Vec2d &q = pts[(i + 1) % pts.size()];
        a += p.x() * q.y() - q.x() * p.y();
    }
    return a / 2.;
}

TEST_CASE("Robot arm workspace validation", "[RobotArmBed]") {
    RobotArmWorkspace ws;
    ws.reach_min = 150.;
    ws.reach_max = 850.;
    ws.sweep_deg = 270.;
    REQUIRE(robot_workspace_valid(ws));

    SECTION("reach_max must exceed reach_min") {
        ws.reach_max = 100.;
        REQUIRE(!robot_workspace_valid(ws));
    }
    SECTION("reach_min must be non-negative") {
        ws.reach_min = -1.;
        REQUIRE(!robot_workspace_valid(ws));
    }
    SECTION("sweep must be in (0, 360]") {
        ws.sweep_deg = 0.;
        REQUIRE(!robot_workspace_valid(ws));
        ws.sweep_deg = 361.;
        REQUIRE(!robot_workspace_valid(ws));
        ws.sweep_deg = 360.;
        REQUIRE(robot_workspace_valid(ws));
    }
}

TEST_CASE("Annular sector polygon", "[RobotArmBed]") {
    RobotArmWorkspace ws;
    ws.reach_min = 150.;
    ws.reach_max = 850.;
    ws.sweep_deg = 270.;

    const Pointfs poly = robot_workspace_printable_area(ws);
    REQUIRE(poly.size() > 4);

    SECTION("CCW orientation and correct area") {
        const double expect = (ws.sweep_deg / 360.) * M_PI * (ws.reach_max * ws.reach_max - ws.reach_min * ws.reach_min);
        const double area   = polygon_area(poly);
        REQUIRE(area > 0.); // CCW
        REQUIRE(area == Catch::Approx(expect).epsilon(0.01));
    }
    SECTION("no vertex inside the dead zone or beyond reach") {
        for (const Vec2d &p : poly) {
            const double r = p.norm();
            REQUIRE(r >= ws.reach_min - 1e-9);
            REQUIRE(r <= ws.reach_max + 1e-9);
        }
    }
    SECTION("no exclude area for a sector") {
        REQUIRE(robot_workspace_exclude_area(ws).empty());
    }
    SECTION("zero dead zone collapses inner arc to apex") {
        ws.reach_min = 0.;
        const Pointfs p2 = robot_workspace_printable_area(ws);
        REQUIRE(std::count_if(p2.begin(), p2.end(), [](const Vec2d &p) { return p.norm() < 1e-9; }) == 1);
    }
}

TEST_CASE("Full circle workspace", "[RobotArmBed]") {
    RobotArmWorkspace ws;
    ws.reach_min = 150.;
    ws.reach_max = 850.;
    ws.sweep_deg = 360.;

    const Pointfs outer = robot_workspace_printable_area(ws);
    REQUIRE(outer.size() == 72);
    REQUIRE(polygon_area(outer) == Catch::Approx(M_PI * ws.reach_max * ws.reach_max).epsilon(0.01));

    const Pointfs inner = robot_workspace_exclude_area(ws);
    REQUIRE(!inner.empty());
    REQUIRE(inner.size() % 4 == 0);
    REQUIRE(polygon_area(inner) == Catch::Approx(M_PI * ws.reach_min * ws.reach_min).epsilon(0.01));

    SECTION("no exclude area when dead zone is zero") {
        ws.reach_min = 0.;
        REQUIRE(robot_workspace_exclude_area(ws).empty());
    }
}

TEST_CASE("Base offset and rotation", "[RobotArmBed]") {
    RobotArmWorkspace ws;
    ws.reach_min   = 100.;
    ws.reach_max   = 500.;
    ws.sweep_deg   = 180.;
    ws.base_offset = Vec2d(250., 250.);
    ws.base_rot_deg = 90.;

    const Pointfs poly = robot_workspace_printable_area(ws);
    REQUIRE(!poly.empty());

    SECTION("all radii measured from the offset base") {
        for (const Vec2d &p : poly) {
            const double r = (p - ws.base_offset).norm();
            REQUIRE(r >= ws.reach_min - 1e-9);
            REQUIRE(r <= ws.reach_max + 1e-9);
        }
    }
    SECTION("sector lies on the bisector side (rot 90 deg -> +Y half-plane)") {
        for (const Vec2d &p : poly)
            REQUIRE(p.y() >= ws.base_offset.y() - 1e-9);
    }
    SECTION("invalid workspace yields empty polygons") {
        ws.reach_max = 0.;
        REQUIRE(robot_workspace_printable_area(ws).empty());
        REQUIRE(robot_workspace_exclude_area(ws).empty());
    }
}
