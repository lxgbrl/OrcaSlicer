#include <catch2/catch_all.hpp>

#include "libslic3r/GCode/ContinuousToolpath.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/libslic3r.h"

using namespace Slic3r;
using namespace Slic3r::ContinuousToolpath;

static Polyline mk(std::initializer_list<std::pair<double, double>> pts)
{
    Polyline pl;
    for (auto& p : pts) pl.points.emplace_back(Point(scale_(p.first), scale_(p.second)));
    return pl;
}

static int bead_moves(const std::vector<Move>& ms)
{
    int n = 0;
    for (const Move& m : ms) if (!m.travel && !m.bridge) ++n;
    return n;
}

TEST_CASE("ContinuousToolpath: each bead traversed once (no reprint)", "[CTP]")
{
    Polylines in;
    in.push_back(mk({{0,0},{20,0},{20,20},{0,20},{0,0}}));  // closed wall loop
    in.push_back(mk({{2,10},{18,10}}));                      // one infill line

    Params p; p.single_path = true;
    auto moves = order(in, p);

    REQUIRE_FALSE(moves.empty());
    // Both original beads emitted exactly once -> no full-bead reprints.
    REQUIRE(bead_moves(moves) == 2);
}

TEST_CASE("ContinuousToolpath: single_path bridges disjoint components", "[CTP]")
{
    Polylines in;
    in.push_back(mk({{0,0},{10,0},{10,10},{0,10},{0,0}}));   // ring A
    in.push_back(mk({{40,0},{50,0},{50,10},{40,10},{40,0}})); // ring B (far away)

    Params p; p.single_path = true; p.sacrificial_max = scale_(3.0);
    auto moves = order(in, p);

    REQUIRE(bead_moves(moves) == 2);
    // the far jump between the two rings must be a travel (gap >> sacrificial_max)
    bool has_travel = false;
    for (const Move& m : moves) if (m.travel) has_travel = true;
    REQUIRE(has_travel);
}

TEST_CASE("ContinuousToolpath: empty input is safe", "[CTP]")
{
    Polylines in;
    REQUIRE(order(in, Params{}).empty());
}
