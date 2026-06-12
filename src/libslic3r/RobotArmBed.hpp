#ifndef slic3r_RobotArmBed_hpp_
#define slic3r_RobotArmBed_hpp_

#include "Point.hpp"

namespace Slic3r {

// Parametric description of a robot arm's reachable workspace at table height:
// an annular sector between reach_min (base dead zone) and reach_max, spanning
// sweep_deg around the base, centered on base_rot_deg.
struct RobotArmWorkspace
{
    double reach_min   = 0.;   // mm, inner dead-zone radius, >= 0
    double reach_max   = 0.;   // mm, maximum reach, > reach_min
    double sweep_deg   = 360.; // base joint sweep, (0, 360]
    Vec2d  base_offset{0., 0.}; // robot base position relative to bed origin, mm
    double base_rot_deg = 0.;  // orientation of the sweep bisector, degrees
};

bool robot_workspace_valid(const RobotArmWorkspace &ws);

// Outer printable_area contour (CCW). Annular sector as a single simple polygon
// when sweep_deg < 360; full circle polygon when sweep_deg >= 360 (the inner
// dead zone is then expressed separately via robot_workspace_exclude_area).
Pointfs robot_workspace_printable_area(const RobotArmWorkspace &ws, int segments_per_full_circle = 72);

// Inner dead-zone polygon for bed_exclude_area. Empty unless sweep_deg >= 360
// and reach_min > 0 (a sector embeds the inner arc in its own contour).
// The point count is kept a multiple of 4 because some legacy consumers group
// exclude-area points into quads.
Pointfs robot_workspace_exclude_area(const RobotArmWorkspace &ws, int segments_per_full_circle = 64);

} // namespace Slic3r

#endif // slic3r_RobotArmBed_hpp_
