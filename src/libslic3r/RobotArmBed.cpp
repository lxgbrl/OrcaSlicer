#include "RobotArmBed.hpp"

#include <cmath>

namespace Slic3r {

static constexpr double ROBOT_SWEEP_FULL = 360.;

bool robot_workspace_valid(const RobotArmWorkspace &ws)
{
    if (!(ws.reach_min >= 0.))
        return false;
    if (!(ws.reach_max > ws.reach_min))
        return false;
    if (!(ws.sweep_deg > 0. && ws.sweep_deg <= ROBOT_SWEEP_FULL))
        return false;
    return true;
}

static Vec2d arc_point(const RobotArmWorkspace &ws, double radius, double angle_deg)
{
    const double a = angle_deg * M_PI / 180.;
    return ws.base_offset + radius * Vec2d(std::cos(a), std::sin(a));
}

// Sample an arc of `radius` from angle a0 to a1 (degrees, inclusive on both ends).
static void append_arc(Pointfs &pts, const RobotArmWorkspace &ws, double radius, double a0, double a1, int segments_per_full_circle)
{
    const double span = a1 - a0;
    // At least 2 points per arc so radial edges stay straight lines.
    const int n = std::max(2, (int) std::ceil(std::abs(span) / ROBOT_SWEEP_FULL * segments_per_full_circle));
    for (int i = 0; i <= n; ++i)
        pts.emplace_back(arc_point(ws, radius, a0 + span * double(i) / double(n)));
}

Pointfs robot_workspace_printable_area(const RobotArmWorkspace &ws, int segments_per_full_circle)
{
    Pointfs pts;
    if (!robot_workspace_valid(ws))
        return pts;

    if (ws.sweep_deg >= ROBOT_SWEEP_FULL) {
        // Full circle of reach_max; the dead zone is handled by the exclude area.
        for (int i = 0; i < segments_per_full_circle; ++i)
            pts.emplace_back(arc_point(ws, ws.reach_max, ROBOT_SWEEP_FULL * double(i) / double(segments_per_full_circle)));
        return pts;
    }

    // Annular sector: outer arc CCW, then inner arc back CW (or the apex when
    // reach_min == 0). Sector is centered on base_rot_deg.
    const double a_start = ws.base_rot_deg - ws.sweep_deg / 2.;
    const double a_end   = ws.base_rot_deg + ws.sweep_deg / 2.;
    append_arc(pts, ws, ws.reach_max, a_start, a_end, segments_per_full_circle);
    if (ws.reach_min > 0.)
        append_arc(pts, ws, ws.reach_min, a_end, a_start, segments_per_full_circle);
    else
        pts.emplace_back(ws.base_offset);
    return pts;
}

Pointfs robot_workspace_exclude_area(const RobotArmWorkspace &ws, int segments_per_full_circle)
{
    Pointfs pts;
    if (!robot_workspace_valid(ws))
        return pts;
    if (ws.sweep_deg < ROBOT_SWEEP_FULL || ws.reach_min <= 0.)
        return pts;

    // Keep the point count a multiple of 4 (legacy consumers chunk exclude
    // areas into quads).
    int n = std::max(8, segments_per_full_circle);
    n -= n % 4;
    for (int i = 0; i < n; ++i)
        pts.emplace_back(arc_point(ws, ws.reach_min, ROBOT_SWEEP_FULL * double(i) / double(n)));
    return pts;
}

} // namespace Slic3r
