#ifndef slic3r_ContinuousToolpath_hpp_
#define slic3r_ContinuousToolpath_hpp_

#include <vector>
#include "../libslic3r.h"
#include "../Polyline.hpp"

namespace Slic3r {
namespace ContinuousToolpath {

// Parameters for the continuous-toolpath reordering (see
// AppEntwicklung/continuous_toolpath_plan.md). C++ port of the validated Python
// prototype: order printable polylines into one near-continuous path per layer by
// Eulerizing the segment graph (odd endpoints joined by short connectors, not bead
// reprints) and chaining disjoint components.
struct Params
{
    // Snap endpoints within this distance (scaled coords) to unify graph nodes.
    coord_t snap_tol      = scale_(0.05);
    // Force one continuous tour per layer (bridge all components). If false, each
    // connected component yields its own continuous stroke.
    bool    single_path   = true;
    // A bridge/connector extrudes only if its gap is at most this (scaled); longer
    // jumps stay travels so material never strings across a hole/void.
    coord_t sacrificial_max = scale_(3.0);
};

// One emitted move in the continuous order.
struct Move
{
    Polyline polyline;     // geometry (>=2 points), already oriented for traversal
    bool     travel = false;   // true => non-extruding move (jump); false => extrude
    bool     bridge = false;   // synthetic connector (short, extruding) vs original bead
};

// Reorder `input` polylines into a continuous traversal. Geometry of original beads
// is preserved verbatim; only ordering + connector moves are added.
std::vector<Move> order(const Polylines& input, const Params& params);

// Convenience: just the ordered polylines (bridges/travels included as segments).
Polylines order_polylines(const Polylines& input, const Params& params);

} // namespace ContinuousToolpath
} // namespace Slic3r

#endif // slic3r_ContinuousToolpath_hpp_
