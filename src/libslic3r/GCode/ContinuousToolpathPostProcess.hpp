#ifndef slic3r_ContinuousToolpathPostProcess_hpp_
#define slic3r_ContinuousToolpathPostProcess_hpp_

#include <string>
#include "ContinuousToolpath.hpp"

namespace Slic3r {
namespace ContinuousToolpath {

// Approach "B": a self-contained post-process pass over a finished G-code file.
// Re-orders each layer's extrusions into a continuous, retraction-free path using
// ContinuousToolpath::order(). Original bead geometry AND their per-segment E are
// preserved (flow calibration kept); only ordering + short bridges are added, and
// retractions are dropped. Header (start g-code) and footer (end g-code + Orca
// CONFIG_BLOCK) are passed through unchanged so the file still previews and prints.
//
// Returns true if the file was rewritten, false on parse failure (file untouched).
bool post_process_file(const std::string& gcode_path,
                       const Params&      params,
                       double             filament_diameter_mm,
                       double             bridge_flow = 0.6);

} // namespace ContinuousToolpath
} // namespace Slic3r

#endif
