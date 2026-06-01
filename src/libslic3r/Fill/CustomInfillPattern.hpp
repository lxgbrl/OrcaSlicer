#ifndef slic3r_CustomInfillPattern_hpp_
#define slic3r_CustomInfillPattern_hpp_

#include <string>
#include <vector>
#include <map>
#include <mutex>

#include <memory>

#include "../libslic3r.h"
#include "../Polyline.hpp"
#include "ImplicitExpr.hpp"

namespace Slic3r {

// A 2D periodic tile loaded from a ".tile" text file.
// Coordinates are stored tile-local, already scaled to Slic3r coord_t.
// The tile origin is (0,0); width/height declare the period in mm.
struct TilePattern
{
    double     width  = 10.0;   // tile period in X (mm)
    double     height = 10.0;   // tile period in Y (mm)
    bool       serpentine = false;
    Polylines  paths;           // scaled coords, tile-local
    bool       valid = false;
};

// Implicit-surface family for the 3D volumetric mode.
//   Gyroid / SchwarzP : built-in TPMS families.
//   Expr              : user formula f(x,y,z,t) (incl. imported MathMod Iso3D).
enum class VolumeType { Gyroid, SchwarzP, Expr };

// A 3D volumetric / TPMS pattern definition loaded from a parametric config
// (JSON or simple INI-like text).
struct VolumePattern
{
    VolumeType type          = VolumeType::Gyroid;
    double     cell_size     = 6.0;   // periodic cell size (mm)
    double     level         = 0.0;   // iso-contour threshold
    double     thickness     = 0.6;   // nominal wall thickness (mm), reserved
    double     density_scale = 1.0;   // multiplier applied to config density

    // Expr-only: compiled formula and the per-cell coordinate domain that the
    // cell (size cell_size) is mapped onto, plus the fixed animation parameter.
    std::shared_ptr<const ImplicitExpr> expr;
    double     domain_min = -3.14159265358979323846;
    double     domain_max =  3.14159265358979323846;
    double     t          = 0.0;

    bool       valid = false;
};

// Locates, parses and caches custom infill pattern definitions for the
// duration of a slicing job. Definitions are searched in:
//   1. <resources_dir>/custom_infill/
//   2. <data_dir>/custom_infill/        (user folder, takes precedence)
// A pattern id may be a bare name (extension inferred) or a filename.
class PatternManager
{
public:
    static PatternManager& instance();

    // Returns a cached tile pattern. default_w/default_h are used only when the
    // .tile file omits a TILE header. Falls back to a built-in rectilinear-like
    // tile if the id cannot be resolved (so slicing never hard-fails).
    const TilePattern& get_tile_pattern(const std::string& id, double default_w, double default_h);

    // Returns a cached volumetric pattern. Falls back to a default gyroid if the
    // id cannot be resolved.
    const VolumePattern& get_volume_pattern(const std::string& id);

    // Drop cached patterns (call between slicing jobs if definitions changed).
    void clear_cache();

    // Parsers exposed for testing.
    static bool parse_tile_file(const std::string& path, TilePattern& out);
    static bool parse_volume_config(const std::string& path, VolumePattern& out);

private:
    PatternManager() = default;
    // Resolve a pattern id to an on-disk path with one of the given extensions.
    static std::string resolve_path(const std::string& id, const std::vector<std::string>& exts);

    std::mutex                          m_mutex;
    std::map<std::string, TilePattern>   m_tile_cache;
    std::map<std::string, VolumePattern> m_volume_cache;
};

} // namespace Slic3r

#endif // slic3r_CustomInfillPattern_hpp_
