#ifndef slic3r_FillCustomScripted_hpp_
#define slic3r_FillCustomScripted_hpp_

#include "../libslic3r.h"

#include "FillBase.hpp"

namespace Slic3r {

// A single scriptable infill engine driven by external pattern definitions.
//
//   * Tile2D   - tiles a 2D polyline pattern (loaded from a .tile file) across
//                the region, with serpentine connection and density culling.
//   * Volume3D - samples a 3D implicit/TPMS field (gyroid, Schwarz P, ...) on the
//                current layer Z and extracts iso-contours via marching squares.
//
// New patterns are added as data (tile files / parametric configs) resolved by
// PatternManager; no recompile is required once this engine exists.
//
// See custom_infill_architecture.md for the full design.
class FillCustomScripted : public Fill
{
public:
    FillCustomScripted() {}
    Fill* clone() const override { return new FillCustomScripted(*this); }

    bool is_self_crossing() override { return false; }

protected:
    void _fill_surface_single(
        const FillParams                &params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction,
        ExPolygon                        expolygon,
        Polylines                       &polylines_out) override;

private:
    // 2D periodic tile pipeline.
    Polylines fill_tile_2d(const FillParams& params, const ExPolygon& expolygon, const BoundingBox& bb);
    // 3D implicit-field iso-contour pipeline.
    Polylines fill_volume_3d(const FillParams& params, const ExPolygon& expolygon, const BoundingBox& bb);
};

} // namespace Slic3r

#endif // slic3r_FillCustomScripted_hpp_
