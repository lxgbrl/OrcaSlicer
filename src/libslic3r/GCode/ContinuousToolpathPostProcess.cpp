#include "ContinuousToolpathPostProcess.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/log/trivial.hpp>

namespace Slic3r {
namespace ContinuousToolpath {

namespace {

struct Bead {
    std::vector<Vec2d>  pts;     // mm
    std::vector<double> seg_e;   // E delta per segment (pts[i-1]->pts[i])
};
struct LayerData {
    double z = 0, height = 0.2;
    std::vector<Bead> beads;
};

inline bool num_after(const std::string& s, char tag, double& out)
{
    for (size_t i = 0; i + 1 < s.size(); ++i)
        if ((s[i] == tag || s[i] == tag + 32 /*lower*/) &&
            (i == 0 || s[i - 1] == ' ')) {
            try { out = std::stod(s.substr(i + 1)); return true; } catch (...) { return false; }
        }
    return false;
}

void tess_arc(double x0, double y0, double x1, double y1, double i, double j, bool cw,
              std::vector<Vec2d>& out)
{
    double cx = x0 + i, cy = y0 + j;
    double r = std::hypot(x0 - cx, y0 - cy);
    if (r < 1e-9) { out.emplace_back(x1, y1); return; }
    double a0 = std::atan2(y0 - cy, x0 - cx), a1 = std::atan2(y1 - cy, x1 - cx);
    if (cw)  { while (a1 >= a0) a1 -= 2 * M_PI; }
    else     { while (a1 <= a0) a1 += 2 * M_PI; }
    int n = std::max(2, int(std::ceil(std::fabs(a1 - a0) * r / 0.4)));
    for (int s = 1; s <= n; ++s) {
        double a = a0 + (a1 - a0) * (double(s) / n);
        out.emplace_back(cx + r * std::cos(a), cy + r * std::sin(a));
    }
    out.back() = Vec2d(x1, y1);
}

} // namespace

bool post_process_file(const std::string& path, const Params& params,
                       double filament_diameter_mm, double bridge_flow)
{
    std::ifstream in(path);
    if (! in.good()) return false;
    std::vector<std::string> lines;
    { std::string l; while (std::getline(in, l)) lines.push_back(l); }
    in.close();
    if (lines.empty()) return false;

    // header = up to first ;LAYER_CHANGE ; footer = after last motion line.
    size_t h_end = lines.size();
    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string& s = lines[i];
        size_t p = s.find_first_not_of(" \t");
        if (p != std::string::npos && s.compare(p, 13, ";LAYER_CHANGE") == 0) { h_end = i; break; }
    }
    auto is_motion = [](const std::string& s) {
        size_t p = s.find_first_not_of(" \t");
        if (p == std::string::npos) return false;
        return (s.compare(p, 2, "G1") == 0 || s.compare(p, 2, "G0") == 0 ||
                s.compare(p, 2, "G2") == 0 || s.compare(p, 2, "G3") == 0) &&
               (s.find('X', p) != std::string::npos || s.find('Y', p) != std::string::npos);
    };
    size_t f_start = lines.size();
    for (size_t i = lines.size(); i-- > 0;) if (is_motion(lines[i])) { f_start = i + 1; break; }

    // ---- parse layers --------------------------------------------------------
    std::vector<LayerData> layers;
    double x = 0, y = 0, z = 0, lastE = 0, curH = 0.2, print_f = 1500, travel_f = 6000;
    bool relE = false, pending_layer = true;   // start a layer at the first ;LAYER_CHANGE / G1 Z
    Bead cur; LayerData* L = nullptr;
    auto flush_bead = [&]() {
        if (cur.pts.size() >= 2 && L) L->beads.push_back(cur);
        cur = Bead();
    };
    auto new_layer = [&](double zv) {
        flush_bead();
        layers.emplace_back(); L = &layers.back(); L->z = zv; L->height = curH;
    };

    for (size_t i = h_end; i < f_start; ++i) {
        const std::string& s = lines[i];
        size_t p = s.find_first_not_of(" \t");
        if (p == std::string::npos) continue;
        if (s[p] == ';') {
            if (s.compare(p, 7, ";WIDTH:") == 0) { flush_bead(); }   // width change = bead boundary
            else if (s.compare(p, 8, ";HEIGHT:") == 0) { try { curH = std::stod(s.substr(p + 8)); } catch (...) {} if (L) L->height = curH; }
            else if (s.compare(p, 6, ";TYPE:") == 0) flush_bead();
            else if (s.compare(p, 13, ";LAYER_CHANGE") == 0) { flush_bead(); pending_layer = true; }
            else if (s.compare(p, 3, ";Z:") == 0) {
                try { double zz = std::stod(s.substr(p + 3));
                      if (pending_layer) { new_layer(zz); pending_layer = false; }
                      else if (L) L->z = zz;
                      z = zz; } catch (...) {}
            }
            continue;
        }
        if (s.compare(p, 3, "M83") == 0) { relE = true; continue; }
        if (s.compare(p, 3, "M82") == 0) { relE = false; continue; }
        if (s.compare(p, 3, "G92") == 0) { double e; if (num_after(s, 'E', e)) lastE = e; continue; }
        bool g2 = s.compare(p, 2, "G2") == 0, g3 = s.compare(p, 2, "G3") == 0;
        bool g01 = s.compare(p, 2, "G1") == 0 || s.compare(p, 2, "G0") == 0;
        if (!(g01 || g2 || g3)) continue;
        double nx = x, ny = y, vz = 0, ve = 0;
        bool hx = num_after(s, 'X', nx), hy = num_after(s, 'Y', ny);
        if (num_after(s, 'Z', vz)) {
            if (pending_layer) { new_layer(vz); pending_layer = false; }   // layer Z from G1 Z when no ;Z: present
            else if (L == nullptr) new_layer(vz);
            z = vz;
        }
        bool he = num_after(s, 'E', ve);
        bool extr = false; double de = 0;
        if (he) { de = relE ? ve : ve - lastE; if (!relE) lastE = ve; extr = de > 1e-6 && (hx || hy); }
        if (extr) {
            double vf2; if (num_after(s, 'F', vf2)) print_f = vf2;
            if (cur.pts.empty()) cur.pts.emplace_back(x, y);
            if (g2 || g3) {
                double ii = 0, jj = 0; num_after(s, 'I', ii); num_after(s, 'J', jj);
                std::vector<Vec2d> arc; tess_arc(x, y, nx, ny, ii, jj, g2, arc);
                // distribute E evenly across tessellated segments
                double per = de / std::max<size_t>(1, arc.size());
                for (auto& pt : arc) { cur.pts.push_back(pt); cur.seg_e.push_back(per); }
            } else {
                cur.pts.emplace_back(nx, ny); cur.seg_e.push_back(de);
            }
        } else {
            double vf2; if (num_after(s, 'F', vf2)) travel_f = vf2;
            flush_bead();
        }
        x = nx; y = ny;
    }
    flush_bead();
    if (layers.empty()) return false;

    // ---- re-emit -------------------------------------------------------------
    const double fil_area = M_PI * (filament_diameter_mm / 2.0) * (filament_diameter_mm / 2.0);
    std::ostringstream out;
    for (size_t i = 0; i < h_end; ++i) out << lines[i] << "\n";
    out << "M83\n; --- continuous toolpath (post-process) ---\n";

    char buf[160];
    for (LayerData& ld : layers) {
        if (ld.beads.empty()) continue;
        Polylines polys; polys.reserve(ld.beads.size());
        for (Bead& b : ld.beads) {
            Polyline pl; pl.points.reserve(b.pts.size());
            for (Vec2d& q : b.pts) pl.points.emplace_back(Point(scale_(q.x()), scale_(q.y())));
            polys.push_back(std::move(pl));
        }
        std::vector<Move> moves = order(polys, params);

        out << ";LAYER_CHANGE\n;Z:" << ld.z << "\n";
        std::snprintf(buf, sizeof(buf), "G1 Z%.3f F%.0f\n", ld.z, travel_f); out << buf;
        bool have_last = false; Vec2d last(0, 0);
        auto travel_to = [&](const Vec2d& q) {
            if (have_last && (q - last).squaredNorm() < 1e-8) return;
            std::snprintf(buf, sizeof(buf), "G1 X%.3f Y%.3f F%.0f\n", q.x(), q.y(), travel_f);
            out << buf; last = q; have_last = true;
        };
        for (Move& m : moves) {
            if (m.bridge && !m.travel) {                       // sacrificial extruding bridge
                Vec2d a(unscale<double>(m.polyline.points.front().x()), unscale<double>(m.polyline.points.front().y()));
                Vec2d b(unscale<double>(m.polyline.points.back().x()),  unscale<double>(m.polyline.points.back().y()));
                travel_to(a);
                double len = (b - a).norm();
                double e = bridge_flow * 0.3 * ld.height * len / fil_area;
                std::snprintf(buf, sizeof(buf), "G1 X%.3f Y%.3f E%.5f F%.0f\n", b.x(), b.y(), e, print_f);
                out << buf; last = b; have_last = true;
            } else if (m.travel) {                             // long jump: travel, no extrude
                Vec2d b(unscale<double>(m.polyline.points.back().x()), unscale<double>(m.polyline.points.back().y()));
                travel_to(b);
            } else if (m.src_index >= 0) {                     // original bead: reuse its E
                Bead& b = ld.beads[m.src_index];
                std::vector<Vec2d> pts = b.pts; std::vector<double> se = b.seg_e;
                if (m.reversed) { std::reverse(pts.begin(), pts.end()); std::reverse(se.begin(), se.end()); }
                travel_to(pts.front());
                for (size_t k = 1; k < pts.size(); ++k) {
                    std::snprintf(buf, sizeof(buf), "G1 X%.3f Y%.3f E%.5f F%.0f\n",
                                  pts[k].x(), pts[k].y(), se[k - 1], print_f);
                    out << buf; last = pts[k]; have_last = true;
                }
            }
        }
    }
    for (size_t i = f_start; i < lines.size(); ++i) out << lines[i] << "\n";

    std::ofstream of(path, std::ios::trunc);
    if (! of.good()) return false;
    of << out.str();
    BOOST_LOG_TRIVIAL(info) << "continuous_toolpath: rewrote " << layers.size() << " layers in " << path;
    return true;
}

} // namespace ContinuousToolpath
} // namespace Slic3r
