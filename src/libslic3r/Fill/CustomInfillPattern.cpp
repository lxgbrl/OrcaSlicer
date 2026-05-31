#include "CustomInfillPattern.hpp"

#include <fstream>
#include <sstream>
#include <cmath>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "nlohmann/json.hpp"

#include "../Utils.hpp"

namespace Slic3r {

PatternManager& PatternManager::instance()
{
    static PatternManager s_inst;
    return s_inst;
}

void PatternManager::clear_cache()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_tile_cache.clear();
    m_volume_cache.clear();
}

std::string PatternManager::resolve_path(const std::string& id, const std::vector<std::string>& exts)
{
    namespace fs = boost::filesystem;
    if (id.empty())
        return std::string();

    // Search order: user data dir first (so users can override bundled patterns),
    // then the bundled resources dir.
    std::vector<fs::path> dirs;
    if (! data_dir().empty())
        dirs.emplace_back(fs::path(data_dir()) / "custom_infill");
    if (! resources_dir().empty())
        dirs.emplace_back(fs::path(resources_dir()) / "custom_infill");

    // If the id already carries an extension or a path separator, try it verbatim.
    std::vector<std::string> candidates;
    candidates.push_back(id);
    for (const std::string& ext : exts)
        candidates.push_back(id + ext);

    for (const fs::path& dir : dirs)
        for (const std::string& cand : candidates) {
            fs::path p = dir / cand;
            boost::system::error_code ec;
            if (fs::exists(p, ec) && ! fs::is_directory(p, ec))
                return p.string();
        }
    // Also accept an absolute / cwd-relative path.
    for (const std::string& cand : candidates) {
        boost::system::error_code ec;
        if (fs::exists(cand, ec) && ! fs::is_directory(fs::path(cand), ec))
            return cand;
    }
    return std::string();
}

bool PatternManager::parse_tile_file(const std::string& path, TilePattern& out)
{
    std::ifstream in(path);
    if (! in.good())
        return false;

    out = TilePattern();
    std::string line;
    bool got_path = false;
    while (std::getline(in, line)) {
        boost::trim(line);
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream ss(line);
        std::string kw;
        ss >> kw;
        boost::to_upper(kw);
        if (kw == "TILE") {
            double w = 0, h = 0;
            if (ss >> w >> h && w > 0 && h > 0) {
                out.width  = w;
                out.height = h;
            }
        } else if (kw == "MODE") {
            std::string m;
            ss >> m;
            boost::to_lower(m);
            if (m == "serpentine")
                out.serpentine = true;
        } else if (kw == "PATH") {
            Polyline pl;
            std::string tok;
            while (ss >> tok) {
                // token is "x,y"
                auto comma = tok.find(',');
                if (comma == std::string::npos)
                    continue;
                try {
                    double x = std::stod(tok.substr(0, comma));
                    double y = std::stod(tok.substr(comma + 1));
                    pl.points.emplace_back(Point(scale_(x), scale_(y)));
                } catch (...) {
                    // skip malformed point
                }
            }
            if (pl.points.size() >= 2) {
                out.paths.emplace_back(std::move(pl));
                got_path = true;
            }
        }
        // unknown directives are ignored (forward-compat)
    }
    out.valid = got_path;
    return out.valid;
}

bool PatternManager::parse_volume_config(const std::string& path, VolumePattern& out)
{
    std::ifstream in(path);
    if (! in.good())
        return false;

    out = VolumePattern();

    namespace fs = boost::filesystem;
    std::string ext = fs::path(path).extension().string();
    boost::to_lower(ext);

    auto set_type = [&out](std::string t) {
        boost::to_lower(t);
        if (t == "schwarzp" || t == "schwarz_p" || t == "schwarz-p" || t == "p")
            out.type = VolumeType::SchwarzP;
        else
            out.type = VolumeType::Gyroid;
    };

    if (ext == ".json") {
        try {
            nlohmann::json j;
            in >> j;
            if (j.contains("type"))          set_type(j["type"].get<std::string>());
            if (j.contains("cell_size"))     out.cell_size     = j["cell_size"].get<double>();
            if (j.contains("level"))         out.level         = j["level"].get<double>();
            if (j.contains("thickness"))     out.thickness     = j["thickness"].get<double>();
            if (j.contains("density_scale")) out.density_scale = j["density_scale"].get<double>();
            out.valid = true;
        } catch (...) {
            return false;
        }
    } else {
        // INI-like: "KEY value" per line.
        std::string line;
        while (std::getline(in, line)) {
            boost::trim(line);
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream ss(line);
            std::string kw;
            ss >> kw;
            boost::to_upper(kw);
            if (kw == "PATTERN" || kw == "TYPE") {
                std::string t; ss >> t; set_type(t);
            } else if (kw == "CELL_SIZE") {
                ss >> out.cell_size;
            } else if (kw == "LEVEL") {
                ss >> out.level;
            } else if (kw == "THICKNESS") {
                ss >> out.thickness;
            } else if (kw == "DENSITY_SCALE") {
                ss >> out.density_scale;
            }
        }
        out.valid = true;
    }
    if (out.cell_size <= 0)
        out.cell_size = 6.0;
    return out.valid;
}

// Built-in fallback: a single horizontal line across the tile. Tiled in rows it
// yields a rectilinear-like pattern, so slicing degrades gracefully when an id
// cannot be resolved.
static TilePattern make_fallback_tile(double w, double h)
{
    TilePattern t;
    t.width  = w > 0 ? w : 10.0;
    t.height = h > 0 ? h : 10.0;
    t.serpentine = true;
    Polyline pl;
    pl.points.emplace_back(Point(scale_(0.0),       scale_(t.height * 0.5)));
    pl.points.emplace_back(Point(scale_(t.width),   scale_(t.height * 0.5)));
    t.paths.emplace_back(std::move(pl));
    t.valid = true;
    return t;
}

const TilePattern& PatternManager::get_tile_pattern(const std::string& id, double default_w, double default_h)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    auto it = m_tile_cache.find(id);
    if (it != m_tile_cache.end())
        return it->second;

    TilePattern tp;
    std::string path = resolve_path(id, {".tile", ".txt"});
    if (path.empty() || ! parse_tile_file(path, tp))
        tp = make_fallback_tile(default_w, default_h);
    else {
        // honor config defaults only when the file omitted them
        if (tp.width  <= 0) tp.width  = default_w;
        if (tp.height <= 0) tp.height = default_h;
    }
    auto res = m_tile_cache.emplace(id, std::move(tp));
    return res.first->second;
}

const VolumePattern& PatternManager::get_volume_pattern(const std::string& id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    auto it = m_volume_cache.find(id);
    if (it != m_volume_cache.end())
        return it->second;

    VolumePattern vp;
    std::string path = resolve_path(id, {".json", ".ini", ".txt", ".cfg"});
    if (path.empty() || ! parse_volume_config(path, vp)) {
        // default gyroid
        vp = VolumePattern();
        vp.valid = true;
    }
    auto res = m_volume_cache.emplace(id, std::move(vp));
    return res.first->second;
}

} // namespace Slic3r
