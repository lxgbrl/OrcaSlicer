#include "ContinuousToolpath.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace Slic3r {
namespace ContinuousToolpath {

namespace {

inline double sq(double v) { return v * v; }
inline double dist2(const Point& a, const Point& b)
{
    return sq(double(a.x()) - double(b.x())) + sq(double(a.y()) - double(b.y()));
}

// Union-find for connected components.
struct DSU
{
    std::vector<int> p;
    void init(int n) { p.resize(n); std::iota(p.begin(), p.end(), 0); }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void uni(int a, int b) { p[find(a)] = find(b); }
};

// A graph edge backed by an original bead polyline (or a synthetic connector).
struct GEdge
{
    int u, v;                 // node ids
    bool bridge = false;      // synthetic connector (no original bead)
    int  src = -1;            // index into the input Polylines (-1 for a bridge)
    Polyline poly;            // geometry u->v (empty for bridge; built on demand)
};

} // namespace

std::vector<Move> order(const Polylines& input, const Params& params)
{
    std::vector<Move> out;

    // ---- 1. snap endpoints -> nodes ------------------------------------------
    const coord_t tol = std::max<coord_t>(params.snap_tol, 1);
    std::map<std::pair<coord_t, coord_t>, int> node_ix;
    std::vector<Point> nodes;
    auto node_of = [&](const Point& p) -> int {
        std::pair<coord_t, coord_t> key{ p.x() / tol, p.y() / tol };
        auto it = node_ix.find(key);
        if (it != node_ix.end()) return it->second;
        int id = int(nodes.size());
        node_ix.emplace(key, id);
        nodes.push_back(p);
        return id;
    };

    std::vector<GEdge> edges;
    for (int i = 0; i < int(input.size()); ++i) {
        const Polyline& pl = input[i];
        if (pl.points.size() < 2) continue;
        GEdge e;
        e.u = node_of(pl.points.front());
        e.v = node_of(pl.points.back());
        e.src = i;
        e.poly = pl;
        edges.push_back(std::move(e));
    }
    if (edges.empty()) return out;
    const int N = int(nodes.size());

    // ---- 2. degrees + connected components -----------------------------------
    DSU dsu; dsu.init(N);
    std::vector<int> degree(N, 0);
    for (const GEdge& e : edges) { degree[e.u]++; degree[e.v]++; dsu.uni(e.u, e.v); }

    auto add_connector = [&](int a, int b) {
        GEdge e; e.u = a; e.v = b; e.bridge = true;
        edges.push_back(e);
        degree[a]++; degree[b]++;
        dsu.uni(a, b);
    };

    // ---- 3. single_path: bridge separate components into one graph FIRST.
    // Must happen before the parity fix: bridges change degrees, so pairing the
    // odd nodes afterwards is what restores an Eulerian (all-even) graph.
    if (params.single_path) {
        // Group nodes by component (before connectors change the DSU), then
        // repeatedly attach the component whose CLOSEST NODE PAIR to the grown
        // cluster is shortest. Pairing arbitrary representatives instead can
        // pick connectors tens of mm long between chains that physically
        // interleave a line-width apart.
        std::map<int, std::vector<int>> comp_nodes;
        for (int n = 0; n < N; ++n) comp_nodes[dsu.find(n)].push_back(n);
        if (comp_nodes.size() > 1) {
            std::vector<std::vector<int>> comps;
            for (auto& kv : comp_nodes) comps.push_back(std::move(kv.second));
            std::vector<char> done(comps.size(), 0);
            std::vector<int>  cluster = comps[0];
            done[0] = 1;
            for (size_t step = 1; step < comps.size(); ++step) {
                int    best_comp = -1, best_a = -1, best_b = -1;
                double bd = 0;
                for (size_t j = 0; j < comps.size(); ++j) {
                    if (done[j]) continue;
                    for (int a : cluster)
                        for (int b : comps[j]) {
                            double d = dist2(nodes[a], nodes[b]);
                            if (best_comp < 0 || d < bd) { bd = d; best_comp = int(j); best_a = a; best_b = b; }
                        }
                }
                if (best_comp < 0) break;
                add_connector(best_a, best_b);
                cluster.insert(cluster.end(), comps[size_t(best_comp)].begin(), comps[size_t(best_comp)].end());
                done[size_t(best_comp)] = 1;
            }
        }
    }

    // ---- 4. Eulerize: pair odd-degree nodes within each (current) component ----
    // Greedy nearest pairing (blossom-optimal is a refinement).
    {
        std::map<int, std::vector<int>> odd_by_comp;
        for (int n = 0; n < N; ++n)
            if (degree[n] & 1) odd_by_comp[dsu.find(n)].push_back(n);
        for (auto& kv : odd_by_comp) {
            std::vector<int>& odd = kv.second;
            std::vector<char> used(odd.size(), 0);
            for (size_t i = 0; i < odd.size(); ++i) {
                if (used[i]) continue;
                int best = -1; double bd = 0;
                for (size_t j = i + 1; j < odd.size(); ++j) {
                    if (used[j]) continue;
                    double d = dist2(nodes[odd[i]], nodes[odd[j]]);
                    if (best < 0 || d < bd) { bd = d; best = int(j); }
                }
                if (best >= 0) { used[i] = used[size_t(best)] = 1; add_connector(odd[i], odd[best]); }
            }
        }
    }

    // ---- 5. Hierholzer Eulerian traversal ------------------------------------
    // adjacency: node -> list of (edge index)
    std::vector<std::vector<int>> adj(N);
    for (int i = 0; i < int(edges.size()); ++i) { adj[edges[i].u].push_back(i); adj[edges[i].v].push_back(i); }
    std::vector<char> used_edge(edges.size(), 0);
    std::vector<size_t> it_pos(N, 0);

    // Cover EVERY component: repeat Hierholzer until no unused edge remains.
    // (A single run only covers one connected component, and only fully when that
    // component is Eulerian. Looping guarantees no geometry is dropped.)
    auto find_start = [&]() -> int {
        for (int n = 0; n < N; ++n) if (degree[n] & 1)
            for (int ei : adj[n]) if (! used_edge[ei]) return n;
        for (int n = 0; n < N; ++n)
            for (int ei : adj[n]) if (! used_edge[ei]) return n;
        return -1;
    };

    std::vector<int> node_path;
    for (;;) {
        int start = find_start();
        if (start < 0) break;
        std::vector<int> node_stack{ start }, seg;
        while (! node_stack.empty()) {
            int v = node_stack.back();
            size_t& it = it_pos[v];
            while (it < adj[v].size() && used_edge[adj[v][it]]) ++it;
            if (it == adj[v].size()) {
                seg.push_back(v);
                node_stack.pop_back();
            } else {
                int ei = adj[v][it++];
                used_edge[ei] = 1;
                int to = (edges[ei].u == v) ? edges[ei].v : edges[ei].u;
                node_stack.push_back(to);
            }
        }
        std::reverse(seg.begin(), seg.end());
        node_path.insert(node_path.end(), seg.begin(), seg.end());
    }

    // ---- 6. rebuild ordered moves from the node walk -------------------------
    // For each consecutive node pair, find the matching unused-in-output edge.
    std::vector<char> emitted(edges.size(), 0);
    auto find_edge = [&](int a, int b) -> int {
        for (int ei : adj[a]) if (!emitted[ei]) {
            const GEdge& e = edges[ei];
            if ((e.u == a && e.v == b) || (e.u == b && e.v == a)) return ei;
        }
        return -1;
    };
    for (size_t i = 1; i < node_path.size(); ++i) {
        int a = node_path[i - 1], b = node_path[i];
        int ei = find_edge(a, b);
        if (ei < 0) continue;
        emitted[ei] = 1;
        const GEdge& e = edges[ei];
        Move m;
        if (e.bridge) {
            m.polyline.points = { nodes[a], nodes[b] };
            double gap2 = dist2(nodes[a], nodes[b]);
            m.travel = gap2 > double(params.sacrificial_max) * double(params.sacrificial_max);
            m.bridge = !m.travel;
        } else {
            m.polyline = e.poly;
            m.src_index = e.src;
            if (m.polyline.points.front() != nodes[a] && m.polyline.points.back() == nodes[a]) {
                m.polyline.reverse();
                m.reversed = true;
            }
            m.travel = false; m.bridge = false;
        }
        out.push_back(std::move(m));
    }
    return out;
}

Polylines order_polylines(const Polylines& input, const Params& params)
{
    Polylines pls;
    for (Move& m : order(input, params)) pls.push_back(std::move(m.polyline));
    return pls;
}

std::vector<std::pair<int, bool>> order_indices(const Polylines& endpoints, const Params& params)
{
    std::vector<std::pair<int, bool>> seq;
    if (endpoints.empty())
        return seq;
    std::vector<Move> moves = order(endpoints, params);
    std::vector<char> placed(endpoints.size(), 0);
    for (const Move& m : moves) {
        if (m.src_index < 0 || m.src_index >= int(endpoints.size()) || placed[m.src_index])
            continue;
        placed[m.src_index] = 1;
        seq.emplace_back(m.src_index, m.reversed);
    }
    for (size_t i = 0; i < endpoints.size(); ++i)   // safety: cover any missed item
        if (! placed[i]) seq.emplace_back(int(i), false);
    return seq;
}

static Polylines endpoint_segments(const ExtrusionEntitiesPtr& entities)
{
    // as_polyline() throws on ExtrusionEntityCollection; endpoints suffice for chaining.
    Polylines polys;
    polys.reserve(entities.size());
    for (const ExtrusionEntity* e : entities) {
        Polyline pl;
        pl.points.push_back(e->first_point());
        pl.points.push_back(e->last_point());
        polys.push_back(std::move(pl));
    }
    return polys;
}

void order_entities(ExtrusionEntitiesPtr& entities, const Point& start, const Params& params)
{
    if (entities.size() < 2)
        return;
    auto seq = order_indices(endpoint_segments(entities), params);
    ExtrusionEntitiesPtr ordered;
    ordered.reserve(entities.size());
    for (const auto& s : seq) {
        ExtrusionEntity* e = entities[s.first];
        if (s.second) e->reverse();
        ordered.push_back(e);
    }
    entities = std::move(ordered);
    (void)start;
}

} // namespace ContinuousToolpath
} // namespace Slic3r
