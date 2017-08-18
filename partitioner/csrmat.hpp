#ifndef _1b833e26_2871_4cce_8c46_dd736f8581e0
#define _1b833e26_2871_4cce_8c46_dd736f8581e0

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <unordered_map>

#include "util.hpp"
#include <metis.h>

template <class NodeW = double, class EdgeW = double>
class CSRMat {
  public:
    // construct from components
    CSRMat(std::unordered_map<int, NodeW> nw, std::unordered_map<std::pair<int, int>, EdgeW> ew)
        : _edges{}, _nodes{}, _node_wgts_map(std::move(nw)), _edge_wgts_map(std::move(ew)) {
        for (const auto& p : _node_wgts_map) {
            _nodes.push_back(p.first);
        }
        std::sort(_nodes.begin(), _nodes.end());

        // use of intermediate edge_sets normalizes all edges to be bi-directional
        std::unordered_map<int, std::unordered_set<int>> edge_sets;
        for (const auto& p : _edge_wgts_map) {
            int src = p.first.first, dst = p.first.second;
            edge_sets[src].insert(dst);
            edge_sets[dst].insert(src);
        }
        for (const auto& p : edge_sets) {
            int src = p.first;
            for (int dst : p.second) {
                _edges[src].push_back(dst);
            }
            std::sort(_edges[src].begin(), _edges[src].end());
        }
    }

    size_t size() const { return _nodes.size(); }
    const std::vector<int>& node_ids() const { return _nodes; }
    int get(int index) const { return _nodes[index]; }
    NodeW node_weight(int id) const { return _node_wgts_map.at(id); }

    std::vector<int> xadj() const {
        std::vector<int> result;
        result.push_back(0);
        for (auto id : _nodes) {
            int edge_n = _edges.count(id) ? _edges.at(id).size() : 0;
            int idx = result.back() + edge_n;
            result.push_back(idx);
        }
        return result;
    }

    std::vector<int> adj() const {
        std::vector<int> result;
        // mapping from node IDs to index [0, num_nodes)
        std::unordered_map<int, int> mapping;
        for (unsigned i = 0; i < _nodes.size(); ++i) {
            mapping[_nodes[i]] = i;
        }
        for (auto src : _nodes) {
            if (_edges.count(src)) {
                for (auto dst : _edges.at(src)) {
                    result.push_back(mapping.at(dst));
                }
            }
        }
        return result;
    }

    // vector of node weights in order
    std::vector<NodeW> node_wgts() const {
        std::vector<NodeW> result;
        result.reserve(size());   // pre-allocate
        for (auto id : _nodes) {  // _nodes is sorted
            result.push_back(_node_wgts_map.at(id));
        }
        return result;
    }

    // returns a vector of edge weights in order
    std::vector<EdgeW> edge_wgts() const {
        std::vector<EdgeW> result;
        for (auto src : _nodes) {                  // _nodes is sorted
            if (_edges.count(src)) {               // node might not be connected to anything
                for (auto dst : _edges.at(src)) {  // _edges.at(src) is sorted
                    EdgeW w{0};
                    auto id = std::make_pair(src, dst);
                    if (_edge_wgts_map.count(id)) {
                        w += _edge_wgts_map.at(id);
                    }
                    id = std::make_pair(dst, src);
                    if (_edge_wgts_map.count(id)) {
                        w += _edge_wgts_map.at(id);
                    }
                    result.push_back(w);
                }
            }
        }
        return result;
    }

    std::vector<int> idxs_to_node_ids(std::vector<int> idxs) const {
        std::vector<int> result;
        for (auto idx : idxs) {
            result.push_back(_nodes[idx]);
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& os, const CSRMat& mat) {
        os << "CSRMat (" << mat.size() << ")\n";
        auto nw = mat.node_wgts();
        auto ew = mat.edge_wgts();
        int nidx = 0, eidx = 0;
        for (auto nid : mat._nodes) {
            os << "(" << nid << "," << nw[nidx++] << ") : ";
            for (auto eid : mat._edges.at(nid)) {
                os << "(" << eid << "," << ew[eidx++] << "),";
            }
            os << '\n';
        }
        return os;
    }

    void csr_info() {
        std::cout << "size() " << _nodes.size() << std::endl;

        std::cout << "First ten elements of _nodes \n";
        for (int j = 0; j < 10; ++j)
            std::cout << _nodes[j] << std::endl;

        std::cout << "First ten elements of _edges \n";
        for (int j = 0; j < 10; ++j) {
            for (uint k = 0; k < _edges[j].size(); ++k)
                std::cout << _edges[j][k] << " ";
            std::cout << "\n";
        }
    }

    std::unordered_map<int, std::vector<int>> get_edges() { return _edges; }

  private:
    std::unordered_map<int, std::vector<int>> _edges;
    std::vector<int> _nodes;  // sorted node IDs
    std::unordered_map<int, NodeW> _node_wgts_map;
    std::unordered_map<std::pair<int, int>, EdgeW> _edge_wgts_map;
};

template <class T>
std::vector<int64_t> to_int64(const std::vector<T>& vec) {
    std::vector<int64_t> result;
    result.reserve(vec.size());
    for (auto val : vec) {
        result.push_back(static_cast<int64_t>(val));
    }
    return result;
}

template <class T>
std::vector<int64_t> scale_to_int64(const std::vector<T>& vec) {
    double max_aval = 0;
    for (auto val : vec) {
        if (fabs(val) > max_aval)
            max_aval = fabs(val);
    }
    double target = std::sqrt(double(std::numeric_limits<int64_t>::max()));
    double scale_factor = target / max_aval;
    std::vector<int64_t> result;
    result.reserve(vec.size());
    for (auto val : vec) {
        result.push_back(static_cast<int64_t>(scale_factor * val));
    }
    return result;
}

template <class T>
void summarize_vec(std::vector<T> vec) {
    mota::Say() << "size=" << vec.size();
    auto min_v = std::numeric_limits<T>::max();
    auto max_v = std::numeric_limits<T>::min();
    for (const auto& v : vec) {
        min_v = std::min(v, min_v);
        max_v = std::max(v, max_v);
    }
    mota::Say() << ", min=" << min_v;
    mota::Say() << ", max=" << max_v;
    mota::Say() << ", (";
    for (size_t i = 0; i < std::min(4ul, vec.size()); ++i) {
        mota::Say() << vec[i] << ", ";
    }
    mota::Say() << "..., ";
    for (size_t i = std::max(0ul, vec.size() - 4); i < vec.size(); ++i) {
        mota::Say() << vec[i] << ", ";
    }
    mota::Say() << "\b\b)";
}

inline std::vector<double> get_node_weights(const std::vector<int>& node_ids,
                                            const std::vector<std::function<double(int)>>& cons) {
    std::vector<double> result;
    result.reserve(node_ids.size() * cons.size());
    for (auto id : node_ids) {
        for (auto lam : cons) {
            result.push_back(lam(id));
        }
    }
    return result;
}

template <class NodeW>
std::vector<int64_t> metis_part(const CSRMat<NodeW>& mat,
                                int64_t nparts,
                                std::vector<std::function<double(int)>> cons,
                                const double imba_ratio) {
    static_assert(sizeof(idx_t) == sizeof(int64_t), "Requires 64-bit METIS idx_t");
    static_assert(sizeof(real_t) == sizeof(double), "Requires 64-bit METIS real_t");

    int64_t nvtxs = mat.size();

    // fast, easy case
    if (nparts == 1 || nvtxs == 1) {
        return std::vector<int64_t>(nvtxs, 0);
    }

    // set up metis parameters
    int64_t ncon = cons.size(), objval;
    std::vector<int64_t> options(METIS_NOPTIONS), part(nvtxs);
    std::vector<double> tpwgts(nparts * ncon, 1.0 / nparts), ubvec(ncon, imba_ratio);
    // set up weights vectorse

    std::vector<int64_t> node_wgts = scale_to_int64(get_node_weights(mat.node_ids(), cons)),
                         edge_wgts = scale_to_int64(mat.edge_wgts());
    // check parameters

    /*mota::Say() << "nparts: " << nparts;
    mota::Say() << "nvtx: " << nvtxs;
    mota::Say() << "xadj: ";
    summarize_vec(mat.xadj());
    mota::Say() << " adj: ";
    summarize_vec(mat.adj());
    mota::Say() << "ncon: " << ncon;*/

    // do partitioning
    METIS_SetDefaultOptions(options.data());
    METIS_PartGraphKway(&nvtxs,
                        &ncon,
                        to_int64(mat.xadj()).data(),
                        to_int64(mat.adj()).data(),
                        node_wgts.data(),
                        nullptr,
                        edge_wgts.data(),
                        &nparts,
                        tpwgts.data(),
                        ubvec.data(),
                        options.data(),
                        &objval,
                        part.data());
    return part;
}

#endif
