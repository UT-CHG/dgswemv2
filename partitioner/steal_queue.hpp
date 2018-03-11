#ifndef STEAL_QUEUE_HPP
#define STEAL_QUEUE_HPP
#include "csrmat.hpp"

#include <deque>

class StealQ {
public:
    StealQ(const CSRMat<>& g, const PartitionType& p, std::deque<int>& high_elements);

    int steal_one();
private:
    int64_t high_partition;
    std::unordered_map<int,double> element_gains;

    const CSRMat<>& g;
    const PartitionType& p;
};

StealQ::StealQ(const CSRMat<>& g, const PartitionType& p, std::deque<int>& high_elements) : g(g), p(p) {
    high_partition = p.vertex2partition.at(high_elements.front());

    for ( int elt : high_elements ) {
        element_gains.insert(std::make_pair(elt,0.));
        for ( int neigh : g.get_neighbors(elt) ) {
            if ( p.vertex2partition.count(neigh) ) {
                int factor = ( high_partition == p.vertex2partition.at(neigh) ) ? -1 : 1;

                std::pair<int,int> edge_name{ std::min(elt, neigh), std::max(elt,neigh) };

                element_gains[elt] += factor * g.edge_weight(edge_name);
            }
        }
    }
}

int StealQ::steal_one() {
    double max_gain{std::numeric_limits<double>::lowest()};
    int elt;

    for ( auto& e_g : element_gains ) {
        if ( e_g.second > max_gain ) {
            elt = e_g.first;
            max_gain = e_g.second;
        }
    }

    //update gains for elements on high partition side
    for ( int neigh : g.get_neighbors(elt) ) {
        if ( p.vertex2partition.count(neigh) ) {
            if ( p.vertex2partition.at(neigh) == high_partition ) {
                std::pair<int,int> edge_name{ std::min(elt, neigh), std::max(elt,neigh) };
                element_gains[neigh] += 2 * g.edge_weight(edge_name);
            }
        }
    }

    //If this is the case, you didn't find a file to steal
    assert( max_gain != std::numeric_limits<double>::lowest() && !element_gains.empty() );

    element_gains.erase(elt);
    return elt;
}

#endif