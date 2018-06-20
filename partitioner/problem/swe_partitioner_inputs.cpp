#include "swe_partitioner_inputs.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"

namespace SWE {
PartitionerInputs::PartitionerInputs(const MeshMetaData& mesh, Inputs inputs) {
    // To determine whether an element is wet or dry for initial load balancing purposes,
    // we simply call the element dry if at least one vertex is below h_o. This may be incorrect
    // for higher orders, but should at least provide a good approximation.
    for (const auto& elt : mesh.elements) {
        bool is_wet{true};

        for (uint i = 0; i < elt.second.node_ID.size(); ++i) {
            uint node_ID                = elt.second.node_ID[i];
            const Point<3>& coordinates = mesh.nodes.at(node_ID).coordinates;

            double h_at_vrtx;
            double bath_at_vrtx = coordinates[2];

            if (inputs.initial_conditions.type == SWE::InitialConditionsType::Constant) {
                h_at_vrtx = inputs.initial_conditions.ze_initial + bath_at_vrtx;

            } else if (inputs.initial_conditions.type == SWE::InitialConditionsType::Function) {
                Point<2> pt{coordinates[0], coordinates[1]};
                h_at_vrtx = ic_ze(0, pt) + bath_at_vrtx;
            } else {  // unknown input conditition; set the whole mesh to wet
                h_at_vrtx = 2 * inputs.wet_dry.h_o;
            }

            is_wet &= (h_at_vrtx >= inputs.wet_dry.h_o);
        }

        weights.insert(std::make_pair(elt.first, std::vector<double>{1., is_wet}));
    }

    {
        uint64_t n_elem = weights.size();
        double n_wet{0};
        for (const auto& id_w : weights) {
            n_wet += id_w.second[1];
        }

        std::cout << "Shallow water specific partitioner information\n"
                  << "  Percentage of wet elements: " << n_wet / n_elem * 100 << " %\n\n";
    }
}
}