#ifndef SWE_PRE_CREATE_INTFACE_HPP
#define SWE_PRE_CREATE_INTFACE_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_interfaces_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    InputParameters<ProblemInputType>& input,
    Writer<SWE::Problem>& writer) {
    // *** //
    using InterfaceTypes = Geometry::InterfaceTypeTuple<SWE::Data, SWE::IS::Interface, SWE::IS::Levee>;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        switch (it->first) {
            case SWE::BoundaryConditions::internal:
                using InterfaceTypeInterface = std::tuple_element<0, InterfaceTypes>::type;

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                    if (it->second.find(key_pre_int_ex) != it->second.end()) {
                        auto& raw_boundary_in = itt->second;
                        auto& raw_boundary_ex = it->second.find(key_pre_int_ex)->second;

                        mesh.template CreateInterface<InterfaceTypeInterface>(raw_boundary_in, raw_boundary_ex);
                    }

                    it->second.erase(itt);
                }

                break;
            case SWE::BoundaryConditions::internal_barrier:
                using InterfaceTypeLevee = std::tuple_element<1, InterfaceTypes>::type;

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                    if (it->second.find(key_pre_int_ex) != it->second.end()) {
                        auto& raw_boundary_in = itt->second;
                        auto& raw_boundary_ex = it->second.find(key_pre_int_ex)->second;

                        auto& levee_data = input.problem_input.levee_is_data;

                        std::vector<double> H_barrier(raw_boundary_in.node_ID.size());
                        std::vector<double> C_subrictical(raw_boundary_in.node_ID.size());
                        std::vector<double> C_superctical(raw_boundary_in.node_ID.size());

                        for (uint node = 0; node < raw_boundary_in.node_ID.size(); node++) {
                            std::pair<uint, uint> key_levee_data = std::pair<uint, uint>{
                                raw_boundary_in.node_ID[node],
                                raw_boundary_ex.node_ID[raw_boundary_in.node_ID.size() - node - 1]};

                            std::pair<uint, uint> key_levee_data_swap;
                            key_levee_data_swap.swap(key_levee_data);

                            if (levee_data.find(key_levee_data) != levee_data.end()) {
                                H_barrier[node]     = input.problem_input.levee_is_data.find(key_levee_data)->second[0];
                                C_subrictical[node] = input.problem_input.levee_is_data.find(key_levee_data)->second[1];
                                C_superctical[node] = input.problem_input.levee_is_data.find(key_levee_data)->second[2];
                            } else if (levee_data.find(key_levee_data_swap) != levee_data.end()) {
                                H_barrier[node] =
                                    input.problem_input.levee_is_data.find(key_levee_data_swap)->second[0];
                                C_subrictical[node] =
                                    input.problem_input.levee_is_data.find(key_levee_data_swap)->second[1];
                                C_superctical[node] =
                                    input.problem_input.levee_is_data.find(key_levee_data_swap)->second[2];
                            } else {
                                throw std::logic_error("Error: Unable to find levee data\n");
                            }
                        }

                        mesh.template CreateInterface<InterfaceTypeLevee>(
                            raw_boundary_in, raw_boundary_ex, SWE::IS::Levee(H_barrier, C_subrictical, C_superctical));
                    }

                    it->second.erase(itt);
                }

                break;
        }
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of interfaces: " << mesh.GetNumberInterfaces() << std::endl;
    }
}
}

#endif