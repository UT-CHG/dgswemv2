#ifndef SWE_PRE_CREATE_INTFACE_HPP
#define SWE_PRE_CREATE_INTFACE_HPP

namespace SWE {
template <typename ProblemType, typename RawBoundaryType>
void create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                       typename ProblemType::ProblemMeshType& mesh,
                       typename ProblemType::ProblemInputType& problem_input,
                       typename ProblemType::ProblemWriterType& writer) {
    // *** //
    using InterfaceTypeInternal = typename std::tuple_element<0, typename ProblemType::ProblemInterfaceTypes>::type;
    using InterfaceTypeLevee    = typename std::tuple_element<1, typename ProblemType::ProblemInterfaceTypes>::type;

    for ( auto it = raw_boundaries.begin(); it != raw_boundaries.end(); ++it ) {
        if ( it->first == SWE::BoundaryTypes::internal ) {
            size_t n_internal = it->second.size();
            assert(n_internal % 2 == 0 );
            mesh.template reserve_interfaces<InterfaceTypeInternal>(n_internal/2);
        } else if ( it->first == SWE::BoundaryTypes::levee ) {
            size_t n_levee = it->second.size();
            assert(n_levee % 2 == 0 );
            mesh.template reserve_interfaces<InterfaceTypeLevee>(n_levee/2);
        }
    }

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); ++it) {
        if (it->first == SWE::BoundaryTypes::internal) {
            uint n_intface_old_internal = mesh.GetNumberInterfaces();

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                if (it->second.find(key_pre_int_ex) != it->second.end()) {
                    auto& raw_boundary_in = itt->second;
                    auto& raw_boundary_ex = it->second.find(key_pre_int_ex)->second;

                    mesh.template CreateInterface<InterfaceTypeInternal>(std::move(raw_boundary_in),
                                                                         std::move(raw_boundary_ex));
                }

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of internal interfaces: "
                                    << mesh.GetNumberInterfaces() - n_intface_old_internal << std::endl;
            }
        } else if (it->first == SWE::BoundaryTypes::levee) {
            uint n_intface_old_levee = mesh.GetNumberInterfaces();

            auto& levee_data = problem_input.levee_is_data;

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                if (it->second.find(key_pre_int_ex) != it->second.end()) {
                    auto& raw_boundary_in = itt->second;
                    auto& raw_boundary_ex = it->second.find(key_pre_int_ex)->second;

                    std::vector<LeveeInput> levee;

                    for (uint node = 0; node < raw_boundary_in.node_ID.size(); ++node) {
                        std::pair<uint, uint> key_levee_data{
                            raw_boundary_in.node_ID[node],
                            raw_boundary_ex.node_ID[raw_boundary_in.node_ID.size() - node - 1]};

                        std::pair<uint, uint> key_levee_data_swap{
                            raw_boundary_ex.node_ID[raw_boundary_in.node_ID.size() - node - 1],
                            raw_boundary_in.node_ID[node]};

                        if (levee_data.find(key_levee_data) != levee_data.end()) {
                            levee.push_back(levee_data[key_levee_data]);
                        } else if (levee_data.find(key_levee_data_swap) != levee_data.end()) {
                            levee.push_back(levee_data[key_levee_data_swap]);
                        } else {
                            throw std::logic_error("Fatal Error: unable to find levee data!\n");
                        }
                    }

                    mesh.template CreateInterface<InterfaceTypeLevee>(
                        std::move(raw_boundary_in), std::move(raw_boundary_ex), levee);
                }

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of levee interfaces: "
                                    << mesh.GetNumberInterfaces() - n_intface_old_levee << std::endl;
            }
        }
    }

    mesh.CallForEachInterface([](auto& intface) { intface.specialization.Initialize(intface); });
}
}

#endif
