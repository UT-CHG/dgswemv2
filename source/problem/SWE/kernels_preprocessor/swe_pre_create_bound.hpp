#ifndef SWE_PRE_CREATE_BOUND_HPP
#define SWE_PRE_CREATE_BOUND_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_boundaries_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    InputParameters<ProblemInputType>& input,
    Writer<SWE::Problem>& writer) {
    // *** //
    using BoundaryTypes = Geometry::BoundaryTypeTuple<SWE::Data, SWE::BC::Land, SWE::BC::Tidal, SWE::BC::Flow>;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        if (it->first == SWE::BoundaryTypes::land) {
            using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;

            uint n_bound_old_land = mesh.GetNumberBoundaries();

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                mesh.template CreateBoundary<BoundaryTypeLand>(raw_boundary);

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of land boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_land
                                    << std::endl;
            }
        } else if (it->first == SWE::BoundaryTypes::tidal) {
            using BoundaryTypeTidal = typename std::tuple_element<1, BoundaryTypes>::type;

            uint n_bound_old_tidal = mesh.GetNumberBoundaries();

            auto& tidal_con_data = input.problem_input.tidal_bc_con_data;
            auto& tidal_data     = input.problem_input.tidal_bc_data;

            Array2D<double> amplitude;
            Array2D<double> phase;

            amplitude.resize(tidal_con_data.size());
            phase.resize(tidal_con_data.size());

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                for (uint con = 0; con < tidal_con_data.size(); con++) {
                    amplitude[con].resize(raw_boundary.node_ID.size());
                    phase[con].resize(raw_boundary.node_ID.size());
                }

                for (uint node = 0; node < raw_boundary.node_ID.size(); node++) {
                    uint node_ID = raw_boundary.node_ID[node];

                    if (tidal_data.find(node_ID) != tidal_data.end()) {
                        const auto& tidal = tidal_data.find(node_ID);

                        for (uint con = 0; con < tidal_con_data.size(); con++) {
                            amplitude[con][node] = tidal->second[con][0];
                            phase[con][node]     = tidal->second[con][1];
                        }
                    } else {
                        throw std::logic_error("Error: Unable to find tidal data\n");
                    }
                }

                mesh.template CreateBoundary<BoundaryTypeTidal>(raw_boundary,
                                                                SWE::BC::Tidal(tidal_con_data, amplitude, phase));

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of tidal boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_tidal
                                    << std::endl;
            }
        } else if (it->first == SWE::BoundaryTypes::flow) {
            using BoundaryTypeFlow = typename std::tuple_element<2, BoundaryTypes>::type;

            uint n_bound_old_flow = mesh.GetNumberBoundaries();

            auto& flow_con_data = input.problem_input.flow_bc_con_data;
            auto& flow_data     = input.problem_input.flow_bc_data;

            Array2D<double> amplitude;
            Array2D<double> phase;

            amplitude.resize(flow_con_data.size());
            phase.resize(flow_con_data.size());

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                for (uint con = 0; con < flow_con_data.size(); con++) {
                    amplitude[con].resize(raw_boundary.node_ID.size());
                    phase[con].resize(raw_boundary.node_ID.size());
                }

                for (uint node = 0; node < raw_boundary.node_ID.size(); node++) {
                    uint node_ID = raw_boundary.node_ID[node];

                    if (flow_data.find(node_ID) != flow_data.end()) {
                        const auto& flow = flow_data.find(node_ID);

                        for (uint con = 0; con < flow_con_data.size(); con++) {
                            amplitude[con][node] = flow->second[con][0];
                            phase[con][node]     = flow->second[con][1];
                        }
                    } else {
                        throw std::logic_error("Error: Unable to find flow data\n");
                    }
                }

                mesh.template CreateBoundary<BoundaryTypeFlow>(raw_boundary,
                                                               SWE::BC::Flow(flow_con_data, amplitude, phase));

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of flow boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_flow
                                    << std::endl;
            }
        }
    }

    mesh.CallForEachBoundary([](auto& bound) { bound.boundary_condition.Initialize(bound); });
}
}

#endif