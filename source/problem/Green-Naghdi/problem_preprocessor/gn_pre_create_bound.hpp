#ifndef GN_PRE_CREATE_BOUND_HPP
#define GN_PRE_CREATE_BOUND_HPP

namespace GN {
template <typename ProblemType, typename RawBoundaryType>
void create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                       typename ProblemType::ProblemMeshType& mesh,
                       typename ProblemType::ProblemInputType& problem_input,
                       typename ProblemType::ProblemWriterType& writer) {
    // *** //
    using BoundaryTypeLand     = typename std::tuple_element<0, typename ProblemType::ProblemBoundaryTypes>::type;
    using BoundaryTypeTide     = typename std::tuple_element<1, typename ProblemType::ProblemBoundaryTypes>::type;
    using BoundaryTypeFlow     = typename std::tuple_element<2, typename ProblemType::ProblemBoundaryTypes>::type;
    using BoundaryTypeFunction = typename std::tuple_element<3, typename ProblemType::ProblemBoundaryTypes>::type;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); ++it) {
        if (it->first == GN::BoundaryTypes::land) {
            uint n_bound_old_land = mesh.GetNumberBoundaries();

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                mesh.template CreateBoundary<BoundaryTypeLand>(std::move(raw_boundary));

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of land boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_land
                                    << std::endl;
            }
        } else if (it->first == GN::BoundaryTypes::tide) {
            uint n_bound_old_tide = mesh.GetNumberBoundaries();

            auto& tide_data = problem_input.tide_bc_data;

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                std::vector<SWE::TideNode> tide;

                bool found_data = false;

                for (auto& tide_bound : tide_data) {
                    found_data = tide_bound.get_tide_data(raw_boundary.node_ID, tide);

                    if (found_data)
                        break;
                }

                if (!found_data)
                    throw std::logic_error("Fatal Error: unable to find tide data!\n");

                mesh.template CreateBoundary<BoundaryTypeTide>(std::move(raw_boundary), tide);

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of tide boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_tide
                                    << std::endl;
            }
        } else if (it->first == GN::BoundaryTypes::flow) {
            uint n_bound_old_flow = mesh.GetNumberBoundaries();

            auto& flow_data = problem_input.flow_bc_data;

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                std::vector<SWE::FlowNode> flow;

                bool found_data = false;

                for (auto& flow_bound : flow_data) {
                    found_data = flow_bound.get_flow_data(raw_boundary.node_ID, flow);

                    if (found_data)
                        break;
                }

                if (!found_data)
                    throw std::logic_error("Fatal Error: unable to find flow data!\n");

                mesh.template CreateBoundary<BoundaryTypeFlow>(std::move(raw_boundary), flow);

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of flow boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_flow
                                    << std::endl;
            }
        } else if (it->first == GN::BoundaryTypes::function) {
            uint n_bound_old_func = mesh.GetNumberBoundaries();

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                mesh.template CreateBoundary<BoundaryTypeFunction>(std::move(raw_boundary));

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of function boundaries: "
                                    << mesh.GetNumberBoundaries() - n_bound_old_func << std::endl;
            }
        }
    }

    mesh.CallForEachBoundary([](auto& bound) { bound.boundary_condition.Initialize(bound); });
}
}

#endif
