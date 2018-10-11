#ifndef EHDG_GN_PRE_CREATE_BOUND_HPP
#define EHDG_GN_PRE_CREATE_BOUND_HPP

namespace GN {
namespace EHDG {
template <typename RawBoundaryType>
void Problem::create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& problem_input,
                                ProblemWriterType& writer) {
    // *** //
    using BoundaryTypes = Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow>;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); ++it) {
        if (it->first == GN::BoundaryTypes::land) {
            using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;

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
            using BoundaryTypeTide = typename std::tuple_element<1, BoundaryTypes>::type;

            uint n_bound_old_tide = mesh.GetNumberBoundaries();

            auto& tide_data = problem_input.tide_bc_data;

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                std::vector<SWE::TideInput> tide;

                for (uint node = 0; node < raw_boundary.node_ID.size(); ++node) {
                    uint node_ID = raw_boundary.node_ID[node];

                    if (tide_data.find(node_ID) != tide_data.end()) {
                        tide.push_back(tide_data[node_ID]);
                    } else {
                        throw std::logic_error("Fatal Error: unable to find tide data!\n");
                    }
                }

                mesh.template CreateBoundary<BoundaryTypeTide>(std::move(raw_boundary), BC::Tide(tide));

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of tide boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_tide
                                    << std::endl;
            }
        } else if (it->first == GN::BoundaryTypes::flow) {
            using BoundaryTypeFlow = typename std::tuple_element<2, BoundaryTypes>::type;

            uint n_bound_old_flow = mesh.GetNumberBoundaries();

            auto& flow_data = problem_input.flow_bc_data;

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                std::vector<SWE::FlowInput> flow;

                for (uint node = 0; node < raw_boundary.node_ID.size(); ++node) {
                    uint node_ID = raw_boundary.node_ID[node];

                    if (flow_data.find(node_ID) != flow_data.end()) {
                        flow.push_back(flow_data[node_ID]);
                    } else {
                        throw std::logic_error("Fatal Error: unable to find flow data!\n");
                    }
                }

                mesh.template CreateBoundary<BoundaryTypeFlow>(std::move(raw_boundary), BC::Flow(flow));

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
}

#endif
