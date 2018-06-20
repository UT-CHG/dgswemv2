#ifndef EHDG_SWE_PRE_INIT_DATA_HPP
#define EHDG_SWE_PRE_INIT_DATA_HPP

#include "utilities//file_exists.hpp"

namespace SWE {
namespace EHDG {
void Problem::initialize_data_kernel(ProblemMeshType& mesh,
                                     const MeshMetaData& mesh_data,
                                     const ProblemInputType& problem_specific_input) {
    mesh.CallForEachElement([](auto& elt) { elt.data.initialize(); });

    std::unordered_map<uint, std::vector<double>> bathymetry;

    std::vector<double> bathymetry_temp;
    for (const auto& elt : mesh_data.elements) {
        auto nodal_coordinates = mesh_data.get_nodal_coordinates(elt.first);

        for (auto& node_coordinate : nodal_coordinates) {
            bathymetry_temp.push_back(node_coordinate[GlobalCoord::z]);
        }

        bathymetry.insert({elt.first, bathymetry_temp});

        bathymetry_temp.clear();
    }

    mesh.CallForEachElement([&bathymetry, &problem_specific_input](auto& elt) {
        uint id = elt.GetID();

        auto& shape = elt.GetShape();

        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        if (!bathymetry.count(id)) {
            throw std::logic_error("Fatal Error: could not find bathymetry for element with id: " + std::to_string(id) +
                                   "!\n");
        }

        elt.L2Projection(bathymetry[id], state.bath);

        elt.ComputeUgp(state.bath, internal.bath_at_gp);

        elt.ComputeDUgp(GlobalCoord::x, state.bath, internal.bath_deriv_wrt_x_at_gp);
        elt.ComputeDUgp(GlobalCoord::y, state.bath, internal.bath_deriv_wrt_y_at_gp);

        if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Constant) {
            uint n_node = elt.GetShape().nodal_coordinates.size();

            std::vector<double> ze_node(n_node, problem_specific_input.initial_conditions.ze_initial);
            elt.L2Projection(ze_node, state.ze);

            std::vector<double> qx_node(n_node, problem_specific_input.initial_conditions.qx_initial);
            elt.L2Projection(qx_node, state.qx);

            std::vector<double> qy_node(n_node, problem_specific_input.initial_conditions.qy_initial);
            elt.L2Projection(qy_node, state.qy);
        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            auto ze_init = [](Point<2>& pt) { return SWE::ic_ze(0, pt); };
            elt.L2Projection(ze_init, state.ze);

            auto qx_init = [](Point<2>& pt) { return SWE::ic_qx(0, pt); };
            elt.L2Projection(qx_init, state.qx);

            auto qy_init = [](Point<2>& pt) { return SWE::ic_qy(0, pt); };
            elt.L2Projection(qy_init, state.qy);
        }
    });

    mesh.CallForEachInterface([&problem_specific_input](auto& intface) {
        auto& state_in = intface.data_in.state[0];
        auto& state_ex = intface.data_ex.state[0];

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        intface.ComputeUgpIN(state_in.bath, boundary_in.bath_at_gp);
        intface.ComputeUgpEX(state_ex.bath, boundary_ex.bath_at_gp);
    });

    mesh.CallForEachBoundary([&problem_specific_input](auto& bound) {
        auto& state    = bound.data.state[0];
        auto& boundary = bound.data.boundary[bound.bound_id];

        bound.ComputeUgp(state.bath, boundary.bath_at_gp);
    });

    // SOURCE TERMS INITIALIZE
    // Manning N bottom friction
    if (problem_specific_input.bottom_friction.type == SWE::BottomFrictionType::Manning) {
        if (!Utilities::file_exists(problem_specific_input.bottom_friction.manning_data_file)) {
            throw std::logic_error("Fatal Error: manning factor data file " +
                                   problem_specific_input.bottom_friction.manning_data_file + " was not found!\n");
        }

        std::ifstream manning_file(problem_specific_input.bottom_friction.manning_data_file);

        uint node_ID;
        double manning_n;
        std::map<uint, double> node_manning_n;

        std::string line;
        while (std::getline(manning_file, line)) {
            std::istringstream input_string(line);

            if (!(input_string >> node_ID >> manning_n))
                break;

            node_manning_n[node_ID] = manning_n;
        }

        mesh.CallForEachElement([&node_manning_n](auto& elt) {
            std::vector<uint>& node_ID = elt.GetNodeID();

            for (uint node = 0; node < elt.data.get_nnode(); node++) {
                elt.data.source.manning_n[node] = node_manning_n[node_ID[node]];
            }

            elt.data.source.manning = true;
            elt.data.source.g_manning_n_sq =
                Global::g *
                std::pow(std::accumulate(elt.data.source.manning_n.begin(), elt.data.source.manning_n.end(), 0.0) /
                             elt.data.source.manning_n.size(),
                         2);
        });
    }

    // Coriolis effect
    if (problem_specific_input.coriolis.type == SWE::CoriolisType::Enable) {
        mesh.CallForEachElement([&problem_specific_input](auto& elt) {
            double y_avg = 0;

            for (uint node = 0; node < elt.data.get_nnode(); node++) {
                y_avg += elt.GetShape().nodal_coordinates[node][GlobalCoord::y];
            }

            y_avg /= elt.data.get_nnode();

            // If spherical projection is not specified we use default value for R which is R earth
            double latitude_avg = y_avg / problem_specific_input.spherical_projection.R;

            elt.data.source.coriolis_f = 2 * 7.2921E-5 * std::sin(latitude_avg);
        });
    }
}

void Problem::initialize_data_parallel_pre_send_kernel(ProblemMeshType& mesh,
                                                       const MeshMetaData& mesh_data,
                                                       const ProblemInputType& problem_specific_input) {
    initialize_data_kernel(mesh, mesh_data, problem_specific_input);

    mesh.CallForEachDistributedBoundary([&problem_specific_input](auto& dbound) {
        auto& state    = dbound.data.state[0];
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        dbound.ComputeUgp(state.bath, boundary.bath_at_gp);
    });
}

void Problem::initialize_data_parallel_post_receive_kernel(ProblemMeshType& mesh) {}
}
}

#endif