#ifndef EHDG_SWE_PRE_INIT_DATA_HPP
#define EHDG_SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

namespace SWE {
namespace EHDG {
void Problem::initialize_data_kernel(ProblemMeshType& mesh,
                                     const MeshMetaData& mesh_data,
                                     const ProblemInputType& problem_specific_input) {
    mesh.CallForEachElement([&problem_specific_input](auto& elt) {
        elt.data.initialize();

        auto& shape = elt.GetShape();

        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        uint ngp = elt.data.get_ngp_internal();

        std::vector<double> bathymetry;

        for (uint node_id = 0; node_id < elt.data.get_nnode(); node_id++) {
            bathymetry.push_back(shape.nodal_coordinates[node_id][GlobalCoord::z]);
        }

        elt.L2Projection(bathymetry, state.bath);

        std::vector<double> bath_at_gp(ngp);
        std::vector<double> dbath_dx_at_gp(ngp);
        std::vector<double> dbath_dy_at_gp(ngp);

        elt.ComputeNodalUgp(bathymetry, bath_at_gp);
        elt.ComputeNodalDUgp(bathymetry, internal.dbath_at_gp);

        for (uint gp = 0; gp < ngp; gp++) {
            internal.aux_at_gp[gp][SWE::Auxiliaries::bath] = bath_at_gp[gp];
        }

        if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Constant) {
            uint n_node = elt.GetShape().nodal_coordinates.size();

            StatVector<double, SWE::n_variables> u_init{problem_specific_input.initial_conditions.ze_initial,
                                                        problem_specific_input.initial_conditions.qx_initial,
                                                        problem_specific_input.initial_conditions.qy_initial};

            std::vector<StatVector<double, SWE::n_variables>> u_node(n_node, u_init);

            elt.L2Projection(u_node, state.q);
        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            auto u_init = [](Point<2>& pt) { return SWE::ic_u(0, pt); };

            elt.L2Projection(u_init, state.q);
        }
    });

    mesh.CallForEachInterface([&problem_specific_input](auto& intface) {
        auto& shape_in = intface.GetShapeIN();

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);

        std::vector<double> bathymetry;

        for (uint node_id = 0; node_id < intface.data_in.get_nnode(); node_id++) {
            bathymetry.push_back(shape_in.nodal_coordinates[node_id][GlobalCoord::z]);
        }

        std::vector<double> bath_at_gp(ngp);

        intface.ComputeNodalUgpIN(bathymetry, bath_at_gp);

        uint gp_ex;
        for (uint gp = 0; gp < ngp; gp++) {
            gp_ex = ngp - gp - 1;

            boundary_in.aux_at_gp[gp][SWE::Auxiliaries::bath]    = bath_at_gp[gp];
            boundary_ex.aux_at_gp[gp_ex][SWE::Auxiliaries::bath] = bath_at_gp[gp];
        }
    });

    mesh.CallForEachBoundary([&problem_specific_input](auto& bound) {
        auto& shape = bound.GetShape();

        auto& boundary = bound.data.boundary[bound.bound_id];

        uint ngp = bound.data.get_ngp_boundary(bound.bound_id);

        std::vector<double> bathymetry;

        for (uint node_id = 0; node_id < bound.data.get_nnode(); node_id++) {
            bathymetry.push_back(shape.nodal_coordinates[node_id][GlobalCoord::z]);
        }

        std::vector<double> bath_at_gp(ngp);

        bound.ComputeNodalUgp(bathymetry, bath_at_gp);

        for (uint gp = 0; gp < ngp; gp++) {
            boundary.aux_at_gp[gp][SWE::Auxiliaries::bath] = bath_at_gp[gp];
        }
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
            const DynVector<uint>& node_ID = elt.GetNodeID();

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
        auto& shape = dbound.GetShape();

        auto& boundary = dbound.data.boundary[dbound.bound_id];

        uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

        std::vector<double> bathymetry;

        for (uint node_id = 0; node_id < dbound.data.get_nnode(); node_id++) {
            bathymetry.push_back(shape.nodal_coordinates[node_id][GlobalCoord::z]);
        }

        std::vector<double> bath_at_gp(ngp);

        dbound.ComputeNodalUgp(bathymetry, bath_at_gp);

        for (uint gp = 0; gp < ngp; gp++) {
            boundary.aux_at_gp[gp][SWE::Auxiliaries::bath] = bath_at_gp[gp];
        }
    });
}

void Problem::initialize_data_parallel_post_receive_kernel(ProblemMeshType& mesh) {}
}
}

#endif