#ifndef SWE_PRE_INIT_DATA_HPP
#define SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"

namespace SWE {
template <typename MeshType>
void initialize_data(MeshType& mesh, const SWE::Inputs& problem_specific_input) {
    mesh.CallForEachElement([&problem_specific_input](auto& elt) {
        elt.data.initialize();

        auto& shape = elt.GetShape();

        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        uint nnode = elt.data.get_nnode();
        uint ngp   = elt.data.get_ngp_internal();

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        state.aux = elt.L2ProjectionNode(bathymetry);

        row(internal.aux_at_gp, SWE::Auxiliaries::bath) = elt.ComputeNodalUgp(bathymetry);
        row(internal.dbath_at_gp, GlobalCoord::x)       = elt.ComputeNodalDUgp(GlobalCoord::x, bathymetry);
        row(internal.dbath_at_gp, GlobalCoord::y)       = elt.ComputeNodalDUgp(GlobalCoord::y, bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < nnode; node_id++) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = elt.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; ++gp) {
                internal.aux_at_gp(SWE::Auxiliaries::sp, gp) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                internal.aux_at_gp(SWE::Auxiliaries::sp, gp) = 1.0;
            }
        }

        if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Constant ||
            problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Default) {
            uint nnode = elt.GetShape().nodal_coordinates.size();

            DynRowVector<double> u_node(nnode);

            set_constant(u_node, problem_specific_input.initial_conditions.ze_initial);
            state.q[SWE::Variables::ze] = elt.L2ProjectionNode(u_node);

            set_constant(u_node, problem_specific_input.initial_conditions.qx_initial);
            state.q[SWE::Variables::qx] = elt.L2ProjectionNode(u_node);

            set_constant(u_node, problem_specific_input.initial_conditions.qy_initial);
            state.q[SWE::Variables::qy] = elt.L2ProjectionNode(u_node);

        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            DynMatrix<double> tmp = elt.L2ProjectionF([](Point<2>& pt) { return SWE::ic_q(0, pt); });
            for ( uint var = 0; var < SWE::n_variables; ++var ) {
                state.q[var] = row(tmp,var);
            }
        }
    });

    mesh.CallForEachInterface([&problem_specific_input](auto& intface) {
        auto& shape_in = intface.GetShapeIN();
        auto& shape_ex = intface.GetShapeEX();

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        uint nnode = intface.data_in.get_nnode();
        uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary_in.aux_at_gp[SWE::Auxiliaries::bath] = intface.ComputeNodalUgpIN(bathymetry);

        // bathymetry is not necessarily continuous across weir nodes
        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape_ex.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary_ex.aux_at_gp[SWE::Auxiliaries::bath] = intface.ComputeNodalUgpEX(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < intface.data_in.get_nnode(); node_id++) {
                y_node[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = intface.ComputeNodalUgpIN(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for ( uint gp = 0; gp < ngp; ++gp ) {
                boundary_in.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / cos_vec(y_at_gp[gp] / R);
                boundary_ex.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / cos_vec(y_at_gp[gp] / R);
            }

        } else {

            set_constant(boundary_in.aux_at_gp[SWE::Auxiliaries::sp], 1.0);
            set_constant(boundary_ex.aux_at_gp[SWE::Auxiliaries::sp], 1.0);

        }
    });

    mesh.CallForEachBoundary([&problem_specific_input](auto& bound) {
        auto& shape = bound.GetShape();

        auto& boundary = bound.data.boundary[bound.bound_id];

        uint nnode = bound.data.get_nnode();
        uint ngp   = bound.data.get_ngp_boundary(bound.bound_id);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary.aux_at_gp[SWE::Auxiliaries::bath] = bound.ComputeNodalUgp(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < bound.data.get_nnode(); node_id++) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = bound.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = 1.0;
            }
        }
    });

    mesh.CallForEachDistributedBoundary([&problem_specific_input](auto& dbound) {
        auto& shape = dbound.GetShape();

        auto& boundary = dbound.data.boundary[dbound.bound_id];

        uint nnode = dbound.data.get_nnode();
        uint ngp   = dbound.data.get_ngp_boundary(dbound.bound_id);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary.aux_at_gp[SWE::Auxiliaries::bath] = dbound.ComputeNodalUgp(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < dbound.data.get_nnode(); node_id++) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = dbound.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = 1.0;
            }
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
            const std::vector<uint>& node_ID = elt.GetNodeID();

            for (uint node = 0; node < elt.data.get_nnode(); ++node) {
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

            for (uint node = 0; node < elt.data.get_nnode(); ++node) {
                y_avg += elt.GetShape().nodal_coordinates[node][GlobalCoord::y];
            }

            y_avg /= elt.data.get_nnode();

            // If spherical projection is not specified we use default value for R which is R earth
            double latitude_avg = y_avg / problem_specific_input.spherical_projection.R;

            elt.data.source.coriolis_f = 2 * 7.2921E-5 * std::sin(latitude_avg);
        });
    }
}
}

#endif