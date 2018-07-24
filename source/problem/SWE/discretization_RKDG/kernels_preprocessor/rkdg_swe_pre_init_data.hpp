#ifndef RKDG_SWE_PRE_INIT_DATA_HPP
#define RKDG_SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

namespace SWE {
namespace RKDG {
void Problem::initialize_data_kernel(ProblemMeshType& mesh,
                                     const MeshMetaData& mesh_data,
                                     const ProblemInputType& problem_specific_input) {
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

        row(state.aux, SWE::Auxiliaries::bath) = elt.L2ProjectionNode(bathymetry);

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

            for (uint gp = 0; gp < ngp; gp++) {
                internal.aux_at_gp(SWE::Auxiliaries::sp, gp) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; gp++) {
                internal.aux_at_gp(SWE::Auxiliaries::sp, gp) = 1.0;
            }
        }

        if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Constant ||
            problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Default) {
            uint nnode = elt.GetShape().nodal_coordinates.size();

            HybMatrix<double, SWE::n_variables> u_node(SWE::n_variables, nnode);

            for (uint node_id = 0; node_id < nnode; ++node_id) {
                u_node(SWE::Variables::ze, node_id) = problem_specific_input.initial_conditions.ze_initial;
                u_node(SWE::Variables::qx, node_id) = problem_specific_input.initial_conditions.qx_initial;
                u_node(SWE::Variables::qy, node_id) = problem_specific_input.initial_conditions.qy_initial;
            }

            state.q = elt.L2ProjectionNode(u_node);
        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            auto u_init = [](Point<2>& pt) { return SWE::ic_u(0, pt); };

            state.q = elt.L2ProjectionF(u_init);
        }
    });

    mesh.CallForEachInterface([&problem_specific_input](auto& intface) {
        auto& shape_in = intface.GetShapeIN();

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        uint nnode = intface.data_in.get_nnode();
        uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::z];
        }

        row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath) = intface.ComputeNodalUgpIN(bathymetry);

        uint gp_ex;
        for (uint gp = 0; gp < ngp; gp++) {
            gp_ex = ngp - gp - 1;

            boundary_ex.aux_at_gp(SWE::Auxiliaries::bath, gp_ex) = boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp);
        }

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < intface.data_in.get_nnode(); node_id++) {
                y_node[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = intface.ComputeNodalUgpIN(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            uint gp_ex;
            for (uint gp = 0; gp < ngp; gp++) {
                gp_ex = ngp - gp - 1;

                boundary_in.aux_at_gp(SWE::Auxiliaries::sp, gp)    = cos_phi_o / std::cos(y_at_gp[gp] / R);
                boundary_ex.aux_at_gp(SWE::Auxiliaries::sp, gp_ex) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            uint gp_ex;
            for (uint gp = 0; gp < ngp; gp++) {
                gp_ex = ngp - gp - 1;

                boundary_in.aux_at_gp(SWE::Auxiliaries::sp, gp)    = 1.0;
                boundary_ex.aux_at_gp(SWE::Auxiliaries::sp, gp_ex) = 1.0;
            }
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

        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) = bound.ComputeNodalUgp(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < bound.data.get_nnode(); node_id++) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = bound.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; gp++) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; gp++) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = 1.0;
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

    // WETTING-DRYING INITIALIZE
    mesh.CallForEachElement([](auto& elt) {
        auto& shape = elt.GetShape();

        auto& state    = elt.data.state[0];
        auto& wd_state = elt.data.wet_dry_state;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.bath_at_vrtx[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
        }

        wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());

        wd_state.q_lin = elt.ProjectBasisToLinear(state.q);

        wd_state.q_at_vrtx = elt.ComputeLinearUvrtx(wd_state.q_lin);

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.h_at_vrtx[vrtx] = wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) + wd_state.bath_at_vrtx[vrtx];
        }

        bool set_wet_element = true;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            if (wd_state.h_at_vrtx[vrtx] <= PostProcessing::h_o + PostProcessing::h_o_threshold) {
                wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) = PostProcessing::h_o - wd_state.bath_at_vrtx[vrtx];

                set_wet_element = false;
            }
        }

        if (set_wet_element) {
            wd_state.wet = true;
        } else {
            wd_state.wet = false;

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                wd_state.q_at_vrtx(SWE::Variables::qx, vrtx) = 0.0;
                wd_state.q_at_vrtx(SWE::Variables::qy, vrtx) = 0.0;
            }

            state.q = elt.ProjectLinearToBasis(elt.data.get_ndof(), wd_state.q_at_vrtx);

            set_constant(state.rhs, 0.0);
        }
    });

    // SLOPE LIMIT INITIALIZE
    mesh.CallForEachElement([](auto& elt) {
        auto& shape = elt.GetShape();

        auto& sl_state = elt.data.slope_limit_state;

        DynRowVector<double> bath_lin(elt.data.get_nvrtx());

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            bath_lin[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
        }

        sl_state.bath_at_baryctr = elt.ComputeLinearUbaryctr(bath_lin);
        sl_state.bath_at_midpts  = elt.ComputeLinearUmidpts(bath_lin);

        sl_state.baryctr_coord = elt.GetShape().GetBarycentricCoordinates();
        sl_state.midpts_coord  = elt.GetShape().GetMidpointCoordinates();

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            sl_state.surface_normal[bound] = elt.GetShape().GetSurfaceNormal(bound, std::vector<Point<2>>(0))[0];
        }
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.baryctr_coord_neigh[intface.bound_id_in] = sl_state_ex.baryctr_coord;
        sl_state_ex.baryctr_coord_neigh[intface.bound_id_ex] = sl_state_in.baryctr_coord;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::x] =
            2.0 * sl_state.midpts_coord[bound.bound_id][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::y] =
            2.0 * sl_state.midpts_coord[bound.bound_id][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2));
        std::vector<double> b(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % elt.data.get_nbound();

            if (!is_distributed(elt.GetBoundaryType()[element_1]) &&
                !is_distributed(elt.GetBoundaryType()[element_2])) {
                A[0][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
                A[0][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

                sl_state.alpha_1[bound] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
                sl_state.alpha_2[bound] = (-A[1][0] * b[0] + A[0][0] * b[1]) / det;

                sl_state.r_sq[bound] = std::pow(b[0], 2.0) + std::pow(b[1], 2.0);
            }
        }
    });
}

void Problem::initialize_data_parallel_pre_send_kernel(ProblemMeshType& mesh,
                                                       const MeshMetaData& mesh_data,
                                                       const ProblemInputType& problem_specific_input) {
    initialize_data_kernel(mesh, mesh_data, problem_specific_input);

    mesh.CallForEachDistributedBoundary([&problem_specific_input](auto& dbound) {
        auto& shape = dbound.GetShape();

        auto& boundary = dbound.data.boundary[dbound.bound_id];

        uint nnode = dbound.data.get_nnode();
        uint ngp   = dbound.data.get_ngp_boundary(dbound.bound_id);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) = dbound.ComputeNodalUgp(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < dbound.data.get_nnode(); node_id++) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = dbound.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; gp++) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; gp++) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = 1.0;
            }
        }
    });

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& sl_state = dbound.data.slope_limit_state;

        dbound.boundary_condition.exchanger.SetPreprocEX(sl_state.baryctr_coord[GlobalCoord::x],
                                                         sl_state.baryctr_coord[GlobalCoord::y]);
    });
}

void Problem::initialize_data_parallel_post_receive_kernel(ProblemMeshType& mesh) {
    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& sl_state = dbound.data.slope_limit_state;

        dbound.boundary_condition.exchanger.GetPreprocEX(sl_state.baryctr_coord_neigh[dbound.bound_id][GlobalCoord::x],
                                                         sl_state.baryctr_coord_neigh[dbound.bound_id][GlobalCoord::y]);
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2));
        std::vector<double> b(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % elt.data.get_nbound();

            if (is_distributed(elt.GetBoundaryType()[element_1]) || is_distributed(elt.GetBoundaryType()[element_2])) {
                A[0][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
                A[0][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

                sl_state.alpha_1[bound] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
                sl_state.alpha_2[bound] = (-A[1][0] * b[0] + A[0][0] * b[1]) / det;

                sl_state.r_sq[bound] = std::pow(b[0], 2.0) + std::pow(b[1], 2.0);
            }
        }
    });
}
}
}

#endif