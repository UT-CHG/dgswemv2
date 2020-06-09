#ifndef SWE_PRE_INIT_DATA_HPP
#define SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"
#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"

namespace SWE {
template <typename MeshType, typename ProblemSpecificInputType>
void initialize_data_serial(MeshType& mesh, const ProblemSpecificInputType& problem_specific_input) {
    mesh.CallForEachElement([&problem_specific_input](auto& elt) {
        elt.data.initialize();

        auto& shape    = elt.GetShape();
        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        const uint nnode = elt.data.get_nnode();
        const uint ngp   = elt.data.get_ngp_internal();

        DynRowVector<double> bathymetry(nnode);
        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }
        row(state.aux, SWE::Auxiliaries::bath) = elt.L2ProjectionNode(bathymetry);

        row(internal.aux_at_gp, SWE::Auxiliaries::bath) = elt.ComputeUgp(row(state.aux, SWE::Auxiliaries::bath));
        row(internal.db_at_gp, GlobalCoord::x) =
            elt.ComputeDUgp(GlobalCoord::x, row(state.aux, SWE::Auxiliaries::bath));
        row(internal.db_at_gp, GlobalCoord::y) =
            elt.ComputeDUgp(GlobalCoord::y, row(state.aux, SWE::Auxiliaries::bath));

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);
            for (uint node_id = 0; node_id < nnode; ++node_id) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }
            DynRowVector<double> y_at_gp(ngp);
            y_at_gp = elt.ComputeNodalUgp(y_node);

            const double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            const double R         = problem_specific_input.spherical_projection.R;
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
            const uint nnode = elt.GetShape().nodal_coordinates.size();
            HybMatrix<double, SWE::n_variables> u_node(SWE::n_variables, nnode);
            for (uint node_id = 0; node_id < nnode; ++node_id) {
                u_node(SWE::Variables::ze, node_id) = problem_specific_input.initial_conditions.ze_initial;
                u_node(SWE::Variables::qx, node_id) = problem_specific_input.initial_conditions.qx_initial;
                u_node(SWE::Variables::qy, node_id) = problem_specific_input.initial_conditions.qy_initial;
                u_node(SWE::Variables::hc, node_id) = problem_specific_input.initial_conditions.hc_initial;
            }
            state.q = elt.L2ProjectionNode(u_node);
        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            const auto q_init = [](Point<2>& pt) { return SWE::ic_q(0, pt); };
            state.q           = elt.L2ProjectionF(q_init);
        }
    });

    mesh.CallForEachInterface([&problem_specific_input](auto& intface) {
        auto& shape_in    = intface.GetShapeIN();
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        const uint nnode = intface.data_in.get_nnode();
        const uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);

        row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath) =
            intface.ComputeUgpIN(row(intface.data_in.state[0].aux, SWE::Auxiliaries::bath));
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath) =
            intface.ComputeUgpEX(row(intface.data_ex.state[0].aux, SWE::Auxiliaries::bath));

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);
            for (uint node_id = 0; node_id < intface.data_in.get_nnode(); ++node_id) {
                y_node[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::y];
            }
            DynRowVector<double> y_at_gp(ngp);
            y_at_gp = intface.ComputeNodalUgpIN(y_node);

            const double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            const double R         = problem_specific_input.spherical_projection.R;
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex                                   = ngp - gp - 1;
                boundary_in.aux_at_gp(SWE::Auxiliaries::sp, gp)    = cos_phi_o / std::cos(y_at_gp[gp] / R);
                boundary_ex.aux_at_gp(SWE::Auxiliaries::sp, gp_ex) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex                                   = ngp - gp - 1;
                boundary_in.aux_at_gp(SWE::Auxiliaries::sp, gp)    = 1.0;
                boundary_ex.aux_at_gp(SWE::Auxiliaries::sp, gp_ex) = 1.0;
            }
        }
    });

    mesh.CallForEachBoundary([&problem_specific_input](auto& bound) {
        auto& shape    = bound.GetShape();
        auto& boundary = bound.data.boundary[bound.bound_id];

        const uint nnode = bound.data.get_nnode();
        const uint ngp   = bound.data.get_ngp_boundary(bound.bound_id);

        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) =
            bound.ComputeUgp(row(bound.data.state[0].aux, SWE::Auxiliaries::bath));

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);
            for (uint node_id = 0; node_id < bound.data.get_nnode(); ++node_id) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }
            DynRowVector<double> y_at_gp(ngp);
            y_at_gp = bound.ComputeNodalUgp(y_node);

            const double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            const double R         = problem_specific_input.spherical_projection.R;
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = 1.0;
            }
        }
    });

    mesh.CallForEachDistributedBoundary([&problem_specific_input](auto& dbound) {
        auto& shape    = dbound.GetShape();
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        const uint nnode = dbound.data.get_nnode();
        const uint ngp   = dbound.data.get_ngp_boundary(dbound.bound_id);

        row(boundary.aux_at_gp, SWE::Auxiliaries::bath) =
            dbound.ComputeUgp(row(dbound.data.state[0].aux, SWE::Auxiliaries::bath));

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);
            for (uint node_id = 0; node_id < dbound.data.get_nnode(); ++node_id) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }
            DynRowVector<double> y_at_gp(ngp);
            y_at_gp = dbound.ComputeNodalUgp(y_node);

            const double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            const double R         = problem_specific_input.spherical_projection.R;
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp(SWE::Auxiliaries::sp, gp) = 1.0;
            }
        }
    });

    // SOURCE TERMS INITIALIZE
    // Manning N bottom friction
    if (problem_specific_input.bottom_friction.type == SWE::BottomFrictionType::Manning) {
        if (!problem_specific_input.bottom_friction.manning_data_file.empty()) {
            if (!Utilities::file_exists(problem_specific_input.bottom_friction.manning_data_file)) {
                throw std::logic_error("Fatal Error: manning factor data file " +
                                       problem_specific_input.bottom_friction.manning_data_file + " was not found!\n");
            }

            std::ifstream manning_file(problem_specific_input.bottom_friction.manning_data_file);
            std::map<uint, double> node_manning_n;
            std::string line;
            while (std::getline(manning_file, line)) {
                std::istringstream input_string(line);
                uint node_ID;
                double manning_n;
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
        } else {
            mesh.CallForEachElement([&problem_specific_input](auto& elt) {
                for (uint node = 0; node < elt.data.get_nnode(); ++node) {
                    elt.data.source.manning_n[node] = problem_specific_input.bottom_friction.manning_n;
                }

                elt.data.source.manning = true;
                elt.data.source.g_manning_n_sq =
                    Global::g *
                    std::pow(std::accumulate(elt.data.source.manning_n.begin(), elt.data.source.manning_n.end(), 0.0) /
                                 elt.data.source.manning_n.size(),
                             2);
            });
        }
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
            double latitude_avg        = y_avg / problem_specific_input.spherical_projection.R;
            elt.data.source.coriolis_f = 2 * 7.2921E-5 * std::sin(latitude_avg);
        });
    }

    // WETTING-DRYING INITIALIZE
    if (SWE::PostProcessing::wetting_drying) {
        mesh.CallForEachElement([](auto& elt) {
            auto& shape    = elt.GetShape();
            auto& state    = elt.data.state[0];
            auto& wd_state = elt.data.wet_dry_state;

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                wd_state.bath_at_vrtx[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
            }
            wd_state.bath_min  = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());
            wd_state.q_lin     = elt.ProjectBasisToLinear(state.q);
            wd_state.q_at_vrtx = elt.ComputeLinearUvrtx(wd_state.q_lin);
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                wd_state.h_at_vrtx[vrtx] = wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) + wd_state.bath_at_vrtx[vrtx];
            }

            bool set_wet_element = true;
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                if (wd_state.h_at_vrtx[vrtx] <= PostProcessing::h_o + PostProcessing::h_o_threshold) {
                    wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) = PostProcessing::h_o - wd_state.bath_at_vrtx[vrtx];
                    set_wet_element                              = false;
                }
            }

            if (set_wet_element) {
                wd_state.wet = true;
            } else {
                wd_state.wet = false;
                for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                    wd_state.q_at_vrtx(SWE::Variables::qx, vrtx) = 0.0;
                    wd_state.q_at_vrtx(SWE::Variables::qy, vrtx) = 0.0;
                    wd_state.q_at_vrtx(SWE::Variables::hc, vrtx) = 0.0;
                }
                state.q = elt.ProjectLinearToBasis(wd_state.q_at_vrtx);
                set_constant(state.rhs, 0.0);
            }
        });
    }

    // SLOPE LIMIT INITIALIZE
    if (SWE::PostProcessing::slope_limiting || SWE::SedimentTransport::bed_slope_limiting) {
        mesh.CallForEachElement([](auto& elt) {
            auto& shape    = elt.GetShape();
            auto& sl_state = elt.data.slope_limit_state;

            DynRowVector<double> bath_lin(elt.data.get_nvrtx());
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                bath_lin[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
            }
            sl_state.bath_at_baryctr = elt.ComputeLinearUbaryctr(bath_lin);
            sl_state.bath_at_vrtx    = bath_lin;
            sl_state.bath_at_midpts  = elt.ComputeLinearUmidpts(bath_lin);
            sl_state.baryctr_coord   = shape.GetBarycentricCoordinates();
            sl_state.midpts_coord    = shape.GetMidpointCoordinates();
            sl_state.lengths         = shape.GetLengths();
            sl_state.radius          = shape.GetRadius();

            StatVector<double, SWE::n_dimensions> median;
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                median                 = sl_state.midpts_coord[bound] - sl_state.baryctr_coord;
                sl_state.median[bound] = median / norm(median);
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

            sl_state.baryctr_coord_neigh[bound.bound_id] = sl_state.midpts_coord[bound.bound_id];
        });

        mesh.CallForEachElement([](auto& elt) {
            auto& sl_state = elt.data.slope_limit_state;

            StatMatrix<double, 2, 2> A;
            StatVector<double, 2> b;
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                for (uint perm = 0; perm < 2; ++perm) {
                    const uint element_1 = bound;
                    const uint element_2 = (bound + perm + 1) % elt.data.get_nbound();
                    if (!is_distributed(elt.GetBoundaryType()[element_1]) &&
                        !is_distributed(elt.GetBoundaryType()[element_2])) {
                        column(A, 0) = sl_state.baryctr_coord_neigh[element_1] - sl_state.baryctr_coord;
                        column(A, 1) = sl_state.baryctr_coord_neigh[element_2] - sl_state.baryctr_coord;
                        b            = sl_state.midpts_coord[bound] - sl_state.baryctr_coord;
                        if (perm == 0) {
                            sl_state.A_inv[bound]       = inverse(transpose(A));
                            sl_state.inv_theta_r[bound] = 1 / (norm(transpose(A)) * norm(sl_state.A_inv[bound]));
                        }
                        const double sqnorm_b = sq_norm(b);
                        solve_sle(A, b);
                        for (uint i = 0; i < 2; ++i)
                            if (Utilities::almost_equal(b[i], 0.0))
                                b[i] = 0.0;
                        if (b[0] >= 0.0 && b[1] >= 0.0) {
                            sl_state.alpha[bound]          = b;
                            sl_state.r_sq[bound]           = sqnorm_b;
                            sl_state.a_elem[2 * bound]     = element_1;
                            sl_state.a_elem[2 * bound + 1] = element_2;
                            break;
                        }
                    }
                }
            }
            const double inv_theta_r_sum =
                std::accumulate(sl_state.inv_theta_r.begin(), sl_state.inv_theta_r.end(), 0.0);
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                sl_state.d_r[bound] = sl_state.inv_theta_r[bound] / inv_theta_r_sum;
            }
        });
    }
}

template <typename MeshType, typename ProblemSpecificInputType>
void initialize_data_parallel_pre_send(MeshType& mesh,
                                       const ProblemSpecificInputType& problem_specific_input,
                                       uint comm_type) {
    initialize_data_serial(mesh, problem_specific_input);

    if (SWE::PostProcessing::slope_limiting || SWE::SedimentTransport::bed_slope_limiting) {
        mesh.CallForEachDistributedBoundary([comm_type](auto& dbound) {
            auto& sl_state = dbound.data.slope_limit_state;

            std::vector<double> message(SWE::n_dimensions);
            for (uint dim = 0; dim < SWE::n_dimensions; ++dim) {
                message[dim] = sl_state.baryctr_coord[dim];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(comm_type, message);
        });
    }
}

template <typename MeshType>
void initialize_data_parallel_post_receive(MeshType& mesh, uint comm_type) {
    if (SWE::PostProcessing::slope_limiting || SWE::SedimentTransport::bed_slope_limiting) {
        mesh.CallForEachDistributedBoundary([comm_type](auto& dbound) {
            auto& sl_state = dbound.data.slope_limit_state;

            std::vector<double> message(SWE::n_dimensions);
            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(comm_type, message);
            for (uint dim = 0; dim < SWE::n_dimensions; ++dim) {
                sl_state.baryctr_coord_neigh[dbound.bound_id][dim] = message[dim];
            }
        });

        mesh.CallForEachElement([](auto& elt) {
            auto& sl_state = elt.data.slope_limit_state;

            StatMatrix<double, 2, 2> A;
            StatVector<double, 2> b;
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                for (uint perm = 0; perm < 2; ++perm) {
                    const uint element_1 = bound;
                    const uint element_2 = (bound + perm + 1) % elt.data.get_nbound();
                    if (is_distributed(elt.GetBoundaryType()[element_1]) ||
                        is_distributed(elt.GetBoundaryType()[element_2])) {
                        column(A, 0) = sl_state.baryctr_coord_neigh[element_1] - sl_state.baryctr_coord;
                        column(A, 1) = sl_state.baryctr_coord_neigh[element_2] - sl_state.baryctr_coord;
                        b            = sl_state.midpts_coord[bound] - sl_state.baryctr_coord;
                        if (perm == 0) {
                            sl_state.A_inv[bound]       = inverse(transpose(A));
                            sl_state.inv_theta_r[bound] = 1 / (norm(transpose(A)) * norm(sl_state.A_inv[bound]));
                        }
                        const double sqnorm_b = sq_norm(b);
                        solve_sle(A, b);
                        for (uint i = 0; i < 2; ++i)
                            if (Utilities::almost_equal(b[i], 0.0))
                                b[i] = 0.0;
                        if (b[0] >= 0.0 && b[1] >= 0.0) {
                            sl_state.alpha[bound]          = b;
                            sl_state.r_sq[bound]           = sqnorm_b;
                            sl_state.a_elem[2 * bound]     = element_1;
                            sl_state.a_elem[2 * bound + 1] = element_2;
                            break;
                        }
                    }
                }
            }
            const double inv_theta_r_sum =
                std::accumulate(sl_state.inv_theta_r.begin(), sl_state.inv_theta_r.end(), 0.0);
            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                sl_state.d_r[bound] = sl_state.inv_theta_r[bound] / inv_theta_r_sum;
            }
        });
    }
}
}

#endif
