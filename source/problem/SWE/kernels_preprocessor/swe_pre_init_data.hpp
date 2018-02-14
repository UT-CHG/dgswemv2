#ifndef SWE_PRE_INIT_DATA_HPP
#define SWE_PRE_INIT_DATA_HPP

namespace SWE {
void Problem::initialize_data_kernel(ProblemMeshType& mesh,
                                     const MeshMetaData& mesh_data,
                                     const ProblemInputType& problem_specific_input) {
    mesh.CallForEachElement([](auto& elt) { elt.data.initialize(); });

    std::unordered_map<uint, std::vector<double>> bathymetry;

    std::vector<double> bathymetry_temp;
    for (const auto& elt : mesh_data.elements) {
        auto nodal_coordinates = mesh_data.GetNodalCoordinates(elt.first);

        for (auto& node_coordinate : nodal_coordinates) {
            bathymetry_temp.push_back(node_coordinate[GlobalCoord::z]);
        }

        bathymetry.insert({elt.first, bathymetry_temp});

        bathymetry_temp.clear();
    }

    mesh.CallForEachElement([&bathymetry, &problem_specific_input](auto& elt) {
        uint id = elt.GetID();

        auto& state = elt.data.state[0];
        auto& internal = elt.data.internal;

        if (!bathymetry.count(id)) {
            throw std::logic_error("Error: could not find bathymetry for element with id: " + id);
        }

        state.bath = elt.L2Projection(bathymetry[id]);

        elt.ComputeUgp(state.bath, internal.bath_at_gp);

        elt.ComputeDUgp(GlobalCoord::x, state.bath, internal.bath_deriv_wrt_x_at_gp);
        elt.ComputeDUgp(GlobalCoord::y, state.bath, internal.bath_deriv_wrt_y_at_gp);

        if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Constant) {
            uint n_node = elt.GetShape().nodal_coordinates.size();

            std::vector<double> ze_node(n_node, problem_specific_input.initial_conditions.ze_initial);
            state.ze = elt.L2Projection(ze_node);

            std::vector<double> qx_node(n_node, problem_specific_input.initial_conditions.qx_initial);
            state.qx = elt.L2Projection(qx_node);

            std::vector<double> qy_node(n_node, problem_specific_input.initial_conditions.qy_initial);
            state.qy = elt.L2Projection(qy_node);
        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            auto ze_init = [](Point<2>& pt) { return SWE::ic_ze(0, pt); };
            state.ze = elt.L2Projection(ze_init);

            auto qx_init = [](Point<2>& pt) { return SWE::ic_qx(0, pt); };
            state.qx = elt.L2Projection(qx_init);

            auto qy_init = [](Point<2>& pt) { return SWE::ic_qy(0, pt); };
            state.qy = elt.L2Projection(qy_init);
        }
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& state_in = intface.data_in->state[0];
        auto& state_ex = intface.data_ex->state[0];

        auto& boundary_in = intface.data_in->boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex->boundary[intface.bound_id_ex];

        intface.ComputeUgpIN(state_in.bath, boundary_in.bath_at_gp);
        intface.ComputeUgpEX(state_ex.bath, boundary_ex.bath_at_gp);
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& state = bound.data->state[0];
        auto& boundary = bound.data->boundary[bound.bound_id];

        bound.ComputeUgp(state.bath, boundary.bath_at_gp);
    });

    // SOURCE TERMS INITIALIZE
    if (problem_specific_input.bottom_friction.type == SWE::BottomFrictionType::Manning) {
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

            //# of node != # of vrtx in case we have an iso-p element with p>1
            // I assume we will have values only at vrtx in files
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                elt.data.source.manning_n[vrtx] = node_manning_n[node_ID[vrtx]];
            }

            elt.data.source.manning = true;
            elt.data.source.g_manning_n_sq =
                Global::g *
                std::pow(std::accumulate(elt.data.source.manning_n.begin(), elt.data.source.manning_n.end(), 0.0) /
                             elt.data.source.manning_n.size(),
                         2);
        });
    }

    /*if (problem_specific_input.coriolis.type != SWE::CoriolisType::None) {
        mesh.CallForEachElement([](auto& elt) {
            //# of node != # of vrtx in case we have an iso-p element with p>1
            // I assume we will have values only at vrtx in files
            std::vector<double>& nodal_coordinates = elt.GetShape().GetNodalCoordinates();
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                elt.GetShape().nodal_coordinates[vrtx][GlobalCoord::x];
                elt.GetShape().nodal_coordinates[vrtx][GlobalCoord::y];
            }
        });
        }*/

    // WETTING-DRYING INITIALIZE

    mesh.CallForEachElement([&bathymetry](auto& elt) {
        auto& state = elt.data.state[0];
        auto& internal = elt.data.internal;
        auto& wd_state = elt.data.wet_dry_state;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.bath_at_vrtx[vrtx] = bathymetry[elt.GetID()][vrtx];
        }

        wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());

        elt.ProjectBasisToLinear(state.ze, wd_state.ze_lin);

        elt.ComputeLinearUvrtx(wd_state.ze_lin, wd_state.ze_at_vrtx);

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.h_at_vrtx[vrtx] = wd_state.ze_at_vrtx[vrtx] + wd_state.bath_at_vrtx[vrtx];
        }

        bool set_wet_element = true;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            if (wd_state.h_at_vrtx[vrtx] <= Global::h_o) {
                wd_state.ze_at_vrtx[vrtx] = Global::h_o - wd_state.bath_at_vrtx[vrtx];

                set_wet_element = false;
            }
        }

        if (set_wet_element) {
            wd_state.wet = true;
        } else {
            wd_state.wet = false;

            elt.ProjectLinearToBasis(wd_state.ze_at_vrtx, state.ze);
            std::fill(state.qx.begin(), state.qx.end(), 0.0);
            std::fill(state.qy.begin(), state.qy.end(), 0.0);

            std::fill(state.rhs_ze.begin(), state.rhs_ze.end(), 0.0);
            std::fill(state.rhs_qx.begin(), state.rhs_qx.end(), 0.0);
            std::fill(state.rhs_qy.begin(), state.rhs_qy.end(), 0.0);
        }

        elt.ComputeUgp(state.ze, internal.ze_at_gp);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            internal.h_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];
        }

        wd_state.water_volume = elt.Integration(internal.h_at_gp);
    });

    // SLOPE LIMIT INITIALIZE

    mesh.CallForEachElement([&bathymetry](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        std::vector<double> bath_lin(elt.data.get_nvrtx());

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            bath_lin[vrtx] = bathymetry[elt.GetID()][vrtx];
        }

        elt.ComputeLinearUbaryctr(bath_lin, sl_state.bath_at_baryctr);
        elt.ComputeLinearUmidpts(bath_lin, sl_state.bath_at_midpts);

        sl_state.baryctr_coord = elt.GetShape().GetBarycentricCoordinates();
        sl_state.midpts_coord = elt.GetShape().GetMidpointCoordinates();

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            sl_state.surface_normal[bound] = elt.GetShape().GetSurfaceNormal(bound, std::vector<Point<2>>(0))[0];
        }
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in->slope_limit_state;
        auto& sl_state_ex = intface.data_ex->slope_limit_state;

        sl_state_in.bath_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.bath_at_baryctr;
        sl_state_ex.bath_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.bath_at_baryctr;

        sl_state_in.baryctr_coord_neigh[intface.bound_id_in] = sl_state_ex.baryctr_coord;
        sl_state_ex.baryctr_coord_neigh[intface.bound_id_ex] = sl_state_in.baryctr_coord;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data->slope_limit_state;

        sl_state.bath_at_baryctr_neigh[bound.bound_id] = sl_state.bath_at_baryctr;

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

            if ((elt.GetBoundaryType()[element_1] != SWE::BoundaryConditions::distributed) &&
                (elt.GetBoundaryType()[element_2] != SWE::BoundaryConditions::distributed)) {
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

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& state = dbound.data->state[0];
        auto& boundary = dbound.data->boundary[dbound.bound_id];

        dbound.ComputeUgp(state.bath, boundary.bath_at_gp);
    });

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& sl_state = dbound.data->slope_limit_state;

        dbound.boundary_condition.SetPreprocEX(sl_state.baryctr_coord[GlobalCoord::x],
                                               sl_state.baryctr_coord[GlobalCoord::y]);
    });
}

void Problem::initialize_data_parallel_post_receive_kernel(ProblemMeshType& mesh,
                                                           const ProblemInputType& problem_specific_input) {
    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& sl_state = dbound.data->slope_limit_state;

        dbound.boundary_condition.GetPreprocEX(sl_state.baryctr_coord_neigh[dbound.bound_id][GlobalCoord::x],
                                               sl_state.baryctr_coord_neigh[dbound.bound_id][GlobalCoord::y]);
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2));
        std::vector<double> b(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % elt.data.get_nbound();

            if ((elt.GetBoundaryType()[element_1] == SWE::BoundaryConditions::distributed) ||
                (elt.GetBoundaryType()[element_2] == SWE::BoundaryConditions::distributed)) {
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

#endif