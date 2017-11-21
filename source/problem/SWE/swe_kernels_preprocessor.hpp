#ifndef SWE_KERNELS_PREPROCESSOR_HPP
#define SWE_KERNELS_PREPROCESSOR_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_boundaries_kernel(ProblemMeshType& mesh,
                                       std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries,
                                       Writer<SWE::Problem>& writer) {
    uint n_bound_old_land = 0;
    uint n_bound_old_tidal = 0;
    uint n_bound_old_flow = 0;

    using BoundaryTypes = Geometry::BoundaryTypeTuple<SWE::Data, SWE::Land, SWE::Tidal, SWE::Flow>;

    for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
        switch (it->first) {
            case SWE::land:
                using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;

                n_bound_old_land = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeLand>(*itt);
                }

                if (writer.WritingLog()) {
                    writer.GetLogFile() << "Number of land boundaries: " << mesh.GetNumberBoundaries() -
                                                                                n_bound_old_land << std::endl;
                }

                break;
            case SWE::tidal:
                using BoundaryTypeTidal = typename std::tuple_element<1, BoundaryTypes>::type;

                n_bound_old_tidal = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeTidal>(*itt);
                }

                if (writer.WritingLog()) {
                    writer.GetLogFile() << "Number of tidal boundaries: " << mesh.GetNumberBoundaries() -
                                                                                 n_bound_old_tidal << std::endl;
                }

                break;
            case SWE::flow:
                using BoundaryTypeFlow = typename std::tuple_element<2, BoundaryTypes>::type;

                n_bound_old_flow = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeFlow>(*itt);
                }

                if (writer.WritingLog()) {
                    writer.GetLogFile() << "Number of flow boundaries: " << mesh.GetNumberBoundaries() -
                                                                                n_bound_old_flow << std::endl;
                }

                break;
        }
    }
}

template <typename RawBoundaryType>
void Problem::create_distributed_boundaries_kernel(ProblemMeshType&,
                                                   std::tuple<>&,
                                                   std::map<uint, std::map<uint, RawBoundaryType>>&,
                                                   Writer<SWE::Problem>& writer) {}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries_kernel(
    ProblemMeshType& mesh,
    Communicator& communicator,
    std::map<uint, std::map<uint, RawBoundaryType>>& pre_distributed_boundaries,
    Writer<SWE::Problem>& writer) {
    using DistributedBoundaryType =
        std::tuple_element<0, Geometry::DistributedBoundaryTypeTuple<SWE::Data, SWE::Distributed>>::type;

    typename DistributedBoundaryType::BoundaryIntegrationType boundary_integration;

    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); rank_boundary_id++) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        std::vector<double>& send_buffer_reference = rank_boundary.send_buffer;
        std::vector<double>& receive_buffer_reference = rank_boundary.receive_buffer;

        uint element_id, bound_id, p, ngp, ze_in_index, qx_in_index, qy_in_index, ze_ex_index, qx_ex_index, qy_ex_index;

        uint begin_index = 0;
        for (uint dboundary_id = 0; dboundary_id < rank_boundary.elements.size(); dboundary_id++) {
            element_id = rank_boundary.elements.at(dboundary_id);
            bound_id = rank_boundary.bound_ids.at(dboundary_id);
            p = rank_boundary.p.at(dboundary_id);
            ngp = boundary_integration.GetNumGP(2 * p);

            ze_in_index = begin_index;
            qx_in_index = begin_index + ngp;
            qy_in_index = begin_index + 2 * ngp;

            ze_ex_index = begin_index + ngp - 1;
            qx_ex_index = begin_index + 2 * ngp - 1;
            qy_ex_index = begin_index + 3 * ngp - 1;

            begin_index += 3 * ngp;

            auto& pre_dboundary = pre_distributed_boundaries.at(element_id).at(bound_id);
            pre_dboundary.p = p;

            mesh.template CreateDistributedBoundary<DistributedBoundaryType>(pre_dboundary,
                                                                             SWE::Distributed(send_buffer_reference,
                                                                                              receive_buffer_reference,
                                                                                              ze_in_index,
                                                                                              qx_in_index,
                                                                                              qy_in_index,
                                                                                              ze_ex_index,
                                                                                              qx_ex_index,
                                                                                              qy_ex_index));
        }
        send_buffer_reference.resize(begin_index);
        receive_buffer_reference.resize(begin_index);
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed boundaries: " << mesh.GetNumberDistributedBoundaries()
                            << std::endl;
    }
}

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

    mesh.CallForEachElement([&bathymetry](auto& elt) {
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

        auto ze_init = [](Point<2>& pt) { return SWE::ic_ze(0, pt); };

        state.ze = elt.L2Projection(ze_init);

        auto qx_init = [](Point<2>& pt) { return SWE::ic_qx(0, pt); };

        state.qx = elt.L2Projection(qx_init);

        auto qy_init = [](Point<2>& pt) { return SWE::ic_qy(0, pt); };

        state.qy = elt.L2Projection(qy_init);
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& state_in = intface.data_in.state[0];
        auto& state_ex = intface.data_ex.state[0];

        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        intface.ComputeUgpIN(state_in.bath, boundary_in.bath_at_gp);
        intface.ComputeUgpEX(state_ex.bath, boundary_ex.bath_at_gp);
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& state = bound.data.state[0];
        auto& boundary = bound.data.boundary[bound.bound_id];

        bound.ComputeUgp(state.bath, boundary.bath_at_gp);
    });

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        dbound.ComputeUgp(dbound.data.state[0].bath, dbound.data.boundary[dbound.bound_id].bath_at_gp);
    });
}

void Problem::initialize_wd_data_kernel(ProblemMeshType& mesh) {
    mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];
        auto& internal = elt.data.internal;
        auto& wd_state = elt.data.wet_dry_state;

        elt.ComputeUvrtx(state.bath, wd_state.bath_at_vrtx);
        wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());

        // for testing initialize wet/dry element
        if (state.bath[0] < 0) {
            wd_state.wet = false;

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                wd_state.ze_at_vrtx[vrtx] = Global::h_o - wd_state.bath_at_vrtx[vrtx];
            }

            state.ze = elt.L2Projection(wd_state.ze_at_vrtx);
            std::fill(state.qx.begin(), state.qx.end(), 0.0);
            std::fill(state.qy.begin(), state.qy.end(), 0.0);

            std::fill(state.rhs_ze.begin(), state.rhs_ze.end(), 0.0);
            std::fill(state.rhs_qx.begin(), state.rhs_qx.end(), 0.0);
            std::fill(state.rhs_qy.begin(), state.rhs_qy.end(), 0.0);
        } else {
            wd_state.wet = true;

            elt.ComputeUgp(state.ze, internal.ze_at_gp);
            elt.ComputeUgp(state.qx, internal.qx_at_gp);
            elt.ComputeUgp(state.qy, internal.qy_at_gp);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                internal.h_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];
            }

            wd_state.water_volume = elt.Integration(internal.h_at_gp);
        }
    });
}

void Problem::initialize_sl_data_kernel(ProblemMeshType& mesh) {
    mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];
        auto& sl_state = elt.data.slope_limit_state;

        elt.ComputeUbaryctr(state.bath, sl_state.bath_at_baryctr);
        elt.ComputeUmidpts(state.bath, sl_state.bath_at_midpts);
        
        sl_state.baryctr_coord = elt.GetShape().GetBarycentricCoordinates();
        sl_state.midpts_coord = elt.GetShape().GetMidpointCoordinates();
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.bath_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.bath_at_baryctr;
        sl_state_ex.bath_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.bath_at_baryctr;

        sl_state_in.baryctr_coord_neigh[intface.bound_id_in] = sl_state_ex.baryctr_coord;
        sl_state_ex.baryctr_coord_neigh[intface.bound_id_ex] = sl_state_in.baryctr_coord;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.bath_at_baryctr_neigh[bound.bound_id] = sl_state.bath_at_baryctr;

        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::x] = 
            2.0*sl_state.midpts_coord[bound.bound_id][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::y] = 
            2.0*sl_state.midpts_coord[bound.bound_id][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2))
        std::vector<double> b = std::vector(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % 3;
            
            A[0][0] = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
            A[1][0] = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
            A[0][1] = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
            A[1][1] = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

            b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
            b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

            double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];

            sl_state.alpha_1[bound] = (A[1][1]*b[0] - A[0][1]*b[1])/det;
            sl_state.alpha_2[bound] = (-A[1][0]*b[0] + A[0][0]*b[1])/det;

            sl_state.r_sq[bound] = std::pow(b[0], 2.0) + std::pow(b[1], 2.0);
        }
    });
}
}

#endif
