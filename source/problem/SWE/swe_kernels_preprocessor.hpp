#ifndef SWE_KERNELS_PREPROCESSOR_HPP
#define SWE_KERNELS_PREPROCESSOR_HPP

#include "swe_true_src_functions.hpp"

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_boundaries_kernel(ProblemMeshType& mesh,
                                       std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries) {
    uint n_bound_old_land = 0;
    uint n_bound_old_tidal = 0;

    std::ofstream log_file("output/" + mesh.GetMeshName() + "_log", std::ofstream::app);

    using BoundaryTypes = Geometry::BoundaryTypeTuple<SWE::Data, SWE::Land, SWE::Tidal>;

    for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
        switch (it->first) {
            case SWE::land:
                using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;

                n_bound_old_land = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeLand>(*itt);
                }

                log_file << "Number of land boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_land << std::endl;

                break;
            case SWE::tidal:
                using BoundaryTypeTidal = typename std::tuple_element<1, BoundaryTypes>::type;

                n_bound_old_tidal = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeTidal>(*itt);
                }

                log_file << "Number of tidal boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_tidal
                         << std::endl;

                break;
        }
    }
}

template <typename RawBoundaryType>
void Problem::create_distributed_boundaries_kernel(ProblemMeshType&,
                                                   std::tuple<>&,
                                                   std::map<uint, std::map<uint, RawBoundaryType>>&) {}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries_kernel(
    ProblemMeshType& mesh,
    Communicator& communicator,
    std::map<uint, std::map<uint, RawBoundaryType>>& pre_distributed_boundaries) {
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

    std::ofstream log_file("output/" + mesh.GetMeshName() + "_log", std::ofstream::app);

    log_file << "Number of distributed boundaries: " << mesh.GetNumberDistributedBoundaries() << std::endl;
}

void Problem::initialize_data_kernel(ProblemMeshType& mesh, const MeshMetaData& mesh_data) {
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

        if (!bathymetry.count(id)) {
            throw std::logic_error("Error: could not find bathymetry for element with id: " + id);
        }

        elt.data.state[0].bath = elt.L2Projection(bathymetry[id]);

        elt.ComputeUgp(elt.data.state[0].bath, elt.data.internal.bath_at_gp);

        elt.ComputeDUgp(GlobalCoord::x, elt.data.state[0].bath, elt.data.internal.bath_deriv_wrt_x_at_gp);
        elt.ComputeDUgp(GlobalCoord::y, elt.data.state[0].bath, elt.data.internal.bath_deriv_wrt_y_at_gp);

        auto ze_init = [](Point<2>& pt) { return SWE::true_ze(0, pt); };

        elt.data.state[0].ze = elt.L2Projection(ze_init);

        auto qx_init = [](Point<2>& pt) { return SWE::true_qx(0, pt); };

        elt.data.state[0].qx = elt.L2Projection(qx_init);

        auto qy_init = [](Point<2>& pt) { return SWE::true_qy(0, pt); };

        elt.data.state[0].qy = elt.L2Projection(qy_init);
    });

    mesh.CallForEachInterface([](auto& intface) {
        intface.ComputeUgpIN(intface.data_in.state[0].bath, intface.data_in.boundary[intface.bound_id_in].bath_at_gp);
        intface.ComputeUgpEX(intface.data_ex.state[0].bath, intface.data_ex.boundary[intface.bound_id_ex].bath_at_gp);
    });

    mesh.CallForEachBoundary([](auto& bound) {
        bound.ComputeUgp(bound.data.state[0].bath, bound.data.boundary[bound.bound_id].bath_at_gp);
    });

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        dbound.ComputeUgp(dbound.data.state[0].bath, dbound.data.boundary[dbound.bound_id].bath_at_gp);
    });
}
}

#endif
