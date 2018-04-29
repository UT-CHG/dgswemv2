#ifndef SWE_PRE_CREATE_DBOUND_HPP
#define SWE_PRE_CREATE_DBOUND_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_distributed_boundaries_kernel(
    ProblemMeshType&,
    std::tuple<>&,
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    Writer<SWE::Problem>& writer) {}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries_kernel(
    ProblemMeshType& mesh,
    Communicator& communicator,
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    Writer<SWE::Problem>& writer) {
    // *** //
    using DistributedBoundaryType =
        std::tuple_element<0, Geometry::DistributedBoundaryTypeTuple<SWE::Data, SWE::BC::Distributed>>::type;

    typename DistributedBoundaryType::BoundaryIntegrationType boundary_integration;

    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); rank_boundary_id++) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        std::vector<double>& send_preproc_buffer_reference    = rank_boundary.send_preproc_buffer;
        std::vector<double>& receive_preproc_buffer_reference = rank_boundary.receive_preproc_buffer;

        std::vector<double>& send_buffer_reference    = rank_boundary.send_buffer;
        std::vector<double>& receive_buffer_reference = rank_boundary.receive_buffer;

        std::vector<double>& send_postproc_buffer_reference    = rank_boundary.send_postproc_buffer;
        std::vector<double>& receive_postproc_buffer_reference = rank_boundary.receive_postproc_buffer;

        uint element_id_in, element_id_ex, bound_id_in, bound_id_ex, p, ngp;
        uint wet_dry_index, ze_in_index, qx_in_index, qy_in_index, ze_ex_index, qx_ex_index, qy_ex_index;
        uint ze_at_baryctr_index, qx_at_baryctr_index, qy_at_baryctr_index, bath_at_baryctr_index;
        uint x_at_baryctr_index, y_at_baryctr_index;

        uint begin_index_preproc  = 0;
        uint begin_index          = 0;
        uint begin_index_postproc = 0;

        for (uint dboundary_id = 0; dboundary_id < rank_boundary.elements_in.size(); dboundary_id++) {
            element_id_in = rank_boundary.elements_in.at(dboundary_id);
            element_id_ex = rank_boundary.elements_ex.at(dboundary_id);
            bound_id_in   = rank_boundary.bound_ids_in.at(dboundary_id);
            bound_id_ex   = rank_boundary.bound_ids_ex.at(dboundary_id);
            p             = rank_boundary.p.at(dboundary_id);
            ngp           = boundary_integration.GetNumGP(2 * p);

            x_at_baryctr_index = begin_index_preproc;
            y_at_baryctr_index = begin_index_preproc + 1;

            begin_index_preproc += 2;

            wet_dry_index = begin_index;

            ze_in_index = begin_index + 1;
            qx_in_index = begin_index + ngp + 1;
            qy_in_index = begin_index + 2 * ngp + 1;

            ze_ex_index = begin_index + ngp;
            qx_ex_index = begin_index + 2 * ngp;
            qy_ex_index = begin_index + 3 * ngp;

            begin_index += 3 * ngp + 1;

            ze_at_baryctr_index   = begin_index_postproc;
            qx_at_baryctr_index   = begin_index_postproc + 1;
            qy_at_baryctr_index   = begin_index_postproc + 2;
            bath_at_baryctr_index = begin_index_postproc + 3;

            begin_index_postproc += 4;

            auto& pre_dboundary = raw_boundaries.at(SWE::BoundaryConditions::distributed)
                                      .at(std::pair<uint, uint>{element_id_in, bound_id_in});
            pre_dboundary.p = p;

            mesh.template CreateDistributedBoundary<DistributedBoundaryType>(
                pre_dboundary,
                SWE::BC::Distributed(send_preproc_buffer_reference,
                                     receive_preproc_buffer_reference,
                                     send_buffer_reference,
                                     receive_buffer_reference,
                                     send_postproc_buffer_reference,
                                     receive_postproc_buffer_reference,
                                     x_at_baryctr_index,
                                     y_at_baryctr_index,
                                     wet_dry_index,
                                     ze_in_index,
                                     qx_in_index,
                                     qy_in_index,
                                     ze_ex_index,
                                     qx_ex_index,
                                     qy_ex_index,
                                     ze_at_baryctr_index,
                                     qx_at_baryctr_index,
                                     qy_at_baryctr_index,
                                     bath_at_baryctr_index));

            raw_boundaries.at(SWE::BoundaryConditions::distributed)
                .erase(std::pair<uint, uint>{element_id_in, bound_id_in});
        }

        send_preproc_buffer_reference.resize(begin_index_preproc);
        receive_preproc_buffer_reference.resize(begin_index_preproc);

        send_buffer_reference.resize(begin_index);
        receive_buffer_reference.resize(begin_index);

        send_postproc_buffer_reference.resize(begin_index_postproc);
        receive_postproc_buffer_reference.resize(begin_index_postproc);
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed boundaries: " << mesh.GetNumberDistributedBoundaries()
                            << std::endl;
    }
}
}

#endif