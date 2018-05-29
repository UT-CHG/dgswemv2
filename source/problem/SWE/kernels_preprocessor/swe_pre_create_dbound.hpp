#ifndef SWE_PRE_CREATE_DBOUND_HPP
#define SWE_PRE_CREATE_DBOUND_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_distributed_boundaries_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType&,
    InputParameters<ProblemInputType>& input,
    std::tuple<>&,
    Writer<SWE::Problem>& writer) {}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    InputParameters<ProblemInputType>& input,
    Communicator& communicator,
    Writer<SWE::Problem>& writer) {
    // *** //
    using DistributedBoundaryType =
        std::tuple_element<0, Geometry::DistributedBoundaryTypeTuple<SWE::Data, SWE::DBC::Distributed>>::type;

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

        uint begin_index_preproc  = 0;
        uint begin_index          = 0;
        uint begin_index_postproc = 0;

        // check if the data in rank_boundary_data matches communicator rank boundary
        RankBoundaryMetaData& rb_meta_data = input.mesh_input.dbmd_data.rank_boundary_data[rank_boundary_id];

        for (uint dboundary_id = 0; dboundary_id < rb_meta_data.elements_in.size(); dboundary_id++) {
            element_id_in = rb_meta_data.elements_in.at(dboundary_id);
            element_id_ex = rb_meta_data.elements_ex.at(dboundary_id);
            bound_id_in   = rb_meta_data.bound_ids_in.at(dboundary_id);
            bound_id_ex   = rb_meta_data.bound_ids_ex.at(dboundary_id);
            p             = rb_meta_data.p.at(dboundary_id);
            ngp           = boundary_integration.GetNumGP(2 * p);

            SWE::DBC::DBIndex index;

            index.x_at_baryctr = begin_index_preproc;
            index.y_at_baryctr = begin_index_preproc + 1;

            begin_index_preproc += 2;

            index.wet_dry = begin_index;

            index.ze_in = begin_index + 1;
            index.qx_in = begin_index + ngp + 1;
            index.qy_in = begin_index + 2 * ngp + 1;

            index.ze_ex = begin_index + ngp;
            index.qx_ex = begin_index + 2 * ngp;
            index.qy_ex = begin_index + 3 * ngp;

            begin_index += 3 * ngp + 1;

            index.ze_at_baryctr   = begin_index_postproc;
            index.qx_at_baryctr   = begin_index_postproc + 1;
            index.qy_at_baryctr   = begin_index_postproc + 2;
            index.bath_at_baryctr = begin_index_postproc + 3;

            begin_index_postproc += 4;

            auto& pre_dboundary = raw_boundaries.at(SWE::BoundaryTypes::distributed)
                                      .at(std::pair<uint, uint>{element_id_in, bound_id_in});
            pre_dboundary.p = p;

            mesh.template CreateDistributedBoundary<DistributedBoundaryType>(
                pre_dboundary,
                SWE::DBC::Distributed(SWE::DBC::DBDataExchanger(index,
                                                                send_preproc_buffer_reference,
                                                                receive_preproc_buffer_reference,
                                                                send_buffer_reference,
                                                                receive_buffer_reference,
                                                                send_postproc_buffer_reference,
                                                                receive_postproc_buffer_reference)));

            raw_boundaries.at(SWE::BoundaryTypes::distributed).erase(std::pair<uint, uint>{element_id_in, bound_id_in});
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