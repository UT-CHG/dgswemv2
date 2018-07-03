#ifndef EHDG_SWE_PRE_CREATE_DBOUND_HPP
#define EHDG_SWE_PRE_CREATE_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename RawBoundaryType>
void Problem::create_distributed_boundaries_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType&,
    ProblemInputType& problem_input,
    std::tuple<>&,
    Writer<Problem>& writer) {}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    ProblemInputType& problem_input,
    Communicator& communicator,
    Writer<Problem>& writer) {
    // *** //
    using DistributedBoundaryTypes = Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>;

    auto& raw_bound_distributed = raw_boundaries[distributed(SWE::BoundaryTypes::internal)];

    uint n_distributed = 0;

    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); rank_boundary_id++) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        std::vector<double>& send_preproc_buffer_reference    = rank_boundary.send_preproc_buffer;
        std::vector<double>& receive_preproc_buffer_reference = rank_boundary.receive_preproc_buffer;

        std::vector<double>& send_buffer_reference    = rank_boundary.send_buffer;
        std::vector<double>& receive_buffer_reference = rank_boundary.receive_buffer;

        std::vector<double>& send_postproc_buffer_reference    = rank_boundary.send_postproc_buffer;
        std::vector<double>& receive_postproc_buffer_reference = rank_boundary.receive_postproc_buffer;

        uint element_id_in, bound_id_in, p, ngp;
        // uint element_id_ex, bound_id_ex;

        uint begin_index_preproc  = 0;
        uint begin_index          = 0;
        uint begin_index_postproc = 0;

        // check if the data in rank_boundary_data matches communicator rank boundary
        const RankBoundaryMetaData& rb_meta_data = rank_boundary.db_data;

        for (uint dboundary_id = 0; dboundary_id < rb_meta_data.elements_in.size(); dboundary_id++) {
            element_id_in = rb_meta_data.elements_in.at(dboundary_id);
            bound_id_in   = rb_meta_data.bound_ids_in.at(dboundary_id);
            p             = rb_meta_data.p.at(dboundary_id);

            // element_id_ex = rb_meta_data.elements_ex.at(dboundary_id);
            // bound_id_ex   = rb_meta_data.bound_ids_ex.at(dboundary_id);

            std::pair<uint, uint> dbound_key = std::pair<uint, uint>{element_id_in, bound_id_in};

            // this finds number of gps used in integrations at the current bound
            // revise for a safer more efficient way to do so, be careful about keeping 2 * p + 1 consistent
            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;

                typename DBTypeDistributed::BoundaryIntegrationType boundary_integration;

                ngp = boundary_integration.GetNumGP(2 * p + 1);
            } else {
                throw std::logic_error("Fatal Error: unable to find raw distributed boundary!\n");
            }

            DBC::DBIndex index;

            index.ze_in            = begin_index;
            index.qx_in            = begin_index + ngp;
            index.qy_in            = begin_index + 2 * ngp;
            index.ze_flux_dot_n_in = begin_index + 3 * ngp;
            index.qx_flux_dot_n_in = begin_index + 4 * ngp;
            index.qy_flux_dot_n_in = begin_index + 5 * ngp;

            index.ze_ex            = begin_index + ngp - 1;
            index.qx_ex            = begin_index + 2 * ngp - 1;
            index.qy_ex            = begin_index + 3 * ngp - 1;
            index.ze_flux_dot_n_ex = begin_index + 4 * ngp - 1;
            index.qx_flux_dot_n_ex = begin_index + 5 * ngp - 1;
            index.qy_flux_dot_n_ex = begin_index + 6 * ngp - 1;

            begin_index += 6 * ngp;

            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;

                auto& raw_boundary = raw_bound_distributed.find(dbound_key)->second;

                raw_boundary.p = p;

                mesh.template CreateDistributedBoundary<DBTypeDistributed>(
                    std::move(raw_boundary),
                    DBC::Distributed(DBC::DBDataExchanger(index,
                                                          send_preproc_buffer_reference,
                                                          receive_preproc_buffer_reference,
                                                          send_buffer_reference,
                                                          receive_buffer_reference,
                                                          send_postproc_buffer_reference,
                                                          receive_postproc_buffer_reference)));

                n_distributed++;

                raw_bound_distributed.erase(dbound_key);
            }
        }

        send_preproc_buffer_reference.resize(begin_index_preproc);
        receive_preproc_buffer_reference.resize(begin_index_preproc);

        send_buffer_reference.resize(begin_index);
        receive_buffer_reference.resize(begin_index);

        send_postproc_buffer_reference.resize(begin_index_postproc);
        receive_postproc_buffer_reference.resize(begin_index_postproc);
    }

    mesh.CallForEachDistributedBoundary([](auto& dbound) { dbound.boundary_condition.Initialize(dbound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed boundaries: " << n_distributed << std::endl;
    }
}
}
}

#endif