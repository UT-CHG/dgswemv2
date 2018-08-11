#ifndef RKDG_SWE_PRE_CREATE_DBOUND_HPP
#define RKDG_SWE_PRE_CREATE_DBOUND_HPP

namespace SWE {
namespace RKDG {
template <typename RawBoundaryType>
void Problem::create_distributed_boundaries_kernel(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType&,
    ProblemInputType& input,
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
    using DistributedBoundaryTypes =
        Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed, DBC::DistributedLevee>;

    auto& raw_bound_distributed = raw_boundaries[distributed(SWE::BoundaryTypes::internal)];
    auto& raw_bound_distr_levee = raw_boundaries[distributed(SWE::BoundaryTypes::levee)];

    uint n_distributed = 0;
    uint n_distr_levee = 0;

    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); rank_boundary_id++) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        rank_boundary.send_buffer.resize(SWE::n_communications);
        rank_boundary.receive_buffer.resize(SWE::n_communications);

        std::vector<double>& send_preproc_buffer    = rank_boundary.send_buffer[SWE::CommTypes::preprocessor];
        std::vector<double>& receive_preproc_buffer = rank_boundary.receive_buffer[SWE::CommTypes::preprocessor];

        std::vector<double>& send_buffer    = rank_boundary.send_buffer[SWE::CommTypes::processor];
        std::vector<double>& receive_buffer = rank_boundary.receive_buffer[SWE::CommTypes::processor];

        std::vector<double>& send_postproc_buffer    = rank_boundary.send_buffer[SWE::CommTypes::postprocessor];
        std::vector<double>& receive_postproc_buffer = rank_boundary.receive_buffer[SWE::CommTypes::postprocessor];

        uint element_id_in, bound_id_in, p, ngp;
        // uint element_id_ex, bound_id_ex;

        uint begin_index_preproc  = 0;
        uint begin_index          = 0;
        uint begin_index_postproc = 0;

        // check if the data in rank_boundary_data matches communicator rank boundary
        const RankBoundaryMetaData& rb_meta_data = rank_boundary.db_data;

        for (uint dboundary_id = 0; dboundary_id < rb_meta_data.elements_in.size(); dboundary_id++) {
            element_id_in = rb_meta_data.elements_in[dboundary_id];
            // element_id_ex = rb_meta_data.elements_ex[dboundary_id];
            bound_id_in = rb_meta_data.bound_ids_in[dboundary_id];
            // bound_id_ex   = rb_meta_data.bound_ids_ex[dboundary_id];
            p = rb_meta_data.p[dboundary_id];

            std::pair<uint, uint> dbound_key = std::pair<uint, uint>{element_id_in, bound_id_in};

            // this finds number of gps used in integrations at the current bound
            // revise for a safer more efficient way to do so, be careful about keeping 2 * p + 1 consistent
            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;

                typename DBTypeDistributed::BoundaryIntegrationType boundary_integration;

                ngp = boundary_integration.GetNumGP(2 * p + 1);
            } else if (raw_bound_distr_levee.find(dbound_key) != raw_bound_distr_levee.end()) {
                using DBTypeDistributedLevee = typename std::tuple_element<1, DistributedBoundaryTypes>::type;

                typename DBTypeDistributedLevee::BoundaryIntegrationType boundary_integration;

                ngp = boundary_integration.GetNumGP(2 * p + 1);
            } else {
                throw std::logic_error("Fatal Error: unable to find raw distributed boundary!\n");
            }

            std::vector<uint> offset(SWE::n_communications);

            offset[SWE::CommTypes::preprocessor]  = begin_index_preproc;
            offset[SWE::CommTypes::processor]     = begin_index;
            offset[SWE::CommTypes::postprocessor] = begin_index_postproc;

            DBC::DBIndex index;

            index.x_at_baryctr = begin_index_preproc;
            index.y_at_baryctr = begin_index_preproc + 1;

            begin_index_preproc += 2;

            index.wet_dry = begin_index;
            index.q_in    = begin_index + 1;
            index.q_ex    = begin_index + SWE::n_variables * ngp;

            begin_index += SWE::n_variables * ngp + 1;

            index.wet_dry_postproc = begin_index_postproc;
            index.q_at_baryctr     = begin_index_postproc + 1;

            begin_index_postproc += SWE::n_variables + 1;

            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;

                auto& raw_boundary = raw_bound_distributed.find(dbound_key)->second;

                raw_boundary.p = p;

                mesh.template CreateDistributedBoundary<DBTypeDistributed>(
                    std::move(raw_boundary),
                    DBC::Distributed(DBC::DBDataExchanger(offset,
                                                          rank_boundary.send_buffer,
                                                          rank_boundary.receive_buffer,
                                                          index,
                                                          send_preproc_buffer,
                                                          receive_preproc_buffer,
                                                          send_buffer,
                                                          receive_buffer,
                                                          send_postproc_buffer,
                                                          receive_postproc_buffer)));

                n_distributed++;

                raw_bound_distributed.erase(dbound_key);
            } else if (raw_bound_distr_levee.find(dbound_key) != raw_bound_distr_levee.end()) {
                using DBTypeDistributedLevee = typename std::tuple_element<1, DistributedBoundaryTypes>::type;

                auto& raw_boundary = raw_bound_distr_levee.find(dbound_key)->second;

                raw_boundary.p = p;

                auto& levee_data = problem_input.levee_is_data;

                std::vector<LeveeInput> levee;

                for (uint node = 0; node < raw_boundary.node_ID.size(); node++) {
                    bool found = false;

                    for (auto& levee_node : levee_data) {
                        if (levee_node.first.first == raw_boundary.node_ID[node] ||
                            levee_node.first.second == raw_boundary.node_ID[node]) {
                            levee.push_back(levee_node.second);

                            found = true;
                        }
                    }

                    if (!found) {
                        throw std::logic_error("Fatal Error: unable to find distributed levee data!\n");
                    }
                }

                mesh.template CreateDistributedBoundary<DBTypeDistributedLevee>(
                    std::move(raw_boundary),
                    DBC::DistributedLevee(DBC::DBDataExchanger(offset,
                                                               rank_boundary.send_buffer,
                                                               rank_boundary.receive_buffer,
                                                               index,
                                                               send_preproc_buffer,
                                                               receive_preproc_buffer,
                                                               send_buffer,
                                                               receive_buffer,
                                                               send_postproc_buffer,
                                                               receive_postproc_buffer),
                                          levee));

                n_distr_levee++;

                raw_bound_distr_levee.erase(dbound_key);
            }
        }

        send_preproc_buffer.resize(begin_index_preproc);
        receive_preproc_buffer.resize(begin_index_preproc);

        send_buffer.resize(begin_index);
        receive_buffer.resize(begin_index);

        send_postproc_buffer.resize(begin_index_postproc);
        receive_postproc_buffer.resize(begin_index_postproc);
    }

    mesh.CallForEachDistributedBoundary([](auto& dbound) { dbound.boundary_condition.Initialize(dbound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed boundaries: " << n_distributed << std::endl
                            << "Number of distributed levee boundaries: " << n_distr_levee << std::endl;
    }
}
}
}

#endif