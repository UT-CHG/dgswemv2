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
    ProblemInputType& input,
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
        RankBoundaryMetaData& rb_meta_data = rank_boundary.dbmd_data;

        for (uint dboundary_id = 0; dboundary_id < rb_meta_data.elements_in.size(); dboundary_id++) {
            element_id_in = rb_meta_data.elements_in.at(dboundary_id);
            element_id_ex = rb_meta_data.elements_ex.at(dboundary_id);
            bound_id_in   = rb_meta_data.bound_ids_in.at(dboundary_id);
            bound_id_ex   = rb_meta_data.bound_ids_ex.at(dboundary_id);
            p             = rb_meta_data.p.at(dboundary_id);

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

            DBC::DBIndex index;

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

            index.wet_dry_postproc = begin_index_postproc;

            index.ze_at_baryctr = begin_index_postproc + 1;
            index.qx_at_baryctr = begin_index_postproc + 2;
            index.qy_at_baryctr = begin_index_postproc + 3;

            begin_index_postproc += 4;

            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;

                auto& raw_boundary = raw_bound_distributed.find(dbound_key)->second;

                raw_boundary.p = p;

                mesh.template CreateDistributedBoundary<DBTypeDistributed>(
                    raw_boundary,
                    DBC::Distributed(DBC::DBDataExchanger(index,
                                                          send_preproc_buffer_reference,
                                                          receive_preproc_buffer_reference,
                                                          send_buffer_reference,
                                                          receive_buffer_reference,
                                                          send_postproc_buffer_reference,
                                                          receive_postproc_buffer_reference)));

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
                    raw_boundary,
                    DBC::DistributedLevee(DBC::DBDataExchanger(index,
                                                               send_preproc_buffer_reference,
                                                               receive_preproc_buffer_reference,
                                                               send_buffer_reference,
                                                               receive_buffer_reference,
                                                               send_postproc_buffer_reference,
                                                               receive_postproc_buffer_reference),
                                          levee));

                n_distr_levee++;

                raw_bound_distr_levee.erase(dbound_key);
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
        writer.GetLogFile() << "Number of distributed boundaries: " << n_distributed << std::endl
                            << "Number of distributed levee boundaries: " << n_distr_levee << std::endl;
    }
}
}
}

#endif