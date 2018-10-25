#ifndef SWE_PRE_CREATE_DBOUND_HPP
#define SWE_PRE_CREATE_DBOUND_HPP

namespace SWE {
template <typename ProblemType, typename RawBoundaryType, typename Communicator>
void create_distributed_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                   typename ProblemType::ProblemMeshType& mesh,
                                   typename ProblemType::ProblemInputType& problem_input,
                                   Communicator& communicator,
                                   typename ProblemType::ProblemWriterType& writer) {
    // *** //
    using DBTypeDistributed =
        typename std::tuple_element<0, typename ProblemType::ProblemDistributedBoundaryTypes>::type;
    using DBTypeDistributedLevee =
        typename std::tuple_element<1, typename ProblemType::ProblemDistributedBoundaryTypes>::type;

    auto& raw_bound_distributed = raw_boundaries[distributed(SWE::BoundaryTypes::internal)];
    auto& raw_bound_distr_levee = raw_boundaries[distributed(SWE::BoundaryTypes::levee)];

    uint n_distributed = 0;
    uint n_distr_levee = 0;

    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); rank_boundary_id++) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        uint element_id_in, bound_id_in, p, ngp;
        // uint element_id_ex, bound_id_ex;

        uint locality_in, submesh_in;
        uint locality_ex, submesh_ex;

        std::vector<uint> begin_index(ProblemType::n_communications, 0);

        // check if the data in rank_boundary_data matches communicator rank boundary
        const RankBoundaryMetaData& rb_meta_data = rank_boundary.db_data;

        locality_in = rb_meta_data.locality_in;
        submesh_in  = rb_meta_data.submesh_in;

        locality_ex = rb_meta_data.locality_ex;
        submesh_ex  = rb_meta_data.submesh_ex;

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
                typename DBTypeDistributed::BoundaryIntegrationType boundary_integration;

                ngp = boundary_integration.GetNumGP(2 * p + 1);
            } else if (raw_bound_distr_levee.find(dbound_key) != raw_bound_distr_levee.end()) {
                typename DBTypeDistributedLevee::BoundaryIntegrationType boundary_integration;

                ngp = boundary_integration.GetNumGP(2 * p + 1);
            } else {
                throw std::logic_error("Fatal Error: unable to find raw distributed boundary!\n");
            }

            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                auto& raw_boundary = raw_bound_distributed.find(dbound_key)->second;

                raw_boundary.p = p;

                mesh.template CreateDistributedBoundary<DBTypeDistributed>(
                    std::move(raw_boundary),
                    DBDataExchanger(locality_in,
                                    submesh_in,
                                    locality_ex,
                                    submesh_ex,
                                    ProblemType::comm_buffer_offsets(begin_index, ngp),
                                    rank_boundary.send_buffer,
                                    rank_boundary.receive_buffer));

                n_distributed++;

                raw_bound_distributed.erase(dbound_key);
            } else if (raw_bound_distr_levee.find(dbound_key) != raw_bound_distr_levee.end()) {
                auto& raw_boundary = raw_bound_distr_levee.find(dbound_key)->second;

                raw_boundary.p = p;

                auto& levee_data = problem_input.levee_is_data;

                std::vector<LeveeInput> levee;

                for (uint node = 0; node < raw_boundary.node_ID.size(); ++node) {
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
                    DBDataExchanger(locality_in,
                                    submesh_in,
                                    locality_ex,
                                    submesh_ex,
                                    ProblemType::comm_buffer_offsets(begin_index, ngp),
                                    rank_boundary.send_buffer,
                                    rank_boundary.receive_buffer),
                    levee);

                n_distr_levee++;

                raw_bound_distr_levee.erase(dbound_key);
            }
        }

        rank_boundary.send_buffer.resize(ProblemType::n_communications);
        rank_boundary.receive_buffer.resize(ProblemType::n_communications);

        for (uint comm = 0; comm < ProblemType::n_communications; ++comm) {
            rank_boundary.send_buffer[comm].resize(begin_index[comm]);
            rank_boundary.receive_buffer[comm].resize(begin_index[comm]);
        }
    }

    mesh.CallForEachDistributedBoundary([](auto& dbound) { dbound.boundary_condition.Initialize(dbound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed boundaries: " << n_distributed << std::endl
                            << "Number of distributed levee boundaries: " << n_distr_levee << std::endl;
    }
}
}

#endif