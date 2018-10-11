#ifndef EHDG_GN_PRE_CREATE_DBOUND_HPP
#define EHDG_GN_PRE_CREATE_DBOUND_HPP

namespace GN {
namespace EHDG {
template <typename RawBoundaryType>
void Problem::create_distributed_boundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType&,
    ProblemInputType& problem_input,
    std::tuple<>&,
    ProblemWriterType& writer) {}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    ProblemInputType& problem_input,
    Communicator& communicator,
    ProblemWriterType& writer) {
    // *** //
    using DistributedBoundaryTypes = Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>;

    auto& raw_bound_distributed = raw_boundaries[distributed(GN::BoundaryTypes::internal)];

    uint n_distributed = 0;

    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); rank_boundary_id++) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        uint element_id_in, bound_id_in, p, ngp;
        // uint element_id_ex, bound_id_ex;

        uint locality_in, submesh_in;
        uint locality_ex, submesh_ex;

        std::vector<uint> begin_index(GN::EHDG::n_communications, 0);

        // check if the data in rank_boundary_data matches communicator rank boundary
        const RankBoundaryMetaData& rb_meta_data = rank_boundary.db_data;

        locality_in = rb_meta_data.locality_in;
        submesh_in  = rb_meta_data.submesh_in;

        locality_ex = rb_meta_data.locality_ex;
        submesh_ex  = rb_meta_data.submesh_ex;

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

            std::vector<uint> offset(GN::EHDG::n_communications);

            offset[SWE::EHDG::CommTypes::bound_state] = begin_index[SWE::EHDG::CommTypes::bound_state];
            offset[CommTypes::dc_global_dof_indx]     = begin_index[CommTypes::dc_global_dof_indx];
            offset[CommTypes::dbath]                  = begin_index[CommTypes::dbath];
            offset[CommTypes::derivatives]            = begin_index[CommTypes::derivatives];

            begin_index[SWE::EHDG::CommTypes::bound_state] += 2 * SWE::n_variables * ngp;
            begin_index[CommTypes::dc_global_dof_indx] += 1;
            begin_index[CommTypes::dbath] += GN::n_ddbath_terms * ngp;
            begin_index[CommTypes::derivatives] += GN::n_du_terms * ngp;

            if (raw_bound_distributed.find(dbound_key) != raw_bound_distributed.end()) {
                using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;

                auto& raw_boundary = raw_bound_distributed.find(dbound_key)->second;

                raw_boundary.p = p;

                mesh.template CreateDistributedBoundary<DBTypeDistributed>(
                    std::move(raw_boundary),
                    DBC::Distributed(DBDataExchanger(locality_in,
                                                     submesh_in,
                                                     locality_ex,
                                                     submesh_ex,
                                                     std::move(offset),
                                                     rank_boundary.send_buffer,
                                                     rank_boundary.receive_buffer)));

                n_distributed++;

                raw_bound_distributed.erase(dbound_key);
            }
        }

        rank_boundary.send_buffer.resize(GN::EHDG::n_communications);
        rank_boundary.receive_buffer.resize(GN::EHDG::n_communications);

        for (uint comm = 0; comm < GN::EHDG::n_communications; ++comm) {
            rank_boundary.send_buffer[comm].resize(begin_index[comm]);
            rank_boundary.receive_buffer[comm].resize(begin_index[comm]);
        }
    }

    mesh.CallForEachDistributedBoundary([](auto& dbound) { dbound.boundary_condition.Initialize(dbound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed boundaries: " << n_distributed << std::endl;
    }
}
}
}

#endif