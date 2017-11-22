#ifndef SWE_PRE_CREATE_DBOUND_HPP
#define SWE_PRE_CREATE_DBOUND_HPP

namespace SWE {
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
}

#endif