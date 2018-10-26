#ifndef EHDG_GN_PRE_CREATE_HPP
#define EHDG_GN_PRE_CREATE_HPP

#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_intface.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_bound.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_dbound.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_edge_intface.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_edge_bound.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_edge_dbound.hpp"

namespace GN {
namespace EHDG {
std::vector<uint> Problem::comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp) {
    std::vector<uint> offset(GN::EHDG::n_communications);

    offset[SWE::EHDG::CommTypes::bound_state] = begin_index[SWE::EHDG::CommTypes::bound_state];
    offset[CommTypes::dc_global_dof_indx]     = begin_index[CommTypes::dc_global_dof_indx];
    offset[CommTypes::dbath]                  = begin_index[CommTypes::dbath];
    offset[CommTypes::derivatives]            = begin_index[CommTypes::derivatives];

    begin_index[SWE::EHDG::CommTypes::bound_state] += 2 * SWE::n_variables * ngp;
    begin_index[CommTypes::dc_global_dof_indx] += 1;
    begin_index[CommTypes::dbath] += GN::n_ddbath_terms * ngp;
    begin_index[CommTypes::derivatives] += GN::n_du_terms * ngp;

    return offset;
}

template <typename RawBoundaryType>
void Problem::create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& input,
                                ProblemWriterType& writer) {
    GN::create_interfaces<GN::EHDG::Problem>(raw_boundaries, mesh, input, writer);
}

template <typename RawBoundaryType>
void Problem::create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& input,
                                ProblemWriterType& writer) {
    GN::create_boundaries<GN::EHDG::Problem>(raw_boundaries, mesh, input, writer);
}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    ProblemInputType& input,
    Communicator& communicator,
    ProblemWriterType& writer) {
    // *** //
    GN::create_distributed_boundaries<GN::EHDG::Problem>(raw_boundaries, mesh, input, communicator, writer);
}

void Problem::create_edge_interfaces(ProblemMeshType& mesh,
                                     ProblemMeshSkeletonType& mesh_skeleton,
                                     ProblemWriterType& writer) {
    GN::create_edge_interfaces<GN::EHDG::Problem>(mesh, mesh_skeleton, writer);
}

void Problem::create_edge_boundaries(ProblemMeshType& mesh,
                                     ProblemMeshSkeletonType& mesh_skeleton,
                                     ProblemWriterType& writer) {
    GN::create_edge_boundaries<GN::EHDG::Problem>(mesh, mesh_skeleton, writer);
}

void Problem::create_edge_distributeds(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
    GN::create_edge_distributeds<GN::EHDG::Problem>(mesh, mesh_skeleton, writer);
}
}
}

#endif