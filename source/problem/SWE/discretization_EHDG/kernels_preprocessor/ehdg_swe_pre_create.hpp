#ifndef EHDG_SWE_PRE_CREATE_HPP
#define EHDG_SWE_PRE_CREATE_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_create_intface.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_bound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_dbound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_edge_intface.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_edge_bound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_edge_dbound.hpp"

namespace SWE {
namespace EHDG {
std::vector<uint> Problem::comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp) {
    std::vector<uint> offset(SWE::EHDG::n_communications);

    offset[CommTypes::baryctr_coord]    = begin_index[CommTypes::baryctr_coord];
    offset[CommTypes::init_global_prob] = begin_index[CommTypes::init_global_prob];
    offset[CommTypes::bound_state]      = begin_index[CommTypes::bound_state];
    offset[CommTypes::baryctr_state]    = begin_index[CommTypes::baryctr_state];

    begin_index[CommTypes::baryctr_coord] += 2;
    begin_index[CommTypes::init_global_prob] += ngp;  // bath_at_gp
    begin_index[CommTypes::bound_state] += SWE::n_variables * ngp;
    begin_index[CommTypes::baryctr_state] += SWE::n_variables + 1;  // + w/d state

    return offset;
}

template <typename RawBoundaryType>
void Problem::create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& input,
                                ProblemWriterType& writer) {
    SWE::create_interfaces<SWE::EHDG::Problem>(raw_boundaries, mesh, input, writer);
}

template <typename RawBoundaryType>
void Problem::create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& input,
                                ProblemWriterType& writer) {
    SWE::create_boundaries<SWE::EHDG::Problem>(raw_boundaries, mesh, input, writer);
}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    ProblemInputType& input,
    Communicator& communicator,
    ProblemWriterType& writer) {
    // *** //
    SWE::create_distributed_boundaries<SWE::EHDG::Problem>(raw_boundaries, mesh, input, communicator, writer);
}

void Problem::create_edge_interfaces(ProblemMeshType& mesh,
                                     ProblemMeshSkeletonType& mesh_skeleton,
                                     ProblemWriterType& writer) {
    SWE::create_edge_interfaces<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);
}

void Problem::create_edge_boundaries(ProblemMeshType& mesh,
                                     ProblemMeshSkeletonType& mesh_skeleton,
                                     ProblemWriterType& writer) {
    SWE::create_edge_boundaries<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);
}

void Problem::create_edge_distributeds(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
    SWE::create_edge_distributeds<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);
}
}
}

#endif