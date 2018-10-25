#ifndef RKDG_SWE_PRE_CREATE_HPP
#define RKDG_SWE_PRE_CREATE_HPP

#include "problem/SWE/problem_preprocessor/swe_preprocessor.hpp"

namespace SWE {
namespace RKDG {
/*std::vector<uint> create_offset() {
    std::vector<uint> offset(SWE::EHDG::n_communications);

    // offset[CommTypes::bound_state] = begin_index[CommTypes::bound_state];

    // begin_index[CommTypes::bound_state] += 2 * SWE::n_variables * ngp;

    return offset;
}*/

std::vector<uint> Problem::comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp) {
    std::vector<uint> offset(SWE::RKDG::n_communications);

    offset[CommTypes::baryctr_coord] = begin_index[CommTypes::baryctr_coord];
    offset[CommTypes::bound_state]   = begin_index[CommTypes::bound_state];
    offset[CommTypes::baryctr_state] = begin_index[CommTypes::baryctr_state];

    begin_index[CommTypes::baryctr_coord] += 2;
    begin_index[CommTypes::bound_state] += SWE::n_variables * ngp + 1;  // + w/d state
    begin_index[CommTypes::baryctr_state] += SWE::n_variables + 1;      // + w/d state

    return offset;
}

template <typename RawBoundaryType>
void Problem::create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& input,
                                ProblemWriterType& writer) {
    SWE::create_interfaces<SWE::RKDG::Problem>(raw_boundaries, mesh, input, writer);
}

template <typename RawBoundaryType>
void Problem::create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& input,
                                ProblemWriterType& writer) {
    SWE::create_boundaries<SWE::RKDG::Problem>(raw_boundaries, mesh, input, writer);
}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_boundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    ProblemMeshType& mesh,
    ProblemInputType& input,
    Communicator& communicator,
    ProblemWriterType& writer) {
    // *** //
    SWE::create_distributed_boundaries<SWE::RKDG::Problem>(raw_boundaries, mesh, input, communicator, writer);
}
}
}

#endif