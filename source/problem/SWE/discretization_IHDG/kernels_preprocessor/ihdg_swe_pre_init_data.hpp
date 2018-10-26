#ifndef IHDG_SWE_PRE_INIT_DATA_HPP
#define IHDG_SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

namespace SWE {
namespace IHDG {
template <typename MeshType>
void Problem::initialize_data_serial(MeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::initialize_data(mesh, problem_specific_input);
}

template <typename MeshType>
void Problem::initialize_data_parallel(MeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::initialize_data(mesh, problem_specific_input);
}
}
}

#endif