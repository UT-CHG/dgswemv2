#ifndef EHDG_SWE_PRE_INIT_DATA_HPP
#define EHDG_SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

namespace SWE {
namespace EHDG {
template <typename MeshType>
void Problem::initialize_data_serial(MeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::initialize_data_serial(mesh, problem_specific_input);
}

template <typename MeshType>
void Problem::initialize_data_parallel(MeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::initialize_data_parallel(mesh, problem_specific_input);
}
}
}

#endif