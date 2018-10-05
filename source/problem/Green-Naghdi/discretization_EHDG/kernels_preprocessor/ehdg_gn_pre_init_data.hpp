#ifndef EHDG_GN_PRE_INIT_DATA_HPP
#define EHDG_GN_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

namespace GN {
namespace EHDG {
void Problem::initialize_data_serial(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::EHDG::Problem::initialize_data_serial(mesh, problem_specific_input);
}

void Problem::initialize_data_parallel(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::EHDG::Problem::initialize_data_parallel(mesh, problem_specific_input);
}
}
}

#endif