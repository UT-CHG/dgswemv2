#ifndef EHDG_GN_PRE_INIT_DATA_HPP
#define EHDG_GN_PRE_INIT_DATA_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_dc_data_serial(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input) {
    /* nothing to init specific to DC problem part yet */
}

void Problem::initialize_dc_data_parallel(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input) {
    /* nothing to init specific to DC problem part yet */
}
}
}

#endif