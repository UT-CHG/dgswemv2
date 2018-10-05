#ifndef EHDG_GN_POST_WRITE_VTU_HPP
#define EHDG_GN_POST_WRITE_VTU_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::write_VTU_data(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
    SWE::EHDG::Problem::write_VTU_data(mesh, raw_data_file);
}
}
}

#endif