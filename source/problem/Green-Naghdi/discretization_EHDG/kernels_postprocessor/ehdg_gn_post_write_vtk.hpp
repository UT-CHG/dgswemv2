#ifndef EHDG_GN_POST_WRITE_VTK_HPP
#define EHDG_GN_POST_WRITE_VTK_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
template <typename MeshType>
void Problem::write_VTK_data(MeshType& mesh, std::ofstream& raw_data_file) {
    SWE::EHDG::Problem::write_VTK_data(mesh, raw_data_file);
}
}
}

#endif