#ifndef EHDG_GN_POST_WRITE_MODAL_HPP
#define EHDG_GN_POST_WRITE_MODAL_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::write_modal_data(const RKStepper& stepper, ProblemMeshType& mesh, const std::string& output_path) {
    SWE::EHDG::Problem::write_modal_data(stepper, mesh, output_path);
}
}
}

#endif