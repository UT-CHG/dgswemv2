#ifndef IHDG_SWE_GLOBAL_DATA_HPP
#define IHDG_SWE_GLOBAL_DATA_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct GlobalData {
    bool destruct = false;

    uint n_global_dofs;
    std::vector<uint> global_dof_indx;

    Mat delta_hat_global;
    Vec rhs_global;
    KSP ksp;

    IS from, to;
    VecScatter scatter;
    Vec sol;
};
}
}

#endif
