#ifndef IHDG_SWE_GLOBAL_DATA_HPP
#define IHDG_SWE_GLOBAL_DATA_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct GlobalData {
#ifndef HAS_PETSC
    SparseMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;
#endif

#ifdef HAS_PETSC
    bool destruct = false;

    Mat delta_hat_global;
    Vec rhs_global;
    KSP ksp;
    PC pc;
    
    IS from, to;
    VecScatter scatter;
    Vec sol;
#endif
};
}
}

#endif
