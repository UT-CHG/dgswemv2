#ifndef IHDG_SWE_GLOBAL_DATA_HPP
#define IHDG_SWE_GLOBAL_DATA_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct GlobalData {
    SparseMatrix<double> delta_local_inv;
    SparseMatrix<double> delta_hat_local;
    DynVector<double> rhs_local;

    SparseMatrix<double> delta_global;
    SparseMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;
};
}
}

#endif
