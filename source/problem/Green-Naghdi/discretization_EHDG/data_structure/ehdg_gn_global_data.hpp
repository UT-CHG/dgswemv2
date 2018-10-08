#ifndef EHDG_GN_GLOBAL_DATA_HPP
#define EHDG_GN_GLOBAL_DATA_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct GlobalData {
#ifndef HAS_PETSC
    SparseMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;
#endif
};
}
}

#endif
