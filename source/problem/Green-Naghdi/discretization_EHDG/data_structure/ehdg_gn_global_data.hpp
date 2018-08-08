#ifndef EHDG_GN_GLOBAL_DATA_HPP
#define EHDG_GN_GLOBAL_DATA_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct GlobalData {
    SparseMatrix<double> w1_w1_hat;
    DynVector<double> w1_rhs;

    SparseMatrix<double> w2_w1;
    SparseMatrix<double> w2_w2_inv;
    SparseMatrix<double> w2_w1_hat;

    SparseMatrix<double> w1_hat_w1;
    SparseMatrix<double> w1_hat_w2;
    SparseMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;
};
}
}

#endif
