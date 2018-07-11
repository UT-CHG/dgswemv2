#ifndef IHDG_SWE_DATA_LOCAL_HPP
#define IHDG_SWE_DATA_LOCAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct Local {
    Local() = default;
    Local(const uint ndof)
        : delta_matrix(SWE::n_variables * ndof, SWE::n_variables * ndof),
          delta_matrix_inv(SWE::n_variables * ndof, SWE::n_variables * ndof),
          rhs(SWE::n_variables * ndof) {}

    DMatrix<double> delta_matrix;
    DMatrix<double> delta_matrix_inv;
    DVector<double> rhs;
};
}
}

#endif
