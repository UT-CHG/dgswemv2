#ifndef IHDG_SWE_DATA_LOCAL_HPP
#define IHDG_SWE_DATA_LOCAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct Local {
    Local() = default;
    Local(const uint ndof, const uint ngp)
        : A(SWE::n_variables * ndof, SWE::n_variables * ndof),
          B(SWE::n_variables * ndof, SWE::n_variables * ndof),
          c(SWE::n_variables * ndof) {}

    DMatrix<double> A;
    DMatrix<double> B;
    DVector<double> c;
};
}
}

#endif
