#ifndef LEASTSQUARES_DBATH_SERIAL_HPP
#define LEASTSQUARES_DBATH_SERIAL_HPP

#include "interpolation_dbath.hpp"

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dbath_ls(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_ddbath_ls(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_dddbath_ls(ProblemDiscretizationType& discretization);

void Problem::compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization) {
    compute_dbath_ls(discretization);
    interpolate_dbath(discretization);

    compute_ddbath_ls(discretization);
    interpolate_ddbath(discretization);

    compute_dddbath_ls(discretization);
    interpolate_dddbath(discretization);
}

template <typename ProblemDiscretizationType>
void compute_dbath_ls(ProblemDiscretizationType& discretization) {
}

template <typename ProblemDiscretizationType>
void compute_ddbath_ls(ProblemDiscretizationType& discretization) {
}

template <typename ProblemDiscretizationType>
void compute_dddbath_ls(ProblemDiscretizationType& discretization) {
}
}
}

#endif
