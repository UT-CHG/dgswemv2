#include "../bases_2D.hpp"

namespace Basis {
template <typename InputArrayType>
inline decltype(auto) Dubiner_2D::ProjectBasisToLinear(const InputArrayType& u) {
    uint nvar = rows(u);

    DynMatrix<double> u_lin(nvar, 3);

    column(u_lin, 0) = column(u, 0) - column(u, 1) - column(u, 2);
    column(u_lin, 1) = column(u, 0) - column(u, 1) + column(u, 2);
    column(u_lin, 2) = column(u, 0) + 2.0 * column(u, 1);

    return u_lin;
}

template <typename InputArrayType>
inline decltype(auto) Dubiner_2D::ProjectLinearToBasis(const uint ndof, const InputArrayType& u_lin) {
    uint nvar = rows(u_lin);

    DynMatrix<double> u(nvar, ndof);

    set_constant(u, 0.0);

    column(u, 0) = (column(u_lin, 0) + column(u_lin, 1) + column(u_lin, 2)) / 3.0;
    column(u, 1) = (-column(u_lin, 0) - column(u_lin, 1) + 2.0 * column(u_lin, 2)) / 6.0;
    column(u, 2) = (-column(u_lin, 0) + column(u_lin, 1)) / 2.0;

    return u;
}
}