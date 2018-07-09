#include "../bases_2D.hpp"

namespace Basis {
template <typename T>
inline void Dubiner_2D::ProjectBasisToLinear(const std::vector<T>& u, std::vector<T>& u_lin) {
    u_lin[0] = u[0] - u[1] - u[2];
    u_lin[1] = u[0] - u[1] + u[2];
    u_lin[2] = u[0] + 2.0 * u[1];
}

template <typename T>
inline void Dubiner_2D::ProjectLinearToBasis(const std::vector<T>& u_lin, std::vector<T>& u) {
    std::fill(u.begin(), u.end(), 0.0);

    u[0] = (u_lin[0] + u_lin[1] + u_lin[2]) / 3.0;
    u[1] = (-u_lin[0] - u_lin[1] + 2.0 * u_lin[2]) / 6.0;
    u[2] = (-u_lin[0] + u_lin[1]) / 2.0;
}
}