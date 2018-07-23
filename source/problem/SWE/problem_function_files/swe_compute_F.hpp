#ifndef SWE_COMPUTE_F_HPP
#define SWE_COMPUTE_F_HPP

#include "general_definitions.hpp"

namespace SWE {
inline void compute_F(DynMatrix<double>& q_at_gp,
                      DynMatrix<double>& aux_at_gp,
                      DynMatrix<double>& Fx_at_gp,
                      DynMatrix<double>& Fy_at_gp) {
    row(aux_at_gp, SWE::Auxiliaries::h) = row(q_at_gp, SWE::Variables::ze) + row(aux_at_gp, SWE::Auxiliaries::bath);

    auto u = cwise_division(row(q_at_gp, SWE::Variables::qx), row(aux_at_gp, SWE::Auxiliaries::h));
    auto v = cwise_division(row(q_at_gp, SWE::Variables::qy), row(aux_at_gp, SWE::Auxiliaries::h));

    auto uuh = cwise_multiplication(u, row(q_at_gp, SWE::Variables::qx));
    auto vvh = cwise_multiplication(v, row(q_at_gp, SWE::Variables::qy));
    auto uvh = cwise_multiplication(u, row(q_at_gp, SWE::Variables::qy));
    auto pe =
        Global::g * (0.5 * pow(row(q_at_gp, SWE::Variables::ze), 2.0) +
                     cwise_multiplication(row(q_at_gp, SWE::Variables::ze), row(aux_at_gp, SWE::Auxiliaries::bath)));

    row(Fx_at_gp, SWE::Variables::ze) =
        cwise_multiplication(row(aux_at_gp, SWE::Auxiliaries::sp), row(q_at_gp, SWE::Variables::qx));
    row(Fx_at_gp, SWE::Variables::qx) = cwise_multiplication(row(aux_at_gp, SWE::Auxiliaries::sp), uuh + pe);
    row(Fx_at_gp, SWE::Variables::qy) = cwise_multiplication(row(aux_at_gp, SWE::Auxiliaries::sp), uvh);

    row(Fy_at_gp, SWE::Variables::ze) = row(q_at_gp, SWE::Variables::qy);
    row(Fy_at_gp, SWE::Variables::qx) = uvh;
    row(Fy_at_gp, SWE::Variables::qy) = vvh + pe;
}
}

#endif