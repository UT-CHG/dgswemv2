#ifndef SWE_STAB_PARAM_LF_HPP
#define SWE_STAB_PARAM_LF_HPP

#include "utilities/sign.hpp"

namespace SWE {
void get_tau_LF(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& tau,
                const double gravity = Global::g) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double c  = std::sqrt(gravity * h);
        double un = u * nx + v * ny;

        tau[gp] = (c + std::abs(un)) * IdentityMatrix<double>(SWE::n_variables);
    }
}

void get_tau_LF(const Column<HybMatrix<double, SWE::n_variables>>& q,
                const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux,
                const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                StatMatrix<double, SWE::n_variables, SWE::n_variables>& tau,
                const double gravity = Global::g) {
    double h = aux[SWE::Auxiliaries::h];
    double u = q[SWE::Variables::qx] / h;
    double v = q[SWE::Variables::qy] / h;

    double nx = surface_normal[GlobalCoord::x];
    double ny = surface_normal[GlobalCoord::y];

    double c  = std::sqrt(gravity * h);
    double un = u * nx + v * ny;

    tau = (c + std::abs(un)) * IdentityMatrix<double>(SWE::n_variables);
}

void get_dtau_dze_LF(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dtau_dze) {  // TODO
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double un      = u * nx + v * ny;
        double dc_dze  = std::sqrt(Global::g / h) / 2.0;
        double dun_dze = -un / h;

        dtau_dze[gp] = (dc_dze + dun_dze * Utilities::sign(un)) * IdentityMatrix<double>(3);  // TODO
    }
}

void get_dtau_dqx_LF(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dtau_dqx) {  // TODO
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double un      = u * nx + v * ny;
        double dun_dqx = nx / h;

        dtau_dqx[gp] = dun_dqx * Utilities::sign(un) * IdentityMatrix<double>(3);  // TODO
    }
}

void get_dtau_dqy_LF(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dtau_dqy) {  // TODO
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        double un      = u * nx + v * ny;
        double dun_dqy = ny / h;

        dtau_dqy[gp] = dun_dqy * Utilities::sign(un) * IdentityMatrix<double>(3);  // TODO
    }
}
}

#endif