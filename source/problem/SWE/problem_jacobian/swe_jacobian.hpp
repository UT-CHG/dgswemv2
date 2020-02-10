#ifndef SWE_JACOBIAN_HPP
#define SWE_JACOBIAN_HPP

#include "utilities/sign.hpp"

namespace SWE {
struct parameters {
    double h, u, v, c, nx, ny;
};

StatMatrix<double, SWE::n_variables, SWE::n_variables> L(const parameters& param);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dze(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqx(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqy(const parameters& param);  // TODO

StatMatrix<double, SWE::n_variables, SWE::n_variables> absL(const parameters& param);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dze(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqx(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqy(const parameters& param);  // TODO

StatMatrix<double, SWE::n_variables, SWE::n_variables> R(const parameters& param);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dze(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqx(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqy(const parameters& param);  // TODO

StatMatrix<double, SWE::n_variables, SWE::n_variables> invR(const parameters& param);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dze(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqx(const parameters& param);  // TODO
StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqy(const parameters& param);  // TODO

void get_A(const HybMatrix<double, SWE::n_variables>& q,
           const HybMatrix<double, SWE::n_auxiliaries>& aux,
           const HybMatrix<double, SWE::n_dimensions>& surface_normal,
           AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& A) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        A[gp] = R(param) * L(param) * invR(param);
    }
}

void get_dA_dze(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dA_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dA_dze[gp] = dR_dze(param) * L(param) * invR(param) + R(param) * dL_dze(param) * invR(param) +
                     R(param) * L(param) * dinvR_dze(param);
    }
}

void get_dA_dqx(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dA_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dA_dqx[gp] = dR_dqx(param) * L(param) * invR(param) + R(param) * dL_dqx(param) * invR(param) +
                     R(param) * L(param) * dinvR_dqx(param);
    }
}

void get_dA_dqy(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dA_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dA_dqy[gp] = dR_dqy(param) * L(param) * invR(param) + R(param) * dL_dqy(param) * invR(param) +
                     R(param) * L(param) * dinvR_dqy(param);
    }
}

void get_absA(const HybMatrix<double, SWE::n_variables>& q,
              const HybMatrix<double, SWE::n_auxiliaries>& aux,
              const HybMatrix<double, SWE::n_dimensions>& surface_normal,
              AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& absA) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        absA[gp] = R(param) * absL(param) * invR(param);
    }
}

void get_dabsA_dze(const HybMatrix<double, SWE::n_variables>& q,
                   const HybMatrix<double, SWE::n_auxiliaries>& aux,
                   const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                   AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dabsA_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dabsA_dze[gp] = dR_dze(param) * absL(param) * invR(param) + R(param) * dabsL_dze(param) * invR(param) +
                        R(param) * absL(param) * dinvR_dze(param);
    }
}

void get_dabsA_dqx(const HybMatrix<double, SWE::n_variables>& q,
                   const HybMatrix<double, SWE::n_auxiliaries>& aux,
                   const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                   AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dabsA_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dabsA_dqx[gp] = dR_dqx(param) * absL(param) * invR(param) + R(param) * dabsL_dqx(param) * invR(param) +
                        R(param) * absL(param) * dinvR_dqx(param);
    }
}

void get_dabsA_dqy(const HybMatrix<double, SWE::n_variables>& q,
                   const HybMatrix<double, SWE::n_auxiliaries>& aux,
                   const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                   AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dabsA_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dabsA_dqy[gp] = dR_dqy(param) * absL(param) * invR(param) + R(param) * dabsL_dqy(param) * invR(param) +
                        R(param) * absL(param) * dinvR_dqy(param);
    }
}

void get_Aplus(const HybMatrix<double, SWE::n_variables>& q,
               const HybMatrix<double, SWE::n_auxiliaries>& aux,
               const HybMatrix<double, SWE::n_dimensions>& surface_normal,
               AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& Aplus) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        Aplus[gp] = 0.5 * R(param) * (L(param) + absL(param)) * invR(param);
    }
}

void get_dAplus_dze(const HybMatrix<double, SWE::n_variables>& q,
                    const HybMatrix<double, SWE::n_auxiliaries>& aux,
                    const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dAplus_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dAplus_dze[gp] = 0.5 * (dR_dze(param) * (L(param) + absL(param)) * invR(param) +
                                R(param) * (dL_dze(param) + dabsL_dze(param)) * invR(param) +
                                R(param) * (L(param) + absL(param)) * dinvR_dze(param));
    }
}

void get_dAplus_dqx(const HybMatrix<double, SWE::n_variables>& q,
                    const HybMatrix<double, SWE::n_auxiliaries>& aux,
                    const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dAplus_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dAplus_dqx[gp] = 0.5 * (dR_dqx(param) * (L(param) + absL(param)) * invR(param) +
                                R(param) * (dL_dqx(param) + dabsL_dqx(param)) * invR(param) +
                                R(param) * (L(param) + absL(param)) * dinvR_dqx(param));
    }
}

void get_dAplus_dqy(const HybMatrix<double, SWE::n_variables>& q,
                    const HybMatrix<double, SWE::n_auxiliaries>& aux,
                    const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dAplus_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dAplus_dqy[gp] = 0.5 * (dR_dqy(param) * (L(param) + absL(param)) * invR(param) +
                                R(param) * (dL_dqy(param) + dabsL_dqy(param)) * invR(param) +
                                R(param) * (L(param) + absL(param)) * dinvR_dqy(param));
    }
}

void get_Aminus(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& Aminus) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        Aminus[gp] = 0.5 * R(param) * (L(param) - absL(param)) * invR(param);
    }
}

void get_dAminus_dze(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dAminus_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dAminus_dze[gp] = 0.5 * (dR_dze(param) * (L(param) - absL(param)) * invR(param) +
                                 R(param) * (dL_dze(param) - dabsL_dze(param)) * invR(param) +
                                 R(param) * (L(param) - absL(param)) * dinvR_dze(param));
    }
}

void get_dAminus_dqx(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dAminus_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dAminus_dqx[gp] = 0.5 * (dR_dqx(param) * (L(param) - absL(param)) * invR(param) +
                                 R(param) * (dL_dqx(param) - dabsL_dqx(param)) * invR(param) +
                                 R(param) * (L(param) - absL(param)) * dinvR_dqx(param));
    }
}

void get_dAminus_dqy(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>>& dAminus_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        const parameters param{aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp),
                               q(SWE::Variables::hc, gp) / aux(SWE::Auxiliaries::h, gp),
                               surface_normal(GlobalCoord::x, gp),
                               surface_normal(GlobalCoord::y, gp)};
        dAminus_dqy[gp] = 0.5 * (dR_dqy(param) * (L(param) - absL(param)) * invR(param) +
                                 R(param) * (dL_dqy(param) - dabsL_dqy(param)) * invR(param) +
                                 R(param) * (L(param) - absL(param)) * dinvR_dqy(param));
    }
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> L(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> L;

    double c  = std::sqrt(Global::g * param.h);
    double un = param.u * param.nx + param.v * param.ny;

    set_constant(L, 0.0);

    L(SWE::Variables::ze, SWE::Variables::ze) = un - c;
    L(SWE::Variables::qx, SWE::Variables::qx) = un;
    L(SWE::Variables::qy, SWE::Variables::qy) = un + c;
    L(SWE::Variables::hc, SWE::Variables::hc) = un;  // TODO

    return L;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dze(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dze;

    double dc_dze  = std::sqrt(Global::g / param.h) / 2.0;
    double dun_dze = -(param.u * param.nx + param.v * param.ny) / param.h;

    set_constant(dL_dze, 0.0);

    dL_dze(SWE::Variables::ze, SWE::Variables::ze) = dun_dze - dc_dze;
    dL_dze(SWE::Variables::qx, SWE::Variables::qx) = dun_dze;
    dL_dze(SWE::Variables::qy, SWE::Variables::qy) = dun_dze + dc_dze;

    return dL_dze;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqx(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqx;

    double dun_dqx = param.nx / param.h;

    set_constant(dL_dqx, 0.0);

    dL_dqx(SWE::Variables::ze, SWE::Variables::ze) = dun_dqx;
    dL_dqx(SWE::Variables::qx, SWE::Variables::qx) = dun_dqx;
    dL_dqx(SWE::Variables::qy, SWE::Variables::qy) = dun_dqx;

    return dL_dqx;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqy(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqy;

    double dun_dqy = param.ny / param.h;

    set_constant(dL_dqy, 0.0);

    dL_dqy(SWE::Variables::ze, SWE::Variables::ze) = dun_dqy;
    dL_dqy(SWE::Variables::qx, SWE::Variables::qx) = dun_dqy;
    dL_dqy(SWE::Variables::qy, SWE::Variables::qy) = dun_dqy;

    return dL_dqy;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> absL(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> absL;

    double c  = std::sqrt(Global::g * param.h);
    double un = param.u * param.nx + param.v * param.ny;

    set_constant(absL, 0.0);

    absL(SWE::Variables::ze, SWE::Variables::ze) = std::abs(un - c);
    absL(SWE::Variables::qx, SWE::Variables::qx) = std::abs(un);
    absL(SWE::Variables::qy, SWE::Variables::qy) = std::abs(un + c);
    absL(SWE::Variables::hc, SWE::Variables::hc) = std::abs(un);

    return absL;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dze(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dze;

    double c  = std::sqrt(Global::g * param.h);
    double un = param.u * param.nx + param.v * param.ny;

    double dc_dze  = std::sqrt(Global::g / param.h) / 2.0;
    double dun_dze = -(param.u * param.nx + param.v * param.ny) / param.h;

    set_constant(dabsL_dze, 0.0);

    dabsL_dze(SWE::Variables::ze, SWE::Variables::ze) = (dun_dze - dc_dze) * Utilities::sign(un - c);
    dabsL_dze(SWE::Variables::qx, SWE::Variables::qx) = dun_dze * Utilities::sign(un);
    dabsL_dze(SWE::Variables::qy, SWE::Variables::qy) = (dun_dze + dc_dze) * Utilities::sign(un + c);

    return dabsL_dze;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqx(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqx;

    double c  = std::sqrt(Global::g * param.h);
    double un = param.u * param.nx + param.v * param.ny;

    double dun_dqx = param.nx / param.h;

    set_constant(dabsL_dqx, 0.0);

    dabsL_dqx(SWE::Variables::ze, SWE::Variables::ze) = dun_dqx * Utilities::sign(un - c);
    dabsL_dqx(SWE::Variables::qx, SWE::Variables::qx) = dun_dqx * Utilities::sign(un);
    dabsL_dqx(SWE::Variables::qy, SWE::Variables::qy) = dun_dqx * Utilities::sign(un + c);

    return dabsL_dqx;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqy(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqy;

    double c  = std::sqrt(Global::g * param.h);
    double un = param.u * param.nx + param.v * param.ny;

    double dun_dqy = param.ny / param.h;

    set_constant(dabsL_dqy, 0.0);

    dabsL_dqy(SWE::Variables::ze, SWE::Variables::ze) = dun_dqy * Utilities::sign(un - c);
    dabsL_dqy(SWE::Variables::qx, SWE::Variables::qx) = dun_dqy * Utilities::sign(un);
    dabsL_dqy(SWE::Variables::qy, SWE::Variables::qy) = dun_dqy * Utilities::sign(un + c);

    return dabsL_dqy;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> R(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> R;

    double c = std::sqrt(Global::g * param.h);

    R(SWE::Variables::ze, SWE::Variables::ze) = 1.0;
    R(SWE::Variables::ze, SWE::Variables::qx) = 0.0;
    R(SWE::Variables::ze, SWE::Variables::qy) = 1.0;
    R(SWE::Variables::ze, SWE::Variables::hc) = 0.0;

    R(SWE::Variables::qx, SWE::Variables::ze) = param.u - c * param.nx;
    R(SWE::Variables::qx, SWE::Variables::qx) = -param.ny;
    R(SWE::Variables::qx, SWE::Variables::qy) = param.u + c * param.nx;
    R(SWE::Variables::qx, SWE::Variables::hc) = 0.0;

    R(SWE::Variables::qy, SWE::Variables::ze) = param.v - c * param.ny;
    R(SWE::Variables::qy, SWE::Variables::qx) = param.nx;
    R(SWE::Variables::qy, SWE::Variables::qy) = param.v + c * param.ny;
    R(SWE::Variables::qy, SWE::Variables::hc) = 0;

    R(SWE::Variables::hc, SWE::Variables::ze) = param.c;  // TODO
    R(SWE::Variables::hc, SWE::Variables::qx) = 0.0;      // TODO
    R(SWE::Variables::hc, SWE::Variables::qy) = param.c;  // TODO
    R(SWE::Variables::hc, SWE::Variables::hc) = 1.0;      // TODO

    return R;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dze(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dze;

    double dc_dze = std::sqrt(Global::g / param.h) / 2.0;

    set_constant(dR_dze, 0.0);

    dR_dze(SWE::Variables::qx, SWE::Variables::ze) = -param.u / param.h - dc_dze * param.nx;
    dR_dze(SWE::Variables::qx, SWE::Variables::qy) = -param.u / param.h + dc_dze * param.nx;
    dR_dze(SWE::Variables::qy, SWE::Variables::ze) = -param.v / param.h - dc_dze * param.ny;
    dR_dze(SWE::Variables::qy, SWE::Variables::qy) = -param.v / param.h + dc_dze * param.ny;

    return dR_dze;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqx(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqx;

    set_constant(dR_dqx, 0.0);

    dR_dqx(SWE::Variables::qx, SWE::Variables::ze) = 1 / param.h;
    dR_dqx(SWE::Variables::qx, SWE::Variables::qy) = 1 / param.h;

    return dR_dqx;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqy(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqy;

    set_constant(dR_dqy, 0.0);

    dR_dqy(SWE::Variables::qy, SWE::Variables::ze) = 1 / param.h;
    dR_dqy(SWE::Variables::qy, SWE::Variables::qy) = 1 / param.h;

    return dR_dqy;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> invR(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> invR;

    double c  = std::sqrt(Global::g * param.h);
    double un = param.u * param.nx + param.v * param.ny;

    invR(SWE::Variables::ze, SWE::Variables::ze) = (c + un) / (2 * c);
    invR(SWE::Variables::ze, SWE::Variables::qx) = -param.nx / (2 * c);
    invR(SWE::Variables::ze, SWE::Variables::qy) = -param.ny / (2 * c);
    invR(SWE::Variables::ze, SWE::Variables::hc) = 0.0;

    invR(SWE::Variables::qx, SWE::Variables::ze) = param.u * param.ny - param.v * param.nx;
    invR(SWE::Variables::qx, SWE::Variables::qx) = -param.ny;
    invR(SWE::Variables::qx, SWE::Variables::qy) = param.nx;
    invR(SWE::Variables::qx, SWE::Variables::hc) = 0.0;

    invR(SWE::Variables::qy, SWE::Variables::ze) = (c - un) / (2 * c);
    invR(SWE::Variables::qy, SWE::Variables::qx) = param.nx / (2 * c);
    invR(SWE::Variables::qy, SWE::Variables::qy) = param.ny / (2 * c);
    invR(SWE::Variables::qy, SWE::Variables::hc) = 0.0;

    invR(SWE::Variables::hc, SWE::Variables::ze) = -param.c;  // TODO
    invR(SWE::Variables::hc, SWE::Variables::qx) = 0.0;       // TODO
    invR(SWE::Variables::hc, SWE::Variables::qy) = 0.0;       // TODO
    invR(SWE::Variables::hc, SWE::Variables::hc) = 1.0;       // TODO

    return invR;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dze(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dze;

    double c    = std::sqrt(Global::g * param.h);
    double c_sq = Global::g * param.h;
    double un   = param.u * param.nx + param.v * param.ny;

    double dc_dze = std::sqrt(Global::g / param.h) / 2.0;

    set_constant(dinvR_dze, 0.0);

    dinvR_dze(SWE::Variables::ze, SWE::Variables::ze) = -3 * un / (4 * param.h * c);
    dinvR_dze(SWE::Variables::ze, SWE::Variables::qx) = dc_dze * param.nx / (2 * c_sq);
    dinvR_dze(SWE::Variables::ze, SWE::Variables::qy) = dc_dze * param.ny / (2 * c_sq);

    dinvR_dze(SWE::Variables::qx, SWE::Variables::ze) = -(param.u * param.ny - param.v * param.nx) / param.h;

    dinvR_dze(SWE::Variables::qy, SWE::Variables::ze) = 3 * un / (4 * param.h * c);
    dinvR_dze(SWE::Variables::qy, SWE::Variables::qx) = -dc_dze * param.nx / (2 * c_sq);
    dinvR_dze(SWE::Variables::qy, SWE::Variables::qy) = -dc_dze * param.ny / (2 * c_sq);

    return dinvR_dze;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqx(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqx;

    double c = std::sqrt(Global::g * param.h);

    set_constant(dinvR_dqx, 0.0);

    dinvR_dqx(SWE::Variables::ze, SWE::Variables::ze) = param.nx / param.h / (2 * c);
    dinvR_dqx(SWE::Variables::qx, SWE::Variables::ze) = param.ny / param.h;
    dinvR_dqx(SWE::Variables::qy, SWE::Variables::ze) = -param.nx / param.h / (2 * c);

    return dinvR_dqx;
}

inline StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqy(const parameters& param) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqy;

    double c = std::sqrt(Global::g * param.h);

    set_constant(dinvR_dqy, 0.0);

    dinvR_dqy(SWE::Variables::ze, SWE::Variables::ze) = param.ny / param.h / (2 * c);
    dinvR_dqy(SWE::Variables::qx, SWE::Variables::ze) = -param.nx / param.h;
    dinvR_dqy(SWE::Variables::qy, SWE::Variables::ze) = -param.ny / param.h / (2 * c);

    return dinvR_dqy;
}
}

#endif