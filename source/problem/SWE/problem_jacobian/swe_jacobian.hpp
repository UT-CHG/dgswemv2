#ifndef SWE_JACOBIAN_HPP
#define SWE_JACOBIAN_HPP

namespace SWE {
StatMatrix<double, SWE::n_variables, SWE::n_variables> L(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dze(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqx(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqy(double h, double u, double v, double nx, double ny);

StatMatrix<double, SWE::n_variables, SWE::n_variables> absL(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dze(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqx(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqy(double h, double u, double v, double nx, double ny);

StatMatrix<double, SWE::n_variables, SWE::n_variables> R(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dze(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqx(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqy(double h, double u, double v, double nx, double ny);

StatMatrix<double, SWE::n_variables, SWE::n_variables> invR(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dze(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqx(double h, double u, double v, double nx, double ny);
StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqy(double h, double u, double v, double nx, double ny);

void get_A(const HybMatrix<double, SWE::n_variables>& q,
           const HybMatrix<double, SWE::n_auxiliaries>& aux,
           const HybMatrix<double, SWE::n_dimensions>& surface_normal,
           HybMatrix<double, SWE::n_variables * SWE::n_variables>& A) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(A, gp) = flatten<double>(R(h, u, v, nx, ny) * L(h, u, v, nx, ny) * invR(h, u, v, nx, ny));
    }
}

void get_dA_dze(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                HybMatrix<double, SWE::n_variables * SWE::n_variables>& dA_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dA_dze, gp) = flatten<double>(dR_dze(h, u, v, nx, ny) * L(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                                             R(h, u, v, nx, ny) * dL_dze(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                                             R(h, u, v, nx, ny) * L(h, u, v, nx, ny) * dinvR_dze(h, u, v, nx, ny));
    }
}

void get_dA_dqx(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                HybMatrix<double, SWE::n_variables * SWE::n_variables>& dA_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dA_dqx, gp) = flatten<double>(dR_dqx(h, u, v, nx, ny) * L(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                                             R(h, u, v, nx, ny) * dL_dqx(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                                             R(h, u, v, nx, ny) * L(h, u, v, nx, ny) * dinvR_dqx(h, u, v, nx, ny));
    }
}

void get_dA_dqy(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                HybMatrix<double, SWE::n_variables * SWE::n_variables>& dA_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dA_dqy, gp) = flatten<double>(dR_dqy(h, u, v, nx, ny) * L(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                                             R(h, u, v, nx, ny) * dL_dqy(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                                             R(h, u, v, nx, ny) * L(h, u, v, nx, ny) * dinvR_dqy(h, u, v, nx, ny));
    }
}

void get_absA(const HybMatrix<double, SWE::n_variables>& q,
              const HybMatrix<double, SWE::n_auxiliaries>& aux,
              const HybMatrix<double, SWE::n_dimensions>& surface_normal,
              HybMatrix<double, SWE::n_variables * SWE::n_variables>& absA) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(absA, gp) = flatten<double>(R(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * invR(h, u, v, nx, ny));
    }
}

void get_dabsA_dze(const HybMatrix<double, SWE::n_variables>& q,
                   const HybMatrix<double, SWE::n_auxiliaries>& aux,
                   const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                   HybMatrix<double, SWE::n_variables * SWE::n_variables>& dabsA_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dabsA_dze, gp) =
            flatten<double>(dR_dze(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                            R(h, u, v, nx, ny) * dabsL_dze(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                            R(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * dinvR_dze(h, u, v, nx, ny));
    }
}

void get_dabsA_dqx(const HybMatrix<double, SWE::n_variables>& q,
                   const HybMatrix<double, SWE::n_auxiliaries>& aux,
                   const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                   HybMatrix<double, SWE::n_variables * SWE::n_variables>& dabsA_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dabsA_dqx, gp) =
            flatten<double>(dR_dqx(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                            R(h, u, v, nx, ny) * dabsL_dqx(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                            R(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * dinvR_dqx(h, u, v, nx, ny));
    }
}

void get_dabsA_dqy(const HybMatrix<double, SWE::n_variables>& q,
                   const HybMatrix<double, SWE::n_auxiliaries>& aux,
                   const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                   HybMatrix<double, SWE::n_variables * SWE::n_variables>& dabsA_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dabsA_dqy, gp) =
            flatten<double>(dR_dqy(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                            R(h, u, v, nx, ny) * dabsL_dqy(h, u, v, nx, ny) * invR(h, u, v, nx, ny) +
                            R(h, u, v, nx, ny) * absL(h, u, v, nx, ny) * dinvR_dqy(h, u, v, nx, ny));
    }
}

void get_Aplus(const HybMatrix<double, SWE::n_variables>& q,
               const HybMatrix<double, SWE::n_auxiliaries>& aux,
               const HybMatrix<double, SWE::n_dimensions>& surface_normal,
               HybMatrix<double, SWE::n_variables * SWE::n_variables>& Aplus) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(Aplus, gp) = 0.5 * flatten<double>(R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) *
                                                  invR(h, u, v, nx, ny));
    }
}

void get_dAplus_dze(const HybMatrix<double, SWE::n_variables>& q,
                    const HybMatrix<double, SWE::n_auxiliaries>& aux,
                    const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                    HybMatrix<double, SWE::n_variables * SWE::n_variables>& dAplus_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dAplus_dze, gp) =
            0.5 *
            flatten<double>(
                dR_dze(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (dL_dze(h, u, v, nx, ny) + dabsL_dze(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) * dinvR_dze(h, u, v, nx, ny));
    }
}

void get_dAplus_dqx(const HybMatrix<double, SWE::n_variables>& q,
                    const HybMatrix<double, SWE::n_auxiliaries>& aux,
                    const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                    HybMatrix<double, SWE::n_variables * SWE::n_variables>& dAplus_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dAplus_dqx, gp) =
            0.5 *
            flatten<double>(
                dR_dqx(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (dL_dqx(h, u, v, nx, ny) + dabsL_dqx(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) * dinvR_dqx(h, u, v, nx, ny));
    }
}

void get_dAplus_dqy(const HybMatrix<double, SWE::n_variables>& q,
                    const HybMatrix<double, SWE::n_auxiliaries>& aux,
                    const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                    HybMatrix<double, SWE::n_variables * SWE::n_variables>& dAplus_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dAplus_dqy, gp) =
            0.5 *
            flatten<double>(
                dR_dqy(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (dL_dqy(h, u, v, nx, ny) + dabsL_dqy(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) + absL(h, u, v, nx, ny)) * dinvR_dqy(h, u, v, nx, ny));
    }
}

void get_Aminus(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                HybMatrix<double, SWE::n_variables * SWE::n_variables>& Aminus) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(Aminus, gp) = 0.5 * flatten<double>(R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) *
                                                   invR(h, u, v, nx, ny));
    }
}

void get_dAminus_dze(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     HybMatrix<double, SWE::n_variables * SWE::n_variables>& dAminus_dze) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dAminus_dze, gp) =
            0.5 *
            flatten<double>(
                dR_dze(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (dL_dze(h, u, v, nx, ny) - dabsL_dze(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) * dinvR_dze(h, u, v, nx, ny));
    }
}

void get_dAminus_dqx(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     HybMatrix<double, SWE::n_variables * SWE::n_variables>& dAminus_dqx) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dAminus_dqx, gp) =
            0.5 *
            flatten<double>(
                dR_dqx(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (dL_dqx(h, u, v, nx, ny) - dabsL_dqx(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) * dinvR_dqx(h, u, v, nx, ny));
    }
}

void get_dAminus_dqy(const HybMatrix<double, SWE::n_variables>& q,
                     const HybMatrix<double, SWE::n_auxiliaries>& aux,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     HybMatrix<double, SWE::n_variables * SWE::n_variables>& dAminus_dqy) {
    for (uint gp = 0; gp < columns(q); ++gp) {
        double h = aux(SWE::Auxiliaries::h, gp);
        double u = q(SWE::Variables::qx, gp) / aux(SWE::Auxiliaries::h, gp);
        double v = q(SWE::Variables::qy, gp) / aux(SWE::Auxiliaries::h, gp);

        double nx = surface_normal(GlobalCoord::x, gp);
        double ny = surface_normal(GlobalCoord::y, gp);

        column(dAminus_dqy, gp) =
            0.5 *
            flatten<double>(
                dR_dqy(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (dL_dqy(h, u, v, nx, ny) - dabsL_dqy(h, u, v, nx, ny)) * invR(h, u, v, nx, ny) +
                R(h, u, v, nx, ny) * (L(h, u, v, nx, ny) - absL(h, u, v, nx, ny)) * dinvR_dqy(h, u, v, nx, ny));
    }
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> L(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> L;

    double c  = std::sqrt(Global::g * h);
    double un = u * nx + v * ny;

    set_constant(L, 0.0);

    L(SWE::Variables::ze, SWE::Variables::ze) = un - c;
    L(SWE::Variables::qx, SWE::Variables::qx) = un;
    L(SWE::Variables::qy, SWE::Variables::qy) = un + c;

    return L;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dze(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dze;

    double dc_dze  = std::sqrt(Global::g / h) / 2.0;
    double dun_dze = -(u * nx + v * ny) / h;

    set_constant(dL_dze, 0.0);

    dL_dze(SWE::Variables::ze, SWE::Variables::ze) = dun_dze - dc_dze;
    dL_dze(SWE::Variables::qx, SWE::Variables::qx) = dun_dze;
    dL_dze(SWE::Variables::qy, SWE::Variables::qy) = dun_dze + dc_dze;

    return dL_dze;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqx(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqx;

    double dun_dqx = nx / h;

    set_constant(dL_dqx, 0.0);

    dL_dqx(SWE::Variables::ze, SWE::Variables::ze) = dun_dqx;
    dL_dqx(SWE::Variables::qx, SWE::Variables::qx) = dun_dqx;
    dL_dqx(SWE::Variables::qy, SWE::Variables::qy) = dun_dqx;

    return dL_dqx;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqy(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dL_dqy;

    double dun_dqy = ny / h;

    set_constant(dL_dqy, 0.0);

    dL_dqy(SWE::Variables::ze, SWE::Variables::ze) = dun_dqy;
    dL_dqy(SWE::Variables::qx, SWE::Variables::qx) = dun_dqy;
    dL_dqy(SWE::Variables::qy, SWE::Variables::qy) = dun_dqy;

    return dL_dqy;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> absL(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> absL;

    double c  = std::sqrt(Global::g * h);
    double un = u * nx + v * ny;

    set_constant(absL, 0.0);

    absL(SWE::Variables::ze, SWE::Variables::ze) = std::abs(un - c);
    absL(SWE::Variables::qx, SWE::Variables::qx) = std::abs(un);
    absL(SWE::Variables::qy, SWE::Variables::qy) = std::abs(un + c);

    return absL;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dze(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dze;

    double c  = std::sqrt(Global::g * h);
    double un = u * nx + v * ny;

    double dc_dze  = std::sqrt(Global::g / h) / 2.0;
    double dun_dze = -(u * nx + v * ny) / h;

    set_constant(dabsL_dze, 0.0);

    dabsL_dze(SWE::Variables::ze, SWE::Variables::ze) = (dun_dze - dc_dze) * std::pow(-1.0, std::signbit(un - c));
    dabsL_dze(SWE::Variables::qx, SWE::Variables::qx) = dun_dze * std::pow(-1.0, std::signbit(un));
    dabsL_dze(SWE::Variables::qy, SWE::Variables::qy) = (dun_dze + dc_dze) * std::pow(-1.0, std::signbit(un + c));

    return dabsL_dze;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqx(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqx;

    double c  = std::sqrt(Global::g * h);
    double un = u * nx + v * ny;

    double dun_dqx = nx / h;

    set_constant(dabsL_dqx, 0.0);

    dabsL_dqx(SWE::Variables::ze, SWE::Variables::ze) = dun_dqx * std::pow(-1.0, std::signbit(un - c));
    dabsL_dqx(SWE::Variables::qx, SWE::Variables::qx) = dun_dqx * std::pow(-1.0, std::signbit(un));
    dabsL_dqx(SWE::Variables::qy, SWE::Variables::qy) = dun_dqx * std::pow(-1.0, std::signbit(un + c));

    return dabsL_dqx;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqy(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dabsL_dqy;

    double c  = std::sqrt(Global::g * h);
    double un = u * nx + v * ny;

    double dun_dqy = ny / h;

    set_constant(dabsL_dqy, 0.0);

    dabsL_dqy(SWE::Variables::ze, SWE::Variables::ze) = dun_dqy * std::pow(-1.0, std::signbit(un - c));
    dabsL_dqy(SWE::Variables::qx, SWE::Variables::qx) = dun_dqy * std::pow(-1.0, std::signbit(un));
    dabsL_dqy(SWE::Variables::qy, SWE::Variables::qy) = dun_dqy * std::pow(-1.0, std::signbit(un + c));

    return dabsL_dqy;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> R(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> R;

    double c = std::sqrt(Global::g * h);

    R(SWE::Variables::ze, SWE::Variables::ze) = 1.0;
    R(SWE::Variables::ze, SWE::Variables::qx) = 0.0;
    R(SWE::Variables::ze, SWE::Variables::qy) = 1.0;

    R(SWE::Variables::qx, SWE::Variables::ze) = u - c * nx;
    R(SWE::Variables::qx, SWE::Variables::qx) = -ny;
    R(SWE::Variables::qx, SWE::Variables::qy) = u + c * nx;

    R(SWE::Variables::qy, SWE::Variables::ze) = v - c * ny;
    R(SWE::Variables::qy, SWE::Variables::qx) = nx;
    R(SWE::Variables::qy, SWE::Variables::qy) = v + c * ny;

    return R;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dze(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dze;

    double dc_dze = std::sqrt(Global::g / h) / 2.0;

    set_constant(dR_dze, 0.0);

    dR_dze(SWE::Variables::qx, SWE::Variables::ze) = -u / h - dc_dze * nx;
    dR_dze(SWE::Variables::qx, SWE::Variables::qy) = -u / h + dc_dze * nx;
    dR_dze(SWE::Variables::qy, SWE::Variables::ze) = -v / h - dc_dze * ny;
    dR_dze(SWE::Variables::qy, SWE::Variables::qy) = -v / h + dc_dze * ny;

    return dR_dze;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqx(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqx;

    set_constant(dR_dqx, 0.0);

    dR_dqx(SWE::Variables::qx, SWE::Variables::ze) = 1 / h;
    dR_dqx(SWE::Variables::qx, SWE::Variables::qy) = 1 / h;

    return dR_dqx;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqy(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dR_dqy;

    set_constant(dR_dqy, 0.0);

    dR_dqy(SWE::Variables::qy, SWE::Variables::ze) = 1 / h;
    dR_dqy(SWE::Variables::qy, SWE::Variables::qy) = 1 / h;

    return dR_dqy;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> invR(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> invR;

    double c  = std::sqrt(Global::g * h);
    double un = u * nx + v * ny;

    invR(SWE::Variables::ze, SWE::Variables::ze) = (c + un) / (2 * c);
    invR(SWE::Variables::ze, SWE::Variables::qx) = -nx / (2 * c);
    invR(SWE::Variables::ze, SWE::Variables::qy) = -ny / (2 * c);

    invR(SWE::Variables::qx, SWE::Variables::ze) = u * ny - v * nx;
    invR(SWE::Variables::qx, SWE::Variables::qx) = -ny;
    invR(SWE::Variables::qx, SWE::Variables::qy) = nx;

    invR(SWE::Variables::qy, SWE::Variables::ze) = (c - un) / (2 * c);
    invR(SWE::Variables::qy, SWE::Variables::qx) = nx / (2 * c);
    invR(SWE::Variables::qy, SWE::Variables::qy) = ny / (2 * c);

    return invR;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dze(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dze;

    double c    = std::sqrt(Global::g * h);
    double c_sq = Global::g * h;
    double un   = u * nx + v * ny;

    double dc_dze = std::sqrt(Global::g / h) / 2.0;

    set_constant(dinvR_dze, 0.0);

    dinvR_dze(SWE::Variables::ze, SWE::Variables::ze) = -3 * un / (4 * h * c);
    dinvR_dze(SWE::Variables::ze, SWE::Variables::qx) = dc_dze * nx / (2 * c_sq);
    dinvR_dze(SWE::Variables::ze, SWE::Variables::qy) = dc_dze * ny / (2 * c_sq);

    dinvR_dze(SWE::Variables::qx, SWE::Variables::ze) = -(u * ny - v * nx) / h;

    dinvR_dze(SWE::Variables::qy, SWE::Variables::ze) = 3 * un / (4 * h * c);
    dinvR_dze(SWE::Variables::qy, SWE::Variables::qx) = -dc_dze * nx / (2 * c_sq);
    dinvR_dze(SWE::Variables::qy, SWE::Variables::qy) = -dc_dze * ny / (2 * c_sq);

    return dinvR_dze;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqx(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqx;

    double c = std::sqrt(Global::g * h);

    set_constant(dinvR_dqx, 0.0);

    dinvR_dqx(SWE::Variables::ze, SWE::Variables::ze) = nx / h / (2 * c);
    dinvR_dqx(SWE::Variables::qx, SWE::Variables::ze) = ny / h;
    dinvR_dqx(SWE::Variables::qy, SWE::Variables::ze) = -nx / h / (2 * c);

    return dinvR_dqx;
}

StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqy(double h, double u, double v, double nx, double ny) {
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dinvR_dqy;

    double c = std::sqrt(Global::g * h);

    set_constant(dinvR_dqy, 0.0);

    dinvR_dqy(SWE::Variables::ze, SWE::Variables::ze) = ny / h / (2 * c);
    dinvR_dqy(SWE::Variables::qx, SWE::Variables::ze) = -nx / h;
    dinvR_dqy(SWE::Variables::qy, SWE::Variables::ze) = -ny / h / (2 * c);

    return dinvR_dqy;
}
}

#endif