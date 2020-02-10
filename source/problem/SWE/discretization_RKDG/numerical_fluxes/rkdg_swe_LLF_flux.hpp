#ifndef RKDG_SWE_LLF_FLUX_HPP
#define RKDG_SWE_LLF_FLUX_HPP

namespace SWE {
namespace RKDG {
// The normal points form the interior side (in) to the exterior side (ex)
void LLF_flux(const double gravity,
              const Column<HybMatrix<double, SWE::n_variables>>& q_in,
              const Column<HybMatrix<double, SWE::n_variables>>& q_ex,
              const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
              const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_ex,
              const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
              Column<HybMatrix<double, SWE::n_variables>>&& F_hat) {
    double sp      = aux_in[SWE::Auxiliaries::sp];
    double bath_in = aux_in[SWE::Auxiliaries::bath];
    double bath_ex = aux_ex[SWE::Auxiliaries::bath];

    double h_in = q_in[SWE::Variables::ze] + bath_in;
    double u_in = q_in[SWE::Variables::qx] / h_in;
    double v_in = q_in[SWE::Variables::qy] / h_in;

    double h_ex = q_ex[SWE::Variables::ze] + bath_ex;
    double u_ex = q_ex[SWE::Variables::qx] / h_ex;
    double v_ex = q_ex[SWE::Variables::qy] / h_ex;

    double un_in = u_in * surface_normal[GlobalCoord::x] + v_in * surface_normal[GlobalCoord::y];
    double un_ex = u_ex * surface_normal[GlobalCoord::x] + v_ex * surface_normal[GlobalCoord::y];

    double sp_correction =
        std::pow(surface_normal[GlobalCoord::x] * sp, 2) + std::pow(surface_normal[GlobalCoord::y], 2);

    double max_eigenvalue = std::max(std::abs(un_in) + std::sqrt(gravity * h_in * sp_correction),
                                     std::abs(un_ex) + std::sqrt(gravity * h_ex * sp_correction));

    StatVector<double, SWE::n_variables> Fn_in;
    StatVector<double, SWE::n_variables> Fn_ex;

    double nx = surface_normal[GlobalCoord::x];
    double ny = surface_normal[GlobalCoord::y];

    // compute internal flux matrix
    double uuh_in = u_in * q_in[SWE::Variables::qx];
    double vvh_in = v_in * q_in[SWE::Variables::qy];
    double uvh_in = u_in * q_in[SWE::Variables::qy];
    double pe_in  = gravity * (std::pow(q_in[SWE::Variables::ze], 2) / 2 + q_in[SWE::Variables::ze] * bath_in);

    Fn_in[SWE::Variables::ze] = sp * q_in[SWE::Variables::qx] * nx + q_in[SWE::Variables::qy] * ny;
    Fn_in[SWE::Variables::qx] = sp * (uuh_in + pe_in) * nx + uvh_in * ny;
    Fn_in[SWE::Variables::qy] = sp * uvh_in * nx + (vvh_in + pe_in) * ny;
    Fn_in[SWE::Variables::hc] = 0.0;  // TODO

    // compute external flux matrix
    double uuh_ex = u_ex * q_ex[SWE::Variables::qx];
    double vvh_ex = v_ex * q_ex[SWE::Variables::qy];
    double uvh_ex = u_ex * q_ex[SWE::Variables::qy];
    double pe_ex  = gravity * (std::pow(q_ex[SWE::Variables::ze], 2) / 2 + q_ex[SWE::Variables::ze] * bath_ex);

    Fn_ex[SWE::Variables::ze] = q_ex[SWE::Variables::qx] * nx + q_ex[SWE::Variables::qy] * ny;
    Fn_ex[SWE::Variables::qx] = (uuh_ex + pe_ex) * nx + uvh_ex * ny;
    Fn_ex[SWE::Variables::qy] = uvh_ex * nx + (vvh_ex + pe_ex) * ny;
    Fn_ex[SWE::Variables::hc] = 0.0;  // TODO

    F_hat = 0.5 * (Fn_in + Fn_ex + max_eigenvalue * (q_in - q_ex));
}
}
}

#endif