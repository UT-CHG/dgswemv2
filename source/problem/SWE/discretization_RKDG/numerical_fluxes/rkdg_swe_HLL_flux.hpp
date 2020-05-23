#ifndef RKDG_SWE_HLL_FLUX_HPP
#define RKDG_SWE_HLL_FLUX_HPP

namespace SWE {
namespace RKDG {
// The normal points form the interior side (in) to the exterior side (ex)
void HLL_flux(const double gravity,
              const Column<HybMatrix<double, SWE::n_variables>>& q_in,
              const Column<HybMatrix<double, SWE::n_variables>>& q_ex,
              const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
              const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_ex,
              const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
              Column<HybMatrix<double, SWE::n_variables>>&& F_hat) {
    const double bath_in = aux_in[SWE::Auxiliaries::bath];
    const double h_in    = aux_in[SWE::Auxiliaries::h];
    const double u_in    = q_in[SWE::Variables::qx] / h_in;
    const double v_in    = q_in[SWE::Variables::qy] / h_in;

    const double bath_ex = aux_ex[SWE::Auxiliaries::bath];
    const double h_ex    = aux_ex[SWE::Auxiliaries::h];
    const double u_ex    = q_ex[SWE::Variables::qx] / h_ex;
    const double v_ex    = q_ex[SWE::Variables::qy] / h_ex;

    const double un_in = u_in * surface_normal[GlobalCoord::x] + v_in * surface_normal[GlobalCoord::y];
    const double un_ex = u_ex * surface_normal[GlobalCoord::x] + v_ex * surface_normal[GlobalCoord::y];

    const double s_in = std::min(un_in - std::sqrt(gravity * h_in), un_ex - std::sqrt(gravity * h_ex));
    const double s_ex = std::max(un_in + std::sqrt(gravity * h_in), un_ex + std::sqrt(gravity * h_ex));

    StatVector<double, SWE::n_variables> Fn_in;
    StatVector<double, SWE::n_variables> Fn_ex;

    const double nx = surface_normal[GlobalCoord::x];
    const double ny = surface_normal[GlobalCoord::y];

    // compute internal flux matrix
    const double uuh_in = u_in * q_in[SWE::Variables::qx];
    const double vvh_in = v_in * q_in[SWE::Variables::qy];
    const double uvh_in = u_in * q_in[SWE::Variables::qy];
    const double pe_in  = gravity * (std::pow(q_in[SWE::Variables::ze], 2) / 2 + q_in[SWE::Variables::ze] * bath_in);
    const double hcu_in = u_in * q_in[SWE::Variables::hc];
    const double hcv_in = v_in * q_in[SWE::Variables::hc];

    Fn_in[SWE::Variables::ze] = q_in[SWE::Variables::qx] * nx + q_in[SWE::Variables::qy] * ny;
    Fn_in[SWE::Variables::qx] = (uuh_in + pe_in) * nx + uvh_in * ny;
    Fn_in[SWE::Variables::qy] = uvh_in * nx + (vvh_in + pe_in) * ny;
    Fn_in[SWE::Variables::hc] = hcu_in * nx + hcv_in * ny;

    // compute external flux matrix
    const double uuh_ex = u_ex * q_ex[SWE::Variables::qx];
    const double vvh_ex = v_ex * q_ex[SWE::Variables::qy];
    const double uvh_ex = u_ex * q_ex[SWE::Variables::qy];
    const double pe_ex  = gravity * (std::pow(q_ex[SWE::Variables::ze], 2) / 2 + q_ex[SWE::Variables::ze] * bath_ex);
    const double hcu_ex = u_ex * q_ex[SWE::Variables::hc];
    const double hcv_ex = v_ex * q_ex[SWE::Variables::hc];

    Fn_ex[SWE::Variables::ze] = q_ex[SWE::Variables::qx] * nx + q_ex[SWE::Variables::qy] * ny;
    Fn_ex[SWE::Variables::qx] = (uuh_ex + pe_ex) * nx + uvh_ex * ny;
    Fn_ex[SWE::Variables::qy] = uvh_ex * nx + (vvh_ex + pe_ex) * ny;
    Fn_ex[SWE::Variables::hc] = hcu_ex * nx + hcv_ex * ny;

    if (s_in >= 0) {
        F_hat = Fn_in;
    } else if (s_ex <= 0) {
        F_hat = Fn_ex;
    } else if (s_in < 0 && s_ex > 0) {
        F_hat = ((s_ex * Fn_in - s_in * Fn_ex) - s_in * s_ex * (q_in - q_ex)) / (s_ex - s_in);
    }
}
}
}

#endif