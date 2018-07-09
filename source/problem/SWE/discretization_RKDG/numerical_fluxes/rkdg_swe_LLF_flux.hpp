#ifndef RKDG_SWE_LLF_FLUX_HPP
#define RKDG_SWE_LLF_FLUX_HPP

#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
namespace RKDG {
// The normal points form the interior side (in) to the exterior side (ex)
inline void LLF_flux(const double gravity,
                     const Vector<double, SWE::n_variables>& q_in,
                     const Vector<double, SWE::n_variables>& q_ex,
                     const Vector<double, SWE::n_auxiliaries>& aux,
                     const std::vector<double>& surface_normal,
                     Vector<double, SWE::n_variables>& F_hat) {
    double bath = aux[SWE::Auxiliaries::bath];
    double sp   = aux[SWE::Auxiliaries::sp];

    double h_in = q_in[SWE::Variables::ze] + bath;
    double u_in = q_in[SWE::Variables::qx] / h_in;
    double v_in = q_in[SWE::Variables::qy] / h_in;

    double h_ex = q_ex[SWE::Variables::ze] + bath;
    double u_ex = q_ex[SWE::Variables::qx] / h_ex;
    double v_ex = q_ex[SWE::Variables::qy] / h_ex;

    double un_in = u_in * surface_normal[GlobalCoord::x] + v_in * surface_normal[GlobalCoord::y];
    double un_ex = u_ex * surface_normal[GlobalCoord::x] + v_ex * surface_normal[GlobalCoord::y];

    double sp_correction =
        std::pow(surface_normal[GlobalCoord::x] * sp, 2) + std::pow(surface_normal[GlobalCoord::y], 2);

    double max_eigenvalue = std::max(std::abs(un_in) + std::sqrt(gravity * h_in * sp_correction),
                                     std::abs(un_ex) + std::sqrt(gravity * h_ex * sp_correction));

    // compute internal flux matrix
    double uuh_in = u_in * q_in[SWE::Variables::qx];
    double vvh_in = v_in * q_in[SWE::Variables::qy];
    double uvh_in = u_in * q_in[SWE::Variables::qy];
    double pe_in  = gravity * (std::pow(q_in[SWE::Variables::ze], 2) / 2 + q_in[SWE::Variables::ze] * bath);

    double ze_flux_x_in = sp * q_in[SWE::Variables::qx];
    double ze_flux_y_in = q_in[SWE::Variables::qy];

    double qx_flux_x_in = sp * (uuh_in + pe_in);
    double qx_flux_y_in = uvh_in;

    double qy_flux_x_in = sp * uvh_in;
    double qy_flux_y_in = vvh_in + pe_in;

    // compute external flux matrix
    double uuh_ex = u_ex * q_ex[SWE::Variables::qx];
    double vvh_ex = v_ex * q_ex[SWE::Variables::qy];
    double uvh_ex = u_ex * q_ex[SWE::Variables::qy];
    double pe_ex  = gravity * (std::pow(q_ex[SWE::Variables::ze], 2) / 2 + q_ex[SWE::Variables::ze] * bath);

    double ze_flux_x_ex = sp * q_ex[SWE::Variables::qx];
    double ze_flux_y_ex = q_ex[SWE::Variables::qy];

    double qx_flux_x_ex = sp * (uuh_ex + pe_ex);
    double qx_flux_y_ex = uvh_ex;

    double qy_flux_x_ex = sp * uvh_ex;
    double qy_flux_y_ex = vvh_ex + pe_ex;

    // compute numerical flux
    double ze_flux_avg = (ze_flux_x_in + ze_flux_x_ex) * surface_normal[GlobalCoord::x] +
                         (ze_flux_y_in + ze_flux_y_ex) * surface_normal[GlobalCoord::y];
    double qx_flux_avg = (qx_flux_x_in + qx_flux_x_ex) * surface_normal[GlobalCoord::x] +
                         (qx_flux_y_in + qx_flux_y_ex) * surface_normal[GlobalCoord::y];
    double qy_flux_avg = (qy_flux_x_in + qy_flux_x_ex) * surface_normal[GlobalCoord::x] +
                         (qy_flux_y_in + qy_flux_y_ex) * surface_normal[GlobalCoord::y];

    F_hat[SWE::Variables::ze] =
        0.5 * (ze_flux_avg + max_eigenvalue * (q_in[SWE::Variables::ze] - q_ex[SWE::Variables::ze]));
    F_hat[SWE::Variables::qx] =
        0.5 * (qx_flux_avg + max_eigenvalue * (q_in[SWE::Variables::qx] - q_ex[SWE::Variables::qx]));
    F_hat[SWE::Variables::qy] =
        0.5 * (qy_flux_avg + max_eigenvalue * (q_in[SWE::Variables::qy] - q_ex[SWE::Variables::qy]));
}
}
}

#endif