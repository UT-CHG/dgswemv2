#ifndef SWE_LLF_FLUX_HPP
#define SWE_LLF_FLUX_HPP

#include "../swe_definitions.hpp"

namespace SWE {
// The normal points form the interior side (in) to the exterior side (ex)
inline void LLF_flux(const double gravity,
                     const double ze_in,
                     const double ze_ex,
                     const double qx_in,
                     const double qx_ex,
                     const double qy_in,
                     const double qy_ex,
                     const double bath,
                     const double sp,
                     const std::vector<double>& normal,
                     double& ze_flux,
                     double& qx_flux,
                     double& qy_flux) {
    double h_in = ze_in + bath;
    double u_in = qx_in / h_in;
    double v_in = qy_in / h_in;

    double h_ex = ze_ex + bath;
    double u_ex = qx_ex / h_ex;
    double v_ex = qy_ex / h_ex;

    double un_in = u_in * normal[GlobalCoord::x] + v_in * normal[GlobalCoord::y];
    double un_ex = u_ex * normal[GlobalCoord::x] + v_ex * normal[GlobalCoord::y];

    double sp_correction = std::pow(normal[GlobalCoord::x] * sp, 2) + std::pow(normal[GlobalCoord::y], 2);

    double max_eigenvalue = std::max(std::abs(un_in) + std::sqrt(gravity * h_in * sp_correction),
                                     std::abs(un_ex) + std::sqrt(gravity * h_ex * sp_correction));

    // compute internal flux matrix
    double uuh_in = u_in * qx_in;
    double vvh_in = v_in * qy_in;
    double uvh_in = u_in * qy_in;
    double pe_in = gravity * (std::pow(ze_in, 2) / 2 + ze_in * bath);

    double ze_flux_x_in = sp * qx_in;
    double ze_flux_y_in = qy_in;

    double qx_flux_x_in = sp * (uuh_in + pe_in);
    double qx_flux_y_in = uvh_in;

    double qy_flux_x_in = sp * uvh_in;
    double qy_flux_y_in = vvh_in + pe_in;

    // compute external flux matrix
    double uuh_ex = u_ex * qx_ex;
    double vvh_ex = v_ex * qy_ex;
    double uvh_ex = u_ex * qy_ex;
    double pe_ex = gravity * (std::pow(ze_ex, 2) / 2 + ze_ex * bath);

    double ze_flux_x_ex = sp * qx_ex;
    double ze_flux_y_ex = qy_ex;

    double qx_flux_x_ex = sp * (uuh_ex + pe_ex);
    double qx_flux_y_ex = uvh_ex;

    double qy_flux_x_ex = sp * uvh_ex;
    double qy_flux_y_ex = vvh_ex + pe_ex;

    // compute numerical flux
    double ze_flux_avg =
        (ze_flux_x_in + ze_flux_x_ex) * normal[GlobalCoord::x] + (ze_flux_y_in + ze_flux_y_ex) * normal[GlobalCoord::y];
    double qx_flux_avg =
        (qx_flux_x_in + qx_flux_x_ex) * normal[GlobalCoord::x] + (qx_flux_y_in + qx_flux_y_ex) * normal[GlobalCoord::y];
    double qy_flux_avg =
        (qy_flux_x_in + qy_flux_x_ex) * normal[GlobalCoord::x] + (qy_flux_y_in + qy_flux_y_ex) * normal[GlobalCoord::y];

    ze_flux = 0.5 * (ze_flux_avg + max_eigenvalue * (ze_in - ze_ex));
    qx_flux = 0.5 * (qx_flux_avg + max_eigenvalue * (qx_in - qx_ex));
    qy_flux = 0.5 * (qy_flux_avg + max_eigenvalue * (qy_in - qy_ex));
}
}

#endif