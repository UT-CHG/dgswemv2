#ifndef RKDG_SWE_LLF_FLUX_HPP
#define RKDG_SWE_LLF_FLUX_HPP

namespace SWE {
namespace RKDG {
// The normal points form the interior side (in) to the exterior side (ex)
template <typename InputArrayType>
inline void LLF_flux(const double gravity,
                     const InputArrayType& ze_in,
                     const InputArrayType& qx_in,
                     const InputArrayType& qy_in,
                     const InputArrayType& ze_ex,
                     const InputArrayType& qx_ex,
                     const InputArrayType& qy_ex,
                     const std::array<InputArrayType, SWE::n_auxiliaries>& aux,
                     const std::array<InputArrayType, SWE::n_dimensions>& surface_normal,
                     std::array<typename Result<InputArrayType>::type, SWE::n_variables>& Flux) {
    auto& bath = aux[SWE::Auxiliaries::bath];
    auto& sp   = aux[SWE::Auxiliaries::sp];

    auto h_in = ze_in + bath;
    auto u_in = qx_in / h_in;
    auto v_in = qy_in / h_in;

    auto h_ex = ze_ex + bath;
    auto u_ex = qx_ex / h_ex;
    auto v_ex = qy_ex / h_ex;

    auto un_in = u_in * surface_normal[GlobalCoord::x] + v_in * surface_normal[GlobalCoord::y];
    auto un_ex = u_ex * surface_normal[GlobalCoord::x] + v_ex * surface_normal[GlobalCoord::y];

    auto sp_correction =
        pow_vec(surface_normal[GlobalCoord::x] * sp, 2) + pow_vec(surface_normal[GlobalCoord::y], 2);

    auto max_eigenvalue = max_vec(abs_vec(un_in) + sqrt_vec(gravity * h_in * sp_correction),
                                  abs_vec(un_ex) + sqrt_vec(gravity * h_ex * sp_correction));

    auto& nx = surface_normal[GlobalCoord::x];
    auto& ny = surface_normal[GlobalCoord::y];

    // compute internal flux matrix
    auto uuh_in = u_in * qx_in;
    auto vvh_in = v_in * qy_in;
    auto uvh_in = u_in * qy_in;
    auto pe_in  = gravity * (ze_in * ze_in / 2 + ze_in * bath);

    auto Fn_ze_in = sp * qx_in * nx + qy_in * ny;
    auto Fn_qx_in = sp * (uuh_in + pe_in) * nx + uvh_in * ny;
    auto Fn_qy_in = sp * uvh_in * nx + (vvh_in + pe_in) * ny;

    // compute external flux matrix
    auto uuh_ex = u_ex * qx_ex;
    auto vvh_ex = v_ex * qy_ex;
    auto uvh_ex = u_ex * qy_ex;
    auto pe_ex  = gravity * (pow_vec(ze_ex, 2) / 2 + ze_ex * bath);

    auto Fn_ze_ex = qx_ex * nx + qy_ex * ny;
    auto Fn_qx_ex = (uuh_ex + pe_ex) * nx + uvh_ex * ny;
    auto Fn_qy_ex = uvh_ex * nx + (vvh_ex + pe_ex) * ny;

    Flux[0] = 0.5 * (Fn_ze_in + Fn_ze_ex + max_eigenvalue * (ze_in - ze_ex));
    Flux[1] = 0.5 * (Fn_qx_in + Fn_qx_ex + max_eigenvalue * (qx_in - qx_ex));
    Flux[2] = 0.5 * (Fn_qy_in + Fn_qy_ex + max_eigenvalue * (qy_in - qy_ex));
}
}
}

#endif