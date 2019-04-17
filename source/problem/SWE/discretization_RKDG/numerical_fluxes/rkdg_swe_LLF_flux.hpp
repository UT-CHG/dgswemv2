#ifndef RKDG_SWE_LLF_FLUX_HPP
#define RKDG_SWE_LLF_FLUX_HPP

namespace SWE {
namespace RKDG {
// The normal points form the interior side (in) to the exterior side (ex)
template <typename InputViewType, typename InputArrayType>
inline void LLF_flux(const double gravity,
                     const InputViewType& ze_in,
                     const InputViewType& qx_in,
                     const InputViewType& qy_in,
                     const InputViewType& ze_ex,
                     const InputViewType& qx_ex,
                     const InputViewType& qy_ex,
                     const std::array<InputViewType, SWE::n_auxiliaries>& aux,
                     const std::array<InputArrayType, SWE::n_dimensions>& surface_normal,
                     std::array<typename Result<InputArrayType>::type, SWE::n_variables>& Flux) {

    using result_t = typename Result<InputArrayType>::type;

    auto& bath = aux[SWE::Auxiliaries::bath];
    auto& sp   = aux[SWE::Auxiliaries::sp];

    auto h_in = ze_in + bath;
    result_t u_in = blaze::evaluate(qx_in / h_in);
    auto v_in = qy_in / h_in;

    auto h_ex = ze_ex + bath;
    result_t u_ex = blaze::evaluate(qx_ex / h_ex);
    auto v_ex = qy_ex / h_ex;

    auto un_in = u_in * surface_normal[GlobalCoord::x] + v_in * surface_normal[GlobalCoord::y];
    auto un_ex = u_ex * surface_normal[GlobalCoord::x] + v_ex * surface_normal[GlobalCoord::y];

    auto& nx = surface_normal[GlobalCoord::x];
    auto& ny = surface_normal[GlobalCoord::y];

    auto sp_correction = (nx * sp ) * (nx * sp) + (ny * ny);

    result_t max_eigenvalue = blaze::evaluate(max_vec(abs_vec(un_in) + sqrt_vec(gravity * h_in * sp_correction),
                                                      abs_vec(un_ex) + sqrt_vec(gravity * h_ex * sp_correction)));

    // compute internal flux matrix
    auto uuh_in = u_in * qx_in;
    auto vvh_in = v_in * qy_in;
    auto uvh_in = u_in * qy_in;
    auto pe_in  = gravity * ze_in *( ze_in / 2 + bath);

    auto Fn_ze_in = sp * qx_in * nx + qy_in * ny;
    auto Fn_qx_in = sp * (uuh_in + pe_in) * nx + uvh_in * ny;
    auto Fn_qy_in = sp * uvh_in * nx + (vvh_in + pe_in) * ny;

    // compute external flux matrix
    auto uuh_ex = u_ex * qx_ex;
    auto vvh_ex = v_ex * qy_ex;
    auto uvh_ex = u_ex * qy_ex;
    auto pe_ex  = gravity * ze_ex * (ze_ex / 2 + bath);

    auto Fn_ze_ex = qx_ex * nx + qy_ex * ny;
    auto Fn_qx_ex = (uuh_ex + pe_ex) * nx + uvh_ex * ny;
    auto Fn_qy_ex = uvh_ex * nx + (vvh_ex + pe_ex) * ny;

    Flux[0] = 0.5 * (Fn_ze_in + Fn_ze_ex + max_eigenvalue * (ze_in - ze_ex));
    Flux[1] = 0.5 * (Fn_qx_in + Fn_qx_ex + max_eigenvalue * (qx_in - qx_ex));
    Flux[2] = 0.5 * (Fn_qy_in + Fn_qy_ex + max_eigenvalue * (qy_in - qy_ex));
}

inline void LLF_flux(const double gravity,
                     const double ze_in,
                     const double qx_in,
                     const double qy_in,
                     const double ze_ex,
                     const double qx_ex,
                     const double qy_ex,
                     const std::array<double, SWE::n_auxiliaries>& aux,
                     const std::array<double, SWE::n_dimensions>& surface_normal,
                     std::array<double, SWE::n_variables>& Flux) {
    const double& bath = aux[SWE::Auxiliaries::bath];
    const double& sp   = aux[SWE::Auxiliaries::sp];

    double h_in = ze_in + bath;
    double u_in = qx_in / h_in;
    double v_in = qy_in / h_in;

    double h_ex = ze_ex + bath;
    double u_ex = qx_ex / h_ex;
    double v_ex = qy_ex / h_ex;

    double un_in = u_in * surface_normal[GlobalCoord::x] + v_in * surface_normal[GlobalCoord::y];
    double un_ex = u_ex * surface_normal[GlobalCoord::x] + v_ex * surface_normal[GlobalCoord::y];

    double sp_correction =
        std::pow(surface_normal[GlobalCoord::x] * sp, 2) + std::pow(surface_normal[GlobalCoord::y], 2);

    double max_eigenvalue = std::max(std::abs(un_in) + std::sqrt(gravity * h_in * sp_correction),
                                     std::abs(un_ex) + std::sqrt(gravity * h_ex * sp_correction));

    const double& nx = surface_normal[GlobalCoord::x];
    const double& ny = surface_normal[GlobalCoord::y];

    // compute internal flux matrix
    double uuh_in = u_in * qx_in;
    double vvh_in = v_in * qy_in;
    double uvh_in = u_in * qy_in;
    double pe_in  = gravity * ze_in *( ze_in / 2 + bath);

    double Fn_ze_in = sp * qx_in * nx + qy_in * ny;
    double Fn_qx_in = sp * (uuh_in + pe_in) * nx + uvh_in * ny;
    double Fn_qy_in = sp * uvh_in * nx + (vvh_in + pe_in) * ny;

    // compute external flux matrix
    double uuh_ex = u_ex * qx_ex;
    double vvh_ex = v_ex * qy_ex;
    double uvh_ex = u_ex * qy_ex;
    double pe_ex  = gravity * (pow_vec(ze_ex, 2) / 2 + ze_ex * bath);

    double Fn_ze_ex = qx_ex * nx + qy_ex * ny;
    double Fn_qx_ex = (uuh_ex + pe_ex) * nx + uvh_ex * ny;
    double Fn_qy_ex = uvh_ex * nx + (vvh_ex + pe_ex) * ny;

    Flux[0] = 0.5 * (Fn_ze_in + Fn_ze_ex + max_eigenvalue * (ze_in - ze_ex));
    Flux[1] = 0.5 * (Fn_qx_in + Fn_qx_ex + max_eigenvalue * (qx_in - qx_ex));
    Flux[2] = 0.5 * (Fn_qy_in + Fn_qy_ex + max_eigenvalue * (qy_in - qy_ex));
}

}
}

#ifdef USE_BLAZE
namespace blaze {

namespace workaround {
template <typename T>
BLAZE_ALWAYS_INLINE const SIMDdouble abs( const SIMDf64<T>& a ) noexcept {
  const SIMDdouble minus_one = blaze::set(-1.);
  return max(a, minus_one*a);
}
}

template <typename Datapar>
BLAZE_STRONG_INLINE void LLF_flux(const Datapar& gravity,
                                 const Datapar& ze_in,
                                 const Datapar& qx_in,
                                 const Datapar& qy_in,
                                 const Datapar& ze_ex,
                                 const Datapar& qx_ex,
                                 const Datapar& qy_ex,
                                 const Datapar& bath,
                                 const Datapar& sp,
                                 const Datapar& nx,
                                 const Datapar& ny,
                                 double* flux_ze_addr,
                                 double* flux_qx_addr,
                                 double* flux_qy_addr) {
    assert(flux_ze_addr);
    assert(flux_qx_addr);
    assert(flux_qy_addr);

    const Datapar half = blaze::set(0.5);

    Datapar h_in = ze_in + bath;
    Datapar u_in = qx_in / h_in;
    Datapar v_in = qy_in / h_in;

    Datapar h_ex = ze_ex + bath;
    Datapar u_ex = qx_ex / h_ex;
    Datapar v_ex = qy_ex / h_ex;

    Datapar un_in = u_in * nx + v_in * ny;
    Datapar un_ex = u_ex * nx + v_ex * ny;

    Datapar nx_sp = nx * sp;
    Datapar sp_correction = nx_sp * nx_sp + ny * ny;

    Datapar max_eigenvalue = max(workaround::abs(un_in) + sqrt(gravity * h_in * sp_correction),
                                 workaround::abs(un_ex) + sqrt(gravity * h_ex * sp_correction));

    // compute internal flux matrix
    Datapar uuh_in = u_in * qx_in;
    Datapar vvh_in = v_in * qy_in;
    Datapar uvh_in = u_in * qy_in;
    Datapar pe_in  = gravity * ze_in * ( half * ze_in + bath); //fma

    Datapar Fn_ze_in = nx_sp * qx_in + qy_in * ny;
    Datapar Fn_qx_in = nx_sp * (uuh_in + pe_in) + uvh_in * ny; //fma
    Datapar Fn_qy_in = nx_sp * uvh_in  + (vvh_in + pe_in) * ny; //fma

    // compute external flux matrix
    Datapar uuh_ex = u_ex * qx_ex;
    Datapar vvh_ex = v_ex * qy_ex;
    Datapar uvh_ex = u_ex * qy_ex;
    Datapar pe_ex  = gravity * ze_ex * (ze_ex * half + bath); //fma

    Datapar Fn_ze_ex = qx_ex * nx + qy_ex * ny;
    Datapar Fn_qx_ex = (uuh_ex + pe_ex) * nx + uvh_ex * ny; //fma
    Datapar Fn_qy_ex = uvh_ex * nx + (vvh_ex + pe_ex) * ny; //fma

    storea(flux_ze_addr, half * (Fn_ze_in + Fn_ze_ex + max_eigenvalue * (ze_in - ze_ex))); //fma
    storea(flux_qx_addr, half * (Fn_qx_in + Fn_qx_ex + max_eigenvalue * (qx_in - qx_ex))); //fma
    storea(flux_qy_addr, half * (Fn_qy_in + Fn_qy_ex + max_eigenvalue * (qy_in - qy_ex))); //fma
}
}

#endif

#endif