#ifndef RKDG_SWE_IS_INTERNAL_HPP
#define RKDG_SWE_IS_INTERNAL_HPP

#include "utilities/is_soa.hpp"

namespace SWE {
namespace RKDG {
namespace ISP {
class Internal {
  private:
    BC::Land land_boundary;

  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {
        land_boundary.Initialize(1);
    }

    template <typename InterfaceType>
    typename std::enable_if<!Utilities::is_SoA<InterfaceType>::value>::type ComputeFlux(InterfaceType& intface) {
        //bool wet_in = intface.data_in.wet_dry_state.wet;
        //bool wet_ex = intface.data_ex.wet_dry_state.wet;

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        // assemble numerical fluxes
        LLF_flux(Global::g,
                 boundary_in.q_at_gp[SWE::Variables::ze],
                 boundary_in.q_at_gp[SWE::Variables::qx],
                 boundary_in.q_at_gp[SWE::Variables::qy],
                 boundary_ex.q_at_gp[SWE::Variables::ze],
                 boundary_ex.q_at_gp[SWE::Variables::qx],
                 boundary_ex.q_at_gp[SWE::Variables::qy],
                 boundary_in.aux_at_gp,
                 intface.surface_normal_in,
                 boundary_in.F_hat_at_gp
            );

        boundary_ex.F_hat_at_gp[SWE::Variables::ze] = -boundary_in.F_hat_at_gp[SWE::Variables::ze];
        boundary_ex.F_hat_at_gp[SWE::Variables::qx] = -boundary_in.F_hat_at_gp[SWE::Variables::qx];
        boundary_ex.F_hat_at_gp[SWE::Variables::qy] = -boundary_in.F_hat_at_gp[SWE::Variables::qy];

        /*if (boundary_in.F_hat_at_gp(Variables::ze, gp) > 1e-12) {
            if (!wet_in) {  // water flowing from dry IN element
                // Zero flux on IN element side
                set_constant(column(boundary_in.F_hat_at_gp, gp), 0.0);

                // Reflective Boundary on EX element side
                this->land_boundary.ComputeFlux(column(intface.surface_normal_ex, gp),
                                                column(boundary_ex.q_at_gp, gp),
                                                column(boundary_ex.aux_at_gp, gp),
                                                column(boundary_ex.F_hat_at_gp, gp));

            } else if (!wet_ex) {  // water flowing to dry EX element
                LLF_flux(0.0,
                         column(boundary_ex.q_at_gp, gp),
                         column(boundary_in.q_at_gp, gp),
                         column(boundary_ex.aux_at_gp, gp),
                         column(intface.surface_normal_ex, gp),
                         column(boundary_ex.F_hat_at_gp, gp));

                // Only remove gravity contributions for the momentum fluxes
                boundary_ex.F_hat_at_gp(Variables::ze, gp) = -boundary_in.F_hat_at_gp(Variables::ze, gp);
            }
        } else if (boundary_in.F_hat_at_gp(Variables::ze, gp) < -1e-12) {
            if (!wet_ex) {  // water flowing from dry EX element
                // Zero flux on EX element side
                set_constant(column(boundary_ex.F_hat_at_gp, gp), 0.0);

                // Reflective Boundary on IN element side
                this->land_boundary.ComputeFlux(column(intface.surface_normal_in, gp),
                                                column(boundary_in.q_at_gp, gp),
                                                column(boundary_in.aux_at_gp, gp),
                                                column(boundary_in.F_hat_at_gp, gp));

            } else if (!wet_in) {  // water flowing to dry IN element
                LLF_flux(0.0,
                         column(boundary_in.q_at_gp, gp),
                         column(boundary_ex.q_at_gp, gp),
                         column(boundary_in.aux_at_gp, gp),
                         column(intface.surface_normal_in, gp),
                         column(boundary_in.F_hat_at_gp, gp));

                boundary_in.F_hat_at_gp(Variables::ze, gp) = -boundary_ex.F_hat_at_gp(Variables::ze, gp);
            }
            }*/

//        assert(!std::isnan(boundary_in.F_hat_at_gp[Variables::ze][gp]));
    //   }
    }

    template <typename SoAType>
    typename std::enable_if<Utilities::is_SoA<SoAType>::value>::type ComputeFlux(SoAType& soa) {
        //matrix typedefs
        using matrix_t = typename std::decay<decltype(soa.data.q_in_at_gp[0])>::type;
        using simd_t = typename matrix_t::SIMDType;
        const size_t SIMDSIZE = simd_t::size;
        const simd_t g_vec = blaze::set(Global::g);
        //check all matrices are the same size
        static_assert(matrix_t::simdEnabled && blaze::usePadding);

        size_t rows_ = rows(soa.data.q_in_at_gp[0]);
        size_t columns_ = columns(soa.data.q_in_at_gp[0]);

	//Note that this is not optional. Converting to a row major matrix requires inverting the loops
	constexpr auto SO = SO::ColumnMajor;

        const DynMatrix<double, SO>& ze_in = soa.data.q_in_at_gp[SWE::Variables::ze];
        const DynMatrix<double, SO>& qx_in = soa.data.q_in_at_gp[SWE::Variables::qx];
        const DynMatrix<double, SO>& qy_in = soa.data.q_in_at_gp[SWE::Variables::qy];
        const DynMatrix<double, SO>& ze_ex = soa.data.q_ex_at_gp[SWE::Variables::ze];
        const DynMatrix<double, SO>& qx_ex = soa.data.q_ex_at_gp[SWE::Variables::qx];
        const DynMatrix<double, SO>& qy_ex = soa.data.q_ex_at_gp[SWE::Variables::qy];
        const DynMatrix<double, SO>& bath  = soa.data.aux_at_gp[SWE::Auxiliaries::bath];
        const DynMatrix<double, SO>& sp    = soa.data.aux_at_gp[SWE::Auxiliaries::sp];
        const DynMatrix<double, SO>& nx    = soa.surface_normal[GlobalCoord::x];
        const DynMatrix<double, SO>& ny    = soa.surface_normal[GlobalCoord::y];
        DynMatrix<double, SO>& flux_ze     = soa.data.F_hat_at_gp[SWE::Variables::ze];
        DynMatrix<double, SO>& flux_qx     = soa.data.F_hat_at_gp[SWE::Variables::qx];
        DynMatrix<double, SO>& flux_qy     = soa.data.F_hat_at_gp[SWE::Variables::qy];

        for ( size_t j=0UL; j < columns_; ++j ) {

            for ( size_t i = 0UL; i < rows_; i += SIMDSIZE) {

                blaze::LLF_flux(g_vec,
                                ze_in.load(i,j),
                                qx_in.load(i,j),
                                qy_in.load(i,j),
                                ze_ex.load(i,j),
                                qx_ex.load(i,j),
                                qy_ex.load(i,j),
                                bath.load(i,j),
                                sp.load(i,j),
                                nx.load(i,j),
                                ny.load(i,j),
                                flux_ze.data(j) + i,
                                flux_qx.data(j) + i,
                                flux_qy.data(j) + i);
            }
        }
    }
};
}
}
}

#endif