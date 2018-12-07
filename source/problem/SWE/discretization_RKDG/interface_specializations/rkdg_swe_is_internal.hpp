#ifndef RKDG_SWE_IS_INTERNAL_HPP
#define RKDG_SWE_IS_INTERNAL_HPP

namespace SWE {
namespace RKDG {
namespace ISP {
class Internal {
  private:
    BC::Land land_boundary;

  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename InterfaceType>
    void ComputeFlux(InterfaceType& intface);
};

template <typename InterfaceType>
void Internal::Initialize(InterfaceType& intface) {
    land_boundary.Initialize(1);
}

template <typename InterfaceType>
void Internal::ComputeFlux(InterfaceType& intface) {
    bool wet_in = intface.data_in.wet_dry_state.wet;
    bool wet_ex = intface.data_ex.wet_dry_state.wet;

    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    // assemble numerical fluxes
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {

        LLF_flux(Global::g,
                 column(boundary_in.q_at_gp, gp),
                 column(boundary_ex.q_at_gp, gp),
                 column(boundary_in.aux_at_gp, gp),
                 column(intface.surface_normal_in, gp),
                 column(boundary_in.F_hat_at_gp, gp));

        column(boundary_ex.F_hat_at_gp, gp) = -column(boundary_in.F_hat_at_gp, gp);

        if (boundary_in.F_hat_at_gp(Variables::ze, gp) > 1e-12) {
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
        }

        assert(!std::isnan(boundary_in.F_hat_at_gp(Variables::ze, gp)));
    }
}
}
}
}

#endif