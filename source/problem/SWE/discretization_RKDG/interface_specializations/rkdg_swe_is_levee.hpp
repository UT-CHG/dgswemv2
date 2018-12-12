#ifndef RKDG_SWE_IS_LEVEE_HPP
#define RKDG_SWE_IS_LEVEE_HPP

namespace SWE {
namespace RKDG {
namespace ISP {
class Levee {
  private:
    double H_tolerance = 0.01;

    HybMatrix<double, SWE::n_variables> q_in_ex;
    HybMatrix<double, SWE::n_variables> q_ex_ex;

    DynRowVector<double> H_barrier;
    DynRowVector<double> C_subcritical;
    DynRowVector<double> C_supercritical;

    DynRowVector<double> H_bar_gp;
    DynRowVector<double> C_subcrit_gp;
    DynRowVector<double> C_supercrit_gp;

    BC::Land land_boundary;

  public:
    Levee() = default;
    Levee(const std::vector<LeveeInput>& levee_input);

    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename InterfaceType>
    void ComputeFlux(InterfaceType& intface);
};

Levee::Levee(const std::vector<LeveeInput>& levee_input) {
    uint n_nodes = levee_input.size();

    this->H_barrier.resize(n_nodes);
    this->C_subcritical.resize(n_nodes);
    this->C_supercritical.resize(n_nodes);

    for (uint node = 0; node < n_nodes; ++node) {
        this->H_barrier[node]       = levee_input[node].H_barrier;
        this->C_subcritical[node]   = levee_input[node].C_subcritical;
        this->C_supercritical[node] = levee_input[node].C_supercritical;
    }
}

template <typename InterfaceType>
void Levee::Initialize(InterfaceType& intface) {
    uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);

    this->q_in_ex.resize(SWE::n_variables, ngp);
    this->q_ex_ex.resize(SWE::n_variables, ngp);

    this->H_bar_gp       = intface.ComputeBoundaryNodalUgpIN(this->H_barrier);
    this->C_subcrit_gp   = intface.ComputeBoundaryNodalUgpIN(this->C_subcritical);
    this->C_supercrit_gp = intface.ComputeBoundaryNodalUgpIN(this->C_supercritical);
}

template <typename InterfaceType>
void Levee::ComputeFlux(InterfaceType& intface) {
    /*bool wet_in = intface.data_in.wet_dry_state.wet;
    bool wet_ex = intface.data_ex.wet_dry_state.wet;

    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);

    double H_levee, C_subcrit, C_supercrit;
    double h_above_levee_in, h_above_levee_ex;
    double gravity_in, gravity_ex;

    for (uint gp = 0; gp < ngp; ++gp) {

        H_levee     = this->H_bar_gp[gp];
        C_subcrit   = this->C_subcrit_gp[gp];
        C_supercrit = this->C_supercrit_gp[gp];

        gravity_in = Global::g;
        gravity_ex = Global::g;

        h_above_levee_in = boundary_in.q_at_gp[SWE::Variables::ze][gp] - H_levee;
        h_above_levee_ex = boundary_ex.q_at_gp[SWE::Variables::ze][gp] - H_levee;

        if ((h_above_levee_in <= H_tolerance && h_above_levee_ex <= H_tolerance) ||  // both side below or
            std::abs(h_above_levee_in - h_above_levee_ex) <= H_tolerance) {          // equal within tolerance

            // reflective boundary in
            this->land_boundary.GetEX(
                column(intface.surface_normal_in, gp), column(boundary_in.q_at_gp, gp), column(this->q_in_ex, gp));

            // reflective boundary ex
            this->land_boundary.GetEX(column(intface.surface_normal_ex, gp_ex),
                                      column(boundary_ex.q_at_gp, gp_ex),
                                      column(this->q_ex_ex, gp_ex));
        } else if (h_above_levee_in > h_above_levee_ex) {  // overtopping from in to ex
            double n_x, n_y, t_x, t_y, qn_in, qt_in;

            n_x = intface.surface_normal_in(GlobalCoord::x, gp);
            n_y = intface.surface_normal_in(GlobalCoord::y, gp);
            t_x = -n_y;
            t_y = n_x;

            if (h_above_levee_ex > 2.0 * h_above_levee_in / 3.0) {  // subcritical flow
                qn_in =
                    C_subcrit * h_above_levee_ex * std::sqrt(2.0 * Global::g * (h_above_levee_in - h_above_levee_ex));
                qt_in = 0.0;
            } else {  // supercritical flow
                qn_in = C_supercrit * (2.0 * h_above_levee_in / 3.0) *
                        std::sqrt(Global::g * (2.0 * h_above_levee_in / 3.0));
                qt_in = 0.0;
            }

            this->q_in_ex(SWE::Variables::ze, gp) = boundary_in.q_at_gp(SWE::Variables::ze, gp);
            this->q_in_ex(SWE::Variables::qx, gp) = qn_in * n_x + qt_in * t_x;
            this->q_in_ex(SWE::Variables::qy, gp) = qn_in * n_y + qt_in * t_y;

            this->q_ex_ex(SWE::Variables::ze, gp_ex) = boundary_ex.q_at_gp(SWE::Variables::ze, gp_ex);
            this->q_ex_ex(SWE::Variables::qx, gp_ex) = this->q_in_ex(SWE::Variables::qx, gp);
            this->q_ex_ex(SWE::Variables::qy, gp_ex) = this->q_in_ex(SWE::Variables::qy, gp);

            if (!wet_ex) {
                gravity_ex = 0.0;
            }
        } else if (h_above_levee_in < h_above_levee_ex) {  // overtopping from ex to in
            double n_x, n_y, t_x, t_y, qn_ex, qt_ex;

            n_x = intface.surface_normal_ex(GlobalCoord::x, gp_ex);
            n_y = intface.surface_normal_ex(GlobalCoord::y, gp_ex);
            t_x = -n_y;
            t_y = n_x;

            if (h_above_levee_in > 2.0 * h_above_levee_ex / 3.0) {  // subcritical flow
                qn_ex =
                    C_subcrit * h_above_levee_in * std::sqrt(2.0 * Global::g * (h_above_levee_ex - h_above_levee_in));
                qt_ex = 0.0;
            } else {  // supercritical flow
                qn_ex = C_supercrit * (2.0 * h_above_levee_ex / 3.0) *
                        std::sqrt(Global::g * (2.0 * h_above_levee_ex / 3.0));
                qt_ex = 0.0;
            }

            this->q_ex_ex(SWE::Variables::ze, gp_ex) = boundary_ex.q_at_gp(SWE::Variables::ze, gp_ex);
            this->q_ex_ex(SWE::Variables::qx, gp_ex) = qn_ex * n_x + qt_ex * t_x;
            this->q_ex_ex(SWE::Variables::qy, gp_ex) = qn_ex * n_y + qt_ex * t_y;

            this->q_in_ex(SWE::Variables::ze, gp) = boundary_in.q_at_gp(SWE::Variables::ze, gp);
            this->q_in_ex(SWE::Variables::qx, gp) = this->q_ex_ex(SWE::Variables::qx, gp_ex);
            this->q_in_ex(SWE::Variables::qy, gp) = this->q_ex_ex(SWE::Variables::qy, gp_ex);

            if (!wet_in) {
                gravity_in = 0.0;
            }
        }

        LLF_flux(gravity_in,
                 column(boundary_in.q_at_gp, gp),
                 column(this->q_in_ex, gp),
                 column(boundary_in.aux_at_gp, gp),
                 column(intface.surface_normal_in, gp),
                 column(boundary_in.F_hat_at_gp, gp));

        LLF_flux(gravity_ex,
                 column(boundary_ex.q_at_gp, gp_ex),
                 column(this->q_ex_ex, gp_ex),
                 column(boundary_ex.aux_at_gp, gp_ex),
                 column(intface.surface_normal_ex, gp_ex),
                 column(boundary_ex.F_hat_at_gp, gp_ex));
                 }*/
}
}
}
    }

#endif