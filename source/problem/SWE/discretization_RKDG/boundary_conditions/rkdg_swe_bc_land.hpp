#ifndef RKDG_SWE_BC_LAND_HPP
#define RKDG_SWE_BC_LAND_HPP

namespace SWE {
namespace RKDG {
namespace BC {
class Land {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;

  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void Initialize(uint ngp);

    template <typename StepperType, typename BoundaryType>
    void ComputeFlux(const StepperType& stepper, BoundaryType& bound);

    void ComputeFlux(const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                     const Column<HybMatrix<double, SWE::n_variables>>& q_in,
                     const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
                     Column<HybMatrix<double, SWE::n_variables>>&& F_hat);

    void GetEX(const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
               const Column<HybMatrix<double, SWE::n_variables>>& q_in,
               Column<HybMatrix<double, SWE::n_variables>>&& q_ex);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);
    this->q_ex.resize(SWE::n_variables, ngp);
}

void Land::Initialize(uint ngp) {
    this->q_ex.resize(SWE::n_variables, ngp);
}

template <typename StepperType, typename BoundaryType>
void Land::ComputeFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];

    uint ngp = boundary.F_hat_at_gp[0].size();
    auto& n_x = bound.surface_normal[GlobalCoord::x];
    auto& n_y = bound.surface_normal[GlobalCoord::y];
    auto& t_x = -n_y;
    auto& t_y = n_x;

    auto qn_ex = -(vec_cw_mult(boundary.q_at_gp[SWE::Variables::qx], n_x) +
                   vec_cw_mult(boundary.q_at_gp[SWE::Variables::qy], n_y));
    auto qt_ex = vec_cw_mult(boundary.q_at_gp[SWE::Variables::qx], t_x) +
                 vec_cw_mult(boundary.q_at_gp[SWE::Variables::qy], t_y);

    row(this->q_ex, SWE::Variables::ze) = boundary.q_at_gp[SWE::Variables::ze];
    row(this->q_ex, SWE::Variables::qx) = vec_cw_mult(qn_ex, n_x) + vec_cw_mult(qt_ex, t_x);
    row(this->q_ex, SWE::Variables::qy) = vec_cw_mult(qn_ex, n_y) + vec_cw_mult(qt_ex, t_y);

    for (uint gp = 0; gp < ngp; ++gp) {
        std::array<double, SWE::n_variables> F_hat_tmp{boundary.F_hat_at_gp[SWE::Variables::ze][gp],
                boundary.F_hat_at_gp[SWE::Variables::qx][gp],
                boundary.F_hat_at_gp[SWE::Variables::qy][gp]};

        LLF_flux(Global::g,
                 boundary.q_at_gp[SWE::Variables::ze][gp],
                 boundary.q_at_gp[SWE::Variables::qx][gp],
                 boundary.q_at_gp[SWE::Variables::qy][gp],
                 this->q_ex(SWE::Variables::ze,gp),
                 this->q_ex(SWE::Variables::qx,gp),
                 this->q_ex(SWE::Variables::qy,gp),
                 std::array<double,SWE::n_auxiliaries>{boundary.aux_at_gp[SWE::Auxiliaries::bath][gp],
                         boundary.aux_at_gp[SWE::Auxiliaries::h][gp],
                         boundary.aux_at_gp[SWE::Auxiliaries::sp][gp]},
                 std::array<double,SWE::n_dimensions>{n_x[gp],n_y[gp]},
                 F_hat_tmp
            );

        for ( uint i = 0; i < SWE::n_variables; ++i ) {
            boundary.F_hat_at_gp[i][gp] = F_hat_tmp[i];
        }
    }
}

void Land::ComputeFlux(const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                       const Column<HybMatrix<double, SWE::n_variables>>& q_in,
                       const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
                       Column<HybMatrix<double, SWE::n_variables>>&& F_hat) {
    // *** //
    double n_x = surface_normal[GlobalCoord::x];
    double n_y = surface_normal[GlobalCoord::y];
    double t_x = -n_y;
    double t_y = n_x;

    double qn_ex = -(q_in[SWE::Variables::qx] * n_x + q_in[SWE::Variables::qy] * n_y);
    double qt_ex = q_in[SWE::Variables::qx] * t_x + q_in[SWE::Variables::qy] * t_y;

    this->q_ex(SWE::Variables::ze, 0) = q_in[SWE::Variables::ze];
    this->q_ex(SWE::Variables::qx, 0) = qn_ex * n_x + qt_ex * t_x;
    this->q_ex(SWE::Variables::qy, 0) = qn_ex * n_y + qt_ex * t_y;

    std::array<double, SWE::n_variables> F_hat_tmp{F_hat[SWE::Variables::ze],
            F_hat[SWE::Variables::qx],
            F_hat[SWE::Variables::qy]};

    LLF_flux(Global::g,
             q_in[SWE::Variables::ze],
             q_in[SWE::Variables::qx],
             q_in[SWE::Variables::qy],
             this->q_ex(SWE::Variables::ze,0),
             this->q_ex(SWE::Variables::qx,0),
             this->q_ex(SWE::Variables::qy,0),
             std::array<double,SWE::n_auxiliaries>{aux_in[SWE::Auxiliaries::bath],
                     aux_in[SWE::Auxiliaries::h],
                     aux_in[SWE::Auxiliaries::sp]},
             std::array<double,SWE::n_dimensions>{n_x,n_y},
             F_hat_tmp
        );

    for ( uint i = 0; i < SWE::n_variables; ++i ) {
        F_hat[i] = F_hat_tmp[i];
    }
}

void Land::GetEX(const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                 const Column<HybMatrix<double, SWE::n_variables>>& q_in,
                 Column<HybMatrix<double, SWE::n_variables>>&& q_ex) {
    double n_x, n_y, t_x, t_y, qn_in, qt_in, qn_ex, qt_ex;

    n_x = surface_normal[GlobalCoord::x];
    n_y = surface_normal[GlobalCoord::y];
    t_x = -n_y;
    t_y = n_x;

    qn_in = q_in[SWE::Variables::qx] * n_x + q_in[SWE::Variables::qy] * n_y;
    qt_in = q_in[SWE::Variables::qx] * t_x + q_in[SWE::Variables::qy] * t_y;

    qn_ex = -qn_in;
    qt_ex = qt_in;

    q_ex[SWE::Variables::ze] = q_in[SWE::Variables::ze];
    q_ex[SWE::Variables::qx] = qn_ex * n_x + qt_ex * t_x;
    q_ex[SWE::Variables::qy] = qn_ex * n_y + qt_ex * t_y;
}
}
}
}

#endif
