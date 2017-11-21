#ifndef SWE_KERNELS_PROCESSOR_HPP
#define SWE_KERNELS_PROCESSOR_HPP

#include "swe_LLF_flux.hpp"

namespace SWE {
void Problem::wetting_drying_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& wd_state = elt.data.wet_dry_state;
        auto& internal = elt.data.internal;

        elt.ComputeUvrtx(state.ze, wd_state.ze_at_vrtx);

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.h_at_vrtx[vrtx] = wd_state.ze_at_vrtx[vrtx] + wd_state.bath_at_vrtx[vrtx];
        }

        double h_avg =
            std::accumulate(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end(), 0.0) / elt.data.get_nvrtx();

        bool set_wet_element = false;
        bool set_dry_element = false;
        bool check_element = false;

        if (h_avg <= Global::h_o) {
            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                wd_state.ze_at_vrtx[vrtx] = h_avg - wd_state.bath_at_vrtx[vrtx];
            }

            state.ze = elt.L2Projection(wd_state.ze_at_vrtx);

            set_dry_element = true;
        } else {
            uint n_dry_vrtx = 0;

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                if (wd_state.h_at_vrtx[vrtx] <= Global::h_o)
                    n_dry_vrtx++;
            }

            if (n_dry_vrtx == 0) {
                if (wd_state.wet) {
                    set_wet_element = true;
                } else {
                    check_element = true;
                }
            } else {
                auto h_min_iter = std::min_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
                uint h_min_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_min_iter);

                auto h_max_iter = std::max_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
                uint h_max_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_max_iter);

                uint h_mid_vrtx = 3 - h_max_vrtx - h_min_vrtx;

                wd_state.h_at_vrtx_temp[h_min_vrtx] = Global::h_o;

                wd_state.h_at_vrtx_temp[h_mid_vrtx] =
                    std::max(Global::h_o,
                             wd_state.h_at_vrtx[h_mid_vrtx] -
                                 (wd_state.h_at_vrtx_temp[h_min_vrtx] - wd_state.h_at_vrtx[h_min_vrtx]) / 2);

                wd_state.h_at_vrtx_temp[h_max_vrtx] =
                    wd_state.h_at_vrtx[h_max_vrtx] -
                    (wd_state.h_at_vrtx_temp[h_min_vrtx] - wd_state.h_at_vrtx[h_min_vrtx]) -
                    (wd_state.h_at_vrtx_temp[h_mid_vrtx] - wd_state.h_at_vrtx[h_mid_vrtx]);

                wd_state.h_at_vrtx = wd_state.h_at_vrtx_temp;

                elt.ComputeUvrtx(state.qx, wd_state.qx_at_vrtx);
                elt.ComputeUvrtx(state.qy, wd_state.qy_at_vrtx);

                double del_qx = 0;
                double del_qy = 0;

                n_dry_vrtx = 0;

                for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                    if (wd_state.h_at_vrtx[vrtx] <= Global::h_o) {
                        n_dry_vrtx++;

                        del_qx += wd_state.qx_at_vrtx[vrtx];
                        del_qy += wd_state.qy_at_vrtx[vrtx];

                        wd_state.qx_at_vrtx[vrtx] = 0;
                        wd_state.qy_at_vrtx[vrtx] = 0;
                    }
                }

                for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
                    wd_state.ze_at_vrtx[vrtx] = wd_state.h_at_vrtx[vrtx] - wd_state.bath_at_vrtx[vrtx];

                    if (wd_state.h_at_vrtx[vrtx] > Global::h_o) {
                        wd_state.qx_at_vrtx[vrtx] += del_qx / (3 - n_dry_vrtx);
                        wd_state.qy_at_vrtx[vrtx] += del_qy / (3 - n_dry_vrtx);
                    }
                }

                state.ze = elt.L2Projection(wd_state.ze_at_vrtx);
                state.qx = elt.L2Projection(wd_state.qx_at_vrtx);
                state.qy = elt.L2Projection(wd_state.qy_at_vrtx);

                if (wd_state.wet) {
                    set_wet_element = true;
                } else {
                    check_element = true;
                }
            }

            if (check_element) {
                auto h_max_iter = std::max_element(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end());
                uint h_max_vrtx = std::distance(wd_state.h_at_vrtx.begin(), h_max_iter);

                double ze_h_max_vrtx = wd_state.ze_at_vrtx[h_max_vrtx];

                if (ze_h_max_vrtx > Global::h_o - wd_state.bath_min) {
                    set_wet_element = true;
                } else {
                    set_dry_element = true;
                }
            };
        }

        if (set_dry_element) {
            wd_state.wet = false;

            std::fill(state.qx.begin(), state.qx.end(), 0.0);
            std::fill(state.qy.begin(), state.qy.end(), 0.0);

            std::fill(state.rhs_ze.begin(), state.rhs_ze.end(), 0.0);
            std::fill(state.rhs_qx.begin(), state.rhs_qx.end(), 0.0);
            std::fill(state.rhs_qy.begin(), state.rhs_qy.end(), 0.0);
        } else if (set_wet_element) {
            wd_state.wet = true;

            elt.ComputeUgp(state.ze, internal.ze_at_gp);
            elt.ComputeUgp(state.qx, internal.qx_at_gp);
            elt.ComputeUgp(state.qy, internal.qy_at_gp);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                internal.h_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];
            }

            wd_state.water_volume = elt.Integration(internal.h_at_gp);
        }
    });
}

void Problem::slope_limiting_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& sl_state = elt.data.slope_limit_state;

        elt.ComputeUbaryctr(state.ze, sl_state.ze_at_baryctr);
        elt.ComputeUbaryctr(state.qx, sl_state.qx_at_baryctr);
        elt.ComputeUbaryctr(state.qy, sl_state.qy_at_baryctr);

        elt.ComputeUmidpts(state.ze, sl_state.ze_at_midpts);
        elt.ComputeUmidpts(state.qx, sl_state.qx_at_midpts);
        elt.ComputeUmidpts(state.qy, sl_state.qy_at_midpts);
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.ze_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.ze_at_baryctr;
        sl_state_in.qx_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.qx_at_baryctr;
        sl_state_in.qy_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.qy_at_baryctr;

        sl_state_ex.ze_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.ze_at_baryctr;
        sl_state_ex.qx_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.qx_at_baryctr;
        sl_state_ex.qy_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.qy_at_baryctr;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.ze_at_baryctr_neigh[bound.bound_id] = sl_state.ze_at_baryctr;
        sl_state.qx_at_baryctr_neigh[bound.bound_id] = sl_state.qx_at_baryctr;
        sl_state.qy_at_baryctr_neigh[bound.bound_id] = sl_state.qy_at_baryctr;
    });

    mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& sl_state = elt.data.slope_limit_state;

        double u = sl_state.qx_at_baryctr / (sl_state.ze_at_baryctr + sl_state.bath_at_baryctr);
        double v = sl_state.qy_at_baryctr / (sl_state.ze_at_baryctr + sl_state.bath_at_baryctr);
        double c = std::sqrt(Global::g * (sl_state.ze_at_baryctr + sl_state.bath_at_baryctr));

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % 3;

            sl_state.R[0][0] = 1.0;
            sl_state.R[1][0] = u + c * sl_state.n[bound][GlobalCoord::x];
            sl_state.R[2][0] = v + c * sl_state.n[bound][GlobalCoord::y];

            sl_state.R[0][1] = 0.0;
            sl_state.R[1][1] = -sl_state.n[bound][GlobalCoord::y];
            sl_state.R[2][1] = sl_state.n[bound][GlobalCoord::x];

            sl_state.R[0][2] = 1.0;
            sl_state.R[1][2] = u - c * sl_state.n[bound][GlobalCoord::x];
            sl_state.R[2][2] = v - c * sl_state.n[bound][GlobalCoord::y];

            double det = R[0][0] * R[1][1] * R[2][2] + R[0][1] * R[1][2] * R[2][0] + R[0][2] * R[1][0] * R[2][1] -
                         R[0][0] * R[1][2] * R[2][1] - R[0][1] * R[1][0] * R[2][2] - R[0][2] * R[1][1] * R[2][0];

            sl_state.L[0][0] = (R[1][1] * R[2][2] - R[1][2] * R[2][1]) / det;
            sl_state.L[1][0] = (R[1][2] * R[2][0] - R[1][0] * R[2][2]) / det;
            sl_state.L[2][0] = (R[1][0] * R[2][1] - R[1][1] * R[2][0]) / det;

            sl_state.L[0][1] = (R[0][2] * R[2][1] - R[0][1] * R[2][2]) / det;
            sl_state.L[1][1] = (R[0][0] * R[2][2] - R[0][2] * R[2][0]) / det;
            sl_state.L[2][1] = (R[0][1] * R[2][0] - R[0][0] * R[2][1]) / det;

            sl_state.L[0][2] = (R[0][1] * R[1][2] - R[0][2] * R[1][1]) / det;
            sl_state.L[1][2] = (R[0][2] * R[1][0] - R[0][0] * R[1][2]) / det;
            sl_state.L[2][2] = (R[0][0] * R[1][1] - R[0][1] * R[1][0]) / det;

            for (uint var = 0; var < 3; var++) {
                sl_state.w_midpt_char[var] = sl_state.L[var][0] * sl_state.ze_at_midpts[bound] +
                                             sl_state.L[var][1] * sl_state.qx_at_midpts[bound] +
                                             sl_state.L[var][2] * sl_state.qy_at_midpts[bound];

                sl_state.w_baryctr_char[var][0] = sl_state.L[var][0] * sl_state.ze_at_baryctr +
                                                  sl_state.L[var][1] * sl_state.qx_at_baryctr +
                                                  sl_state.L[var][2] * sl_state.qy_at_baryctr;

                sl_state.w_baryctr_char[var][1] = sl_state.L[var][0] * sl_state.ze_at_baryctr_neigh[element_1] +
                                                  sl_state.L[var][1] * sl_state.qx_at_baryctr_neigh[element_1] +
                                                  sl_state.L[var][2] * sl_state.qy_at_baryctr_neigh[element_1];

                sl_state.w_baryctr_char[var][2] = sl_state.L[var][0] * sl_state.ze_at_baryctr_neigh[element_2] +
                                                  sl_state.L[var][1] * sl_state.qx_at_baryctr_neigh[element_2] +
                                                  sl_state.L[var][2] * sl_state.qy_at_baryctr_neigh[element_2];
            }

            double w_tilda;
            double w_delta;

            double M = 50;
            double nu = 1.5;

            for (uint var = 0; var < 3; var++) {
                w_tilda = sl_state.w_midpt_char[var] - sl_state.w_baryctr_char[var][0];

                w_delta =
                    sl_state.alpha_1[bound] * (sl_state.w_baryctr_char[var][1] - sl_state.w_baryctr_char[var][0]) +
                    sl_state.alpha_2[bound] * (sl_state.w_baryctr_char[var][2] - sl_state.w_baryctr_char[var][0]);

                // TVB modified minmod
                if (std::abs(w_tilda) <= M * sl_state.r_sq[bound]) {
                    sl_state.delta_char[var] = w_tilda;
                } else if (std::signbit(w_tilda) == std::signbit(w_delta)) {
                    sl_state.delta_char[var] =
                        std::copysign(1.0, w_tilda) * std::min(std::abs(w_tilda), std::abs(nu * w_delta));
                } else {
                    sl_state.delta_char[var] = 0.0;
                }
            }

            for (uint var = 0; var < 3; var++) {
                sl_state.delta[var][bound] = sl_state.R[var][0] * sl_state.delta_char[0] +
                                             sl_state.R[var][1] * sl_state.delta_char[1] +
                                             sl_state.R[var][2] * sl_state.delta_char[2];
            }
        }

        for (uint var = 0; var < 3; var++) {
            double delta_sum = 0.0;

            for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
                delta_sum += sl_state.delta[var][bound];
            }

            if (delta_sum != 0.0) {
                double positive = 0.0;
                double negative = 0.0;

                for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
                    positive += std::max(0.0, sl_state.delta[var][bound]);
                    negative += std::max(0.0, -sl_state.delta[var][bound]);
                }

                double theta_positive = std::min(1.0, negative / positive);
                double theta_negative = std::min(1.0, positive / negative);

                for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
                    sl_state.delta[var][bound] = theta_positive * std::max(0.0, sl_state.delta[var][bound]) -
                                                 theta_negative * std::max(0.0, -sl_state.delta[var][bound]);
                }
            }
        }

        Array2D<double> T{{-1.0, 1.0, 1.0}, {1.0, -1.0, 1.0}, {1.0, 1.0, -1.0}};

        for (uint vrtx = 0; vrtx < 3; vrtx++) {
            sl_state.ze_at_vrtx[vrtx] = sl_state.ze_at_baryctr;
            sl_state.qx_at_vrtx[vrtx] = sl_state.qx_at_baryctr;
            sl_state.qy_at_vrtx[vrtx] = sl_state.qy_at_baryctr;

            for (uint bound = 0; elt.data.get_nbound(); bound++) {
                sl_state.ze_at_vrtx[vrtx] += T[vrtx][bound] * sl_state.delta[0][bound];
                sl_state.qx_at_vrtx[vrtx] += T[vrtx][bound] * sl_state.delta[1][bound];
                sl_state.qy_at_vrtx[vrtx] += T[vrtx][bound] * sl_state.delta[2][bound];
            }
        }

        state.ze = elt.L2Projection(sl_state.ze_at_vrtx);
        state.qx = elt.L2Projection(sl_state.qx_at_vrtx);
        state.qy = elt.L2Projection(sl_state.qy_at_vrtx);
    });
}

template <typename ElementType>
void Problem::volume_kernel(const Stepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& internal = elt.data.internal;

        // assemble flux
        // note we assume that the values at gauss points have already been computed
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
            internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];

            internal.qx_flux_at_gp[GlobalCoord::x][gp] = std::pow(internal.qx_at_gp[gp], 2) / internal.h_at_gp[gp] +
                                                         Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) +
                                                                      internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
            internal.qx_flux_at_gp[GlobalCoord::y][gp] =
                internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.h_at_gp[gp];

            internal.qy_flux_at_gp[GlobalCoord::x][gp] =
                internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.h_at_gp[gp];
            internal.qy_flux_at_gp[GlobalCoord::y][gp] = std::pow(internal.qy_at_gp[gp], 2) / internal.h_at_gp[gp] +
                                                         Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) +
                                                                      internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
        }

        // skip dof = 0, which is a constant and thus trivially 0 NOT ALWAYS!
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.ze_flux_at_gp[GlobalCoord::x]) +
                                elt.IntegrationDPhi(GlobalCoord::y, dof, internal.ze_flux_at_gp[GlobalCoord::y]);

            state.rhs_qx[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qx_flux_at_gp[GlobalCoord::x]) +
                                elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qx_flux_at_gp[GlobalCoord::y]);

            state.rhs_qy[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qy_flux_at_gp[GlobalCoord::x]) +
                                elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qy_flux_at_gp[GlobalCoord::y]);
        }
    }
}

template <typename ElementType>
void Problem::source_kernel(const Stepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& internal = elt.data.internal;

        double t = stepper.get_t_at_curr_stage();

        auto source_ze = [t](Point<2>& pt) { return SWE::source_ze(t, pt); };

        auto source_qx = [t](Point<2>& pt) { return SWE::source_qx(t, pt); };

        auto source_qy = [t](Point<2>& pt) { return SWE::source_qy(t, pt); };

        elt.ComputeFgp(source_ze, internal.ze_source_term_at_gp);
        elt.ComputeFgp(source_qx, internal.qx_source_term_at_gp);
        elt.ComputeFgp(source_qy, internal.qy_source_term_at_gp);

        // note we assume that the values at gauss points have already been computed
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute contribution of hydrostatic pressure
            internal.qx_source_term_at_gp[gp] +=
                Global::g * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
            internal.qy_source_term_at_gp[gp] +=
                Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];

            double u_at_gp = internal.qx_at_gp[gp] / internal.h_at_gp[gp];
            double v_at_gp = internal.qy_at_gp[gp] / internal.h_at_gp[gp];

            // compute bottom friction contribution
            double bottom_friction_stress = Global::Cf * std::hypot(u_at_gp, v_at_gp) / internal.h_at_gp[gp];

            internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
            internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
        }

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
            state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
            state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
        }
    }
}

template <typename InterfaceType>
void Problem::interface_kernel(const Stepper& stepper, InterfaceType& intface) {
    auto& wd_state_in = intface.data_in.wet_dry_state;
    auto& wd_state_ex = intface.data_ex.wet_dry_state;

    if (wd_state_in.wet || wd_state_ex.wet) {
        const uint stage = stepper.get_stage();
        const double dt = stepper.get_dt();

        auto& state_in = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        intface.ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
        intface.ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
        intface.ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

        intface.ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
        intface.ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
        intface.ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);

        // assemble numerical fluxes
        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex = 0;
        for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
            gp_ex = ngp - gp - 1;

            LLF_flux(boundary_in.ze_at_gp[gp],
                     boundary_ex.ze_at_gp[gp_ex],
                     boundary_in.qx_at_gp[gp],
                     boundary_ex.qx_at_gp[gp_ex],
                     boundary_in.qy_at_gp[gp],
                     boundary_ex.qy_at_gp[gp_ex],
                     boundary_in.bath_at_gp[gp],
                     intface.surface_normal_in[gp],
                     boundary_in.ze_numerical_flux_at_gp[gp],
                     boundary_in.qx_numerical_flux_at_gp[gp],
                     boundary_in.qy_numerical_flux_at_gp[gp]);

            boundary_ex.ze_numerical_flux_at_gp[gp_ex] = -boundary_in.ze_numerical_flux_at_gp[gp];
            boundary_ex.qx_numerical_flux_at_gp[gp_ex] = -boundary_in.qx_numerical_flux_at_gp[gp];
            boundary_ex.qy_numerical_flux_at_gp[gp_ex] = -boundary_in.qy_numerical_flux_at_gp[gp];
        }

        // compute net volume flux out of IN/EX elements
        double net_volume_flux_in = 0;
        double net_volume_flux_ex = 0;

        net_volume_flux_in = intface.IntegrationIN(boundary_in.ze_numerical_flux_at_gp);
        net_volume_flux_ex = -net_volume_flux_in;

        if (net_volume_flux_in > 0) {
            if (!wd_state_in.wet) {  // water flowing from dry IN element
                // Zero flux on IN element side
                std::fill(boundary_in.ze_numerical_flux_at_gp.begin(), boundary_in.ze_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_in.qx_numerical_flux_at_gp.begin(), boundary_in.qx_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_in.qy_numerical_flux_at_gp.begin(), boundary_in.qy_numerical_flux_at_gp.end(), 0.0);

                // Reflective Boundary on EX element side
                SWE::Land land_boundary;
                double ze_in, qx_in, qy_in;

                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    land_boundary.GetEX(stepper,
                                        gp,
                                        intface.surface_normal_ex,
                                        boundary_ex.ze_at_gp,
                                        boundary_ex.qx_at_gp,
                                        boundary_ex.qy_at_gp,
                                        ze_in,
                                        qx_in,
                                        qy_in);

                    LLF_flux(boundary_ex.ze_at_gp[gp],
                             ze_in,
                             boundary_ex.qx_at_gp[gp],
                             qx_in,
                             boundary_ex.qy_at_gp[gp],
                             qy_in,
                             boundary_ex.bath_at_gp[gp],
                             intface.surface_normal_ex[gp],
                             boundary_ex.ze_numerical_flux_at_gp[gp],
                             boundary_ex.qx_numerical_flux_at_gp[gp],
                             boundary_ex.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = 0;
                net_volume_flux_ex = 0;
            } else if (!wd_state_ex.wet) {  // water flowing to dry EX element
                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    gp_ex = ngp - gp - 1;

                    LLF_flux_zero_g(boundary_ex.ze_at_gp[gp_ex],
                                    boundary_in.ze_at_gp[gp],
                                    boundary_ex.qx_at_gp[gp_ex],
                                    boundary_in.qx_at_gp[gp],
                                    boundary_ex.qy_at_gp[gp_ex],
                                    boundary_in.qy_at_gp[gp],
                                    boundary_ex.bath_at_gp[gp_ex],
                                    intface.surface_normal_ex[gp_ex],
                                    boundary_ex.ze_numerical_flux_at_gp[gp_ex],
                                    boundary_ex.qx_numerical_flux_at_gp[gp_ex],
                                    boundary_ex.qy_numerical_flux_at_gp[gp_ex]);
                }

                net_volume_flux_in = intface.IntegrationIN(boundary_in.ze_numerical_flux_at_gp);
                net_volume_flux_ex = intface.IntegrationEX(boundary_ex.ze_numerical_flux_at_gp);
            }
        } else if (net_volume_flux_in < 0) {
            if (!wd_state_ex.wet) {  // water flowing from dry EX element
                // Zero flux on EX element side
                std::fill(boundary_ex.ze_numerical_flux_at_gp.begin(), boundary_ex.ze_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_ex.qx_numerical_flux_at_gp.begin(), boundary_ex.qx_numerical_flux_at_gp.end(), 0.0);
                std::fill(boundary_ex.qy_numerical_flux_at_gp.begin(), boundary_ex.qy_numerical_flux_at_gp.end(), 0.0);

                // Reflective Boundary on IN element side
                SWE::Land land_boundary;
                double ze_ex, qx_ex, qy_ex;

                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    land_boundary.GetEX(stepper,
                                        gp,
                                        intface.surface_normal_in,
                                        boundary_in.ze_at_gp,
                                        boundary_in.qx_at_gp,
                                        boundary_in.qy_at_gp,
                                        ze_ex,
                                        qx_ex,
                                        qy_ex);

                    LLF_flux(boundary_in.ze_at_gp[gp],
                             ze_ex,
                             boundary_in.qx_at_gp[gp],
                             qx_ex,
                             boundary_in.qy_at_gp[gp],
                             qy_ex,
                             boundary_in.bath_at_gp[gp],
                             intface.surface_normal_in[gp],
                             boundary_in.ze_numerical_flux_at_gp[gp],
                             boundary_in.qx_numerical_flux_at_gp[gp],
                             boundary_in.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = 0;
                net_volume_flux_ex = 0;
            } else if (!wd_state_in.wet) {  // water flowing to dry IN element
                for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
                    gp_ex = ngp - gp - 1;

                    LLF_flux_zero_g(boundary_in.ze_at_gp[gp],
                                    boundary_ex.ze_at_gp[gp_ex],
                                    boundary_in.qx_at_gp[gp],
                                    boundary_ex.qx_at_gp[gp_ex],
                                    boundary_in.qy_at_gp[gp],
                                    boundary_ex.qy_at_gp[gp_ex],
                                    boundary_in.bath_at_gp[gp],
                                    intface.surface_normal_in[gp],
                                    boundary_in.ze_numerical_flux_at_gp[gp],
                                    boundary_in.qx_numerical_flux_at_gp[gp],
                                    boundary_in.qy_numerical_flux_at_gp[gp]);
                }

                net_volume_flux_in = intface.IntegrationIN(boundary_in.ze_numerical_flux_at_gp);
                net_volume_flux_ex = intface.IntegrationEX(boundary_ex.ze_numerical_flux_at_gp);
            }
        }

        // compute water level reduction
        wd_state_in.water_volume -= net_volume_flux_in * dt;
        wd_state_ex.water_volume -= net_volume_flux_ex * dt;

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
            state_in.rhs_ze[dof] -= intface.IntegrationPhiIN(dof, boundary_in.ze_numerical_flux_at_gp);
            state_in.rhs_qx[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qx_numerical_flux_at_gp);
            state_in.rhs_qy[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qy_numerical_flux_at_gp);
        }

        for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
            state_ex.rhs_ze[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.ze_numerical_flux_at_gp);
            state_ex.rhs_qx[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qx_numerical_flux_at_gp);
            state_ex.rhs_qy[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qy_numerical_flux_at_gp);
        }
    }
}

template <typename BoundaryType>
void Problem::boundary_kernel(const Stepper& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.get_stage();

        auto& state = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        bound.ComputeUgp(state.ze, boundary.ze_at_gp);
        bound.ComputeUgp(state.qx, boundary.qx_at_gp);
        bound.ComputeUgp(state.qy, boundary.qy_at_gp);

        double ze_ex, qx_ex, qy_ex;
        for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
            bound.boundary_condition.GetEX(stepper,
                                           gp,
                                           bound.surface_normal,
                                           boundary.ze_at_gp,
                                           boundary.qx_at_gp,
                                           boundary.qy_at_gp,
                                           ze_ex,
                                           qx_ex,
                                           qy_ex);

            LLF_flux(boundary.ze_at_gp[gp],
                     ze_ex,
                     boundary.qx_at_gp[gp],
                     qx_ex,
                     boundary.qy_at_gp[gp],
                     qy_ex,
                     boundary.bath_at_gp[gp],
                     bound.surface_normal[gp],
                     boundary.ze_numerical_flux_at_gp[gp],
                     boundary.qx_numerical_flux_at_gp[gp],
                     boundary.qy_numerical_flux_at_gp[gp]);
        }

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] -= bound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
            state.rhs_qx[dof] -= bound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
            state.rhs_qy[dof] -= bound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
        }
    }
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_send_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.get_stage();

    auto& state = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.ze, boundary.ze_at_gp);
    dbound.ComputeUgp(state.qx, boundary.qx_at_gp);
    dbound.ComputeUgp(state.qy, boundary.qy_at_gp);

    dbound.boundary_condition.SetEX(boundary.ze_at_gp, boundary.qx_at_gp, boundary.qy_at_gp);
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.get_stage();

    auto& state = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        dbound.boundary_condition.GetEX(stepper,
                                        gp,
                                        dbound.surface_normal,
                                        boundary.ze_at_gp,
                                        boundary.qx_at_gp,
                                        boundary.qy_at_gp,
                                        ze_ex,
                                        qx_ex,
                                        qy_ex);

        LLF_flux(boundary.ze_at_gp[gp],
                 ze_ex,
                 boundary.qx_at_gp[gp],
                 qx_ex,
                 boundary.qy_at_gp[gp],
                 qy_ex,
                 boundary.bath_at_gp[gp],
                 dbound.surface_normal[gp],
                 boundary.ze_numerical_flux_at_gp[gp],
                 boundary.qx_numerical_flux_at_gp[gp],
                 boundary.qy_numerical_flux_at_gp[gp]);
    }

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < dbound.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] -= dbound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
        state.rhs_qx[dof] -= dbound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
        state.rhs_qy[dof] -= dbound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
    }
}

template <typename ElementType>
void Problem::update_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();
    const double dt = stepper.get_dt();

    auto& state = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    curr_state.rhs_ze = elt.ApplyMinv(curr_state.rhs_ze);
    curr_state.rhs_qx = elt.ApplyMinv(curr_state.rhs_qx);
    curr_state.rhs_qy = elt.ApplyMinv(curr_state.rhs_qy);

    std::fill(next_state.ze.begin(), next_state.ze.end(), 0);
    std::fill(next_state.qx.begin(), next_state.qx.end(), 0);
    std::fill(next_state.qy.begin(), next_state.qy.end(), 0);

    for (uint s = 0; s <= stage; ++s) {
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            next_state.ze[dof] +=
                stepper.ark[stage][s] * state[s].ze[dof] + dt * stepper.brk[stage][s] * state[s].rhs_ze[dof];

            next_state.qx[dof] +=
                stepper.ark[stage][s] * state[s].qx[dof] + dt * stepper.brk[stage][s] * state[s].rhs_qx[dof];

            next_state.qy[dof] +=
                stepper.ark[stage][s] * state[s].qy[dof] + dt * stepper.brk[stage][s] * state[s].rhs_qy[dof];
        }
    }
}

template <typename ElementType>
void Problem::swap_states_kernel(const Stepper& stepper, ElementType& elt) {
    uint n_stages = stepper.get_num_stages();
    auto& state = elt.data.state;

    std::swap(state[0].ze, state[n_stages].ze);
    std::swap(state[0].qx, state[n_stages].qx);
    std::swap(state[0].qy, state[n_stages].qy);
}

template <typename ElementType>
void Problem::scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt) {
    uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];

    for (auto& ze_mode : state.ze) {
        if (isnan(ze_mode)) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qx_mode : state.qx) {
        if (isnan(qx_mode)) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qy_mode : state.qy) {
        if (isnan(qy_mode)) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_ze_mode : state.rhs_ze) {
        if (isnan(rhs_ze_mode)) {
            std::cerr << "Error: found isnan rhs_ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qx_mode : state.rhs_qx) {
        if (isnan(rhs_qx_mode)) {
            std::cerr << "Error: found isnan rhs_qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qy_mode : state.rhs_qy) {
        if (isnan(rhs_qy_mode)) {
            std::cerr << "Error: found isnan rhs_qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }
}
}

#endif
