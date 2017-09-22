#ifndef SWE_KERNELS_PROCESSOR_HPP
#define SWE_KERNELS_PROCESSOR_HPP

#include "swe_LLF_flux.hpp"

namespace SWE {
template <typename ElementType>
void Problem::volume_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];
    auto& internal = elt.data.internal;

    // get state at Gauss points
    elt.ComputeUgp(state.ze, internal.ze_at_gp);
    elt.ComputeUgp(state.qx, internal.qx_at_gp);
    elt.ComputeUgp(state.qy, internal.qy_at_gp);

    // assemble flux
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.water_column_hgt_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];

        internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
        internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];

        internal.qx_flux_at_gp[GlobalCoord::x][gp] =
            std::pow(internal.qx_at_gp[gp], 2) / internal.water_column_hgt_at_gp[gp] +
            Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
        internal.qx_flux_at_gp[GlobalCoord::y][gp] =
            internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

        internal.qy_flux_at_gp[GlobalCoord::x][gp] =
            internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
        internal.qy_flux_at_gp[GlobalCoord::y][gp] =
            std::pow(internal.qy_at_gp[gp], 2) / internal.water_column_hgt_at_gp[gp] +
            Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
    }

    // skip dof = 0, which is a constant and thus trivially 0 NOT ALWAYS!
    for (uint dof = 1; dof < elt.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.ze_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.ze_flux_at_gp[GlobalCoord::y]);

        state.rhs_qx[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qx_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qx_flux_at_gp[GlobalCoord::y]);

        state.rhs_qy[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qy_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qy_flux_at_gp[GlobalCoord::y]);
    }
}

template <typename ElementType>
void Problem::source_kernel(const Stepper& stepper, ElementType& elt) {
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
        internal.qx_source_term_at_gp[gp] += Global::g * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
        internal.qy_source_term_at_gp[gp] += Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];

        double u_at_gp = internal.qx_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
        double v_at_gp = internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

        // compute bottom friction contribution
        double bottom_friction_stress = Global::Cf * std::hypot(u_at_gp, v_at_gp) / internal.water_column_hgt_at_gp[gp];

        internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
        internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
    }

    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
        state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
        state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
    }
}

template <typename InterfaceType>
void Problem::interface_kernel(const Stepper& stepper, InterfaceType& intface) {
    const uint stage = stepper.get_stage();

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
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.bound_id_in); ++gp) {
        uint gp_ex = intface.data_in.get_ngp_boundary(intface.bound_id_in) - gp - 1;

        LLF_flux(boundary_in.ze_at_gp[gp],
                 boundary_ex.ze_at_gp[gp_ex],
                 boundary_in.qx_at_gp[gp],
                 boundary_ex.qx_at_gp[gp_ex],
                 boundary_in.qy_at_gp[gp],
                 boundary_ex.qy_at_gp[gp_ex],
                 boundary_in.bath_at_gp[gp],
                 intface.surface_normal[gp],
                 boundary_in.ze_numerical_flux_at_gp[gp],
                 boundary_in.qx_numerical_flux_at_gp[gp],
                 boundary_in.qy_numerical_flux_at_gp[gp]);

        boundary_ex.ze_numerical_flux_at_gp[gp_ex] = -boundary_in.ze_numerical_flux_at_gp[gp];
        boundary_ex.qx_numerical_flux_at_gp[gp_ex] = -boundary_in.qx_numerical_flux_at_gp[gp];
        boundary_ex.qy_numerical_flux_at_gp[gp_ex] = -boundary_in.qy_numerical_flux_at_gp[gp];
    }

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

template <typename BoundaryType>
void Problem::boundary_kernel(const Stepper& stepper, BoundaryType& bound) {
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
    double dt = stepper.get_dt();

    auto& state = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    curr_state.rhs_ze = elt.SolveLSE(curr_state.rhs_ze);
    curr_state.rhs_qx = elt.SolveLSE(curr_state.rhs_qx);
    curr_state.rhs_qy = elt.SolveLSE(curr_state.rhs_qy);

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
        if (std::isnan(ze_mode)) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qx_mode : state.qx) {
        if (std::isnan(qx_mode)) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qy_mode : state.qy) {
        if (std::isnan(qy_mode)) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_ze_mode : state.rhs_ze) {
        if (std::isnan(rhs_ze_mode)) {
            std::cerr << "Error: found isnan rhs_ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qx_mode : state.rhs_qx) {
        if (std::isnan(rhs_qx_mode)) {
            std::cerr << "Error: found isnan rhs_qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qy_mode : state.rhs_qy) {
        if (std::isnan(rhs_qy_mode)) {
            std::cerr << "Error: found isnan rhs_qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }
}
}

#endif
