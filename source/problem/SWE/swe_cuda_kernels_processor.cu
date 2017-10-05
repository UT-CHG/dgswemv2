#include "swe_cuda_kernels.cuh"
#include "swe_definitions.hpp"
#include "general_definitions.hpp"

template <typename ElementType>
void cuda_volume_kernel1(uint block_dim, uint thread_dim, ElementType& elt, uint n_gp) {
    std::cout << "running" << std::endl;
    cuda_volume_kernel1<ElementType><<<block_dim, thread_dim>>>(elt, n_gp);
}

template <typename ElementType>
__global__ void cuda_volume_kernel1(ElementType& elt, uint n_gp) {
    auto& internal = elt.data.internal;
    uint gp = blockDim.x * blockIdx.x + threadIdx.x;
    if (gp < n_gp) {
        internal.water_column_hgt_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];
        internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
        internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];
        internal.qx_flux_at_gp[GlobalCoord::x][gp] = pow(internal.qx_at_gp[gp], 2) /
                internal.water_column_hgt_at_gp[gp] +
        SWE::Global::g * (0.5 * pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] *
                internal.bath_at_gp[gp]);
        internal.qx_flux_at_gp[GlobalCoord::y][gp] = internal.qx_at_gp[gp] * internal.qy_at_gp[gp] /
                internal.water_column_hgt_at_gp[gp];
        internal.qy_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp] * internal.qy_at_gp[gp] /
                internal.water_column_hgt_at_gp[gp];
        internal.qy_flux_at_gp[GlobalCoord::y][gp] = pow(internal.qy_at_gp[gp], 2) /
                internal.water_column_hgt_at_gp[gp] + SWE::Global::g * (0.5 *
                pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
    }
}

template <typename ElementType>
void cuda_volume_kernel2(uint block_dim, uint thread_dim, ElementType& elt, uint stage,
        uint n_dof) {
    std::cout << "running" << std::endl;
    cuda_volume_kernel2<ElementType><<<block_dim, thread_dim>>>(elt, stage, n_dof);
}

template <typename ElementType>
__global__ void cuda_volume_kernel2(ElementType& elt, uint stage, uint n_dof) {
    auto& internal = elt.data.internal;
    auto& state = elt.data.state[stage];
    uint dof = blockDim.x * blockIdx.x + threadIdx.x;
    if (dof < n_dof && dof > 0) {
        // skip dof = 0, which is a constant and thus trivially 0 NOT ALWAYS!
        state.rhs_ze[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof,
                internal.ze_flux_at_gp[GlobalCoord::x]) + elt.IntegrationDPhi(GlobalCoord::y, dof,
                internal.ze_flux_at_gp[GlobalCoord::y]);
        state.rhs_qx[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof,
                internal.qx_flux_at_gp[GlobalCoord::x]) + elt.IntegrationDPhi(GlobalCoord::y, dof,
                internal.qx_flux_at_gp[GlobalCoord::y]);
        state.rhs_qy[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof,
                internal.qy_flux_at_gp[GlobalCoord::x]) + elt.IntegrationDPhi(GlobalCoord::y, dof,
                internal.qy_flux_at_gp[GlobalCoord::y]);
    }
}

template <typename ElementType>
void cuda_source_kernel1(uint block_dim, uint thread_dim, ElementType& elt, uint n_gp) {
    std::cout << "running" << std::endl;
    cuda_source_kernel1<ElementType><<<block_dim, thread_dim>>>(elt, n_gp);
}

template <typename ElementType>
__global__ void cuda_source_kernel1(ElementType& elt, uint n_gp) {
    auto& internal = elt.data.internal;
    uint gp = blockDim.x * blockIdx.x + threadIdx.x;
    if (gp < n_gp) {
            // compute contribution of hydrostatic pressure
            internal.qx_source_term_at_gp[gp] += SWE::Global::g * internal.bath_deriv_wrt_x_at_gp[gp] *
                    internal.ze_at_gp[gp];
            internal.qy_source_term_at_gp[gp] += SWE::Global::g * internal.bath_deriv_wrt_y_at_gp[gp] *
                    internal.ze_at_gp[gp];
            double u_at_gp = internal.qx_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
            double v_at_gp = internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

            // compute bottom friction contribution
            double bottom_friction_stress = SWE::Global::Cf * std::hypot(u_at_gp, v_at_gp) /
                    internal.water_column_hgt_at_gp[gp];
            internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
            internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
    }
}

template <typename ElementType>
void cuda_source_kernel2(uint block_dim, uint thread_dim, ElementType& elt, uint stage,
        uint n_dof) {
    std::cout << "running" << std::endl;
    cuda_source_kernel2<ElementType><<<block_dim, thread_dim>>>(elt, stage, n_dof);
}

template <typename ElementType>
__global__ void cuda_source_kernel2(ElementType& elt, uint stage, uint n_dof) {
    auto& internal = elt.data.internal;
    auto& state = elt.data.state[stage];
    uint dof = blockDim.x * blockIdx.x + threadIdx.x;
    if (dof < n_dof) {
        state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
        state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
        state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
    }
}

template <typename InterfaceType>
void cuda_interface_kernel1(uint block_dim, uint thread_dim, InterfaceType& interface, uint n_gp) {
    std::cout << "running" << std::endl;
    cuda_interface_kernel1<InterfaceType><<<block_dim, thread_dim>>>(interface, n_gp);
}

template <typename InterfaceType>
__global__ void cuda_interface_kernel1(InterfaceType& interface, uint n_gp) {
  uint gp = blockDim.x * blockIdx.x + threadIdx.x;
  auto& boundary_ex = interface.data_ex.boundary[interface.bound_id_ex];
  auto& boundary_in = interface.data_in.boundary[interface.bound_id_in];
  if (gp < n_gp) {
        uint gp_ex = interface.data_in.get_ngp_boundary(interface.bound_id_in) - gp - 1;

        LLF_flux(boundary_in.ze_at_gp[gp],
                 boundary_ex.ze_at_gp[gp_ex],
                 boundary_in.qx_at_gp[gp],
                 boundary_ex.qx_at_gp[gp_ex],
                 boundary_in.qy_at_gp[gp],
                 boundary_ex.qy_at_gp[gp_ex],
                 boundary_in.bath_at_gp[gp],
                 interface.surface_normal[gp],
                 boundary_in.ze_numerical_flux_at_gp[gp],
                 boundary_in.qx_numerical_flux_at_gp[gp],
                 boundary_in.qy_numerical_flux_at_gp[gp]);

        boundary_ex.ze_numerical_flux_at_gp[gp_ex] = -boundary_in.ze_numerical_flux_at_gp[gp];
        boundary_ex.qx_numerical_flux_at_gp[gp_ex] = -boundary_in.qx_numerical_flux_at_gp[gp];
        boundary_ex.qy_numerical_flux_at_gp[gp_ex] = -boundary_in.qy_numerical_flux_at_gp[gp];
    }
}

template <typename InterfaceType>
void cuda_interface_kernel2(uint block_dim, uint thread_dim, InterfaceType& interface, uint stage,
       uint n_dof) {
    std::cout << "running" << std::endl;
    cuda_interface_kernel2<InterfaceType><<<block_dim, thread_dim>>>(interface, stage, n_dof);
}

template <typename InterfaceType>
__global__ void cuda_interface_kernel2(InterfaceType& interface, uint stage, uint n_dof) {
    uint dof = blockIdx.x * blockDim.x + threadIdx.x;
    auto& state_in = interface.data_in.state[stage];
    auto& boundary_in = interface.data_in.boundary[interface.bound_id_in];
    if (dof < n_dof) {
        state_in.rhs_ze[dof] -= interface.IntegrationPhiIN(dof, boundary_in.ze_numerical_flux_at_gp);
        state_in.rhs_qx[dof] -= interface.IntegrationPhiIN(dof, boundary_in.qx_numerical_flux_at_gp);
        state_in.rhs_qy[dof] -= interface.IntegrationPhiIN(dof, boundary_in.qy_numerical_flux_at_gp);
    }
}

template <typename InterfaceType>
void cuda_interface_kernel3(uint block_dim, uint thread_dim, InterfaceType& interface, uint stage,
        uint n_dof) {
    std::cout << "running" << std::endl;
  cuda_interface_kernel3<InterfaceType><<<block_dim, thread_dim>>>(interface, stage, n_dof);
}

template <typename InterfaceType>
__global__ void cuda_interface_kernel3(InterfaceType& interface, uint stage, uint n_dof) {
    uint dof = blockIdx.x * blockDim.x + threadIdx.x;
    auto& state_ex = interface.data_ex.state[stage];
    auto& boundary_ex = interface.data_ex.boundary[interface.bound_id_ex];
    if (dof < n_dof) {
        state_ex.rhs_ze[dof] -= interface.IntegrationPhiEX(dof, boundary_ex.ze_numerical_flux_at_gp);
        state_ex.rhs_qx[dof] -= interface.IntegrationPhiEX(dof, boundary_ex.qx_numerical_flux_at_gp);
        state_ex.rhs_qy[dof] -= interface.IntegrationPhiEX(dof, boundary_ex.qy_numerical_flux_at_gp);
    }
}

template <typename BoundaryType>
void cuda_boundary_kernel1(uint block_dim, uint thread_dim, const Stepper& stepper,
        BoundaryType& bound, uint n_gp) {
    std::cout << "running" << std::endl;
    cuda_boundary_kernel1<BoundaryType><<<block_dim, thread_dim>>>(stepper, bound, n_gp);
}

template <typename BoundaryType>
__global__ void cuda_boundary_kernel1(const Stepper& stepper, BoundaryType& bound, uint n_gp) {
    uint gp = blockIdx.x * blockDim.x + threadIdx.x;
    auto& boundary = bound.data.boundary[bound.bound_id];
    if (gp < n_gp) {
        double ze_ex, qx_ex, qy_ex;
        bound.boundary_condition.GetEX(stepper, gp, bound.surface_normal, boundary.ze_at_gp,
                boundary.qx_at_gp, boundary.qy_at_gp, ze_ex, qx_ex, qy_ex);
        LLF_flux(boundary.ze_at_gp[gp], ze_ex, boundary.qx_at_gp[gp], qx_ex, boundary.qy_at_gp[gp],
                qy_ex, boundary.bath_at_gp[gp], bound.surface_normal[gp],
                boundary.ze_numerical_flux_at_gp[gp], boundary.qx_numerical_flux_at_gp[gp],
                boundary.qy_numerical_flux_at_gp[gp]);
    }
}

template <typename BoundaryType>
void cuda_boundary_kernel2(uint block_dim, uint thread_dim, const Stepper& stepper,
        BoundaryType& bound, uint n_dof) {
    std::cout << "running" << std::endl;
    cuda_boundary_kernel2<BoundaryType><<<block_dim, thread_dim>>>(stepper, bound, n_dof);
}

template <typename BoundaryType>
__global__ void cuda_boundary_kernel2(const Stepper& stepper, BoundaryType& bound, uint n_dof) {
    uint dof = blockIdx.x * blockDim.x + threadIdx.x;
    auto& boundary = bound.data.boundary[bound.bound_id];
    const uint stage = stepper.get_stage();
    auto& state = bound.data.state[stage];
    if (dof < n_dof) {
        state.rhs_ze[dof] -= bound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
        state.rhs_qx[dof] -= bound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
        state.rhs_qy[dof] -= bound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
    }
}

template <typename ElementType>
void cuda_update_kernel(uint block_dim, uint thread_dim, const Stepper& stepper, ElementType& elt,
        uint stage, uint n_dof) {
    std::cout << "running" << std::endl;
    cuda_update_kernel<ElementType><<<block_dim, thread_dim>>>(stepper, elt, stage, n_dof);
}

template <typename ElementType>
__global__ void cuda_update_kernel(const Stepper& stepper, ElementType& elt, uint stage,
        uint n_dof) {
    uint dof = blockIdx.x * blockDim.x + threadIdx.x;
    auto& state = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];
    double dt = stepper.get_dt();
    if (dof < n_dof) {
        for (uint s = 0; s <= stage; ++s) {
            next_state.ze[dof] += stepper.ark[stage][s] * state[s].ze[dof] + dt *
                    stepper.brk[stage][s] * state[s].rhs_ze[dof];
            next_state.qx[dof] += stepper.ark[stage][s] * state[s].qx[dof] + dt *
                    stepper.brk[stage][s] * state[s].rhs_qx[dof];
            next_state.qy[dof] += stepper.ark[stage][s] * state[s].qy[dof] + dt *
                    stepper.brk[stage][s] * state[s].rhs_qy[dof];
        }
    }
}

