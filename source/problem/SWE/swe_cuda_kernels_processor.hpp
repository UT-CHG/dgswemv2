/*
 * swe_cuda_kernels_processor.hpp
 *
 *  Created on: Sep 28, 2017
 *      Author: vance
 */

#ifndef SWE_CUDA_KERNELS_PROCESSOR_HPP_
#define SWE_CUDA_KERNELS_PROCESSOR_HPP_

#include "swe_initial_conditions_function.hpp"
#include "swe_cuda_kernels.cuh"
#include "swe_cuda_problem.hpp"
#include "swe_LLF_flux.hpp"

#define THREADS_PER_BLOCK 1024
namespace SWE {

template <typename ElementType>
void CUDAProblem::volume_kernel(const Stepper& stepper, ElementType& elt) {
    std::cout << "volume" << std::endl;
    const uint stage = stepper.get_stage();
    auto& state = elt.data.state[stage];
    auto& internal = elt.data.internal;

    // get state at Gauss points
    elt.ComputeUgp(state.ze, internal.ze_at_gp);
    elt.ComputeUgp(state.qx, internal.qx_at_gp);
    elt.ComputeUgp(state.qy, internal.qy_at_gp);

    // assemble flux
    uint blocksPerGrid = (elt.data.get_ngp_internal() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    cuda_volume_kernel1<ElementType><<<blocksPerGrid, THREADS_PER_BLOCK>>>(
            internal_args(elt.data.internal), elt.data.get_ngp_internal());
    blocksPerGrid = (elt.data.get_ndof() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    cuda_volume_kernel2<ElementType><<<blocksPerGrid, THREADS_PER_BLOCK>>>(
            internal_args(elt.data.internal), state_args(elt.data.state[stage]), NULL, //elt.int_fact_dphi,
            elt.data.get_ndof());
}

template <typename ElementType>
void CUDAProblem::source_kernel(const Stepper& stepper, ElementType& elt) {
    std::cout << "source" << std::endl;
    const uint stage = stepper.get_stage();
    auto& state = elt.data.state[stage];
    auto& internal = elt.data.internal;
    double t = stepper.get_t_at_curr_stage();
    auto _source_ze = [&, t](Point<2>& pt) { return source_ze(t, pt); };
    auto _source_qx = [&, t](Point<2>& pt) { return source_qx(t, pt); };
    auto _source_qy = [&, t](Point<2>& pt) { return source_qy(t, pt); };

    elt.ComputeFgp(_source_ze, internal.ze_source_term_at_gp);
    elt.ComputeFgp(_source_qx, internal.qx_source_term_at_gp);
    elt.ComputeFgp(_source_qy, internal.qy_source_term_at_gp);

    // note we assume that the values at gauss points have already been computed

    uint blocksPerGrid = (elt.data.get_ngp_internal() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_source_kernel1<ElementType>(blocksPerGrid, THREADS_PER_BLOCK, elt,
            elt.data.get_ngp_internal());*/
    blocksPerGrid = (elt.data.get_ndof() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_source_kernel2<ElementType>(blocksPerGrid, THREADS_PER_BLOCK, elt, stage,
            elt.data.get_ndof());*/
}

template <typename InterfaceType>
void CUDAProblem::interface_kernel(const Stepper& stepper, InterfaceType& interface) {
    std::cout << "interface" << std::endl;
    const uint stage = stepper.get_stage();
    auto& state_in = interface.data_in.state[stage];
    auto& state_ex = interface.data_ex.state[stage];
    auto& boundary_in = interface.data_in.boundary[interface.bound_id_in];
    auto& boundary_ex = interface.data_ex.boundary[interface.bound_id_ex];
    interface.ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
    interface.ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
    interface.ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

    interface.ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
    interface.ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
    interface.ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);

    // assemble numerical fluxes
    uint blocksPerGrid = (interface.data_in.get_ngp_boundary(interface.bound_id_in) +
            THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_interface_kernel1<InterfaceType>(blocksPerGrid, THREADS_PER_BLOCK, interface,
            interface.data_in.get_ngp_boundary(interface.bound_id_in));*/

    // now compute contributions to the righthand side
    blocksPerGrid = (interface.data_in.get_ndof() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_interface_kernel2<InterfaceType>(blocksPerGrid, THREADS_PER_BLOCK, interface, stage,
            interface.data_in.get_ndof());*/
    blocksPerGrid = (interface.data_ex.get_ndof() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_interface_kernel3<InterfaceType>(blocksPerGrid, THREADS_PER_BLOCK, interface, stage,
            interface.data_ex.get_ndof());*/
}

template <typename BoundaryType>
void CUDAProblem::boundary_kernel(const Stepper& stepper, BoundaryType& bound) {
    std::cout << "boundary" << std::endl;
    const uint stage = stepper.get_stage();

    auto& state = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    bound.ComputeUgp(state.ze, boundary.ze_at_gp);
    bound.ComputeUgp(state.qx, boundary.qx_at_gp);
    bound.ComputeUgp(state.qy, boundary.qy_at_gp);

    uint blocksPerGrid = (bound.data.get_ngp_boundary(bound.bound_id) + THREADS_PER_BLOCK - 1) /
            THREADS_PER_BLOCK;
    /*cuda_boundary_kernel1<BoundaryType>(blocksPerGrid, THREADS_PER_BLOCK, stepper, bound,
            bound.data.get_ngp_boundary(bound.bound_id));*/

    // now compute contributions to the righthand side
    blocksPerGrid = (bound.data.get_ndof() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_boundary_kernel2<BoundaryType>(blocksPerGrid, THREADS_PER_BLOCK, stepper, bound,
            bound.data.get_ndof());*/
}

template <typename ElementType>
void CUDAProblem::update_kernel(const Stepper& stepper, ElementType& elt) {
    std::cout << "update" << std::endl;
    const uint stage = stepper.get_stage();

    auto& state = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    curr_state.rhs_ze = elt.SolveLSE(curr_state.rhs_ze);
    curr_state.rhs_qx = elt.SolveLSE(curr_state.rhs_qx);
    curr_state.rhs_qy = elt.SolveLSE(curr_state.rhs_qy);

    std::fill(next_state.ze.begin(), next_state.ze.end(), 0);
    std::fill(next_state.qx.begin(), next_state.qx.end(), 0);
    std::fill(next_state.qy.begin(), next_state.qy.end(), 0);

    uint blocksPerGrid = (elt.data.get_ndof() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    /*cuda_update_kernel<ElementType>(blocksPerGrid, THREADS_PER_BLOCK, stepper, elt, stage,
            elt.data.get_ndof());*/
}

template <typename ElementType>
void CUDAProblem::swap_states_kernel(const Stepper& stepper, ElementType& elt) {
    std::cout << "swap" << std::endl;
    uint n_stages = stepper.get_num_stages();
    auto& state = elt.data.state;

    std::swap(state[0].ze, state[n_stages].ze);
    std::swap(state[0].qx, state[n_stages].qx);
    std::swap(state[0].qy, state[n_stages].qy);
}

template <typename ElementType>
void CUDAProblem::scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt) {
    std::cout << "scrutinize" << std::endl;
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



#endif /* SWE_CUDA_KERNELS_PROCESSOR_HPP_ */
