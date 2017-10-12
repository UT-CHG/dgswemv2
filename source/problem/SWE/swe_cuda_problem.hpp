/*
 * swe_cuda_problem.hpp
 *
 *  Created on: October 8, 2017
 *      Author: vance
 */

#ifndef SWE_CUDA_PROBLEM_HPP_
#define SWE_CUDA_PROBLEM_HPP_

#include "problem/SWE/swe_problem.hpp"

namespace SWE {

struct CUDAProblem : public Problem {
    // processor kernels
    template <typename ElementType>
    static void volume_kernel(const Stepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void source_kernel(const Stepper& stepper, ElementType& elt);

    template <typename InterfaceType>
    static void interface_kernel(const Stepper& stepper, InterfaceType& intface);

    template <typename BoundaryType>
    static void boundary_kernel(const Stepper& stepper, BoundaryType& bound);

    template <typename DistributedBoundaryType>
    static void distributed_boundary_send_kernel(const Stepper& stepper, DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    static void distributed_boundary_kernel(const Stepper& stepper, DistributedBoundaryType& dbound);

    template <typename ElementType>
    static void update_kernel(const Stepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void swap_states_kernel(const Stepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt);
};

}

#endif
