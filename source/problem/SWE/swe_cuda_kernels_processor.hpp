/*
 * swe_cuda_kernels_processor.hpp
 *
 *  Created on: Sep 28, 2017
 *      Author: vance
 */

#ifndef SWE_CUDA_KERNELS_PROCESSOR_HPP_
#define SWE_CUDA_KERNELS_PROCESSOR_HPP_

#include "swe_problem.hpp"
#include "swe_LLF_flux.hpp"

#define THREADS_PER_BLOCK 1024
namespace SWE {

struct CUDAProblem : public SWE::Problem {
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

    // postprocessor kernels
    template <typename ElementType>
    static void extract_VTK_data_kernel(ElementType& elt, Array2D<double>& cell_data, Array2D<double>& point_data);

    template <typename MeshType>
    static void write_VTK_data_kernel(const Stepper& stepper, MeshType& mesh);

    template <typename ElementType>
    static void extract_modal_data_kernel(ElementType& elt, std::vector<std::pair<uint, Array2D<double>>>& modal_data);

    template <typename MeshType>
    static void write_modal_data_kernel(const Stepper& stepper, MeshType& mesh);
};
}



#endif /* SWE_CUDA_KERNELS_PROCESSOR_HPP_ */
