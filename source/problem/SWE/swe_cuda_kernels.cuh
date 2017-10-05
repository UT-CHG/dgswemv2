/*
 * swe_cuda_kernels.cuh
 *
 *  Created on: Sep 29, 2017
 *      Author: vance
 */

#ifndef SWE_CUDA_KERNELS_CUH_
#define SWE_CUDA_KERNELS_CUH_
// Forward declarations
class Stepper;

template <typename ElementType>
void cuda_volume_kernel1(uint block_dim, uint thread_dim, ElementType& elt, uint n_gp);

template <typename ElementType>
void cuda_volume_kernel2(uint block_dim, uint thread_dim, ElementType& elt, uint stage, uint n_dof);

template <typename ElementType>
void cuda_source_kernel1(uint block_dim, uint thread_dim, ElementType& elt, uint n_gp);

template <typename ElementType>
void cuda_source_kernel2(uint block_dim, uint thread_dim, ElementType& elt, uint stage, uint n_dof);

template <typename InterfaceType>
void cuda_interface_kernel1(uint block_dim, uint thread_dim, InterfaceType& interface, uint n_gp);

template <typename InterfaceType>
void cuda_interface_kernel2(uint block_dim, uint thread_dim, InterfaceType& interface, uint stage,
       uint n_dof);

template <typename InterfaceType>
void cuda_interface_kernel3(uint block_dim, uint thread_dim, InterfaceType& interface, uint stage,
        uint n_dof);

template <typename BoundaryType>
void cuda_boundary_kernel1(uint block_dim, uint thread_dim, const Stepper& stepper,
        BoundaryType& bound, uint n_gp);

template <typename BoundaryType>
void cuda_boundary_kernel2(uint block_dim, uint thread_dim, const Stepper& stepper,
        BoundaryType& bound, uint n_dof);

template <typename ElementType>
void cuda_update_kernel(uint block_dim, uint thread_dim, const Stepper& stepper, ElementType& elt,
        uint stage, uint n_dof);

#endif /* SWE_CUDA_KERNELS_CUH_ */
