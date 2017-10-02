#ifndef SWE_MANUFACTURED_PROBLEM_HPP
#define SWE_MANUFACTURED_PROBLEM_HPP

#include "general_definitions.hpp"

#include "problem/SWE/swe_problem.hpp"
#include "problem/SWE/swe_kernels_preprocessor.hpp"
#include "problem/SWE/swe_kernels_processor.hpp"
#include "problem/SWE/swe_kernels_postprocessor.hpp"

namespace SWE {
struct ManufacturedProblem : Problem {
    /*  typedef SWE::Data ProblemDataType;

    typedef Problem::ProblemMeshType ProblemMeshType;

      // preprocessor kernels
      template <typename RawBoundaryType>
      static void create_boundaries_kernel(ProblemMeshType& mesh,
                                           std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries);

      template <typename RawBoundaryType>
      static void create_distributed_boundaries_kernel(ProblemMeshType&,
                                                       std::tuple<>&,
                                                       std::map<uint, std::map<uint, RawBoundaryType>>&);

      template <typename RawBoundaryType, typename Communicator>
      static void create_distributed_boundaries_kernel(
          ProblemMeshType& mesh,
          Communicator& communicator,
          std::map<uint, std::map<uint, RawBoundaryType>>& pre_distributed_boundaries);

      static void initialize_data_kernel(ProblemMeshType& mesh, const MeshMetaData& mesh_data);

      // processor kernels
      template <typename ElementType>
      static void volume_kernel(const Stepper& stepper, ElementType& elt);
    */
    template <typename ElementType>
    static void source_kernel(const Stepper& stepper, ElementType& elt);
    /*
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
      static void extract_modal_data_kernel(ElementType& elt, std::vector<std::pair<uint, Array2D<double>>>&
      modal_data);

      template <typename MeshType>
      static void write_modal_data_kernel(const Stepper& stepper, MeshType& mesh);*/
};

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
}

#endif
