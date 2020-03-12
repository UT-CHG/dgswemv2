#ifndef EHDG_GN_PROBLEM_HPP
#define EHDG_GN_PROBLEM_HPP

#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"
#include "simulation/stepper/second_strang_stepper.hpp"
#include "simulation/writer.hpp"
#include "simulation/discretization.hpp"

#include "problem/Green-Naghdi/gn_definitions.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/stabilization_parameters/ehdg_gn_stabilization_params.hpp"

#include "boundary_conditions/ehdg_gn_boundary_conditions.hpp"
#include "dist_boundary_conditions/ehdg_gn_distributed_boundary_conditions.hpp"
#include "interface_specializations/ehdg_gn_interface_specializations.hpp"

#include "problem/Green-Naghdi/problem_data_structure/gn_data.hpp"
#include "problem/Green-Naghdi/problem_data_structure/gn_edge_data.hpp"
#include "problem/Green-Naghdi/problem_data_structure/gn_global_data.hpp"

#include "problem/Green-Naghdi/problem_input/gn_inputs.hpp"
#include "problem/Green-Naghdi/problem_parser/gn_parser.hpp"

#include "geometry/mesh_definitions.hpp"
#include "geometry/mesh_skeleton_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_intface.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_bound.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_dbound.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_edge_intface.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_edge_bound.hpp"
#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_create_edge_dbound.hpp"

namespace GN {
namespace EHDG {
struct Problem {
    using ProblemInputType   = GN::Inputs;
    using ProblemStepperType = SecondStrangStepper<SWE_SIM::Problem::ProblemStepperType, ESSPRKStepper>;
    using ProblemWriterType  = Writer<Problem>;
    using ProblemParserType  = GN::Parser;

    using ProblemDataType       = GN::Data;
    using ProblemEdgeDataType   = GN::EdgeData;
    using ProblemGlobalDataType = GN::GlobalData;

    using ProblemInterfaceTypes = Geometry::InterfaceTypeTuple<Data, ISP::Internal, ISP::Levee>;
    using ProblemBoundaryTypes  = Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow, BC::Function>;
    using ProblemDistributedBoundaryTypes =
        Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed, DBC::DistributedLevee>;

    using ProblemEdgeInterfaceTypes = Geometry::EdgeInterfaceTypeTuple<EdgeData, ProblemInterfaceTypes>::Type;
    using ProblemEdgeBoundaryTypes  = Geometry::EdgeBoundaryTypeTuple<EdgeData, ProblemBoundaryTypes>::Type;
    using ProblemEdgeDistributedTypes =
        Geometry::EdgeDistributedTypeTuple<EdgeData, ProblemDistributedBoundaryTypes>::Type;

    using ProblemMeshType = Geometry::MeshType<Data,
                                               std::tuple<ISP::Internal, ISP::Levee>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow, BC::Function>,
                                               std::tuple<DBC::Distributed, DBC::DistributedLevee>>::Type;

    using ProblemMeshSkeletonType = Geometry::MeshSkeletonType<
        EdgeData,
        Geometry::InterfaceTypeTuple<Data, ISP::Internal, ISP::Levee>,
        Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow, BC::Function>,
        Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed, DBC::DistributedLevee>>::Type;

    using ProblemDiscretizationType = HDGDiscretization<Problem>;

    // preprocessor kernels
    static void initialize_problem_parameters(const ProblemInputType& problem_specific_input) {
        SWE::initialize_problem_parameters(problem_specific_input);

        GN::Global::g         = problem_specific_input.g;
        GN::Global::rho_air   = problem_specific_input.rho_air;
        GN::Global::rho_water = problem_specific_input.rho_water;
    }

    static void preprocess_mesh_data(InputParameters<ProblemInputType>& input) { SWE::preprocess_mesh_data(input); }

    // helpers to create communication
    static constexpr uint n_communications = GN::EHDG::n_communications;

    static std::vector<uint> comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp) {
        std::vector<uint> offset = SWE_SIM::Problem::comm_buffer_offsets(begin_index, ngp);

        offset.resize(Problem::n_communications);

        offset[CommTypes::dc_global_dof_indx] = begin_index[CommTypes::dc_global_dof_indx];
        offset[CommTypes::dbath]              = begin_index[CommTypes::dbath];
        offset[CommTypes::derivatives]        = begin_index[CommTypes::derivatives];

        begin_index[CommTypes::dc_global_dof_indx] += 1;
        begin_index[CommTypes::dbath] += GN::n_dddbath_terms * ngp;
        begin_index[CommTypes::derivatives] += GN::n_ddu_terms * ngp;

        return offset;
    }

    template <typename RawBoundaryType>
    static void create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& input,
                                  ProblemWriterType& writer) {
        GN::create_interfaces<GN::EHDG::Problem>(raw_boundaries, mesh, input, writer);
    }

    template <typename RawBoundaryType>
    static void create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& input,
                                  ProblemWriterType& writer) {
        GN::create_boundaries<GN::EHDG::Problem>(raw_boundaries, mesh, input, writer);
    }

    template <typename RawBoundaryType>
    static void create_distributed_boundaries(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType&,
        ProblemInputType& problem_input,
        std::tuple<>&,
        ProblemWriterType&) {}

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        ProblemInputType& input,
        Communicator& communicator,
        ProblemWriterType& writer) {
        // *** //
        GN::create_distributed_boundaries<GN::EHDG::Problem>(raw_boundaries, mesh, input, communicator, writer);
    }

    static void create_edge_interfaces(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
        GN::create_edge_interfaces<GN::EHDG::Problem>(mesh, mesh_skeleton, writer);
    }

    static void create_edge_boundaries(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
        GN::create_edge_boundaries<GN::EHDG::Problem>(mesh, mesh_skeleton, writer);
    }

    static void create_edge_distributeds(ProblemMeshType& mesh,
                                         ProblemMeshSkeletonType& mesh_skeleton,
                                         ProblemWriterType& writer) {
        GN::create_edge_distributeds<GN::EHDG::Problem>(mesh, mesh_skeleton, writer);
    }

    static void preprocessor_serial(ProblemDiscretizationType& discretization,
                                    ProblemGlobalDataType& global_data,
                                    const ProblemStepperType& stepper,
                                    const ProblemInputType& problem_specific_input);

    template <typename OMPISimUnitType>
    static void preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                  ProblemGlobalDataType& global_data,
                                  const ProblemStepperType& stepper,
                                  const uint begin_sim_id,
                                  const uint end_sim_id);

    static void initialize_global_dc_problem_serial(ProblemDiscretizationType& discretization,
                                                    uint& dc_global_dof_offset);

    static void initialize_global_dc_problem_parallel_pre_send(ProblemDiscretizationType& discretization,
                                                               uint& dc_global_dof_offset);

    static void initialize_global_dc_problem_parallel_finalize_pre_send(ProblemDiscretizationType& discretization,
                                                                        uint dc_global_dof_offset);

    static void initialize_global_dc_problem_parallel_post_receive(ProblemDiscretizationType& discretization,
                                                                   std::vector<uint>& dc_global_dof_indx);

    static void compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization,
                                                      ProblemGlobalDataType& global_data);

    template <typename OMPISimUnitType>
    static void compute_bathymetry_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                                    ProblemGlobalDataType& global_data,
                                                    const uint begin_sim_id,
                                                    const uint end_sim_id);

    // processor kernels
    static void step_serial(ProblemDiscretizationType& discretization,
                            ProblemGlobalDataType& global_data,
                            ProblemStepperType& stepper,
                            ProblemWriterType& writer,
                            ProblemParserType& parser);

    template <typename OMPISimUnitType>
    static void step_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                          ProblemGlobalDataType& global_data,
                          ProblemStepperType& stepper,
                          const uint begin_sim_id,
                          const uint end_sim_id);

    /* Dispersive correction part */

    static void dispersive_correction_serial(ProblemDiscretizationType& discretization,
                                             ProblemGlobalDataType& global_data,
                                             ESSPRKStepper& stepper);

    static void compute_derivatives_serial(ProblemDiscretizationType& discretization,
                                           ProblemGlobalDataType& global_data,
                                           const ESSPRKStepper& stepper);

    template <typename OMPISimUnitType>
    static void dispersive_correction_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                           ProblemGlobalDataType& global_data,
                                           ESSPRKStepper& stepper,
                                           const uint begin_sim_id,
                                           const uint end_sim_id);

    template <typename OMPISimUnitType>
    static void compute_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                         ProblemGlobalDataType& global_data,
                                         const ESSPRKStepper& stepper,
                                         const uint begin_sim_id,
                                         const uint end_sim_id);

    /* local step */

    template <typename ElementType>
    static void local_dc_volume_kernel(const ESSPRKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void local_dc_source_kernel(const ESSPRKStepper& stepper, ElementType& elt);

    template <typename EdgeInterfaceType>
    static void local_dc_edge_interface_kernel(const ESSPRKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeBoundaryType>
    static void local_dc_edge_boundary_kernel(const ESSPRKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeDistributedType>
    static void local_dc_edge_distributed_kernel(const ESSPRKStepper& stepper, EdgeDistributedType& edge_dbound);

    /* global step */

    template <typename EdgeInterfaceType>
    static void global_dc_edge_interface_kernel(const ESSPRKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeBoundaryType>
    static void global_dc_edge_boundary_kernel(const ESSPRKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeDistributedType>
    static void global_dc_edge_distributed_kernel(const ESSPRKStepper& stepper, EdgeDistributedType& edge_dbound);

    static void serial_solve_global_dc_problem(ProblemDiscretizationType& discretization,
                                               ProblemGlobalDataType& global_data,
                                               const ESSPRKStepper& stepper);

    template <typename OMPISimUnitType>
    static void ompi_solve_global_dc_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                             ProblemGlobalDataType& global_data,
                                             const ESSPRKStepper& stepper,
                                             const uint begin_sim_id,
                                             const uint end_sim_id);
    template <typename OMPISimUnitType>
    static void reset_PETSC_solver(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                   ProblemGlobalDataType& global_data,
                                   const ESSPRKStepper& stepper,
                                   const uint begin_sim_id,
                                   const uint end_sim_id);

    template <typename ElementType>
    static void dispersive_correction_kernel(const ESSPRKStepper& stepper, ElementType& elt);

    // writing output kernels
    static void write_VTK_data(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
        return SWE::write_VTK_data(mesh, raw_data_file);
    }

    static void write_VTU_data(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
        return SWE::write_VTU_data(mesh, raw_data_file);
    }

    static void write_modal_data(const ProblemStepperType& stepper,
                                 ProblemMeshType& mesh,
                                 const std::string& output_path) {
        return SWE::write_modal_data(stepper, mesh, output_path);
    }

    template <typename ElementType>
    static double compute_residual_L2(const ProblemStepperType& stepper, ElementType& elt) {
        return SWE::compute_residual_L2(stepper, elt);
    }

    template <typename GlobalDataType>
    static void finalize_simulation(GlobalDataType& global_data) {
#ifdef HAS_PETSC
        global_data.destroy();
#endif
    }
};
}
}

#endif
