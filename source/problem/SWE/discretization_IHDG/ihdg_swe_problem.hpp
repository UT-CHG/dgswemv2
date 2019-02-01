#ifndef IHDG_SWE_PROBLEM_HPP
#define IHDG_SWE_PROBLEM_HPP

#include "problem/SWE/swe_definitions.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

#include "simulation/stepper/implicit_stepper.hpp"
#include "simulation/writer.hpp"
#include "simulation/discretization.hpp"

#include "problem/SWE/problem_stabilization_parameters/swe_stabilization_params.hpp"

#include "boundary_conditions/ihdg_swe_boundary_conditions.hpp"
#include "dist_boundary_conditions/ihdg_swe_distributed_boundary_conditions.hpp"
#include "interface_specializations/ihdg_swe_interface_specializations.hpp"

#include "data_structure/ihdg_swe_data.hpp"
#include "data_structure/ihdg_swe_edge_data.hpp"
#include "data_structure/ihdg_swe_global_data.hpp"

#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"
#include "problem/SWE/problem_preprocessor/swe_preprocessor.hpp"
#include "problem/SWE/problem_postprocessor/swe_postprocessor.hpp"

#include "geometry/mesh_definitions.hpp"
#include "geometry/mesh_skeleton_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

namespace SWE {
namespace IHDG {
struct Problem {
    using ProblemInputType   = SWE::Inputs;
    using ProblemStepperType = ImplicitStepper;
    using ProblemWriterType  = Writer<Problem>;
    using ProblemParserType  = SWE::Parser;

    using ProblemDataType       = Data;
    using ProblemEdgeDataType   = EdgeData;
    using ProblemGlobalDataType = GlobalData;

    using ProblemInterfaceTypes = Geometry::InterfaceTypeTuple<Data, ISP::Internal, ISP::Levee>;
    using ProblemBoundaryTypes =
        Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow, BC::Function, BC::Outflow>;
    using ProblemDistributedBoundaryTypes =
        Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed, DBC::DistributedLevee>;

    using ProblemEdgeInterfaceTypes = Geometry::EdgeInterfaceTypeTuple<EdgeData, ProblemInterfaceTypes>::Type;
    using ProblemEdgeBoundaryTypes  = Geometry::EdgeBoundaryTypeTuple<EdgeData, ProblemBoundaryTypes>::Type;
    using ProblemEdgeDistributedTypes =
        Geometry::EdgeDistributedTypeTuple<EdgeData, ProblemDistributedBoundaryTypes>::Type;

    using ProblemMeshType = Geometry::MeshType<Data,
                                               std::tuple<ISP::Internal, ISP::Levee>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow, BC::Function, BC::Outflow>,
                                               std::tuple<DBC::Distributed, DBC::DistributedLevee>>::Type;

    using ProblemMeshSkeletonType = Geometry::MeshSkeletonType<
        EdgeData,
        Geometry::InterfaceTypeTuple<Data, ISP::Internal, ISP::Levee>,
        Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow, BC::Function, BC::Outflow>,
        Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed, DBC::DistributedLevee>>::Type;

    using ProblemDiscretizationType = HDGDiscretization<Problem>;

    // preprocessor kernels
    static void initialize_problem_parameters(ProblemInputType& problem_specific_input) {
        SWE::initialize_problem_parameters(problem_specific_input);
    }

    static void preprocess_mesh_data(InputParameters<ProblemInputType>& input) { SWE::preprocess_mesh_data(input); }

    template <typename RawBoundaryType>
    static void create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& problem_input,
                                  ProblemWriterType& writer);

    template <typename RawBoundaryType>
    static void create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& problem_input,
                                  ProblemWriterType& writer);

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
        ProblemWriterType& writer);

    // helpers to create communication
    static constexpr uint n_communications = SWE::IHDG::n_communications;
    static std::vector<uint> comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp);

    static void create_edge_interfaces(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer);

    static void create_edge_boundaries(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer);

    static void create_edge_distributeds(ProblemMeshType& mesh,
                                         ProblemMeshSkeletonType& mesh_skeleton,
                                         ProblemWriterType& writer);

    template <typename ProblemType>
    static void preprocessor_serial(HDGDiscretization<ProblemType>& discretization,
                                    const ProblemInputType& problem_specific_input);

    template <typename OMPISimUnitType>
    static void preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                  uint begin_sim_id,
                                  uint end_sim_id);

    template <typename ProblemType>
    static void initialize_global_problem_serial(HDGDiscretization<ProblemType>& discretization,
                                                 uint& global_dof_offset);

    template <typename ProblemType>
    static void initialize_global_problem_parallel_pre_send(HDGDiscretization<ProblemType>& discretization,
                                                            uint& global_dof_offset);

    template <typename ProblemType>
    static void initialize_global_problem_parallel_finalize_pre_send(HDGDiscretization<ProblemType>& discretization,
                                                                     uint global_dof_offset);

    template <typename ProblemType>
    static void initialize_global_problem_parallel_post_receive(HDGDiscretization<ProblemType>& discretization,
                                                                std::vector<uint>& global_dof_indx);

    // processor kernels
    template <typename SerialSimType>
    static void step_serial(SerialSimType* sim);

    template <typename StepperType, typename ProblemType>
    static void stage_serial(StepperType& stepper, HDGDiscretization<ProblemType>& discretization);

    template <typename OMPISimType>
    static void step_ompi(OMPISimType* sim, uint begin_sim_id, uint end_sim_id);

    template <typename OMPISimUnitType>
    static void stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                           uint begin_sim_id,
                           uint end_sim_id);

    /* init interation begin */

    template <typename StepperType, typename ProblemType>
    static void init_iteration(const StepperType& stepper, HDGDiscretization<ProblemType>& discretization);

    template <typename StepperType, typename ElementType>
    static void init_volume_kernel(const StepperType& stepper, ElementType& elt);

    template <typename StepperType, typename ElementType>
    static void init_source_kernel(const StepperType& stepper, ElementType& elt);

    template <typename StepperType, typename InterfaceType>
    static void init_interface_kernel(const StepperType& stepper, InterfaceType& intface);

    template <typename StepperType, typename BoundaryType>
    static void init_boundary_kernel(const StepperType& stepper, BoundaryType& bound);

    template <typename StepperType, typename DistributedBoundaryType>
    static void init_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound);

    template <typename StepperType, typename EdgeInterfaceType>
    static void init_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int);

    template <typename StepperType, typename EdgeBoundaryType>
    static void init_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeDistributedType>
    static void init_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound);

    /* init interation end */

    /* local step begin */

    template <typename StepperType, typename ElementType>
    static void local_volume_kernel(const StepperType& stepper, ElementType& elt);

    template <typename StepperType, typename ElementType>
    static void local_source_kernel(const StepperType& stepper, ElementType& elt);

    template <typename StepperType, typename InterfaceType>
    static void local_interface_kernel(const StepperType& stepper, InterfaceType& intface);

    template <typename StepperType, typename BoundaryType>
    static void local_boundary_kernel(const StepperType& stepper, BoundaryType& bound);

    template <typename StepperType, typename DistributedBoundaryType>
    static void local_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound);

    template <typename StepperType, typename EdgeInterfaceType>
    static void local_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int);

    template <typename StepperType, typename EdgeBoundaryType>
    static void local_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeDistributedType>
    static void local_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound);

    /* local step end */

    /* global step begin */

    template <typename StepperType, typename EdgeInterfaceType>
    static void global_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int);

    template <typename StepperType, typename EdgeBoundaryType>
    static void global_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeDistributedType>
    static void global_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_bound);

    template <typename StepperType, typename ProblemType>
    static bool serial_solve_global_problem(const StepperType& stepper, HDGDiscretization<ProblemType>& discretization);

    template <typename OMPISimUnitType>
    static bool ompi_solve_global_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                          uint begin_sim_id,
                                          uint end_sim_id);

    /* global step end */

    // writing output kernels
    template <typename MeshType>
    static void write_VTK_data(MeshType& mesh, std::ofstream& raw_data_file) {
        SWE::write_VTK_data(mesh, raw_data_file);
    }

    template <typename MeshType>
    static void write_VTU_data(MeshType& mesh, std::ofstream& raw_data_file) {
        SWE::write_VTU_data(mesh, raw_data_file);
    }

    template <typename StepperType, typename MeshType>
    static void write_modal_data(const StepperType& stepper, MeshType& mesh, const std::string& output_path) {
        SWE::write_modal_data(stepper, mesh, output_path);
    }

    template <typename StepperType, typename ElementType>
    static double compute_residual_L2(const StepperType& stepper, ElementType& elt) {
        return SWE::compute_residual_L2(stepper, elt);
    }

    template <typename SimType>
    static void finalize_simulation(SimType* sim) {
#ifdef HAS_PETSC
        // Here one assumes that there is at lease one sim unit present
        // This is of course not always true
        if (!sim->sim_units.empty())
            sim->sim_units[0]->discretization.global_data.destroy();
#endif
    }
};
}
}

#endif
