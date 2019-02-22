#ifndef EHDG_SWE_PROBLEM_HPP
#define EHDG_SWE_PROBLEM_HPP

#include "problem/SWE/swe_definitions.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"
#include "simulation/writer.hpp"
#include "simulation/discretization.hpp"

#include "problem/SWE/problem_stabilization_parameters/swe_stabilization_params.hpp"

#include "boundary_conditions/ehdg_swe_boundary_conditions.hpp"
#include "dist_boundary_conditions/ehdg_swe_distributed_boundary_conditions.hpp"
#include "interface_specializations/ehdg_swe_interface_specializations.hpp"

#include "data_structure/ehdg_swe_data.hpp"
#include "data_structure/ehdg_swe_edge_data.hpp"
#include "data_structure/ehdg_swe_global_data.hpp"

#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"
#include "problem/SWE/problem_preprocessor/swe_preprocessor.hpp"
#include "problem/SWE/problem_postprocessor/swe_postprocessor.hpp"

#include "geometry/mesh_definitions.hpp"
#include "geometry/mesh_skeleton_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "problem/SWE/problem_preprocessor/swe_pre_create_intface.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_bound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_dbound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_edge_intface.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_edge_bound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_edge_dbound.hpp"

namespace SWE {
namespace EHDG {
struct Problem {
    using ProblemInputType   = SWE::Inputs;
    using ProblemStepperType = ESSPRKStepper;
    using ProblemWriterType  = Writer<Problem>;
    using ProblemParserType  = SWE::Parser;

    using ProblemDataType       = Data;
    using ProblemEdgeDataType   = EdgeData;
    using ProblemGlobalDataType = std::tuple<>;

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

    // helpers to create communication
    static constexpr uint n_communications = SWE::EHDG::n_communications;

    static std::vector<uint> comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp) {
        std::vector<uint> offset(SWE::EHDG::n_communications);

        offset[CommTypes::baryctr_coord]    = begin_index[CommTypes::baryctr_coord];
        offset[CommTypes::init_global_prob] = begin_index[CommTypes::init_global_prob];
        offset[CommTypes::bound_state]      = begin_index[CommTypes::bound_state];
        offset[CommTypes::baryctr_state]    = begin_index[CommTypes::baryctr_state];

        begin_index[CommTypes::baryctr_coord] += 2;
        begin_index[CommTypes::init_global_prob] += ngp;  // bath_at_gp
        begin_index[CommTypes::bound_state] += SWE::n_variables * ngp;
        begin_index[CommTypes::baryctr_state] += SWE::n_variables + 1;  // + w/d state

        return offset;
    }

    template <typename RawBoundaryType>
    static void create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& input,
                                  ProblemWriterType& writer) {
        SWE::create_interfaces<SWE::EHDG::Problem>(raw_boundaries, mesh, input, writer);
    }

    template <typename RawBoundaryType>
    static void create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& input,
                                  ProblemWriterType& writer) {
        SWE::create_boundaries<SWE::EHDG::Problem>(raw_boundaries, mesh, input, writer);
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
        SWE::create_distributed_boundaries<SWE::EHDG::Problem>(raw_boundaries, mesh, input, communicator, writer);
    }

    static void create_edge_interfaces(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
        SWE::create_edge_interfaces<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);
    }

    static void create_edge_boundaries(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
        SWE::create_edge_boundaries<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);
    }

    static void create_edge_distributeds(ProblemMeshType& mesh,
                                         ProblemMeshSkeletonType& mesh_skeleton,
                                         ProblemWriterType& writer) {
        SWE::create_edge_distributeds<SWE::EHDG::Problem>(mesh, mesh_skeleton, writer);
    }

    template <typename ProblemType>
    static void preprocessor_serial(HDGDiscretization<ProblemType>& discretization,
                                    typename ProblemType::ProblemGlobalDataType& global_data,
                                    const typename ProblemType::ProblemStepperType& stepper,
                                    const typename ProblemType::ProblemInputType& problem_specific_input);

    template <template <typename> typename OMPISimUnitType, typename ProblemType>
    static void preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                                  typename ProblemType::ProblemGlobalDataType& global_data,
                                  const typename ProblemType::ProblemStepperType& stepper,
                                  const uint begin_sim_id,
                                  const uint end_sim_id);

    template <typename HPXSimUnitType>
    static auto preprocessor_hpx(HPXSimUnitType* sim_unit);

    template <typename ProblemType>
    static void initialize_global_problem_serial(HDGDiscretization<ProblemType>& discretization);

    template <typename ProblemType>
    static void initialize_global_problem_parallel_pre_send(HDGDiscretization<ProblemType>& discretization);

    template <typename ProblemType>
    static void initialize_global_problem_parallel_post_receive(HDGDiscretization<ProblemType>& discretization);

    // processor kernels
    template <typename ProblemType>
    static void step_serial(HDGDiscretization<ProblemType>& discretization,
                            typename ProblemType::ProblemGlobalDataType& global_data,
                            typename ProblemType::ProblemStepperType& stepper,
                            typename ProblemType::ProblemWriterType& writer,
                            typename ProblemType::ProblemParserType& parser);

    template <typename ProblemType>
    static void stage_serial(HDGDiscretization<ProblemType>& discretization,
                             typename ProblemType::ProblemGlobalDataType& global_data,
                             typename ProblemType::ProblemStepperType& stepper);

    template <template <typename> typename OMPISimUnitType, typename ProblemType>
    static void step_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                          typename ProblemType::ProblemGlobalDataType& global_data,
                          typename ProblemType::ProblemStepperType& stepper,
                          const uint begin_sim_id,
                          const uint end_sim_id);

    template <template <typename> typename OMPISimUnitType, typename ProblemType>
    static void stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                           typename ProblemType::ProblemGlobalDataType& global_data,
                           typename ProblemType::ProblemStepperType& stepper,
                           const uint begin_sim_id,
                           const uint end_sim_id);

    template <typename HPXSimUnitType>
    static auto stage_hpx(HPXSimUnitType* sim_unit);

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

    /* local step end */

    /* global step begin */

    template <typename StepperType, typename InterfaceType>
    static void global_interface_kernel(const StepperType& stepper, InterfaceType& intface);

    template <typename StepperType, typename EdgeInterfaceType>
    static void global_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int);

    template <typename StepperType, typename EdgeInterfaceType>
    static void global_edge_interface_iteration(const StepperType& stepper, EdgeInterfaceType& edge_int);

    template <typename StepperType, typename BoundaryType>
    static void global_boundary_kernel(const StepperType& stepper, BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    static void global_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    static void global_edge_boundary_iteration(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename DistributedBoundaryType>
    static void global_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound);

    template <typename StepperType, typename EdgeDistributedType>
    static void global_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound);

    template <typename StepperType, typename EdgeDistributedType>
    static void global_edge_distributed_iteration(const StepperType& stepper, EdgeDistributedType& edge_dbound);

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

    template <typename GlobalDataType>
    static void finalize_simulation(GlobalDataType& global_data) {}
};
}
}

#endif
