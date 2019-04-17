#ifndef RKDG_SWE_PROBLEM_HPP
#define RKDG_SWE_PROBLEM_HPP

#include "problem/SWE/swe_definitions.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

#include "simulation/stepper/explicit_ssp_rk_stepper.hpp"
#include "simulation/writer.hpp"
#include "simulation/discretization.hpp"

#include "boundary_conditions/rkdg_swe_boundary_conditions.hpp"
#include "dist_boundary_conditions/rkdg_swe_distributed_boundary_conditions.hpp"
#include "interface_specializations/rkdg_swe_interface_specializations.hpp"

#include "problem/SWE/problem_data_structure/swe_data.hpp"
#include "problem/SWE/problem_data_structure/swe_data_interface.hpp"
#include "problem/SWE/problem_data_structure/swe_edge_data.hpp"
#include "problem/SWE/problem_data_structure/swe_global_data.hpp"

#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"
#include "problem/SWE/problem_preprocessor/swe_preprocessor.hpp"
#include "problem/SWE/problem_postprocessor/swe_postprocessor.hpp"

#include "geometry/mesh_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "problem/SWE/problem_preprocessor/swe_pre_create_intface.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_bound.hpp"
#include "problem/SWE/problem_preprocessor/swe_pre_create_dbound.hpp"

namespace SWE {
namespace RKDG {
struct Problem {
    using ProblemInputType   = SWE::Inputs;
    using ProblemStepperType = ESSPRKStepper;
    using ProblemWriterType  = Writer<Problem>;
    using ProblemParserType  = SWE::Parser;

    using ProblemAccessorType = SWE::Accessor;
    using ProblemSoAContainerType  = SWE::SoAContainer;
    using ProblemEdgeDataType   = SWE::EdgeData;
    using ProblemGlobalDataType = SWE::GlobalData;

    using ProblemInterfaceSoAType = SWE::InterfaceData;

    using ProblemInterfaceTypes = Geometry::InterfaceTypeTuple<SWE::Accessor, ISP::Internal, ISP::Levee>;
    using ProblemBoundaryTypes  = Geometry::BoundaryTypeTuple<SWE::Accessor, BC::Land, BC::Tide, BC::Flow, BC::Function, BC::Outflow>;
    using ProblemDistributedBoundaryTypes =
        Geometry::DistributedBoundaryTypeTuple<SWE::Accessor, DBC::Distributed, DBC::DistributedLevee>;

    using ProblemMeshType = Geometry::MeshType<SWE::Accessor,
                                               SWE::SoAContainer,
                                               ProblemInterfaceSoAType,
                                               std::tuple<ISP::Internal, ISP::Levee>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow, BC::Function, BC::Outflow>,
                                               std::tuple<DBC::Distributed, DBC::DistributedLevee>>::Type;

    using ProblemDiscretizationType = DGDiscretization<Problem>;

    // preprocessor kernels
    static void initialize_problem_parameters(const ProblemInputType& problem_specific_input) {
        SWE::initialize_problem_parameters(problem_specific_input);
    }

    static void preprocess_mesh_data(InputParameters<ProblemInputType>& input) { SWE::preprocess_mesh_data(input); }

    // helpers to create communication
    static constexpr uint n_communications = SWE::RKDG::n_communications;

    static std::vector<uint> comm_buffer_offsets(std::vector<uint>& begin_index, uint ngp) {
        std::vector<uint> offset(SWE::RKDG::n_communications);

        offset[CommTypes::baryctr_coord] = begin_index[CommTypes::baryctr_coord];
        offset[CommTypes::bound_state]   = begin_index[CommTypes::bound_state];
        offset[CommTypes::baryctr_state] = begin_index[CommTypes::baryctr_state];

        begin_index[CommTypes::baryctr_coord] += 2;
        begin_index[CommTypes::bound_state] += SWE::n_variables * ngp + 1;  // + w/d state
        begin_index[CommTypes::baryctr_state] += SWE::n_variables + 1;      // + w/d state

        return offset;
    }

    template <typename RawBoundaryType>
    static void create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& input,
                                  ProblemWriterType& writer) {
        SWE::create_interfaces<SWE::RKDG::Problem>(raw_boundaries, mesh, input, writer);
    }

    template <typename RawBoundaryType>
    static void create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& input,
                                  ProblemWriterType& writer) {
        SWE::create_boundaries<SWE::RKDG::Problem>(raw_boundaries, mesh, input, writer);
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
        SWE::create_distributed_boundaries<SWE::RKDG::Problem>(raw_boundaries, mesh, input, communicator, writer);
    }

    template <template <typename> typename DiscretizationType, typename ProblemType>
    static void preprocessor_serial(DiscretizationType<ProblemType>& discretization,
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

    // processor kernels
    template <template <typename> typename DiscretizationType, typename ProblemType>
    static void step_serial(DiscretizationType<ProblemType>& discretization,
                            typename ProblemType::ProblemGlobalDataType& global_data,
                            typename ProblemType::ProblemStepperType& stepper,
                            typename ProblemType::ProblemWriterType& writer,
                            typename ProblemType::ProblemParserType& parser);

    template <template <typename> typename DiscretizationType, typename ProblemType>
    static void stage_serial(DiscretizationType<ProblemType>& discretization,
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

    template <typename StepperType, typename ElementType>
    static void source_kernel(const StepperType& stepper, ElementType& elt);

    template <typename InterfaceType>
    static void interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface);

    template <typename StepperType, typename BoundaryType>
    static void boundary_kernel(const StepperType& stepper, BoundaryType& bound);

    template <typename StepperType, typename DistributedBoundaryType>
    static void distributed_boundary_send_kernel(const StepperType& stepper, DistributedBoundaryType& dbound);

    template <typename StepperType, typename DistributedBoundaryType>
    static void distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound);

    template <typename StepperType, typename ElementType>
    static void wetting_drying_kernel(const StepperType& stepper, ElementType& elt);

    // postprocessor kernels
    static void write_VTK_data(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
        SWE::write_VTK_data(mesh, raw_data_file);
    }

    static void write_VTU_data(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
        SWE::write_VTU_data(mesh, raw_data_file);
    }

    static void write_modal_data(const ProblemStepperType& stepper,
                                 ProblemMeshType& mesh,
                                 const std::string& output_path) {
        SWE::write_modal_data(stepper, mesh, output_path);
    }

    template <typename ElementType>
    static double compute_residual_L2(const ProblemStepperType& stepper, ElementType& elt) {
        return SWE::compute_residual_L2(stepper, elt);
    }

    static void finalize_simulation(ProblemGlobalDataType& global_data) {}
};
}
}

#endif
