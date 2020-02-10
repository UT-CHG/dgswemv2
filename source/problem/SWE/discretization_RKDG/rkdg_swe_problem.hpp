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

    using ProblemDataType       = SWE::Data;
    using ProblemEdgeDataType   = SWE::EdgeData;
    using ProblemGlobalDataType = SWE::GlobalData;

    using ProblemInterfaceTypes = Geometry::InterfaceTypeTuple<Data, ISP::Internal, ISP::Levee>;
    using ProblemBoundaryTypes  = Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow, BC::Function>;
    using ProblemDistributedBoundaryTypes =
        Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed, DBC::DistributedLevee>;

    using ProblemMeshType = Geometry::MeshType<Data,
                                               std::tuple<ISP::Internal, ISP::Levee>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow, BC::Function>,
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
        begin_index[CommTypes::bound_state] += (SWE::n_variables + 1) * ngp + 1;  // + bath + w/d state
        begin_index[CommTypes::baryctr_state] += SWE::n_variables + ngp + 1;      // + w/d state

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

    template <template <typename> class DiscretizationType, typename ProblemType>
    static void preprocessor_serial(DiscretizationType<ProblemType>& discretization,
                                    typename ProblemType::ProblemGlobalDataType& global_data,
                                    const ProblemStepperType& stepper,
                                    const typename ProblemType::ProblemInputType& problem_specific_input);

    template <template <typename> class OMPISimUnitType, typename ProblemType>
    static void preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                                  typename ProblemType::ProblemGlobalDataType& global_data,
                                  const ProblemStepperType& stepper,
                                  const uint begin_sim_id,
                                  const uint end_sim_id);

    // processor kernels
    template <template <typename> class DiscretizationType, typename ProblemType>
    static void step_serial(DiscretizationType<ProblemType>& discretization,
                            typename ProblemType::ProblemGlobalDataType& global_data,
                            ProblemStepperType& stepper,
                            typename ProblemType::ProblemWriterType& writer,
                            typename ProblemType::ProblemParserType& parser);

    template <template <typename> class DiscretizationType, typename ProblemType>
    static void stage_serial(DiscretizationType<ProblemType>& discretization,
                             typename ProblemType::ProblemGlobalDataType& global_data,
                             ProblemStepperType& stepper);

    template <template <typename> class OMPISimUnitType, typename ProblemType>
    static void step_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                          typename ProblemType::ProblemGlobalDataType& global_data,
                          ProblemStepperType& stepper,
                          const uint begin_sim_id,
                          const uint end_sim_id);

    template <template <typename> class OMPISimUnitType, typename ProblemType>
    static void stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                           typename ProblemType::ProblemGlobalDataType& global_data,
                           ProblemStepperType& stepper,
                           const uint begin_sim_id,
                           const uint end_sim_id);

    template <typename ElementType>
    static void volume_kernel(const ProblemStepperType& stepper, ElementType& elt);

    template <typename ElementType>
    static void source_kernel(const ProblemStepperType& stepper, ElementType& elt);

    template <typename InterfaceType>
    static void interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface);

    template <typename BoundaryType>
    static void boundary_kernel(const ProblemStepperType& stepper, BoundaryType& bound);

    template <typename DistributedBoundaryType>
    static void distributed_boundary_send_kernel(const ProblemStepperType& stepper, DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    static void distributed_boundary_kernel(const ProblemStepperType& stepper, DistributedBoundaryType& dbound);

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
