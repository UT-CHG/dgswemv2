#ifndef SWE_PROBLEM_HPP
#define SWE_PROBLEM_HPP

#include "simulation/stepper/rk_stepper.hpp"
#include "simulation/writer.hpp"

#include "problem/SWE/swe_definitions.hpp"
#include "boundary_conditions/swe_boundary_conditions.hpp"
#include "dist_boundary_conditions/swe_distributed_boundary_conditions.hpp"
#include "interface_specializations/swe_interface_specializations.hpp"
#include "data_structure/swe_data.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"

#include "geometry/mesh_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

namespace SWE {
namespace RKDG {
struct Problem {
    typedef SWE::Inputs ProblemInputType;

    typedef SWE::RKDG::Data ProblemDataType;

    typedef Geometry::MeshType<SWE::RKDG::Data,
                               std::tuple<SWE::RKDG::IS::Internal, SWE::RKDG::IS::Levee>,
                               std::tuple<SWE::RKDG::BC::Land, SWE::RKDG::BC::Tide, SWE::RKDG::BC::Flow>,
                               std::tuple<SWE::RKDG::DBC::Distributed, SWE::RKDG::DBC::DistributedLevee>>::Type
        ProblemMeshType;

    typedef SWE::Parser ProblemParserType;

    // preprocessor kernels
    static void initialize_problem_parameters(const ProblemInputType& problem_specific_input);

    static void preprocess_mesh_data(InputParameters<ProblemInputType>& input);

    template <typename RawBoundaryType>
    static void create_interfaces_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        InputParameters<ProblemInputType>& input,
        Writer<SWE::RKDG::Problem>& writer);

    template <typename RawBoundaryType>
    static void create_boundaries_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        InputParameters<ProblemInputType>& input,
        Writer<SWE::RKDG::Problem>& writer);

    template <typename RawBoundaryType>
    static void create_distributed_boundaries_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType&,
        InputParameters<ProblemInputType>& input,
        std::tuple<>&,
        Writer<SWE::RKDG::Problem>&);

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        InputParameters<ProblemInputType>& input,
        Communicator& communicator,
        Writer<SWE::RKDG::Problem>& writer);

    static void initialize_data_kernel(ProblemMeshType& mesh,
                                       const MeshMetaData& mesh_data,
                                       const ProblemInputType& problem_specific_input);

    static void initialize_data_parallel_pre_send_kernel(ProblemMeshType& mesh,
                                                         const MeshMetaData& mesh_data,
                                                         const ProblemInputType& problem_specific_input);

    static void initialize_data_parallel_post_receive_kernel(ProblemMeshType& mesh);

    // processor kernels
    template <typename ElementType>
    static void volume_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void source_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename InterfaceType>
    static void interface_kernel(const RKStepper& stepper, InterfaceType& intface);

    template <typename BoundaryType>
    static void boundary_kernel(const RKStepper& stepper, BoundaryType& bound);

    template <typename DistributedBoundaryType>
    static void distributed_boundary_send_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    static void distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound);

    template <typename ElementType>
    static void update_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static bool scrutinize_solution_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void swap_states_kernel(const RKStepper& stepper, ElementType& elt);

    // postprocessor kernels
    static void postprocessor_serial_kernel(const RKStepper& stepper, ProblemMeshType& mesh);

    static void postprocessor_parallel_pre_send_kernel(const RKStepper& stepper, ProblemMeshType& mesh);

    static void postprocessor_parallel_pre_receive_kernel(const RKStepper& stepper, ProblemMeshType& mesh);

    static void postprocessor_parallel_post_receive_kernel(const RKStepper& stepper, ProblemMeshType& mesh);

    template <typename ElementType>
    static void wetting_drying_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void slope_limiting_prepare_element_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename InterfaceType>
    static void slope_limiting_prepare_interface_kernel(const RKStepper& stepper, InterfaceType& intface);

    template <typename BoundaryType>
    static void slope_limiting_prepare_boundary_kernel(const RKStepper& stepper, BoundaryType& bound);

    template <typename DistributedBoundaryType>
    static void slope_limiting_distributed_boundary_send_kernel(const RKStepper& stepper,
                                                                DistributedBoundaryType& dbound);

    template <typename DistributedBoundaryType>
    static void slope_limiting_prepare_distributed_boundary_kernel(const RKStepper& stepper,
                                                                   DistributedBoundaryType& dbound);

    template <typename ElementType>
    static void slope_limiting_kernel(const RKStepper& stepper, ElementType& elt);

    // writing output kernels
    static void write_VTK_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file);

    static void write_VTU_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file);

    static void write_modal_data_kernel(const RKStepper& stepper,
                                        ProblemMeshType& mesh,
                                        const std::string& output_path);

    template <typename ElementType>
    static double compute_residual_L2_kernel(const RKStepper& stepper, ElementType& elt);
};
}
}

#endif
