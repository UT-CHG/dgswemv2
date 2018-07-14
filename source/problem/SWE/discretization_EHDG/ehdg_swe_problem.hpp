#ifndef EHDG_SWE_PROBLEM_HPP
#define EHDG_SWE_PROBLEM_HPP

#include "simulation/stepper/rk_stepper.hpp"
#include "simulation/writer.hpp"

#include "problem/SWE/swe_definitions.hpp"
#include "boundary_conditions/ehdg_swe_boundary_conditions.hpp"
#include "dist_boundary_conditions/ehdg_swe_distributed_boundary_conditions.hpp"
#include "interface_specializations/ehdg_swe_interface_specializations.hpp"
#include "data_structure/ehdg_swe_data.hpp"
#include "data_structure/ehdg_swe_edge_data.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"

#include "geometry/mesh_definitions.hpp"
#include "geometry/mesh_skeleton_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

namespace SWE {
namespace EHDG {
struct Problem {
    using ProblemInputType = SWE::Inputs;

    using ProblemParserType = SWE::Parser;

    using ProblemDataType = Data;

    using ProblemEdgeDataType = EdgeData;

    using ProblemMeshType = Geometry::MeshType<Data,
                                               std::tuple<IS::Internal>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow>,
                                               std::tuple<DBC::Distributed>>::Type;

    using ProblemMeshSkeletonType =
        Geometry::MeshSkeletonType<EdgeData,
                                   Geometry::InterfaceTypeTuple<Data, IS::Internal>,
                                   Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow>,
                                   Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>>::Type;

    // preprocessor kernels
    static void initialize_problem_parameters(const ProblemInputType& problem_specific_input);

    static void preprocess_mesh_data(InputParameters<ProblemInputType>& input);

    template <typename RawBoundaryType>
    static void create_interfaces_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        ProblemInputType& problem_input,
        Writer<Problem>& writer);

    template <typename RawBoundaryType>
    static void create_boundaries_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        ProblemInputType& problem_input,
        Writer<Problem>& writer);

    template <typename RawBoundaryType>
    static void create_distributed_boundaries_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType&,
        ProblemInputType& problem_input,
        std::tuple<>&,
        Writer<Problem>&);

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries_kernel(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        ProblemInputType& input,
        Communicator& communicator,
        Writer<Problem>& writer);

    static void create_edge_interfaces_kernel(ProblemMeshType& mesh,
                                              ProblemMeshSkeletonType& mesh_skeleton,
                                              Writer<Problem>& writer);

    static void create_edge_boundaries_kernel(ProblemMeshType& mesh,
                                              ProblemMeshSkeletonType& mesh_skeleton,
                                              Writer<Problem>& writer);

    static void create_edge_distributeds_kernel(ProblemMeshType& mesh,
                                                ProblemMeshSkeletonType& mesh_skeleton,
                                                Writer<Problem>& writer);

    static void initialize_data_kernel(ProblemMeshType& mesh,
                                       const MeshMetaData& mesh_data,
                                       const ProblemInputType& problem_specific_input);

    static void initialize_data_parallel_pre_send_kernel(ProblemMeshType& mesh,
                                                         const MeshMetaData& mesh_data,
                                                         const ProblemInputType& problem_specific_input);

    static void initialize_data_parallel_post_receive_kernel(ProblemMeshType& mesh);

    template <typename SimulationType>
    static void initialize_global_problem(SimulationType* simulation);

    // processor kernels

    /* global step */

    template <typename InterfaceType>
    static void global_interface_kernel(const RKStepper& stepper, InterfaceType& intface);

    template <typename EdgeInterfaceType>
    static void global_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    static void global_edge_interface_iteration(const RKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename BoundaryType>
    static void global_boundary_kernel(const RKStepper& stepper, BoundaryType& bound);

    template <typename EdgeBoundaryType>
    static void global_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeBoundaryType>
    static void global_edge_boundary_iteration(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename DistributedBoundaryType>
    static void global_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    static void global_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound);

    template <typename EdgeDistributedType>
    static void global_edge_distributed_iteration(const RKStepper& stepper, EdgeDistributedType& edge_dbound);

    /* local step */

    template <typename ElementType>
    static void local_volume_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void local_source_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename InterfaceType>
    static void local_interface_kernel(const RKStepper& stepper, InterfaceType& intface);

    template <typename BoundaryType>
    static void local_boundary_kernel(const RKStepper& stepper, BoundaryType& bound);

    template <typename DistributedBoundaryType>
    static void local_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound);

    /* local step */

    template <typename ElementType>
    static void update_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static bool scrutinize_solution_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void swap_states_kernel(const RKStepper& stepper, ElementType& elt);

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
