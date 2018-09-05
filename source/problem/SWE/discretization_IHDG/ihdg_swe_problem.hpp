#ifndef IHDG_SWE_PROBLEM_HPP
#define IHDG_SWE_PROBLEM_HPP

#include "simulation/stepper/rk_stepper.hpp"
#include "simulation/writer.hpp"
#include "simulation/discretization.hpp"

#include "problem/SWE/swe_definitions.hpp"

#include "boundary_conditions/ihdg_swe_boundary_conditions.hpp"
#include "dist_boundary_conditions/ihdg_swe_distributed_boundary_conditions.hpp"
#include "interface_specializations/ihdg_swe_interface_specializations.hpp"

#include "data_structure/ihdg_swe_data.hpp"
#include "data_structure/ihdg_swe_edge_data.hpp"
#include "data_structure/ihdg_swe_global_data.hpp"

#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/SWE/problem_parser/swe_parser.hpp"

#include "geometry/mesh_definitions.hpp"
#include "geometry/mesh_skeleton_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

namespace SWE {
namespace IHDG {
struct Problem {
    using ProblemInputType = SWE::Inputs;

    using ProblemParserType = SWE::Parser;

    using ProblemDataType = Data;

    using ProblemEdgeDataType = EdgeData;

    using ProblemGlobalDataType = GlobalData;

    using ProblemMeshType = Geometry::MeshType<Data,
                                               std::tuple<IS::Internal>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow>,
                                               std::tuple<DBC::Distributed>>::Type;

    using ProblemMeshSkeletonType =
        Geometry::MeshSkeletonType<EdgeData,
                                   Geometry::InterfaceTypeTuple<Data, IS::Internal>,
                                   Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow>,
                                   Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>>::Type;

    using ProblemDiscretizationType = HDGDiscretization<Problem>;

    // preprocessor kernels
    static void initialize_problem_parameters(const ProblemInputType& problem_specific_input);

    static void preprocess_mesh_data(InputParameters<ProblemInputType>& input);

    template <typename RawBoundaryType>
    static void create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& problem_input,
                                  Writer<Problem>& writer);

    template <typename RawBoundaryType>
    static void create_boundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                  ProblemMeshType& mesh,
                                  ProblemInputType& problem_input,
                                  Writer<Problem>& writer);

    template <typename RawBoundaryType>
    static void create_distributed_boundaries(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType&,
        ProblemInputType& problem_input,
        std::tuple<>&,
        Writer<Problem>&);

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        ProblemInputType& input,
        Communicator& communicator,
        Writer<Problem>& writer);

    static void create_edge_interfaces(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       Writer<Problem>& writer);

    static void create_edge_boundaries(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       Writer<Problem>& writer);

    static void create_edge_distributeds(ProblemMeshType& mesh,
                                         ProblemMeshSkeletonType& mesh_skeleton,
                                         Writer<Problem>& writer);

    static void preprocessor_serial(ProblemDiscretizationType& discretization,
                                    const ProblemInputType& problem_specific_input);

    template <typename OMPISimUnitType>
    static void preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    static void initialize_data_serial(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input);

    static void initialize_data_parallel(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input);

    static void initialize_global_problem_serial(HDGDiscretization<Problem>& discretization);

    template <typename Communicator>
    static void initialize_global_problem_parallel_pre_send(HDGDiscretization<Problem>& discretization,
                                                            Communicator& communicator,
                                                            uint& local_dof_offset,
                                                            uint& global_dof_offset);

    static void initialize_global_problem_parallel_finalize_pre_send(HDGDiscretization<Problem>& discretization,
                                                                     uint local_dof_offset,
                                                                     uint global_dof_offset);

    template <typename Communicator>
    static void initialize_global_problem_parallel_post_receive(HDGDiscretization<Problem>& discretization,
                                                                Communicator& communicator);

    // processor kernels
    static void stage_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization);

    template <typename OMPISimUnitType>
    static void stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    /* local step begin */

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

    template <typename EdgeInterfaceType>
    static void local_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeBoundaryType>
    static void local_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeDistributedType>
    static void local_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound);

    /* local step end */

    /* global step begin */

    template <typename EdgeInterfaceType>
    static void global_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeBoundaryType>
    static void global_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeDistributedType>
    static void global_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_bound);

    static bool solve_global_problem(const RKStepper& stepper, HDGDiscretization<Problem>& discretization);

    /* global step end */

    template <typename ElementType>
    static void update_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static bool scrutinize_solution_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void swap_states_kernel(const RKStepper& stepper, ElementType& elt);

    // writing output kernels
    static void write_VTK_data(ProblemMeshType& mesh, std::ofstream& raw_data_file);

    static void write_VTU_data(ProblemMeshType& mesh, std::ofstream& raw_data_file);

    static void write_modal_data(const RKStepper& stepper, ProblemMeshType& mesh, const std::string& output_path);

    template <typename ElementType>
    static double compute_residual_L2(const RKStepper& stepper, ElementType& elt);
};
}
}

#endif
