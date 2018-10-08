#ifndef EHDG_GN_PROBLEM_HPP
#define EHDG_GN_PROBLEM_HPP

#include "simulation/stepper/rk_stepper.hpp"
#include "simulation/writer.hpp"
#include "simulation/discretization.hpp"

#include "problem/Green-Naghdi/gn_definitions.hpp"

#include "boundary_conditions/ehdg_gn_boundary_conditions.hpp"
#include "dist_boundary_conditions/ehdg_gn_distributed_boundary_conditions.hpp"
#include "interface_specializations/ehdg_gn_interface_specializations.hpp"

#include "data_structure/ehdg_gn_data.hpp"
#include "data_structure/ehdg_gn_edge_data.hpp"
#include "data_structure/ehdg_gn_global_data.hpp"

#include "problem/Green-Naghdi/problem_input/gn_inputs.hpp"
#include "problem/Green-Naghdi/problem_parser/gn_parser.hpp"

#include "geometry/mesh_definitions.hpp"
#include "geometry/mesh_skeleton_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

namespace GN {
namespace EHDG {
struct Problem {
    using ProblemInputType = GN::Inputs;

    using ProblemWriterType = Writer<Problem>;

    using ProblemParserType = GN::Parser;

    using ProblemDataType = Data;

    using ProblemEdgeDataType = EdgeData;

    using ProblemGlobalDataType = GlobalData;

    using ProblemMeshType = Geometry::MeshType<Data,
                                               std::tuple<ISP::Internal>,
                                               std::tuple<BC::Land, BC::Tide, BC::Flow>,
                                               std::tuple<DBC::Distributed>>::Type;

    using ProblemMeshSkeletonType =
        Geometry::MeshSkeletonType<EdgeData,
                                   Geometry::InterfaceTypeTuple<Data, ISP::Internal>,
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
        ProblemWriterType&);

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries(
        std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
        ProblemMeshType& mesh,
        ProblemInputType& input,
        Communicator& communicator,
        ProblemWriterType& writer);

    static void create_edge_interfaces(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer);

    static void create_edge_boundaries(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer);

    static void create_edge_distributeds(ProblemMeshType& mesh,
                                         ProblemMeshSkeletonType& mesh_skeleton,
                                         ProblemWriterType& writer);

    static void preprocessor_serial(ProblemDiscretizationType& discretization,
                                    const ProblemInputType& problem_specific_input);

    template <typename OMPISimUnitType>
    static void preprocessor_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    static void initialize_dc_data_serial(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input);

    static void initialize_dc_data_parallel(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input);

    static void initialize_global_dc_problem(ProblemDiscretizationType& discretization, uint& dc_global_dof_offset);

    template <typename Communicator>
    static void initialize_global_dc_problem_parallel_pre_send(ProblemDiscretizationType& discretization,
                                                               Communicator& communicator,
                                                               uint& dc_global_dof_offset);

    static void initialize_global_dc_problem_parallel_finalize_pre_send(ProblemDiscretizationType& discretization,
                                                                        uint dc_global_dof_offset);

    template <typename Communicator>
    static void initialize_global_dc_problem_parallel_post_receive(ProblemDiscretizationType& discretization,
                                                                   Communicator& communicator,
                                                                   std::vector<uint>& dc_global_dof_indx);

    static void compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization);

    template <typename OMPISimUnitType>
    static void compute_bathymetry_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    // processor kernels
    static void stage_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization);

    template <typename OMPISimUnitType>
    static void stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    /* Dispersive correction part */

    static void dispersive_correction_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization);

    static void compute_derivatives_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization);

    template <typename OMPISimUnitType>
    static void dispersive_correction_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    template <typename OMPISimUnitType>
    static void compute_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    /* local step */

    template <typename ElementType>
    static void local_dc_volume_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void local_dc_source_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename EdgeInterfaceType>
    static void local_dc_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeBoundaryType>
    static void local_dc_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    /* global step */

    template <typename EdgeInterfaceType>
    static void global_dc_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int);

    template <typename EdgeBoundaryType>
    static void global_dc_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    static void solve_global_dc_problem(const RKStepper& stepper, ProblemDiscretizationType& discretization);

    /* SWE part */

    static void swe_stage_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization);

    template <typename OMPISimUnitType>
    static void swe_stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units);

    /* SWE part */

    template <typename ElementType>
    static void dc_update_kernel(const RKStepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void dispersive_correction_kernel(const RKStepper& stepper, ElementType& elt);

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
