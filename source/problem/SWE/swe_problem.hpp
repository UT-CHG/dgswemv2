#ifndef SWE_PROBLEM_HPP
#define SWE_PROBLEM_HPP

#include "swe_definitions.hpp"
#include "swe_boundary_conditions.hpp"
#include "swe_data.hpp"
#include "../../geometry/mesh_definitions.hpp"
#include "../../preprocessor/mesh_metadata.hpp"

namespace SWE {
struct Problem {
    typedef SWE::Data ProblemDataType;

    typedef Geometry::MeshType<SWE::Data, SWE::Distributed, SWE::Land, SWE::Tidal> ProblemMeshType;

    // preprocessor kernels
    template <typename RawBoundaryType>
    static void create_boundaries_kernel(ProblemMeshType& mesh,
                                         std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries);

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries_kernel(
        ProblemMeshType& mesh,
        Communicator& communicator,
        std::map<uint, std::map<uint, RawBoundaryType>>& pre_distributed_boundaries);

    static void initialize_data_kernel(ProblemMeshType& mesh, const MeshMetaData& mesh_data);

    // processor kernels
    template <typename ElementType>
    static void volume_kernel(const Stepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void source_kernel(const Stepper& stepper, ElementType& elt);

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
    static void extract_modal_data_kernel(ElementType& elt, std::vector<std::pair<uint, Array2D<double>>>& modal_data);

    template <typename MeshType>
    static void write_modal_data_kernel(const Stepper& stepper, MeshType& mesh);
};
}

#endif
