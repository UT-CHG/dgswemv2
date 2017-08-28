#ifndef SWE_PROBLEM_HPP
#define SWE_PROBLEM_HPP

#include "swe_definitions.hpp"
#include "swe_boundary_conditions.hpp"
#include "swe_data.hpp"
#include "../../geometry/mesh_definitions.hpp"
#include "../../preprocessor/mesh_metadata.hpp"

namespace SWE {
struct Problem {
    typedef SWE::Data data_type;

    typedef Geometry::MeshType<SWE::Data, SWE::Land, SWE::Tidal> mesh_type;

    // preprocessor kernels
    template <typename RawBoundaryType>
    static void create_boundaries_kernel(mesh_type&, std::map<uchar, std::vector<RawBoundaryType>>&);

    static void initialize_data_kernel(mesh_type&, const MeshMetaData&);

    // processor kernels
    template <typename ElementType>
    static void volume_kernel(const Stepper&, ElementType&);

    template <typename ElementType>
    static void source_kernel(const Stepper&, ElementType&);

    template <typename InterfaceType>
    static void interface_kernel(const Stepper&, InterfaceType&);

    template <typename BoundaryType>
    static void boundary_kernel(const Stepper&, BoundaryType&);

    template <typename ElementType>
    static void update_kernel(const Stepper&, ElementType&);

    template <typename ElementType>
    static void swap_states_kernel(const Stepper&, ElementType&);

    template <typename ElementType>
    static void scrutinize_solution_kernel(const Stepper&, ElementType&);

    // postprocessor kernels
    template<typename ElementType>
	static double compute_residual_L2_kernel(const Stepper&, ElementType&);

    template <typename ElementType>
    static void extract_VTK_data_kernel(ElementType&, Array2D<double>&, Array2D<double>&);

    template <typename MeshType>
    static void write_VTK_data_kernel(const Stepper&, MeshType&);

    template <typename ElementType>
    static void extract_modal_data_kernel(ElementType&, std::vector<std::pair<uint, Array2D<double>>>&);

    template <typename MeshType>
    static void write_modal_data_kernel(const Stepper&, MeshType&);
};
}

#endif
