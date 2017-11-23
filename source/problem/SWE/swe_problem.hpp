#ifndef SWE_PROBLEM_HPP
#define SWE_PROBLEM_HPP

#include "../../simulation/stepper.hpp"
#include "../../simulation/writer.hpp"

#include "swe_definitions.hpp"
#include "boundary_conditions/swe_boundary_conditions.hpp"
#include "data/swe_data.hpp"
#include "input/swe_inputs.hpp"

#include "../../geometry/mesh_definitions.hpp"
#include "../../preprocessor/mesh_metadata.hpp"

namespace SWE {
struct Problem {
    typedef SWE::Inputs ProblemInputType;

    typedef SWE::Data ProblemDataType;

    typedef Geometry::MeshType<SWE::Data, SWE::Distributed, SWE::Land, SWE::Tidal, SWE::Flow> ProblemMeshType;

    // preprocessor kernels
    template <typename RawBoundaryType>
    static void create_boundaries_kernel(ProblemMeshType& mesh,
                                         std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries,
                                         Writer<SWE::Problem>& writer);

    template <typename RawBoundaryType>
    static void create_distributed_boundaries_kernel(ProblemMeshType&,
                                                     std::tuple<>&,
                                                     std::map<uint, std::map<uint, RawBoundaryType>>&,
                                                     Writer<SWE::Problem>&);

    template <typename RawBoundaryType, typename Communicator>
    static void create_distributed_boundaries_kernel(
        ProblemMeshType& mesh,
        Communicator& communicator,
        std::map<uint, std::map<uint, RawBoundaryType>>& pre_distributed_boundaries,
        Writer<SWE::Problem>& writer);

    static void initialize_data_kernel(ProblemMeshType& mesh,
                                       const MeshMetaData& mesh_data,
                                       const ProblemInputType& problem_specific_input);

    static void initialize_wd_data_kernel(ProblemMeshType& mesh);

    static void initialize_sl_data_kernel(ProblemMeshType& mesh);

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
    static void scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void wetting_drying_kernel(const Stepper& stepper, ElementType& elt);

    template <typename ElementType>
    static void slope_limiting_prepare_element_kernel(const Stepper& stepper,
                                                      ElementType& elt);
                                                      
                                                       template <typename InterfaceType>
    static void slope_limiting_prepare_interface_kernel(const Stepper& stepper,
                                                        InterfaceType& intface);
                                                        
                                                         template <typename BoundaryType>
    static void slope_limiting_prepare_boundary_kernel(const Stepper& stepper,
                                                       BoundaryType& bound);
                                                       
                                                        template <typename ElementType>
    static void slope_limiting_kernel(const Stepper& stepper, ElementType& elt);
    
     template <typename ElementType>
    static void swap_states_kernel(const Stepper& stepper, ElementType& elt);

    // postprocessor kernels
    template <typename ElementType>
    static void extract_VTK_data_kernel(ElementType& elt, Array2D<double>& cell_data, Array2D<double>& point_data);

    static void write_VTK_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file);

    template <typename ElementType>
    static void extract_modal_data_kernel(ElementType& elt, std::vector<std::pair<uint, Array2D<double>>>& modal_data);

    static void write_modal_data_kernel(const Stepper& stepper, ProblemMeshType& mesh, const std::string& output_path);

    template <typename ElementType>
    static double compute_residual_L2_kernel(const Stepper& stepper, ElementType& elt);
};
}

#endif
