#include "preprocessor/input_parameters.hpp"
#include "geometry/mesh_definitions.hpp"

#include "problem/SWE/swe_definitions.hpp"
#include "problem/SWE/boundary_conditions/swe_boundary_conditions.hpp"
#include "problem/SWE/data_structure/swe_data.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"

#include "preprocessor/initialize_mesh.hpp"

#include "simulation/writer.hpp"

#include "utilities/almost_equal.hpp"

using InputType = InputParameters<SWE::Inputs>;
using MeshType=Geometry::MeshType<SWE::Data,SWE::Distributed,SWE::Land,SWE::Tidal,SWE::Flow>;

struct DummyProblem {
    using ProblemInputType = SWE::Inputs;
    using ProblemDataType = SWE::Data;
    using ProblemMeshType = MeshType;

    template <typename RawBoundaryType>
    static void create_boundaries_kernel(ProblemMeshType& mesh,
                                         std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries,
                                         Writer<DummyProblem>& writer) {}

    template <typename RawBoundaryType>
    static void create_distributed_boundaries_kernel(ProblemMeshType&,
                                                     std::tuple<>&,
                                                     std::map<uint, std::map<uint, RawBoundaryType>>&,
                                                     Writer<DummyProblem>&) {}

};

using WriterType = Writer<DummyProblem>;


bool equal(MeshType& o_mesh, MeshType& i_mesh) {
    bool is_equal{true};
    is_equal &= (o_mesh.GetNumberElements() == i_mesh.GetNumberElements());
    is_equal &= (o_mesh.GetMeshName() == i_mesh.GetMeshName());

    int64_t o_val{0};
    int64_t toggle{1};

    o_mesh.CallForEachElement([&o_val,&toggle](auto& element) {
            std::cout << o_val << '\n';
            o_val += toggle * element.GetID();
            toggle *= -1;
        });

    toggle = 1;
    int64_t i_val{0};

    i_mesh.CallForEachElement([&i_val,&toggle](auto& element) {
            i_val += toggle * element.GetID();
            toggle *= -1;
        });

    std::cout << "  o_val: " << o_val << '\n'
              << "  i_val: " << i_val << '\n';

    if ( i_val != o_val ) {
        std::cerr << "Error: toggle times element id not equal for o_mesh and i_mesh\n";
        is_equal = false;
    }
    std::cout << '\n';

    const auto sin_fcn = [](Point<2>& x)->double { return x[0]*x[0] + x[1]*x[1]; };

    double o_integral{0};
    o_mesh.CallForEachElement([&sin_fcn, &o_integral] (auto& element) {
            o_integral += element.IntegrationPhi(0, element.L2Projection(sin_fcn));
        });

    double i_integral{0};
    i_mesh.CallForEachElement([&sin_fcn, &i_integral] (auto& element) {
            i_integral += element.IntegrationPhi(0, element.L2Projection(sin_fcn));
        });

    std::cout << "  o_integral: " << o_integral << '\n'
              << "  i_integral: " << i_integral << '\n';

    if ( !Utilities::almost_equal(o_integral,i_integral) ) {
        std::cerr << "Error: Integration of quadric function not equal\n";
        is_equal = false;
    }
    std::cout << '\n';

    return is_equal;
};

int main(int argc, char** argv) {

    InputType input(argv[1]);
    input.ReadMesh();

    MeshType o_mesh(3);
    o_mesh.SetMeshName(input.mesh_data.mesh_name);

    WriterType writer(input);

    std::tuple<> empty_comm;

    initialize_mesh<DummyProblem>(o_mesh, input.mesh_data, empty_comm, input.problem_input, writer);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_mesh;

    hpx::serialization::input_archive i_archive(buffer);
    MeshType i_mesh;
    i_archive >> i_mesh;

    if ( !equal(o_mesh, i_mesh) ) {
        return 1;
    }
    return 0;
}
