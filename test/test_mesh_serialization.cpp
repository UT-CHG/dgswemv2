#include "preprocessor/input_parameters.hpp"
#include "geometry/mesh_definitions.hpp"

#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"

#include "preprocessor/initialize_mesh.hpp"

#include "simulation/writer.hpp"

#include "utilities/almost_equal.hpp"

using InputType = InputParameters<SWE::Inputs>;
using MeshType  = Geometry::MeshType<SWE::RKDG::Data,
                                    std::tuple<SWE::RKDG::IS::Internal, SWE::RKDG::IS::Levee>,
                                    std::tuple<SWE::RKDG::BC::Land, SWE::RKDG::BC::Tide, SWE::RKDG::BC::Flow>,
                                    std::tuple<SWE::RKDG::DBC::Distributed, SWE::RKDG::DBC::DistributedLevee>>::Type;

struct DummyProblem {
    using ProblemInputType = SWE::Inputs;
    using ProblemDataType  = SWE::RKDG::Data;
    using ProblemMeshType  = MeshType;
};

using WriterType = Writer<DummyProblem>;

// specialization which doesn't connect interfaces
template <>
void initialize_mesh_interfaces_boundaries<DummyProblem, std::tuple<>>(typename DummyProblem::ProblemMeshType& mesh,
                                                                       typename DummyProblem::ProblemInputType& input,
                                                                       std::tuple<>& comm,
                                                                       Writer<DummyProblem>& writer) {}

bool equal(MeshType& o_mesh, MeshType& i_mesh) {
    bool is_equal{true};
    is_equal &= (o_mesh.GetNumberElements() == i_mesh.GetNumberElements());
    is_equal &= (o_mesh.GetMeshName() == i_mesh.GetMeshName());

    int64_t o_val{0};
    int64_t toggle{1};

    o_mesh.CallForEachElement([&o_val, &toggle](auto& element) {
        std::cout << o_val << '\n';
        o_val += toggle * element.GetID();
        toggle *= -1;
    });

    toggle = 1;
    int64_t i_val{0};

    i_mesh.CallForEachElement([&i_val, &toggle](auto& element) {
        i_val += toggle * element.GetID();
        toggle *= -1;
    });

    std::cout << "  o_val: " << o_val << '\n' << "  i_val: " << i_val << '\n';

    if (i_val != o_val) {
        std::cerr << "Error: toggle times element id not equal for o_mesh and i_mesh\n";
        is_equal = false;
    }
    std::cout << '\n';

    const auto sin_fcn = [](Point<2>& x) -> double { return x[0] * x[0] + x[1] * x[1]; };

    double o_integral{0};
    o_mesh.CallForEachElement([&sin_fcn, &o_integral](auto& element) {
        uint p = element.GetMaster().p;
        std::vector<double> modal_values((p + 2) * (p + 1) / 2);
        element.L2Projection(sin_fcn, modal_values);

        o_integral += element.IntegrationPhi(0, modal_values);
    });

    double i_integral{0};
    i_mesh.CallForEachElement([&sin_fcn, &i_integral](auto& element) {
        uint p = element.GetMaster().p;
        std::vector<double> modal_values((p + 2) * (p + 1) / 2);
        element.L2Projection(sin_fcn, modal_values);

        i_integral += element.IntegrationPhi(0, modal_values);
    });

    std::cout << "  o_integral: " << o_integral << '\n' << "  i_integral: " << i_integral << '\n';

    if (!Utilities::almost_equal(o_integral, i_integral)) {
        std::cerr << "Error: Integration of quadric function not equal\n";
        is_equal = false;
    }
    std::cout << '\n';

    return is_equal;
};

int main(int argc, char** argv) {
    InputType input(argv[1]);
    input.read_mesh();

    MeshType o_mesh(3);

    WriterType writer(input.writer_input);

    std::tuple<> empty_comm;

    initialize_mesh<DummyProblem>(o_mesh, input, empty_comm, writer);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_mesh;

    hpx::serialization::input_archive i_archive(buffer);
    MeshType i_mesh;
    i_archive >> i_mesh;
    i_mesh.SetMasters();
    i_mesh.CallForEachElement([&i_mesh](auto& elt) {
        using MasterType = typename std::remove_reference<decltype(elt)>::type::ElementMasterType;

        using MasterElementTypes = typename decltype(i_mesh)::MasterElementTypes;

        MasterType& master_elt =
            std::get<Utilities::index<MasterType, MasterElementTypes>::value>(i_mesh.GetMasters());

        element.SetMaster(master_elt);

        element.Initialize();
    });

    if (!equal(o_mesh, i_mesh)) {
        return 1;
    }
    return 0;
}
