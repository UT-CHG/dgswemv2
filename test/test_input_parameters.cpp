#include "preprocessor/input_parameters.hpp"
#include "utilities/almost_equal.hpp"

const static auto equal = [](const InputParameters<typename SWE::Inputs> & ipa,
                             const InputParameters<typename SWE::Inputs> & ipb) -> bool {

    const SWE::Inputs& ia = ipa.problem_input;
    const SWE::Inputs& ib = ipb.problem_input;

    bool inputs_are_equal =
        (ia.g == ib.g) && (ia.bottom_friction.type == ib.bottom_friction.type) &&
        Utilities::almost_equal(ia.bottom_friction.coefficient, ib.bottom_friction.coefficient) &&
        (ia.initial_conditions.type == ib.initial_conditions.type) &&
        Utilities::almost_equal(ia.initial_conditions.ze_initial, ib.initial_conditions.ze_initial) &&
        Utilities::almost_equal(ia.initial_conditions.qx_initial, ib.initial_conditions.qx_initial) &&
        Utilities::almost_equal(ia.initial_conditions.qy_initial, ib.initial_conditions.qy_initial);

    return (ipa.mesh_file_name == ipb.mesh_file_name) && (ipa.mesh_format == ipb.mesh_format) &&
           (ipa.rk.nstages == ipb.rk.nstages) && (ipa.rk.order == ipb.rk.order) &&
           Utilities::almost_equal(ipa.dt, ipb.dt) && Utilities::almost_equal(ipa.T_end, ipb.T_end) &&
           (ipa.polynomial_order == ipb.polynomial_order) && inputs_are_equal;
};

const static auto equal2 = [](const InputParameters<> & ipa, const InputParameters<> & ipb) -> bool {
    const YAML::Node& na = ipa.problem_input.node;
    const YAML::Node& nb = ipb.problem_input.node;

    bool nodes_are_equal =
        (na["name"].as<std::string>() == nb["name"].as<std::string>()) &&
        (na["bottom_friction"]["type"].as<std::string>() == nb["bottom_friction"]["type"].as<std::string>()) &&
        Utilities::almost_equal(na["bottom_friction"]["coefficient"].as<double>(),
                                nb["bottom_friction"]["coefficient"].as<double>()) &&
        Utilities::almost_equal(na["gravity"].as<double>(), nb["gravity"].as<double>()) &&
        (na["initial_conditions"]["type"].as<std::string>() == nb["initial_conditions"]["type"].as<std::string>()) &&
        Utilities::almost_equal(na["initial_conditions"]["initial_surface_height"].as<double>(),
                                nb["initial_conditions"]["initial_surface_height"].as<double>()) &&
        Utilities::almost_equal(na["initial_conditions"]["initial_momentum_x"].as<double>(),
                                nb["initial_conditions"]["initial_momentum_x"].as<double>()) &&
        Utilities::almost_equal(na["initial_conditions"]["initial_momentum_y"].as<double>(),
                                nb["initial_conditions"]["initial_momentum_y"].as<double>());

    return (ipa.mesh_file_name == ipb.mesh_file_name) && (ipa.mesh_format == ipb.mesh_format) &&
           (ipa.rk.nstages == ipb.rk.nstages) && (ipa.rk.order == ipb.rk.order) &&
           Utilities::almost_equal(ipa.dt, ipb.dt) && Utilities::almost_equal(ipa.T_end, ipb.T_end) &&
           (ipa.polynomial_order == ipb.polynomial_order) && (ipa.problem_input.node == ipb.problem_input.node) &&
           nodes_are_equal;
};

int main(int argc, char** argv) {
    bool error_found{false};

    // try reading a well formatted input file
    std::cout << "Try a correct input file\n";
    {
        bool local_error{false};
        try {
            InputParameters<typename SWE::Inputs> input(argv[1]);
            std::string output_file_name = std::string(argv[1]) + ".emitted";
            std::cout << output_file_name << '\n';
            input.WriteTo(output_file_name);

            InputParameters<typename SWE::Inputs> input2(output_file_name);
            local_error = !equal(input, input2);

            if (local_error) {
                std::cerr << "Error found in correct file with SWE::Inputs type\n";
            }
        }
        catch (const std::exception& e) {
            local_error = true;
            std::cout << "Bad News: Exception was thrown ( " << e.what() << " )\n";
        }

        error_found = error_found || local_error;
    }

    std::cout << "Try a correct input file with anonymous input type\n";
    {
        bool local_error{false};
        try {
            InputParameters<> input(argv[1]);
            std::string output_file_name = std::string(argv[1]) + ".emitted";
            std::cout << output_file_name << '\n';
            input.WriteTo(output_file_name);

            InputParameters<> input2(output_file_name);
            local_error = !equal2(input, input);

            if (local_error) {
                std::cerr << "Error found in correct file with anonymous input type\n";
            }
        }
        catch (const std::exception& e) {
            local_error = true;
            std::cout << "Bad News: Exception was thrown ( " << e.what() << " )\n";
        }

        error_found = error_found || local_error;
    }

    std::cout << "Now we will pass in two bad meshes and see if exceptions get "
                 "thrown\n";
    // try reading a file that doesn't exist
    {
        bool local_error{true};
        try {
            InputParameters<typename SWE::Inputs> input(argv[2]);
        }
        catch (const std::exception& e) {
            local_error = false;
            std::cout << "Good News: Exception was thrown ( " << e.what() << " )\n";
        }

        error_found = error_found || local_error;
        if (local_error) {
            std::cout << "Bad News: Exception should have been thrown\n";
        }
    }

    // try reading a file that exists but that's missing a field
    {
        bool local_error{true};
        try {
            InputParameters<typename SWE::Inputs> input(argv[3]);
        }
        catch (const std::exception& e) {
            local_error = false;
            std::cout << "Good News: Exception was thrown ( " << e.what() << " )\n";
        }

        if (local_error) {
            std::cout << "Bad News: Exception should have been thrown\n";
        }

        error_found = error_found || local_error;
    }

    if (error_found) {
        return -1;
    }

    return 0;
}
