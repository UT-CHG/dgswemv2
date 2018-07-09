#include "preprocessor/input_parameters.hpp"
#include "utilities/almost_equal.hpp"

const static auto equal_writer = [](const WriterInput& wa, const WriterInput& wb) -> bool {
    bool output_nodes_are_equal =
        (wa.writing_output == wb.writing_output) && (wa.output_path == wb.output_path) &&
        (wa.log_file_name == wb.log_file_name) && (wa.writing_vtk_output == wb.writing_vtk_output) &&
        (wa.vtk_output_frequency == wb.vtk_output_frequency) && (wa.writing_modal_output == wb.writing_modal_output) &&
        (wa.modal_output_frequency == wb.modal_output_frequency);

    if (!output_nodes_are_equal) {
        std::cerr << "writing_output: " << wa.writing_output << " : " << wb.writing_output << '\n'
                  << "              " << wa.output_path << " : " << wb.output_path << '\n'
                  << "              " << wa.log_file_name << " : " << wb.log_file_name << '\n'
                  << "              " << wa.writing_vtk_output << " : " << wb.writing_vtk_output << '\n'
                  << "              " << wa.vtk_output_frequency << " : " << wb.vtk_output_frequency << '\n'
                  << "              " << wa.writing_modal_output << " : " << wb.writing_modal_output << '\n'
                  << "              " << wa.modal_output_frequency << " : " << wb.modal_output_frequency << '\n';

        std::cerr << "Error: Writer Inputs not equal\n";
    }

    return output_nodes_are_equal;
};

const static auto equal_load_balancer = [](const LoadBalancerInput& lba, const LoadBalancerInput& lbb) -> bool {
    bool load_balancer_nodes_are_equal = (lba.use_load_balancer == lbb.use_load_balancer) && (lba.name == lbb.name) &&
                                         Utilities::almost_equal(lba.rebalance_frequency, lbb.rebalance_frequency);

    if (!load_balancer_nodes_are_equal) {
        std::cerr << "Load Balancer nodes are not equal" << std::endl;
        std::cerr << "        " << std::boolalpha << lba.use_load_balancer << " : " << lbb.use_load_balancer << '\n'
                  << "        " << lba.name << " : " << lbb.name << '\n'
                  << "        " << lba.rebalance_frequency << " : " << lbb.rebalance_frequency << '\n';
    }

    return load_balancer_nodes_are_equal;
};

const static auto equal = [](const InputParameters<typename SWE::Inputs>& ipa,
                             const InputParameters<typename SWE::Inputs>& ipb) -> bool {
    const SWE::Inputs& ia = ipa.problem_input;
    const SWE::Inputs& ib = ipb.problem_input;

    bool inputs_are_equal =
        (ia.g == ib.g) && (ia.bottom_friction.type == ib.bottom_friction.type) &&
        Utilities::almost_equal(ia.bottom_friction.coefficient, ib.bottom_friction.coefficient) &&
        (ia.initial_conditions.type == ib.initial_conditions.type) &&
        Utilities::almost_equal(ia.initial_conditions.ze_initial, ib.initial_conditions.ze_initial) &&
        Utilities::almost_equal(ia.initial_conditions.qx_initial, ib.initial_conditions.qx_initial) &&
        Utilities::almost_equal(ia.initial_conditions.qy_initial, ib.initial_conditions.qy_initial);

    if (!inputs_are_equal) {
        std::cerr << "Error: Problem Specific inputs not equal\n";
    }

    return (ipa.mesh_input.mesh_file_name == ipb.mesh_input.mesh_file_name) &&
           (ipa.mesh_input.mesh_format == ipb.mesh_input.mesh_format) &&
           (ipa.stepper_input.nstages == ipb.stepper_input.nstages) &&
           (ipa.stepper_input.order == ipb.stepper_input.order) &&
           Utilities::almost_equal(ipa.stepper_input.dt, ipb.stepper_input.dt) &&
           Utilities::almost_equal(ipa.stepper_input.run_time, ipb.stepper_input.run_time) &&
           (ipa.polynomial_order == ipb.polynomial_order) && inputs_are_equal &&
           equal_writer(ipa.writer_input, ipb.writer_input) &&
           equal_load_balancer(ipa.load_balancer_input, ipb.load_balancer_input);
    ;
};

const static auto equal2 = [](const InputParameters<>& ipa, const InputParameters<>& ipb) -> bool {
    const YAML::Node& na = ipa.problem_input.node;
    const YAML::Node& nb = ipb.problem_input.node;

    bool swe_nodes_are_equal =
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

    if (!swe_nodes_are_equal) {
        std::cerr << "Error: Problem Specific inputs not equal\n";
    }

    return (ipa.mesh_input.mesh_file_name == ipb.mesh_input.mesh_file_name) &&
           (ipa.mesh_input.mesh_format == ipb.mesh_input.mesh_format) &&
           (ipa.stepper_input.nstages == ipb.stepper_input.nstages) &&
           (ipa.stepper_input.order == ipb.stepper_input.order) &&
           Utilities::almost_equal(ipa.stepper_input.dt, ipb.stepper_input.dt) &&
           Utilities::almost_equal(ipa.stepper_input.run_time, ipb.stepper_input.run_time) &&
           (ipa.polynomial_order == ipb.polynomial_order) && swe_nodes_are_equal &&
           equal_writer(ipa.writer_input, ipb.writer_input) &&
           equal_load_balancer(ipa.load_balancer_input, ipb.load_balancer_input);
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
            std::cout << "Emitted filename: " << output_file_name << '\n';
            input.write_to(output_file_name);

            InputParameters<typename SWE::Inputs> input2(output_file_name);
            local_error = !equal(input, input2);

            if (local_error) {
                std::cerr << "Error found in correct file with SWE::Inputs type\n";
            }
        } catch (const std::exception& e) {
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
            std::cout << "Emitted filename: " << output_file_name << '\n';
            input.write_to(output_file_name);

            InputParameters<> input2(output_file_name);
            local_error = !equal2(input, input2);

            if (local_error) {
                std::cerr << "Error found in correct file with anonymous input type\n";
            }
        } catch (const std::exception& e) {
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
        } catch (const std::exception& e) {
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
        } catch (const std::exception& e) {
            local_error = false;
            std::cout << "Good News: Exception was thrown ( " << e.what() << " )\n";
        }

        if (local_error) {
            std::cout << "Bad News: Exception should have been thrown\n";
        }

        error_found = error_found || local_error;
    }

    // try a correct input file with no output node specified
    std::cout << "Try a correct input file\n";
    {
        bool local_error{false};
        try {
            InputParameters<typename SWE::Inputs> input(argv[4]);
            std::string output_file_name = std::string(argv[4]) + ".emitted";
            std::cout << "Emitted filename: " << output_file_name << '\n';
            input.write_to(output_file_name);

            InputParameters<typename SWE::Inputs> input2(output_file_name);
            local_error = !equal(input, input2);

            if (local_error) {
                std::cerr << "Error found in correct file with SWE::Inputs type\n";
            }
        } catch (const std::exception& e) {
            local_error = true;
            std::cout << "Bad News: Exception was thrown ( " << e.what() << " )\n";
        }

        error_found = error_found || local_error;
    }

    if (error_found) {
        return -1;
    }

    return 0;
}
