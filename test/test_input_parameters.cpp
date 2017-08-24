#include "preprocessor/input_parameters.hpp"
#include "utilities/almost_equal.hpp"

const static auto equal = [](const InputParameters & ipa, const InputParameters & ipb) -> bool {
    return (ipa.mesh_file_name == ipb.mesh_file_name) && (ipa.mesh_format == ipb.mesh_format) &&
           (ipa.rk.nstages == ipb.rk.nstages) && (ipa.rk.order == ipb.rk.order) &&
           Utilities::almost_equal(ipa.dt, ipb.dt) && Utilities::almost_equal(ipa.T_end, ipb.T_end) &&
           (ipa.polynomial_order == ipb.polynomial_order);
};
int main(int argc, char** argv) {
    bool error_found{false};

    // try reading a well formatted input file
    std::cout << "Try a correct input file\n";
    {
        bool local_error{false};
        try {
            InputParameters input(argv[1]);
            std::string output_file_name = std::string(argv[1]) + ".emitted";
            std::cout << output_file_name << '\n';
            input.WriteTo(output_file_name);

            InputParameters input2(output_file_name);
            local_error = !equal(input, input2);
        }
        catch (const std::exception& e) {
            local_error = true;
            std::cout << "Bad News: Exception was thrown ( " << e.what() << " )\n";
        }

        error_found = error_found || local_error;
    }

    std::cout << "Now we will pass in two bad meshes and see if exceptions get thrown\n";
    // try reading a file that doesn't exist
    {
        bool local_error{true};
        try {
            InputParameters input(argv[2]);
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
            InputParameters input(argv[3]);
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
