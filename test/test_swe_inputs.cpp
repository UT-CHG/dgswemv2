#include "problem/SWE/problem_input/swe_inputs.hpp"

#include "utilities/almost_equal.hpp"

int main() {
    bool error_found{false};

    // The first test checks that the default behavior is correct for all fields
    {
        std::cout << "Beginning test 1\n";

        YAML::Node  test;
        SWE::Inputs result(test);

        if (!Utilities::almost_equal(9.81, result.g)) {
            std::cerr << "Gravity was set incorrectly\n";
            error_found = true;
        }

        const SWE::BottomFriction& bf = result.bottom_friction;
        if (!(bf.type == SWE::BottomFrictionType::None && Utilities::almost_equal(bf.coefficient, 0.0))) {
            std::cerr << "Default bottom friction is incorrectly set\n";
            error_found = true;
        }

        const SWE::InitialConditions& ics = result.initial_conditions;
        if (!(ics.type == SWE::InitialConditionsType::Default && Utilities::almost_equal(0., ics.ze_initial) &&
              Utilities::almost_equal(0., ics.qx_initial) && Utilities::almost_equal(0., ics.qy_initial))) {
            std::cerr << "Default initial conditions are incorrectly set\n";
            error_found = true;
        }
    }

    // The second test examines setting all of the values
    {
        std::cout << "\nBeginning test 2\n";
        YAML::Node test;
        test["gravity"] = 1.0;

        YAML::Node bf_node;
        bf_node["type"]        = std::string("Chezy");
        bf_node["coefficient"] = 0.01;

        YAML::Node ic_node;
        ic_node["type"]                   = "Constant";
        ic_node["initial_surface_height"] = 1.0;
        ic_node["initial_momentum_x"]     = 2.0;
        ic_node["initial_momentum_y"]     = 3.0;

        test["bottom_friction"]    = bf_node;
        test["initial_conditions"] = ic_node;

        SWE::Inputs result(test);
        if (!Utilities::almost_equal(1.0, result.g)) {
            std::cerr << "Gravity was set incorrectly\n";
            error_found = true;
        }

        const SWE::BottomFriction& bf = result.bottom_friction;
        if (!(bf.type == SWE::BottomFrictionType::Chezy && Utilities::almost_equal(bf.coefficient, 0.01))) {
            std::cerr << "Bottom friction is incorrectly set\n";
            error_found = true;
        }

        const SWE::InitialConditions& ics = result.initial_conditions;
        if (!(ics.type == SWE::InitialConditionsType::Constant && Utilities::almost_equal(1., ics.ze_initial) &&
              Utilities::almost_equal(2., ics.qx_initial) && Utilities::almost_equal(3., ics.qy_initial))) {
            std::cerr << "Initial conditions are incorrectly set\n";
            error_found = true;
        }
    }

    // The  test checks that an error gets throw for negative friction coefficients
    {
        std::cout << "\nBeginning test 3\n";
        YAML::Node bf_node;
        bf_node["type"]        = "Chezy";
        bf_node["coefficient"] = -0.01;

        YAML::Node test;
        test["bottom_friction"] = bf_node;

        bool local_error{true};
        try {
            SWE::Inputs results(test);
        } catch (std::exception& e) {
            std::cout << "Good news (this error should have been thrown)\n"
                      << "    " << e.what() << '\n';
            local_error = false;
        }
        if (local_error) {
            std::cerr << "Error: No exception was thrown for a Chezy"
                      << " friction law with a negative friction coefficient\n";
            error_found = true;
        }
    }

    // Check function gets set
    {
        std::cout << "\nBeginning test 4\n";

        YAML::Node ic_node;
        ic_node["type"] = "Function";
        YAML::Node test;
        test["initial_conditions"] = ic_node;

        SWE::Inputs result(test);
        if (result.initial_conditions.type != SWE::InitialConditionsType::Function) {
            std::cerr << "Error: Unable to correctly set initial conditions via a function\n";
            error_found = true;
        }
    }

    if (error_found) {
        return 1;
    }
    return 0;
}