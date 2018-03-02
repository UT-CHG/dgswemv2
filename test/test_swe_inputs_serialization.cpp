#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "utilities/almost_equal.hpp"

bool inputs_are_equal(const SWE::Inputs& a, const SWE::Inputs& b) {
    using Utilities::almost_equal;

    bool are_same{true};
    are_same &= almost_equal(a.g, b.g);
    are_same &= almost_equal(a.h_o, b.h_o);
    are_same &= (a.parse_input == b.parse_input);
    are_same &= ((a.initial_conditions.type == b.initial_conditions.type) &&
                 almost_equal(a.initial_conditions.ze_initial, b.initial_conditions.ze_initial) &&
                 almost_equal(a.initial_conditions.qx_initial, b.initial_conditions.qx_initial) &&
                 almost_equal(a.initial_conditions.qy_initial, b.initial_conditions.qy_initial));

    are_same &= (a.function_source.type == b.function_source.type);

    are_same &= ((a.bottom_friction.type == b.bottom_friction.type) &&
                 almost_equal(a.bottom_friction.coefficient, b.bottom_friction.coefficient) &&
                 (a.bottom_friction.manning_data_file == b.bottom_friction.manning_data_file));

    are_same &= ((a.meteo_forcing.type == b.meteo_forcing.type) &&
                 (a.meteo_forcing.meteo_data_file == b.meteo_forcing.meteo_data_file) &&
                 almost_equal(a.meteo_forcing.frequency, b.meteo_forcing.frequency));

    are_same &= (a.tidal_potential.type == b.tidal_potential.type);
    are_same &= (a.coriolis.type == b.coriolis.type);
    return are_same;
}

int main() {

    SWE::Inputs o_input;
    o_input.g = 4.2;
    o_input.h_o = 10.2;
    o_input.parse_input = true;

    SWE::InitialConditions ic_{SWE::InitialConditionsType::Function, -1.0, -2.0, -3.0};
    o_input.initial_conditions = ic_;

    SWE::BottomFriction bf_{SWE::BottomFrictionType::Chezy, 0.2, "foo.13"};
    o_input.bottom_friction = bf_;

    SWE::MeteoForcing mf_{SWE::MeteoForcingType::Test, "bar.13", 3600.};
    o_input.meteo_forcing = mf_;

    o_input.tidal_potential.type = SWE::TidalPotentialType::Test;
    o_input.coriolis.type = SWE::CoriolisType::Test;

    std::vector<char> buffer;
    hpx::serialization::output_archive o_ar(buffer);
    o_ar << o_input;

    hpx::serialization::input_archive i_ar(buffer);
    SWE::Inputs i_input;
    i_ar >> i_input;

    if (!inputs_are_equal(o_input, i_input)) {
        return 1;
    }
    return 0;
}