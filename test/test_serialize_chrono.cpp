#include <hpx/runtime/serialization/serialize.hpp>
#include "utilities/serialize_chrono.hpp"

#include <iostream>
#include <iomanip>

int main() {

    using clock_t = std::chrono::system_clock;
    using time_point_t = std::chrono::time_point<clock_t>;
    using duration = clock_t::duration;

    time_point_t o_start = clock_t::now();
    duration o_duration = clock_t::now() - o_start;

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_start << o_duration;

    hpx::serialization::input_archive i_archive(buffer);
    time_point_t i_start;
    duration i_duration;
    i_archive >> i_start >> i_duration;

    bool error_found{false};
    if ( o_duration.count() != i_duration.count() ) {
        error_found = true;
        std::cerr << "Error in serializing std::chrono::duration<>\n"
                  << "o_duration.count(): " << o_duration.count() << '\n'
                  << "i_duration.count(): " << i_duration.count() << '\n';
    }
    if ( o_start.time_since_epoch().count() != i_start.time_since_epoch().count() ) {
        error_found = true;
        std::cerr << "Error in serializing std::chrono::time_point<>\n"
                  << "o_start.time_since_epoch().count(): " << o_start.time_since_epoch().count() << '\n'
                  << "i_start.time_since_epoch().count(): " << i_start.time_since_epoch().count() << '\n';
    }
    return error_found;
}