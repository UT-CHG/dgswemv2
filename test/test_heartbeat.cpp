#include "utilities/heartbeat.hpp"

#include <iostream>

int main() {

    using clock_t = std::chrono::system_clock;
    using time_point_t = std::chrono::time_point<clock_t>;

    time_point_t start = clock_t::now();
    Utilities::HeartBeat heartbeat(std::chrono::milliseconds(500));

    int counter{0};

    while ( std::chrono::duration_cast<std::chrono::seconds>(clock_t::now() - start) < std::chrono::seconds(4) ) {
        if (heartbeat.Thump() ) {
            std::cout << "Thump @ " << std::chrono::duration_cast<std::chrono::milliseconds>(clock_t::now() - start).count() << '\n';
            ++counter;
        }
    }

    return counter != 5;
}