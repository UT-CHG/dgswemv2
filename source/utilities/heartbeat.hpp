#ifndef HEARTBEAT_HPP
#define HEARTBEAT_HPP

#include "general_definitions.hpp"

namespace Utilities {
class HeartBeat {
  public:
    using clock_t      = std::chrono::system_clock;
    using time_point_t = std::chrono::time_point<clock_t>;

    HeartBeat() = default;
    HeartBeat(const std::chrono::duration<double>& period);

    bool Thump();

  private:
    clock_t::duration period;
    time_point_t t_next;

  public:
};

// Note that by default a duration of double is interpreted as seconds
HeartBeat::HeartBeat(const std::chrono::duration<double>& period_)
    : period(std::chrono::duration_cast<clock_t::duration>(period_)),
      t_next((clock_t::now().time_since_epoch() / period + 2) * period) {}

bool HeartBeat::Thump() {
    time_point_t now = clock_t::now();
    if (now > this->t_next) {
        while (now > this->t_next) {
            this->t_next += this->period;
        }
        return true;
    }
    return false;
}
}
#endif