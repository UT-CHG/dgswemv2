#ifndef HEARTBEAT_HPP
#define HEARTBEAT_HPP

#include <chrono>

#ifdef HAS_HPX
#include "serialize_chrono.hpp"
#endif

namespace Utilities {
class HeartBeat {
  public:
    using clock_t = std::chrono::system_clock;
    using time_point_t = std::chrono::time_point<clock_t>;

    HeartBeat()=default;
    HeartBeat(const std::chrono::duration<double>& period);

    bool Thump();

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif

  private:
    clock_t::duration period;
    time_point_t t_next;
};

HeartBeat::HeartBeat(const std::chrono::duration<double>& period_)
    : period(std::chrono::duration_cast<clock_t::duration>(period_)),
      t_next((clock_t::now().time_since_epoch()/period+2)*period) {}

bool HeartBeat::Thump() {
    time_point_t now = clock_t::now();
    if ( now > this->t_next ) {
        while ( now > this->t_next ) {
            this->t_next += this->period;
        }
        return true;
    }
    return false;
}

#ifdef HAS_HPX
template <typename Archive>
void HeartBeat::serialize(Archive& ar, unsigned) {
    ar & period & t_next;
}
#endif
}
#endif