#ifndef SERIALIZE_CHRONO_HPP
#define SERIALIZE_CHRONO_HPP

#include <hpx/config.hpp>
#include <hpx/runtime/serialization/input_archive.hpp>
#include <hpx/runtime/serialization/output_archive.hpp>

#include <chrono>

namespace hpx { namespace serialization {
template <typename Rep, typename Period>
void serialize(input_archive& ar, std::chrono::duration<Rep,Period>& t, unsigned) {
    Rep rep;
    ar >> rep;
    t = std::chrono::duration<Rep,Period>(rep);
}

template <typename Rep, typename Period>
void serialize(output_archive& ar, const std::chrono::duration<Rep,Period>& t, unsigned) {
    ar << t.count();
}

template <typename Clock, typename Duration>
void serialize(input_archive& ar, std::chrono::time_point<Clock,Duration>& t, unsigned) {
    Duration d;
    ar >> d;
    t = std::chrono::time_point<Clock,Duration>(d);
}

template <typename Clock, typename Duration>
void serialize(output_archive& ar, const std::chrono::time_point<Clock,Duration>& t, unsigned) {
    ar << t.time_since_epoch();
}
}
}
#endif