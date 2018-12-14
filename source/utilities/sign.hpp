#ifndef SIGN_HPP
#define SIGN_HPP

#include "general_definitions.hpp"

namespace Utilities {
// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}
}

#endif