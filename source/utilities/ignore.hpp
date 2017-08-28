#ifndef IGNORE_HPP
#define IGNORE_HPP

namespace Utilities {
// using Herb Sutter's approach to swallowing warnings
// https://herbsutter.com/2009/10/18/mailbag-shutting-up-compiler-warnings/
template <typename T>
void ignore(const T&) {}
}

#endif