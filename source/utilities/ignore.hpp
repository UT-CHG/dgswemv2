#ifndef IGNORE_HPP
#define IGNORE_HPP

namespace Utilities {
/**
 * Ignore the argument.
 * Ignore the argument to hide unused variable warnings. This file is implemented following
 * <a href="https://herbsutter.com/2009/10/18/mailbag-shutting-up-compiler-warnings/">Herb Sutter's approach</a>.
 */
template <typename T>
void ignore(const T&) {}
}

#endif