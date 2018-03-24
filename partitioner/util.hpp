#ifndef _cb3576bc_af75_4da3_b0c5_a69b2f21f2c2
#define _cb3576bc_af75_4da3_b0c5_a69b2f21f2c2

#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>
#include <tuple>
#include <unordered_map>

namespace stdx {
template <int i, int n, class tup>
struct hash_tuple {
    std::size_t operator()(size_t h, const tup& x) const {
        h ^= std::hash<typename std::tuple_element<i, tup>::type>()(std::get<i>(x));
        h *= 0x9e3779b97f4a7c15u;
        return hash_tuple<i + 1, n, tup>()(h, x);
    }
};

template <int n, class tup>
struct hash_tuple<n, n, tup> {
    std::size_t operator()(size_t h, const tup& x) const { return h; }
};
}

namespace std {
template <class... Ts>
std::size_t hash_of(const Ts... xs) {
    return stdx::hash_tuple<0, sizeof...(Ts), std::tuple<Ts...>>()(0, std::tuple<Ts...>(xs...));
}

/*template<class T>
struct hash {
  std::size_t operator()(const T &x) const {
    return hash_of(x);
  }
};*/

template <class... Ts>
struct hash<std::tuple<Ts...>> {
    std::size_t operator()(std::tuple<Ts...> x) const {
        return stdx::hash_tuple<0, sizeof...(Ts), std::tuple<Ts...>>()(0, x);
    }
};

// grabbed from jdbachan's typeclass.hxx
template <class A, class B>
struct hash<pair<A, B>> {
    inline size_t operator()(const pair<A, B>& x) const {
        size_t h = hash<A>()(x.first);
        h ^= h >> 13;
        h *= 41;
        h += hash<B>()(x.second);
        return h;
    }
};
}

namespace stdx {
template <int i, int n, class tup>
struct print_tuple {
    void operator()(std::ostream& o, const tup& x) const {
        if (i != 0)
            o << ',';
        o << std::get<i>(x);
        print_tuple<i + 1, n, tup>()(o, x);
    }
};

template <int n, class tup>
struct print_tuple<n, n, tup> {
    void operator()(std::ostream& o, const tup& x) const {}
};
}

template <class... Ts>
std::ostream& operator<<(std::ostream& o, const std::tuple<Ts...>& tup) {
    o << '{';
    stdx::print_tuple<0, sizeof...(Ts), std::tuple<Ts...>>()(o, tup);
    o << '}';
    return o;
}

// mhb: copied from progAMR
namespace mota {
// simplified Say grabbed from jdbachan's diagnostic.hxx
class Say {
    std::stringstream ss;

  public:
    Say() {}
    ~Say() {
        ss << '\n';
        std::cerr << ss.str();
        std::cerr.flush();
    }
    template <class T>
    Say& operator<<(const T& x) {
        ss << x;
        return *this;
    }
};
}

template <class key, class T>
std::unordered_map<key, T> intersect(const std::unordered_map<key, T>& um1, const std::unordered_map<key, T>& um2) {
    const std::unordered_map<key, T>& smaller = (um1.size() < um2.size()) ? um1 : um2;
    const std::unordered_map<key, T>& larger  = (um1.size() >= um2.size()) ? um1 : um2;
    std::unordered_map<key, T>        r_um;

    for (auto key_T : smaller) {
        if (larger.count(key_T.first)) {
            assert(key_T.second == larger.at(key_T.first));
            r_um.insert(key_T);
        }
    }

    return r_um;
}

#endif
