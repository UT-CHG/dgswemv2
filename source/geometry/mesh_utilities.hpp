#ifndef MESH_UTILITIES_HPP
#define MESH_UTILITIES_HPP

namespace Geometry {
template <typename E>
struct make_master_type;

template <typename E>
struct make_master_type<std::tuple<E>> {
    typedef std::tuple<typename E::master_element_type> type;
};

template <typename E, typename... Es>
struct make_master_type<std::tuple<E, Es...>> {
    typedef typename Utilities::tuple_join<std::tuple<typename E::master_element_type>,
                                           typename make_master_type<std::tuple<Es...>>::type>::type type;
};

template <typename... M>
struct master_maker;

template <typename... M>
struct master_maker<std::tuple<M...>> {
    static std::tuple<M...> construct_masters(uint p) {
        return std::make_tuple(M(p)...);
    };
};
}

#endif