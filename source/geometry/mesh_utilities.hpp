#ifndef MESH_UTILITIES_HPP
#define MESH_UTILITIES_HPP

namespace Geometry {
template <typename E>
struct make_master_type;

template <typename E>
struct make_master_type<std::tuple<E>> {
    typedef std::tuple<typename E::ElementMasterType> type;
};

template <typename E, typename... Es>
struct make_master_type<std::tuple<E, Es...>> {
    typedef typename Utilities::tuple_join<std::tuple<typename E::ElementMasterType>,
                                           typename make_master_type<std::tuple<Es...>>::type>::type type;
};

template <typename... Ms>
struct master_maker;

template <typename... Ms>
struct master_maker<std::tuple<Ms...>> {
    static std::tuple<Ms...> construct_masters(uint p) { return std::make_tuple(Ms(p)...); };
};
}

#endif