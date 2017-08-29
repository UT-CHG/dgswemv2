#ifndef MESH_UTILITIES_HPP
#define MESH_UTILITIES_HPP

namespace Geometry {
template <typename Element>
struct make_master_type;

template <typename Element>
struct make_master_type<std::tuple<Element>> {
    typedef std::tuple<typename Element::ElementMasterType> type;
};

template <typename Element, typename... Es>
struct make_master_type<std::tuple<Element, Es...>> {
    typedef typename Utilities::tuple_join<std::tuple<typename Element::ElementMasterType>,
                                           typename make_master_type<std::tuple<Es...>>::type>::type type;
};

template <typename... Master>
struct master_maker;

template <typename... Master>
struct master_maker<std::tuple<Master...>> {
    static std::tuple<Master...> construct_masters(uint p) {
        return std::make_tuple(Master(p)...);
    };
};
}

#endif