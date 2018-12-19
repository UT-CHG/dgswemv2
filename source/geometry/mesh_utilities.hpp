#ifndef MESH_UTILITIES_HPP
#define MESH_UTILITIES_HPP

#include "general_definitions.hpp"

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

template <typename ContainerTuple>
struct container_maker;

template <typename... Cs>
struct container_maker<std::tuple<Cs...>> {
    template <typename... Args>
    static std::tuple<Cs...> construct_containers(Args&&... args) {
        return std::make_tuple(Cs(std::forward<Args>(args)...)...);
    }
};
}

#endif