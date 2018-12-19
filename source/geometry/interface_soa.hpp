#ifndef GEOMETRY_INTERFACE_SOA_HPP
#define GEOMETRY_INTERFACE_SOA_HPP

namespace Geometry {

template <typename InterfaceType, typename ElementContainers>
class InterfaceSoA;

template <typename InterfaceType, typename... ElementContainer>
class InterfaceSoA<InterfaceType,
                   std::tuple<ElementContainer...>> {
public:
    using ElementContainers = std::tuple<ElementContainer...>;

    InterfaceSoA() = default;
    InterfaceSoA(ElementContainers& element_data_) : element_data(&element_data_) {}

/*    void reserve(uint n_interfaces, uint ngp) {
        data = ProblemBoundarySoA(n_interfaces, ngp);
    }

    template <typename... Args>
    AccessorType at(const size_t index,
                    Args&&... args) {
        assert(element_data);
        return AccessorType(boundary_data.at(index),
                            std::forward<Args>(args)...);
                            }*/

private:
    std::tuple<ElementContainer...> * element_data = nullptr;
};
}

#endif