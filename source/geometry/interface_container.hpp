#ifndef INTERFACE_CONTAINER_HPP
#define INTERFACE_CONTAINER_HPP

#include "geometry/interface.hpp"

namespace Geometry {

template <typename Interface, typename BoundarySoAType>
class InterfaceContainer;

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType, typename BoundarySoAType>
class InterfaceContainer<Interface<dimension, IntegrationType, DataType, SpecializationType>, BoundarySoAType> {
public:
    using InterfaceType = Interface<dimension, IntegrationType, DataType, SpecializationType>;

    InterfaceContainer() : capacity(0UL), size_(0UL) {}

    void reserve(uint ninterfaces) {
        this->capacity = ninterfaces;

        interface_accessors.reserve(ninterfaces);
    }

    template <typename... Args>
    void CreateInterface(Args&&... args) {
        if ( size_++ >= capacity ) {
            throw std::runtime_error{"Not enough interfaces reserved. Emplacing new interfaces will cause accessors to be invalidated\n"};
        }

        interface_accessors.emplace_back(std::forward<Args>(args)...);
    }

    template <typename F>
    void CallForEachInterface(const F& f, std::false_type /*is not vectorized*/) {
        std::for_each(interface_accessors.begin(), interface_accessors.end(), [&f](InterfaceType& intfc) {
                f(intfc);
            });
    }

    size_t size() const { return size_; }
private:
    AlignedVector<InterfaceType> interface_accessors;
    //InterfaceSoA<InterfaceType, BoundarySoAType> interface_data;

    size_t capacity;
    size_t size_;
};


}

#endif