#ifndef GEOMETRY_ELEMENT_CONTAINER_HPP
#define GEOMETRY_ELEMENT_CONTAINER_HPP

#include "geometry/element.hpp"
#include "geometry/element_soa.hpp"

#include "utilities/tuple_helpers.hpp"

namespace Geometry {

template <typename Element>
class ElementContainer;

template <uint dimension, typename MasterType, typename ShapeType, typename AccessorType>
class ElementContainer<Element<dimension,MasterType,ShapeType, AccessorType>> {
public:
    using ElementType = Element<dimension,MasterType,ShapeType,AccessorType>;

    ElementContainer() = default;
    ElementContainer(uint p) : master_element(p), size_(0u), capacity(0u) {}

    template <typename...Args>
    void CreateElement(const uint ID, Args&&... args) {
        ++size_;
        element_accessors.emplace_back(ID, master_element, std::forward<Args>(args)...);
    }

    template <typename F>
    void CallForEachElement(const F& f) {
        std::for_each(element_accessors.begin(), element_accessors.end(), [&f](ElementType& elt) { f(elt); });
    }

    size_t size() { return size_; }

private:
    MasterType master_element;

    AlignedVector<ElementType> element_accessors;
    ElementSoA<ElementType> element_data;

    size_t size_;
    size_t capacity;
};
}

#endif