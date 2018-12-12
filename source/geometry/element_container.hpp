#ifndef GEOMETRY_ELEMENT_CONTAINER_HPP
#define GEOMETRY_ELEMENT_CONTAINER_HPP

#include "geometry/element.hpp"
#include "geometry/element_soa.hpp"

#include "utilities/tuple_helpers.hpp"

namespace Geometry {

template <typename Element, typename DataSoAType>
class ElementContainer;

template <uint dimension, typename MasterType, typename ShapeType, typename AccessorType, typename DataSoAType>
class ElementContainer<Element<dimension,MasterType,ShapeType, AccessorType>, DataSoAType> {
public:
    using ElementType = Element<dimension,MasterType,ShapeType,AccessorType>;

    ElementContainer() = default;
    ElementContainer(uint p) : master_element(p), size_(0u), capacity(0u) {}

    void reserve(uint nstages, uint nelements) {

        std::cout << "master address in element ctor " << &master_element << std::endl
                  << " nvrtxs: " << master_element.nvrtx << std::endl;

        this->capacity = nelements;

        //Element data gets constructed here since we seem to create to copy the values in construct_containers(p)
        //which invalidates pointers to this->master_element.
        element_data = ElementSoA<ElementType, DataSoAType>(this->master_element);

        element_accessors.reserve(nelements);
        element_data.reserve(master_element.ndof, nstages, nelements);
    }

    template <typename...Args>
    void CreateElement(const uint ID, Args&&... args) {
        if ( size_ >= capacity ) {
            throw std::runtime_error{"Not enough elements reserved. Emplacing new elements will cause accessors to be invalidated\n"};
        }

        element_accessors.emplace_back(std::move(element_data.at(size_++,
                                                                 ID,
                                                                 std::forward<Args>(args)...)));
    }

    template <typename F>/*,
                           typename = std::enable_if<!Utilities::is_vectorized<F>::value>::type>*/
    void CallForEachElement(const F& f) {
        std::for_each(element_accessors.begin(), element_accessors.end(), [&f](ElementType& elt) { f(elt); });
    }

/*    template <typename F,
              typename = std::enable_if<Utilities::is_vectorized<F>::value>::type>
    void CallForEachElement(const F& f) {
        f(element_data);
        }*/

    size_t size() { return size_; }

private:
    MasterType master_element;

    AlignedVector<ElementType> element_accessors;
    ElementSoA<ElementType, DataSoAType> element_data;

    size_t size_;
    size_t capacity;
};
}

#endif