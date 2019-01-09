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

    void reserve(uint nstages, uint nelements, uint ngp_edge) {

        this->capacity = nelements;

        //Element data gets constructed here since we seem to create to copy the values in construct_containers(p)
        //which invalidates pointers to this->master_element.
        element_data = ElementSoA<ElementType, DataSoAType>(this->master_element);

        element_accessors.reserve(nelements);
        element_data.reserve(nstages, nelements, ngp_edge);
    }

    void finalize_initialization() {
        for ( uint elt_id = 0; elt_id < element_accessors.size(); ++elt_id ) {
            element_data.set_abs_J(elt_id, element_accessors.at(elt_id).GetAbsJ());
            element_data.set_J_inv(elt_id, element_accessors.at(elt_id).GetJinv());
        }
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

    ptrdiff_t GetLocalIndex(const AccessorType& acc) const {
        const AccessorType* problem_data_start = &element_accessors.front().data;

        ptrdiff_t  indx = ( &acc - problem_data_start ) / sizeof(AccessorType);

        //return a negative number if the accessor is not contained in this element container
        return indx < 0                                        ? indx :
               (unsigned long) indx < element_accessors.size() ? indx : -1L;
    }

    template <typename F>
    void CallForEachElement(const F& f, std::false_type /*is not vectorized*/) {
        std::for_each(element_accessors.begin(), element_accessors.end(), [&f](ElementType& elt) { f(elt); });
    }

    template <typename F>
    void CallForEachElement(const F& f, std::true_type /*is vectorized*/) {
        f(element_data);
    }

    size_t size() { return size_; }

    const MasterType& GetMaster() const {
        return master_element;
    }

private:
    MasterType master_element;

    AlignedVector<ElementType> element_accessors;
    ElementSoA<ElementType, DataSoAType> element_data;

    size_t size_;
    size_t capacity;
};
}

#endif