#ifndef BOUNDARY_CONTAINER_HPP
#define BOUNDARY_CONTAINER_HPP

#include "geometry/boundary.hpp"

namespace Geometry {

template <typename Boundary>
class BoundaryContainer;

template <uint dimension, typename IntegrationType, typename DataType, typename ConditionType>
class BoundaryContainer<Boundary<dimension, IntegrationType, DataType, ConditionType>> {
public:
    using BoundaryType = Boundary<dimension, IntegrationType, DataType, ConditionType>;

    BoundaryContainer() : capacity(0UL), size_(0UL) {}

    void reserve(uint nboundaries) {
        this->capacity = nboundaries;

        boundary_accessors.reserve(nboundaries);
    }

    template <typename... Args>
    void CreateBoundary(Args&&... args) {
        if ( size_++ >= capacity ) {
            throw std::runtime_error{"Not enough boundaries reserved. Emplacing new boundaries will cause accessors to be invalidated\n"};
        }

        boundary_accessors.emplace_back(std::forward<Args>(args)...);
    }

    template <typename F>
    void CallForEachBoundary(const F& f, std::false_type /*is not vectorized*/) {
        std::for_each(boundary_accessors.begin(), boundary_accessors.end(), [&f](BoundaryType& bdry) {
                f(bdry);
            });
    }

    size_t size() const { return size_; }

private:
    AlignedVector<BoundaryType> boundary_accessors;
    //BoundarySoA<BoundaryType, BoundarySoAType> interface_data;

    size_t capacity;
    size_t size_;

};

}

#endif