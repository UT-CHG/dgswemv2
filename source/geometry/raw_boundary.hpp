#ifndef CLASS_RAW_BOUNDARY_HPP
#define CLASS_RAW_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, class DataType>
class RawBoundary {
  public:
    uint p;
    uint bound_id;

    DataType& data;

    Basis::Basis<dimension + 1>& basis;
    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

    RawBoundary(uint p,
                uint bound_id,
                DataType& data,
                Basis::Basis<dimension + 1>& basis,
                Master::Master<dimension + 1>& master,
                Shape::Shape<dimension + 1>& shape)
        : p(p), bound_id(bound_id), data(data), basis(basis), master(master), shape(shape) {}
};
}

#endif