#ifndef CLASS_RAW_BOUNDARY_HPP
#define CLASS_RAW_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename DataType>
class RawBoundary {
  public:
    uint p;
    uint bound_id;
    std::vector<uint> node_ID;

    DataType& data;

    Basis::Basis<dimension + 1>& basis;
    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

  public:
    RawBoundary(uint p,
                uint bound_id,
                std::vector<uint>& node_ID,
                DataType& data,
                Basis::Basis<dimension + 1>& basis,
                Master::Master<dimension + 1>& master,
                Shape::Shape<dimension + 1>& shape)
        : p(p), bound_id(bound_id), node_ID(node_ID), data(data), basis(basis), master(master), shape(shape) {}
};
}

#endif