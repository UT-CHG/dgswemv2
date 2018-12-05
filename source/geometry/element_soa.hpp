#ifndef GEOMETRY_ELEMENT_SOA_HPP
#define GEOMETRY_ELEMENT_SOA_HPP

namespace Geometry {

template <typename Element>
class ElementSoA;

template <uint dimension, typename MasterType, typename ShapeType, typename AccessorType>
class ElementSoA<Element<dimension, MasterType, ShapeType, AccessorType>> {
public:
private:
};
}

#endif