#ifndef GEOMETRY_ELEMENT_SOA_HPP
#define GEOMETRY_ELEMENT_SOA_HPP

namespace Geometry {

template <typename ElementType, typename ProblemSoA>
class ElementSoA {
public:
    void reserve(uint ndof, uint nstages , uint nelements) {
        std::cout << "Reserving " << nelements << " in element_soa.hpp\n";
        this->problem_data = ProblemSoA(ndof,nstages,nelements);
    }

    //fixme: this signature needs to change down the road to an element accessor type (which is element)
    typename ProblemSoA::AccessorType at(const size_t index) {
        return problem_data.at(index);
    }
private:
    ProblemSoA problem_data;
};
}

#endif