#ifndef GEOMETRY_ELEMENT_SOA_HPP
#define GEOMETRY_ELEMENT_SOA_HPP

namespace Geometry {

template <typename ElementType, typename ProblemSoA>
class ElementSoA {
private:
    using MasterType = typename ElementType::ElementMasterType;
public:
    ElementSoA()=default;
    ElementSoA(MasterType& master_) : master(&master_) {}

    void reserve(uint ndof, uint nstages , uint nelements) {
        std::cout << "Reserving " << nelements << " in element_soa.hpp\n";
        this->problem_data = ProblemSoA(ndof,nstages,nelements);
        this->abs_J = DynVector<double>(nelements);
    }

    //fixme: this signature needs to change down the road to an element accessor type (which is element)
    template <typename... Args>
    ElementType at(const size_t index,
                   const uint ID,
                   Args&&... args) {
        assert(this->master);
        return ElementType(ID,
                           *master,
                           this->problem_data.at(index),
                           this->abs_J[index],
                           std::forward<Args>(args)...);
    }
private:
    ProblemSoA problem_data;
    MasterType* master = nullptr;

    DynVector<double> abs_J;
};
}

#endif