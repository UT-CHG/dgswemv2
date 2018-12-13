#ifndef GEOMETRY_ELEMENT_SOA_HPP
#define GEOMETRY_ELEMENT_SOA_HPP

namespace Geometry {

template <typename ElementType, typename ProblemSoA>
class ElementSoA {
public:
    using AccessorType = ElementType;
private:
    using MasterType = typename ElementType::ElementMasterType;
public:
    ElementSoA()=default;
    ElementSoA(MasterType& master_) : master(&master_) {}

    void set_abs_J(uint index, double abs_J_) {
        this->abs_J(index,index) = abs_J_;
        this->inv_abs_J(index,index) = 1/abs_J_;
    }

    void reserve(uint ndof, uint nstages , uint nelements) {
        std::cout << "Reserving " << nelements << " in element_soa.hpp\n";
        this->data = ProblemSoA(ndof,nstages,nelements);
        this->abs_J = DiagonalMatrix<double>(nelements);
        this->inv_abs_J = DiagonalMatrix<double>(nelements);
    }

    //fixme: this signature needs to change down the road to an element accessor type (which is element)
    template <typename... Args>
    AccessorType at(const size_t index,
                   Args&&... args) {
        assert(this->master);
        return AccessorType(*master,
                            this->data.at(index),
                            std::forward<Args>(args)...);
    }

    ProblemSoA data;


    template <typename ArrayType>
    decltype(auto) ApplyMinv(const ArrayType& rhs) {
        return inv_abs_J*(rhs * master->m_inv);
    }

private:
    MasterType* master = nullptr;

    DiagonalMatrix<double> abs_J;
    DiagonalMatrix<double> inv_abs_J;
};
}

#endif