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

    void set_J_inv(uint index, const StatMatrix<double,2,2>& J_inv_) {
        for ( uint x = 0; x < 2; ++x ) {
            for ( uint y = 0; y < 2; ++y ) {
                this->J_inv(x,y)(index,index) = J_inv_(x,y);
            }
        }
    }

    void reserve(uint ndof, uint nstages , uint nelements) {
        std::cout << "Reserving " << nelements << " in element_soa.hpp\n";
        this->data = ProblemSoA(ndof, master->ngp, nstages, nelements);
        this->abs_J = DiagonalMatrix<double>(nelements);
        this->inv_abs_J = DiagonalMatrix<double>(nelements);
        this->J_inv = DiagonalMatrix<double, SO::ColumnMajor>(nelements);
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
    decltype(auto) ComputeUgp(const ArrayType& u) {
        return u * this->master->phi_gp;
    }

  decltype(auto) IntegrationDPhi(const StatVector<DynMatrix<double, SO::ColumnMajor>,2>& u_gp) {
//        DynMatrix<double> temp2 =  blaze::evaluate(
//            ;
//        auto temp3 = blaze::evaluate( * temp2);
        
        return blaze::evaluate(abs_J* (transpose(this->J_inv * u_gp) * this->master->int_dphi_fact));
    }

    template <typename ArrayType>
    decltype(auto) ApplyMinv(const ArrayType& rhs) {
        return inv_abs_J*(rhs * master->m_inv);
    }

private:
    MasterType* master = nullptr;

    StatMatrix<DiagonalMatrix<double, SO::ColumnMajor>,2,2> J_inv;
    DiagonalMatrix<double> abs_J;
    DiagonalMatrix<double> inv_abs_J;
};
}

#endif