#ifndef RKDG_SWE_DATA_BOUNDARY_HPP
#define RKDG_SWE_DATA_BOUNDARY_HPP

namespace SWE {
namespace RKDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp) {
        for ( uint var = 0; var < SWE::n_variables; ++var) {
            q_at_gp[var] = DynRowVector<double>(ngp);
            F_hat_at_gp[var] = DynRowVector<double>(ngp);
        }

        for ( uint aux = 0; aux < SWE::n_auxiliaries; ++aux ) {
            aux_at_gp[aux] = DynRowVector<double>(ngp);
        }
    }

    std::array<DynRowVector<double>, SWE::n_variables> q_at_gp;
    std::array<DynRowVector<double>, SWE::n_auxiliaries> aux_at_gp;

    std::array<DynRowVector<double>, SWE::n_variables> F_hat_at_gp;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & q_at_gp
            & aux_at_gp
            & F_hat_at_gp;
        // clang-format on
    }
#endif
};
}
}

#endif
