#ifndef SWE_DATA_SOURCE_HPP
#define SWE_DATA_SOURCE_HPP

namespace SWE {
struct Source {
    Source() = default;
    Source(const uint nnode)
        : parsed_meteo_data(nnode), tau_s(nnode), p_atm(nnode), tide_pot(nnode), manning_n(nnode) {}

    double coriolis_f = 0.0;

    bool manning          = false;
    double g_manning_n_sq = 0.0;

    std::vector<std::vector<double>*> parsed_meteo_data;

    AlignedVector<StatVector<double, SWE::n_dimensions>> tau_s;
    std::vector<double> p_atm;

    std::vector<double> tide_pot;
    std::vector<double> manning_n;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & coriolis_f
            & manning
            & g_manning_n_sq
            & tau_s
            & p_atm
            & tide_pot
            & manning_n;
        // clang-format on
    }
#endif
};
}

#endif