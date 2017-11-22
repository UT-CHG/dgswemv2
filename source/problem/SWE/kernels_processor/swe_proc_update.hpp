#ifndef SWE_PROC_UPDATE_HPP
#define SWE_PROC_UPDATE_HPP

namespace SWE {
template <typename ElementType>
void Problem::update_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();
    const double dt = stepper.get_dt();

    auto& state = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    curr_state.rhs_ze = elt.ApplyMinv(curr_state.rhs_ze);
    curr_state.rhs_qx = elt.ApplyMinv(curr_state.rhs_qx);
    curr_state.rhs_qy = elt.ApplyMinv(curr_state.rhs_qy);

    std::fill(next_state.ze.begin(), next_state.ze.end(), 0);
    std::fill(next_state.qx.begin(), next_state.qx.end(), 0);
    std::fill(next_state.qy.begin(), next_state.qy.end(), 0);

    for (uint s = 0; s <= stage; ++s) {
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            next_state.ze[dof] +=
                stepper.ark[stage][s] * state[s].ze[dof] + dt * stepper.brk[stage][s] * state[s].rhs_ze[dof];

            next_state.qx[dof] +=
                stepper.ark[stage][s] * state[s].qx[dof] + dt * stepper.brk[stage][s] * state[s].rhs_qx[dof];

            next_state.qy[dof] +=
                stepper.ark[stage][s] * state[s].qy[dof] + dt * stepper.brk[stage][s] * state[s].rhs_qy[dof];
        }
    }
}

template <typename ElementType>
void Problem::swap_states_kernel(const Stepper& stepper, ElementType& elt) {
    uint n_stages = stepper.get_num_stages();
    auto& state = elt.data.state;

    std::swap(state[0].ze, state[n_stages].ze);
    std::swap(state[0].qx, state[n_stages].qx);
    std::swap(state[0].qy, state[n_stages].qy);
}

template <typename ElementType>
void Problem::scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt) {
    uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];

    for (auto& ze_mode : state.ze) {
        if (isnan(ze_mode)) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qx_mode : state.qx) {
        if (isnan(qx_mode)) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qy_mode : state.qy) {
        if (isnan(qy_mode)) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_ze_mode : state.rhs_ze) {
        if (isnan(rhs_ze_mode)) {
            std::cerr << "Error: found isnan rhs_ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qx_mode : state.rhs_qx) {
        if (isnan(rhs_qx_mode)) {
            std::cerr << "Error: found isnan rhs_qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qy_mode : state.rhs_qy) {
        if (isnan(rhs_qy_mode)) {
            std::cerr << "Error: found isnan rhs_qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }
}
}

#endif
