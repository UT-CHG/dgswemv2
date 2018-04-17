#ifndef SIMULATION_LOAD_BALANCER_HPP
#define SIMULATION_LOAD_BALANCER_HPP

#include "base_models.hpp"
#include "random.hpp"

namespace LoadBalancer {
template <typename ProblemType>
hpx::future<void> initialize_locality_and_world_models(const uint locality_id) {
    return Random<ProblemType>::initialize_locality_and_world_models(locality_id);
}

template <typename ProblemType>
void reset_locality_and_world_models() {
    Random<ProblemType>::reset_locality_and_world_models();
}
}

#define DGSWEMV2_REGISTER_LOAD_BALANCERS(ProblemType)\
    DGSWEMV2_REGISTER_RANDOM_LOAD_BALANCER(ProblemType);
/**/
#endif