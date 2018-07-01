#ifndef SIMULATION_LOAD_BALANCER_HPP
#define SIMULATION_LOAD_BALANCER_HPP

#include "base_model.hpp"
#include "random.hpp"
#include "abstract_load_balancer_factory.hpp" //added for clarity
#include "abstract_load_balancer_factory_impl.hpp"

#define DGSWEMV2_REGISTER_LOAD_BALANCERS(ProblemType)\
    DGSWEMV2_REGISTER_RANDOM_LOAD_BALANCER(ProblemType);
/**/
#endif