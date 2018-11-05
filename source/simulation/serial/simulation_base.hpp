#ifndef SERIAL_SIMULATION_BASE_HPP
#define SERIAL_SIMULATION_BASE_HPP

#include <memory>

#include "utilities/is_defined.hpp"

namespace Serial {
class SimulationBase {
public:
    virtual ~SimulationBase()=default;

    virtual void Run() = 0;
    virtual void ComputeL2Residual() = 0;
};

template<typename ProblemType>
class Simulation;

struct SimulationFactory {
public:
    static std::unique_ptr<SimulationBase> Create(const std::string& input_string);
private:
    template <typename ProblemType>
    static std::unique_ptr<SimulationBase> CreateSimulation(const std::string& input_string) {
        return SimulationFactory::CreateSimulationImpl<ProblemType>(input_string, Utilities::is_defined<ProblemType>{});
    }

    template<typename ProblemType>
    static std::unique_ptr<SimulationBase> CreateSimulationImpl(const std::string& input_string, std::true_type) {
        return std::make_unique<Simulation<ProblemType>>(input_string);
    }

    template<typename ProblemType>
    static std::unique_ptr<SimulationBase> CreateSimulationImpl(const std::string& input_string, std::false_type) {
        throw std::runtime_error( "Problem class not supported, please check cmake configuration for proper support\n");
        return nullptr;
}
};
}
#endif