#ifndef SERIAL_SIMULATION_BASE_HPP
#define SERIAL_SIMULATION_BASE_HPP

#include <memory>

namespace Serial {
class SimulationBase {
public:
    virtual ~SimulationBase()=default;

    virtual void Run() = 0;
    virtual void ComputeL2Residual() = 0;
};

struct SimulationFactory {
    static std::unique_ptr<SimulationBase> Create(const std::string& input_string);
};
}
#endif