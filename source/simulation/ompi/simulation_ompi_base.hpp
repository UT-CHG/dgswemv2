#ifndef SIMULATION_OMPI_BASE_HPP
#define SIMULATION_OMPI_BASE_HPP

#include "general_definitions.hpp"
#include "utilities/is_defined.hpp"

class OMPISimulationBase {
public:
    virtual ~OMPISimulationBase()=default;

    virtual void Run() = 0;
    virtual void ComputeL2Residual() = 0;
    virtual void DestroyPETSc() {}
};

template <typename ProblemType>
class OMPISimulation;

struct OMPISimulationFactory {
    static std::unique_ptr<OMPISimulationBase> Create(const std::string& input_string);
private:
    template <typename ProblemType>
    static std::unique_ptr<OMPISimulationBase> CreateSimulation(const std::string& input_string) {
        return OMPISimulationFactory::CreateSimulationImpl<ProblemType>(input_string, Utilities::is_defined<ProblemType>{});
    }

    template<typename ProblemType>
    static std::unique_ptr<OMPISimulationBase> CreateSimulationImpl(const std::string& input_string, std::true_type) {
        return std::make_unique<OMPISimulation<ProblemType>>(input_string);
    }

    template<typename ProblemType>
    static std::unique_ptr<OMPISimulationBase> CreateSimulationImpl(const std::string& input_string, std::false_type) {
        throw std::runtime_error( "Problem class not supported, please check cmake configuration for proper support\n");
        return nullptr;
}
};

#endif