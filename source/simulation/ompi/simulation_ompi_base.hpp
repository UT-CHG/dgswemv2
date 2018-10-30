#ifndef SIMULATION_OMPI_BASE_HPP
#define SIMULATION_OMPI_BASE_HPP

#include "general_definitions.hpp"

class OMPISimulationBase {
public:
    virtual ~OMPISimulationBase()=default;

    virtual void Run() = 0;
    virtual void ComputeL2Residual() = 0;
    virtual void DestroyPETSc() {}
};

struct OMPISimulationFactory {
    static std::unique_ptr<OMPISimulationBase> Create(const std::string& input_string);
};

#endif