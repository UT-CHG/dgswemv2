#ifndef EHDG_GN_PROC_OMPI_STAGE_HPP
#define EHDG_GN_PROC_OMPI_STAGE_HPP

#include "general_definitions.hpp"

#include "ehdg_gn_kernels_processor.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {}
}
}

#endif