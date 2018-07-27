#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include "preprocessor/initialize_mesh.hpp"

template <typename ProblemType>
struct DGDiscretization {
    typename ProblemType::ProblemMeshType mesh;

    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input, Writer<ProblemType>& writer) {
        std::tuple<> empty_comm;

        initialize_mesh<ProblemType>(this->mesh, input, empty_comm, writer);
    }
};

#endif