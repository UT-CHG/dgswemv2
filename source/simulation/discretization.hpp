#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include "general_definitions.hpp"
#include "preprocessor/initialize_mesh.hpp"
#include "preprocessor/initialize_mesh_skeleton.hpp"

template <typename ProblemType>
struct DGDiscretization {
    typename ProblemType::ProblemMeshType mesh;

    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input,
                    typename ProblemType::ProblemWriterType& writer) {
        std::tuple<> empty_comm;

        initialize_mesh<ProblemType>(this->mesh, input, empty_comm, writer);
    }

    template <typename CommunicatorType>
    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input,
                    CommunicatorType& communicator,
                    typename ProblemType::ProblemWriterType& writer) {
        initialize_mesh<ProblemType>(this->mesh, input, communicator, writer);
    }

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        ar& mesh;
    }
#endif
};

template <typename ProblemType>
struct HDGDiscretization {
    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemMeshSkeletonType mesh_skeleton;
    typename ProblemType::ProblemGlobalDataType global_data;

    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input,
                    typename ProblemType::ProblemWriterType& writer) {
        std::tuple<> empty_comm;

        initialize_mesh<ProblemType>(this->mesh, input, empty_comm, writer);
        initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, writer);
    }

    template <typename CommunicatorType>
    void initialize(InputParameters<typename ProblemType::ProblemInputType>& input,
                    CommunicatorType& communicator,
                    typename ProblemType::ProblemWriterType& writer) {
        initialize_mesh<ProblemType>(this->mesh, input, communicator, writer);
        initialize_mesh_skeleton<ProblemType>(this->mesh, this->mesh_skeleton, writer);
    }

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive&, unsigned) {
        throw std::logic_error("Error: Serialization of HDGDiscretization not supported");
    }
#endif
};

#endif