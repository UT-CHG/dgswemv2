#ifndef EHDG_GN_POST_WRITE_MODAL_HPP
#define EHDG_GN_POST_WRITE_MODAL_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::write_modal_data(const RKStepper& stepper, ProblemMeshType& mesh, const std::string& output_path) {
    std::vector<std::pair<uint, HybMatrix<double, GN::n_variables>>> modal_q;
    std::vector<std::pair<uint, HybMatrix<double, 1>>> modal_aux;

    mesh.CallForEachElement([&modal_q, &modal_aux](auto& elt) {
        modal_q.push_back(std::make_pair(elt.GetID(), elt.data.state[0].q));
        modal_aux.push_back(std::make_pair(elt.GetID(), elt.data.state[0].aux));
    });

    std::ofstream file;

    std::string file_name = output_path + mesh.GetMeshName() + "_modal_ze.txt";
    if (stepper.GetStep() == 0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.GetTimeAtCurrentStage()) << std::endl;
    for (auto it = modal_q.begin(); it != modal_q.end(); ++it) {
        uint ndof = columns((*it).second);

        for (uint dof = 0; dof < ndof; ++dof) {
            file << (*it).first << ' ' << std::scientific << (*it).second(GN::Variables::ze, dof) << std::endl;
        }
    }

    file.close();

    file_name = output_path + mesh.GetMeshName() + "_modal_qx.txt";
    if (stepper.GetStep() == 0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.GetTimeAtCurrentStage()) << std::endl;
    for (auto it = modal_q.begin(); it != modal_q.end(); ++it) {
        uint ndof = columns((*it).second);

        for (uint dof = 0; dof < ndof; ++dof) {
            file << (*it).first << ' ' << std::scientific << (*it).second(GN::Variables::qx, dof) << std::endl;
        }
    }

    file.close();

    file_name = output_path + mesh.GetMeshName() + "_modal_qy.txt";
    if (stepper.GetStep() == 0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.GetTimeAtCurrentStage()) << std::endl;
    for (auto it = modal_q.begin(); it != modal_q.end(); ++it) {
        uint ndof = columns((*it).second);

        for (uint dof = 0; dof < ndof; ++dof) {
            file << (*it).first << ' ' << std::scientific << (*it).second(GN::Variables::qy, dof) << std::endl;
        }
    }

    file.close();

    file_name = output_path + mesh.GetMeshName() + "_modal_bath.txt";
    if (stepper.GetStep() == 0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.GetTimeAtCurrentStage()) << std::endl;
    for (auto it = modal_aux.begin(); it != modal_aux.end(); ++it) {
        uint ndof = columns((*it).second);

        for (uint dof = 0; dof < ndof; ++dof) {
            file << (*it).first << ' ' << std::scientific << (*it).second(GN::Auxiliaries::bath, dof) << std::endl;
        }
    }

    file.close();
}
}
}

#endif