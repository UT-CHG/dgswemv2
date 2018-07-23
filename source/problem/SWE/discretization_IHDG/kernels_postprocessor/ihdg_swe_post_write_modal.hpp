#ifndef IHDG_SWE_POST_WRITE_MODAL_HPP
#define IHDG_SWE_POST_WRITE_MODAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
void Problem::write_modal_data_kernel(const RKStepper& stepper, ProblemMeshType& mesh, const std::string& output_path) {
    std::vector<std::pair<uint, HybMatrix<double, SWE::n_variables>>> modal_q;
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
    for (auto it = modal_q.begin(); it != modal_q.end(); it++) {
        for (auto itt = column((*it).second, SWE::Variables::ze).begin();
             itt != column((*it).second, SWE::Variables::ze).end();
             itt++) {
            file << (*it).first << ' ' << std::scientific << *itt << std::endl;
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
    for (auto it = modal_q.begin(); it != modal_q.end(); it++) {
        for (auto itt = column((*it).second, SWE::Variables::qx).begin();
             itt != column((*it).second, SWE::Variables::qx).end();
             itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << std::endl;
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
    for (auto it = modal_q.begin(); it != modal_q.end(); it++) {
        for (auto itt = column((*it).second, SWE::Variables::qy).begin();
             itt != column((*it).second, SWE::Variables::qy).end();
             itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << std::endl;
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
    for (auto it = modal_aux.begin(); it != modal_aux.end(); it++) {
        for (auto itt = column((*it).second, SWE::Auxiliaries::bath).begin();
             itt != column((*it).second, SWE::Auxiliaries::bath).end();
             itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << std::endl;
        }
    }

    file.close();
}
}
}

#endif