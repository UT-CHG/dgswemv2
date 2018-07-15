#ifndef RKDG_SWE_POST_WRITE_MODAL_HPP
#define RKDG_SWE_POST_WRITE_MODAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
void Problem::write_modal_data_kernel(const RKStepper& stepper, ProblemMeshType& mesh, const std::string& output_path) {
    std::vector<std::pair<uint, std::vector<StatVector<double, SWE::n_variables>>>> modal_u;
    std::vector<std::pair<uint, std::vector<double>>> modal_bath;

    mesh.CallForEachElement([&modal_u, &modal_bath](auto& elt) {
        modal_u.push_back(std::make_pair(elt.GetID(), elt.data.state[0].q));

        modal_bath.push_back(std::make_pair(elt.GetID(), elt.data.state[0].bath));
    });

    std::ofstream file;

    std::string file_name = output_path + mesh.GetMeshName() + "_modal_ze.txt";
    if (stepper.GetStep() == 0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.GetTimeAtCurrentStage()) << std::endl;
    for (auto it = modal_u.begin(); it != modal_u.end(); it++) {
        for (auto itt = (*it).second.begin(); itt != (*it).second.end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt)[SWE::Variables::ze] << std::endl;
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
    for (auto it = modal_u.begin(); it != modal_u.end(); it++) {
        for (auto itt = (*it).second.begin(); itt != (*it).second.end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt)[SWE::Variables::qx] << std::endl;
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
    for (auto it = modal_u.begin(); it != modal_u.end(); it++) {
        for (auto itt = (*it).second.begin(); itt != (*it).second.end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt)[SWE::Variables::qy] << std::endl;
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
    for (auto it = modal_bath.begin(); it != modal_bath.end(); it++) {
        for (auto itt = (*it).second.begin(); itt != (*it).second.end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << std::endl;
        }
    }

    file.close();
}
}
}

#endif