#ifndef IHDG_SWE_POST_WRITE_MODAL_HPP
#define IHDG_SWE_POST_WRITE_MODAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
void Problem::write_modal_data_kernel(const RKStepper& stepper, ProblemMeshType& mesh, const std::string& output_path) {
    std::vector<std::pair<uint, Array2D<double>>> modal_data;

    mesh.CallForEachElement([&modal_data](auto& elt) {
        modal_data.push_back(std::make_pair(
            elt.GetID(),
            Array2D<double>{elt.data.state[0].ze, elt.data.state[0].qx, elt.data.state[0].qy, elt.data.state[0].bath}));
    });

    std::ofstream file;

    std::string file_name = output_path + mesh.GetMeshName() + "_modal_ze.txt";
    if (stepper.GetStep() == 0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.GetTimeAtCurrentStage()) << std::endl;
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[0].begin(); itt != (*it).second[0].end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << std::endl;
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
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[1].begin(); itt != (*it).second[1].end(); itt++) {
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
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[2].begin(); itt != (*it).second[2].end(); itt++) {
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
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[3].begin(); itt != (*it).second[3].end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << std::endl;
        }
    }

    file.close();
}
}
}

#endif