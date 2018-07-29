#ifndef EHDG_GN_POST_WRITE_VTK_HPP
#define EHDG_GN_POST_WRITE_VTK_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::write_VTK_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
    std::vector<StatVector<double, GN::n_variables>> q_point_data;
    std::vector<StatVector<double, GN::n_variables>> q_cell_data;

    std::vector<StatVector<double, 1>> aux_point_data;
    std::vector<StatVector<double, 1>> aux_cell_data;

    mesh.CallForEachElement([&q_point_data, &q_cell_data, &aux_point_data, &aux_cell_data](auto& elt) {
        elt.WritePointDataVTK(elt.data.state[0].q, q_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].q, q_cell_data);

        elt.WritePointDataVTK(row(elt.data.state[0].aux, GN::Auxiliaries::bath), aux_point_data);
        elt.WriteCellDataVTK(row(elt.data.state[0].aux, GN::Auxiliaries::bath), aux_cell_data);
    });

    std::vector<uint> elt_id_data;

    mesh.CallForEachElement([&elt_id_data](auto& elt) {
        for (uint cell = 0; cell < N_DIV * N_DIV; cell++) {
            elt_id_data.push_back(elt.GetID());
        }
    });

    raw_data_file << "CELL_DATA " << q_cell_data.size() << std::endl;

    raw_data_file << "SCALARS ze_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = q_cell_data.begin(); it != q_cell_data.end(); it++)
        raw_data_file << (*it)[GN::Variables::ze] << std::endl;

    raw_data_file << "SCALARS qx_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = q_cell_data.begin(); it != q_cell_data.end(); it++)
        raw_data_file << (*it)[GN::Variables::qx] << std::endl;

    raw_data_file << "SCALARS qy_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = q_cell_data.begin(); it != q_cell_data.end(); it++)
        raw_data_file << (*it)[GN::Variables::qy] << std::endl;

    raw_data_file << "SCALARS bath_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = aux_cell_data.begin(); it != aux_cell_data.end(); it++)
        raw_data_file << (*it)[GN::Auxiliaries::bath] << std::endl;

    raw_data_file << "SCALARS ID unsigned_int 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = elt_id_data.begin(); it != elt_id_data.end(); it++)
        raw_data_file << *it << std::endl;

    raw_data_file << "POINT_DATA " << q_point_data.size() << std::endl;

    raw_data_file << "SCALARS ze_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = q_point_data.begin(); it != q_point_data.end(); it++)
        raw_data_file << (*it)[GN::Variables::ze] << std::endl;

    raw_data_file << "SCALARS qx_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = q_point_data.begin(); it != q_point_data.end(); it++)
        raw_data_file << (*it)[GN::Variables::qx] << std::endl;

    raw_data_file << "SCALARS qy_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = q_point_data.begin(); it != q_point_data.end(); it++)
        raw_data_file << (*it)[GN::Variables::qx] << std::endl;

    raw_data_file << "SCALARS bath_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = aux_point_data.begin(); it != aux_point_data.end(); it++)
        raw_data_file << (*it)[GN::Auxiliaries::bath] << std::endl;
}
}
}

#endif