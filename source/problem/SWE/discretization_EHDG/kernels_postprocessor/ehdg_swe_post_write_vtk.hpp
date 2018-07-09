#ifndef EHDG_SWE_POST_WRITE_VTK_HPP
#define EHDG_SWE_POST_WRITE_VTK_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
void Problem::write_VTK_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
    std::vector<Vector<double, SWE::n_variables>> u_point_data;
    std::vector<Vector<double, SWE::n_variables>> u_cell_data;

    std::vector<double> bath_point_data;
    std::vector<double> bath_cell_data;

    mesh.CallForEachElement([&u_point_data, &u_cell_data, &bath_point_data, &bath_cell_data](auto& elt) {
        elt.WritePointDataVTK(elt.data.state[0].q, u_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].q, u_cell_data);

        elt.WritePointDataVTK(elt.data.state[0].bath, bath_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].bath, bath_cell_data);
    });

    std::vector<uint> elt_id_data;

    mesh.CallForEachElement([&elt_id_data](auto& elt) {
        for (uint cell = 0; cell < N_DIV * N_DIV; cell++) {
            elt_id_data.push_back(elt.GetID());
        }
    });

    raw_data_file << "CELL_DATA " << u_cell_data.size() << std::endl;

    raw_data_file << "SCALARS ze_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = u_cell_data.begin(); it != u_cell_data.end(); it++)
        raw_data_file << (*it)[SWE::Variables::ze] << std::endl;

    raw_data_file << "SCALARS qx_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = u_cell_data.begin(); it != u_cell_data.end(); it++)
        raw_data_file << (*it)[SWE::Variables::qx] << std::endl;

    raw_data_file << "SCALARS qy_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = u_cell_data.begin(); it != u_cell_data.end(); it++)
        raw_data_file << (*it)[SWE::Variables::qy] << std::endl;

    raw_data_file << "SCALARS bath_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = bath_cell_data.begin(); it != bath_cell_data.end(); it++)
        raw_data_file << *it << std::endl;

    raw_data_file << "SCALARS ID unsigned_int 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = elt_id_data.begin(); it != elt_id_data.end(); it++)
        raw_data_file << *it << std::endl;

    raw_data_file << "POINT_DATA " << u_point_data.size() << std::endl;

    raw_data_file << "SCALARS ze_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = u_point_data.begin(); it != u_point_data.end(); it++)
        raw_data_file << (*it)[SWE::Variables::ze] << std::endl;

    raw_data_file << "SCALARS qx_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = u_point_data.begin(); it != u_point_data.end(); it++)
        raw_data_file << (*it)[SWE::Variables::qx] << std::endl;

    raw_data_file << "SCALARS qy_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = u_point_data.begin(); it != u_point_data.end(); it++)
        raw_data_file << (*it)[SWE::Variables::qx] << std::endl;

    raw_data_file << "SCALARS bath_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto it = bath_point_data.begin(); it != bath_point_data.end(); it++)
        raw_data_file << *it << std::endl;
}
}
}

#endif