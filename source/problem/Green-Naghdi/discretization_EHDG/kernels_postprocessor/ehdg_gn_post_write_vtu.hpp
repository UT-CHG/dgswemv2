#ifndef EHDG_GN_POST_WRITE_VTU_HPP
#define EHDG_GN_POST_WRITE_VTU_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::write_VTU_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
    std::vector<StatVector<double, GN::n_variables>> q_point_data;
    std::vector<StatVector<double, GN::n_variables>> q_cell_data;

    std::vector<StatVector<double, 1>> aux_point_data;
    std::vector<StatVector<double, 1>> aux_cell_data;

    mesh.CallForEachElement([&q_point_data, &q_cell_data, &aux_point_data, &aux_cell_data](auto& elt) {
        elt.WritePointDataVTK(elt.data.state[0].q, q_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].q, q_cell_data);

        elt.WritePointDataVTK(elt.data.state[0].aux, aux_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].aux, aux_cell_data);
    });

    std::vector<uint> elt_id_data;

    mesh.CallForEachElement([&elt_id_data](auto& elt) {
        for (uint cell = 0; cell < N_DIV * N_DIV; cell++) {
            elt_id_data.push_back(elt.GetID());
        }
    });

    raw_data_file << "\t\t\t<PointData>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"ze_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = q_point_data.begin(); it != q_point_data.end(); ++it)
        raw_data_file << (*it)[GN::Variables::ze] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qx_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = q_point_data.begin(); it != q_point_data.end(); ++it)
        raw_data_file << (*it)[GN::Variables::qx] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qy_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = q_point_data.begin(); it != q_point_data.end(); ++it)
        raw_data_file << (*it)[GN::Variables::qy] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"bath_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = aux_point_data.begin(); it != aux_point_data.end(); ++it)
        raw_data_file << (*it)[GN::Auxiliaries::bath] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t</PointData>\n";

    raw_data_file << "\t\t\t<CellData>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"ze_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = q_cell_data.begin(); it != q_cell_data.end(); ++it)
        raw_data_file << (*it)[GN::Variables::ze] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qx_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = q_cell_data.begin(); it != q_cell_data.end(); ++it)
        raw_data_file << (*it)[GN::Variables::qx] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qy_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = q_cell_data.begin(); it != q_cell_data.end(); ++it)
        raw_data_file << (*it)[GN::Variables::qy] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"bath_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = aux_cell_data.begin(); it != aux_cell_data.end(); ++it)
        raw_data_file << (*it)[GN::Auxiliaries::bath] << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"UInt32\" Name=\"ID\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = elt_id_data.begin(); it != elt_id_data.end(); ++it)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t</CellData>\n";
}
}
}

#endif