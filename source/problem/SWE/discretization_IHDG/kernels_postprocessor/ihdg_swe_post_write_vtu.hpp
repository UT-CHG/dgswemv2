#ifndef IHDG_SWE_POST_WRITE_VTU_HPP
#define IHDG_SWE_POST_WRITE_VTU_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
void Problem::write_VTU_data_kernel(ProblemMeshType& mesh, std::ofstream& raw_data_file) {
    Array2D<double> cell_data;
    Array2D<double> point_data;

    cell_data.resize(4);
    point_data.resize(4);

    mesh.CallForEachElement([&cell_data, &point_data](auto& elt) {
        elt.WriteCellDataVTK(elt.data.state[0].ze, cell_data[0]);
        elt.WriteCellDataVTK(elt.data.state[0].qx, cell_data[1]);
        elt.WriteCellDataVTK(elt.data.state[0].qy, cell_data[2]);
        elt.WriteCellDataVTK(elt.data.state[0].bath, cell_data[3]);

        elt.WritePointDataVTK(elt.data.state[0].ze, point_data[0]);
        elt.WritePointDataVTK(elt.data.state[0].qx, point_data[1]);
        elt.WritePointDataVTK(elt.data.state[0].qy, point_data[2]);
        elt.WritePointDataVTK(elt.data.state[0].bath, point_data[3]);
    });

    std::vector<uint> elt_id_data;

    mesh.CallForEachElement([&elt_id_data](auto& elt) {
        for (uint cell = 0; cell < N_DIV * N_DIV; cell++) {
            elt_id_data.push_back(elt.GetID());
        }
    });

    raw_data_file << "\t\t\t<PointData>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"ze_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = point_data[0].begin(); it != point_data[0].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qx_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = point_data[1].begin(); it != point_data[1].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qy_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = point_data[2].begin(); it != point_data[2].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"bath_point\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = point_data[3].begin(); it != point_data[3].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t</PointData>\n";

    raw_data_file << "\t\t\t<CellData>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"ze_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = cell_data[0].begin(); it != cell_data[0].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qx_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = cell_data[1].begin(); it != cell_data[1].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"qy_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = cell_data[2].begin(); it != cell_data[2].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"bath_cell\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = cell_data[3].begin(); it != cell_data[3].end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t\t<DataArray type=\"UInt32\" Name=\"ID\" format=\"ascii\">\n\t\t\t\t\t";
    for (auto it = elt_id_data.begin(); it != elt_id_data.end(); it++)
        raw_data_file << *it << ' ';
    raw_data_file << "\n\t\t\t\t</DataArray>\n";

    raw_data_file << "\t\t\t</CellData>\n";
}
}
}

#endif