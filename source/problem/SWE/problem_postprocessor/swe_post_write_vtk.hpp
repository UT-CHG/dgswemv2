#ifndef SWE_POST_WRITE_VTK_HPP
#define SWE_POST_WRITE_VTK_HPP

namespace SWE {
template <typename MeshType>
void write_VTK_data(MeshType& mesh, std::ofstream& raw_data_file) {
    AlignedVector<StatVector<double, SWE::n_variables>> q_point_data;
    AlignedVector<StatVector<double, SWE::n_variables>> q_cell_data;

    AlignedVector<StatVector<double, 1>> aux_point_data;
    AlignedVector<StatVector<double, 1>> aux_cell_data;

    mesh.CallForEachElement([&q_point_data, &q_cell_data, &aux_point_data, &aux_cell_data](auto& elt) {
        elt.WritePointDataVTK(elt.data.state[0].q, q_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].q, q_cell_data);

        elt.WritePointDataVTK(elt.data.state[0].aux, aux_point_data);
        elt.WriteCellDataVTK(elt.data.state[0].aux, aux_cell_data);
    });

    std::vector<uint> elt_id_data;
    std::vector<std::array<bool, 2>> wd_data;

    mesh.CallForEachElement([&elt_id_data, &wd_data](auto& elt) {
        for (uint cell = 0; cell < N_DIV * N_DIV; ++cell) {
            elt_id_data.push_back(elt.GetID());
            wd_data.push_back({elt.data.wet_dry_state.wet, elt.data.wet_dry_state.went_completely_dry});
        }
    });

    raw_data_file << "CELL_DATA " << q_cell_data.size() << std::endl;

    raw_data_file << "SCALARS ze_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : q_cell_data)
        raw_data_file << (float)it[SWE::Variables::ze] << std::endl;

    raw_data_file << "SCALARS qx_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : q_cell_data)
        raw_data_file << (float)it[SWE::Variables::qx] << std::endl;

    raw_data_file << "SCALARS qy_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : q_cell_data)
        raw_data_file << (float)it[SWE::Variables::qy] << std::endl;

    raw_data_file << "SCALARS bath_cell double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : aux_cell_data)
        raw_data_file << (float)it[SWE::Auxiliaries::bath] << std::endl;

    raw_data_file << "SCALARS ID unsigned_int 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (uint& it : elt_id_data)
        raw_data_file << it << std::endl;

    raw_data_file << "SCALARS wet_dry bit 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : wd_data)
        raw_data_file << it[0] << std::endl;

    raw_data_file << "SCALARS went_completely_dry bit 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : wd_data)
        raw_data_file << it[1] << std::endl;

    raw_data_file << "POINT_DATA " << q_point_data.size() << std::endl;

    raw_data_file << "SCALARS ze_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : q_point_data)
        raw_data_file << (float)it[SWE::Variables::ze] << std::endl;

    raw_data_file << "SCALARS qx_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : q_point_data)
        raw_data_file << (float)it[SWE::Variables::qx] << std::endl;

    raw_data_file << "SCALARS qy_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : q_point_data)
        raw_data_file << (float)it[SWE::Variables::qy] << std::endl;

    raw_data_file << "SCALARS bath_point double 1" << std::endl;
    raw_data_file << "LOOKUP_TABLE default" << std::endl;
    for (auto& it : aux_point_data)
        raw_data_file << (float)it[SWE::Auxiliaries::bath] << std::endl;
}
}

#endif