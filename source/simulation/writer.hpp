#ifndef WRITER_HPP
#define WRITER_HPP

#include "preprocessor/input_parameters.hpp"

template <typename ProblemType>
class Writer {
  public:
    Writer() = default;
    Writer(const InputParameters<typename ProblemType::InputType>& input,
           const uint locality_id = std::numeric_limits<uint>::max(),
           const uint submesh_id = std::numeric_limits<uint>::max());

    void WriteFirstStep(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh);

    void WriteOutput(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh);

    std::ofstream& get_log_file() { return log_file; }

  private:
    std::string output_path;

    std::string log_file_name;
    std::ofstream log_file;

    bool writing_vtk;
    int vtk_output_frequency;
    std::string vtk_file_name_geom;
    std::string vtk_file_name_raw;

    bool writing_modal_output;
    int modal_output_frequency;

    void initialize_mesh_VTK_geometry(typename ProblemType::ProblemMeshType& mesh);
};

template <typename ProblemType>
Writer<ProblemType>::Writer(const InputParameters<typename ProblemType::InputType>& input,
                            const uint locality_id,
                            const uint submesh_id)
    : output_path(input.writer_input.output_path),
      writing_vtk(input.writer_input.writing_vtk),
      vtk_output_frequency(input.writer_input.vtk_output_frequency),
      writing_modal_output(input.writer_input.writing_modal_output),
      modal_output_frequency(input.writer_input.modal_output_frequency) {
    if (!input.writer_input.log_file_name.empty()) {
        log_file_name = output_path + input.writer_input.log_file_name;
        if (locality_id != std::numeric_limits<uint>::max() && submesh_id != std::numeric_limits<uint>::max()) {
            log_file_name = log_file_name + '_' + std::to_string(locality_id) + '_' + std::to_string(submesh_id);
        }

        log_file = std::ofstream(log_file_name, std::ofstream::out);

        if (!log_file) {
            std::cerr << "Error in opening log file, presumably the output directory does not exists.\n";
        }

        log_file << "Starting simulation with p=" << input.polynomial_order << " for " << input.mesh_data.mesh_name
                 << " mesh" << std::endl << std::endl;
    }
}

template <typename ProblemType>
void Writer<ProblemType>::WriteFirstStep(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh) {

    vtk_file_name_geom = output_path + mesh.GetMeshName() + "_geometry.vtk";
    vtk_file_name_raw = output_path + mesh.GetMeshName() + "_raw_data.vtk";

    this->initialize_mesh_VTK_geometry(mesh);
    this->WriteOutput(stepper, mesh);
}

template <typename ProblemType>
void Writer<ProblemType>::WriteOutput(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh) {

    if (writing_vtk && (stepper.get_step() % vtk_output_frequency == 0)) {

        std::string file_name = output_path + mesh.GetMeshName() + "_raw_data.vtk";
        std::ofstream raw_data_file(vtk_file_name_raw);

        ProblemType::write_VTK_data_kernel(mesh, raw_data_file);

        std::ifstream file_geom(vtk_file_name_geom, std::ios_base::binary);
        std::ifstream file_data(vtk_file_name_raw, std::ios_base::binary);

        uint step = stepper.get_step();
        std::string file_name_merge = output_path + mesh.GetMeshName() + "_data_" + std::to_string(step) + ".vtk";
        std::ofstream file_merge(file_name_merge, std::ios_base::binary);

        file_merge << file_geom.rdbuf() << file_data.rdbuf();
        file_merge.close();
    }

    if (writing_modal_output && (stepper.get_step() % modal_output_frequency == 0)) {
        ProblemType::write_modal_data_kernel(stepper, mesh, output_path);
    }
}

template <typename ProblemType>
void Writer<ProblemType>::initialize_mesh_VTK_geometry(typename ProblemType::ProblemMeshType& mesh) {
    std::vector<Point<3>> points;
    Array2D<uint> cells;

    mesh.CallForEachElement([&points, &cells](auto& elem) { elem.InitializeVTK(points, cells); });

    std::ofstream file(vtk_file_name_geom);

    file << "# vtk DataFile Version 3.0\n";
    file << "OUTPUT DATA\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << points.size() << " double\n";

    for (auto it = points.begin(); it != points.end(); it++) {
        file << (*it)[0] << '\t' << (*it)[1] << '\t' << (*it)[2] << '\n';
    }

    uint n_cell_entries = 0;
    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
            case VTKElementTypes::straight_triangle:
                n_cell_entries += 4;
                break;
            default:
                printf("\n");
                printf("MESH InitializeVTK - Fatal error!\n");
                printf("Undefined cell type = %d\n", (*it)[0]);
                exit(1);
        }
    }

    file << "CELLS " << cells.size() << ' ' << n_cell_entries << '\n';

    uint n_nodes;

    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
            case VTKElementTypes::straight_triangle:
                file << 3 << '\t';
                n_nodes = 3;
                break;
            default:
                printf("\n");
                printf("MESH InitializeVTK - Fatal error!\n");
                printf("Undefined cell type = %d\n", (*it)[0]);
                exit(1);
        }

        for (uint i = 1; i <= n_nodes; i++) {
            file << (*it)[i] << '\t';
        }
        file << '\n';
    }

    file << "CELL_TYPES " << cells.size() << '\n';

    for (auto it = cells.begin(); it != cells.end(); it++) {
        file << (*it)[0] << '\n';
    }
}
#endif