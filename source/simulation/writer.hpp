#ifndef WRITER_HPP
#define WRITER_HPP

#include "../preprocessor/input_parameters.hpp"

template <typename ProblemType>
class Writer {
  private:
    bool        writing_output;
    std::string output_path;

    bool          writing_log_file;
    bool          verbose_log_file;
    std::string   log_file_name;
    std::ofstream log_file;

    bool        writing_vtk_output;
    int         vtk_output_frequency;
    std::string vtk_file_name_geom;
    std::string vtk_file_name_raw;

    bool writing_modal_output;
    int  modal_output_frequency;

  public:
    Writer() = default;
    Writer(const InputParameters<typename ProblemType::ProblemInputType>& input);
    Writer(const InputParameters<typename ProblemType::ProblemInputType>& input,
           const uint                                                     locality_id,
           const uint                                                     submesh_id);

    bool           WritingLog() { return this->writing_log_file; }
    bool           WritingVerboseLog() { return (this->writing_log_file && this->verbose_log_file); }
    std::ofstream& GetLogFile() { return this->log_file; }
    void           StartLog();

    bool WritingOutput() { return this->writing_output; }
    void WriteFirstStep(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh);
    void WriteOutput(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh);

  private:
    void InitializeMeshGeometryVTK(typename ProblemType::ProblemMeshType& mesh);
};

template <typename ProblemType>
Writer<ProblemType>::Writer(const InputParameters<typename ProblemType::ProblemInputType>& input)
    : writing_output(input.writer_input.writing_output),
      output_path(input.writer_input.output_path),
      writing_log_file(input.writer_input.writing_log_file),
      verbose_log_file(input.writer_input.verbose_log_file),
      writing_vtk_output(input.writer_input.writing_vtk_output),
      vtk_output_frequency(input.writer_input.vtk_output_frequency),
      writing_modal_output(input.writer_input.writing_modal_output),
      modal_output_frequency(input.writer_input.modal_output_frequency) {
    if (this->writing_log_file) {
        this->log_file_name = this->output_path + input.writer_input.log_file_name;
    }
}

template <typename ProblemType>
Writer<ProblemType>::Writer(const InputParameters<typename ProblemType::ProblemInputType>& input,
                            const uint                                                     locality_id,
                            const uint                                                     submesh_id)
    : Writer(input) {
    if (this->writing_log_file) {
        this->log_file_name = this->output_path + input.writer_input.log_file_name + '_' + std::to_string(locality_id) +
                              '_' + std::to_string(submesh_id);
    }
}

template <typename ProblemType>
void Writer<ProblemType>::StartLog() {
    this->log_file = std::ofstream(this->log_file_name, std::ofstream::out);

    if (!this->log_file) {
        std::cerr << "Error in opening log file, presumably the output directory does not exists.\n";
    }
}

template <typename ProblemType>
void Writer<ProblemType>::WriteFirstStep(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh) {
    if (this->writing_vtk_output) {
        this->vtk_file_name_geom = this->output_path + mesh.GetMeshName() + "_geometry.vtk";
        this->vtk_file_name_raw  = this->output_path + mesh.GetMeshName() + "_raw_data.vtk";

        this->InitializeMeshGeometryVTK(mesh);
    }

    this->WriteOutput(stepper, mesh);
}

template <typename ProblemType>
void Writer<ProblemType>::WriteOutput(const Stepper& stepper, typename ProblemType::ProblemMeshType& mesh) {
    if (this->writing_vtk_output && (stepper.GetStep() % this->vtk_output_frequency == 0)) {
        std::ofstream raw_data_file(this->vtk_file_name_raw);

        ProblemType::write_VTK_data_kernel(mesh, raw_data_file);

        std::ifstream file_geom(this->vtk_file_name_geom, std::ios_base::binary);
        std::ifstream file_data(this->vtk_file_name_raw, std::ios_base::binary);

        uint        step            = stepper.GetStep();
        std::string file_name_merge = this->output_path + mesh.GetMeshName() + "_data_" + std::to_string(step) + ".vtk";
        std::ofstream file_merge(file_name_merge, std::ios_base::binary);

        file_merge << file_geom.rdbuf() << file_data.rdbuf();
        file_merge.close();
    }

    if (this->writing_modal_output && (stepper.GetStep() % this->modal_output_frequency == 0)) {
        ProblemType::write_modal_data_kernel(stepper, mesh, this->output_path);
    }
}

template <typename ProblemType>
void Writer<ProblemType>::InitializeMeshGeometryVTK(typename ProblemType::ProblemMeshType& mesh) {
    std::vector<Point<3>> points;
    Array2D<uint>         cells;

    mesh.CallForEachElement([&points, &cells](auto& elem) { elem.InitializeVTK(points, cells); });

    std::ofstream file(this->vtk_file_name_geom);

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