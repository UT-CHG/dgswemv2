#ifndef WRITER_HPP
#define WRITER_HPP

#include "general_definitions.hpp"
#include "preprocessor/input_parameters.hpp"

template <typename ProblemType>
class Writer {
  private:
    bool writing_output;
    std::string output_path;

    bool writing_log_file;
    bool verbose_log_file;
    std::string log_file_name;
    mutable std::ofstream log_file;

    bool writing_vtk_output;
    uint vtk_output_frequency;
    std::string vtk_file_name_geom;
    std::string vtk_file_name_raw;

    bool writing_vtu_output;
    uint vtu_output_frequency;
    std::string vtu_file_name_geom_head;
    std::string vtu_file_name_geom_foot;
    std::string vtu_file_name_raw;

    bool writing_modal_output;
    uint modal_output_frequency;

    uint version;

  public:
    Writer() = default;
    Writer(const WriterInput& writer_input);
    Writer(const WriterInput& writer_input, const uint locality_id, const uint submesh_id);

    Writer(Writer&& rhs) = default;

    Writer& operator=(Writer&& rhs) = default;

    bool WritingLog() { return this->writing_log_file; }
    bool WritingVerboseLog() const { return (this->writing_log_file && this->verbose_log_file); }
    std::ofstream& GetLogFile() const { return this->log_file; }
    void StartLog();

    bool WritingOutput() { return this->writing_output; }
    void WriteFirstStep(const typename ProblemType::ProblemStepperType& stepper,
                        typename ProblemType::ProblemMeshType& mesh);
    void WriteOutput(const typename ProblemType::ProblemStepperType& stepper,
                     typename ProblemType::ProblemMeshType& mesh);

  private:
    void InitializeMeshGeometryVTK(typename ProblemType::ProblemMeshType& mesh);
    void InitializeMeshGeometryVTU(typename ProblemType::ProblemMeshType& mesh);

  public:
#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & writing_output
            & output_path
            & writing_log_file
            & verbose_log_file
            & log_file_name
            & writing_vtk_output
            &vtk_output_frequency
            & vtk_file_name_geom
            & vtk_file_name_raw
            & writing_modal_output
            & modal_output_frequency
            & version;
        // clang-format on
    }
#endif
};

template <typename ProblemType>
Writer<ProblemType>::Writer(const WriterInput& writer_input)
    : writing_output(writer_input.writing_output),
      output_path(writer_input.output_path),
      writing_log_file(writer_input.writing_log_file),
      verbose_log_file(writer_input.verbose_log_file),
      writing_vtk_output(writer_input.writing_vtk_output),
      vtk_output_frequency(writer_input.vtk_output_freq_step),
      writing_vtu_output(writer_input.writing_vtu_output),
      vtu_output_frequency(writer_input.vtu_output_freq_step),
      writing_modal_output(writer_input.writing_modal_output),
      modal_output_frequency(writer_input.modal_output_freq_step),
      version(0) {
    if (this->writing_log_file) {
        this->log_file_name = this->output_path + writer_input.log_file_name;
    }
}

template <typename ProblemType>
Writer<ProblemType>::Writer(const WriterInput& writer_input, const uint locality_id, const uint submesh_id)
    : Writer(writer_input) {
    if (this->writing_log_file) {
        this->log_file_name = this->output_path + writer_input.log_file_name + '_' + std::to_string(locality_id) + '_' +
                              std::to_string(submesh_id);
    }
}

template <typename ProblemType>
void Writer<ProblemType>::StartLog() {
    this->log_file = std::ofstream(this->log_file_name + '_' + std::to_string(version++));

    if (!this->log_file) {
        std::cerr << "Error in opening log file, presumably the output directory does not exists.\n";
    }
}

template <typename ProblemType>
void Writer<ProblemType>::WriteFirstStep(const typename ProblemType::ProblemStepperType& stepper,
                                         typename ProblemType::ProblemMeshType& mesh) {
    if (this->writing_vtk_output) {
        this->vtk_file_name_geom = this->output_path + mesh.GetMeshName() + "_geometry.vtk";
        this->vtk_file_name_raw  = this->output_path + mesh.GetMeshName() + "_raw_data.vtk";

        this->InitializeMeshGeometryVTK(mesh);
    }

    if (this->writing_vtu_output) {
        this->vtu_file_name_geom_head = this->output_path + mesh.GetMeshName() + "_geom_head.vtu";
        this->vtu_file_name_geom_foot = this->output_path + mesh.GetMeshName() + "_geom_foot.vtu";
        this->vtu_file_name_raw       = this->output_path + mesh.GetMeshName() + "_raw_data.vtu";

        this->InitializeMeshGeometryVTU(mesh);
    }

    this->WriteOutput(stepper, mesh);
}

template <typename ProblemType>
void Writer<ProblemType>::WriteOutput(const typename ProblemType::ProblemStepperType& stepper,
                                      typename ProblemType::ProblemMeshType& mesh) {
    if (this->writing_vtk_output && (stepper.GetStep() % this->vtk_output_frequency == 0)) {
        std::ofstream raw_data_file(this->vtk_file_name_raw);

        ProblemType::write_VTK_data(mesh, raw_data_file);

        raw_data_file.close();

        if (!Utilities::file_exists(this->vtk_file_name_geom)) {
            throw std::logic_error("Fatal Error: vtk geometry data file " + this->vtk_file_name_geom +
                                   " was not found!\n");
        }

        if (!Utilities::file_exists(this->vtk_file_name_raw)) {
            throw std::logic_error("Fatal Error: vtk raw data file " + this->vtk_file_name_raw + " was not found!\n");
        }

        std::ifstream file_geom(this->vtk_file_name_geom, std::ios_base::binary);
        std::ifstream file_data(this->vtk_file_name_raw, std::ios_base::binary);

        uint step                   = stepper.GetStep();
        std::string file_name_merge = this->output_path + mesh.GetMeshName() + "_data_" + std::to_string(step) + ".vtk";
        std::ofstream file_merge(file_name_merge, std::ios_base::binary);

        file_merge << file_geom.rdbuf() << file_data.rdbuf();

        file_merge.close();
        file_geom.close();
        file_data.close();
    }

    if (this->writing_vtu_output && (stepper.GetStep() % this->vtu_output_frequency == 0)) {
        std::ofstream raw_data_file(this->vtu_file_name_raw);

        ProblemType::write_VTU_data(mesh, raw_data_file);

        raw_data_file.close();

        if (!Utilities::file_exists(this->vtu_file_name_geom_head)) {
            throw std::logic_error("Fatal Error: vtu geometry head data file " + this->vtu_file_name_geom_head +
                                   " was not found!\n");
        }

        if (!Utilities::file_exists(this->vtu_file_name_geom_foot)) {
            throw std::logic_error("Fatal Error: vtu geometry foot data file " + this->vtu_file_name_geom_foot +
                                   " was not found!\n");
        }

        if (!Utilities::file_exists(this->vtu_file_name_raw)) {
            throw std::logic_error("Fatal Error: vtu raw data file " + this->vtu_file_name_raw + " was not found!\n");
        }

        std::ifstream file_geom_head(this->vtu_file_name_geom_head, std::ios_base::binary);
        std::ifstream file_geom_foot(this->vtu_file_name_geom_foot, std::ios_base::binary);
        std::ifstream file_data(this->vtu_file_name_raw, std::ios_base::binary);

        uint step                   = stepper.GetStep();
        std::string file_name_merge = this->output_path + mesh.GetMeshName() + "_data_" + std::to_string(step) + ".vtu";
        std::ofstream file_merge(file_name_merge, std::ios_base::binary);

        file_merge << file_geom_head.rdbuf() << file_data.rdbuf() << file_geom_foot.rdbuf();

        file_merge.close();
        file_geom_head.close();
        file_geom_foot.close();
        file_data.close();
    }

    if (this->writing_modal_output && (stepper.GetStep() % this->modal_output_frequency == 0)) {
        ProblemType::write_modal_data(stepper, mesh, this->output_path);
    }
}

template <typename ProblemType>
void Writer<ProblemType>::InitializeMeshGeometryVTK(typename ProblemType::ProblemMeshType& mesh) {
    AlignedVector<Point<3>> points;
    Array2D<uint> cells;

    mesh.CallForEachElement([&points, &cells](auto& elem) { elem.InitializeVTK(points, cells); });

    std::ofstream file(this->vtk_file_name_geom);

    file << "# vtk DataFile Version 3.0\n";
    file << "OUTPUT DATA\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << points.size() << " double\n";

    for (auto it = points.begin(); it != points.end(); ++it) {
        file << (*it)[0] << '\t' << (*it)[1] << '\t' << (*it)[2] << '\n';
    }

    uint n_cell_entries = 0;
    for (auto it = cells.begin(); it != cells.end(); ++it) {
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

    for (auto it = cells.begin(); it != cells.end(); ++it) {
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

        for (uint i = 1; i <= n_nodes; ++i) {
            file << (*it)[i] << '\t';
        }
        file << '\n';
    }

    file << "CELL_TYPES " << cells.size() << '\n';

    for (auto it = cells.begin(); it != cells.end(); ++it) {
        file << (*it)[0] << '\n';
    }

    file.close();
}

template <typename ProblemType>
void Writer<ProblemType>::InitializeMeshGeometryVTU(typename ProblemType::ProblemMeshType& mesh) {
    AlignedVector<Point<3>> points;
    Array2D<uint> cells;

    mesh.CallForEachElement([&points, &cells](auto& elem) { elem.InitializeVTK(points, cells); });

    std::ofstream file(this->vtu_file_name_geom_head);

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\">\n";
    file << "\t<UnstructuredGrid>\n";
    file << "\t\t<Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";

    file.close();

    file = std::ofstream(this->vtu_file_name_geom_foot);

    file << "\t\t\t<Points>\n";
    file << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n\t\t\t\t\t";

    for (auto it = points.begin(); it != points.end(); ++it) {
        file << (*it)[0] << ' ' << (*it)[1] << ' ' << (*it)[2] << ' ';
    }

    file << "\n\t\t\t\t</DataArray>\n";
    file << "\t\t\t</Points>\n";

    file << "\t\t\t<Cells>\n";

    file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n\t\t\t\t\t";

    uint n_nodes;

    for (auto it = cells.begin(); it != cells.end(); ++it) {
        switch ((*it)[0]) {
            case VTKElementTypes::straight_triangle:
                n_nodes = 3;
                break;
            default:
                printf("\n");
                printf("MESH InitializeVTK - Fatal error!\n");
                printf("Undefined cell type = %d\n", (*it)[0]);
                exit(1);
        }

        for (uint i = 1; i <= n_nodes; ++i) {
            file << (*it)[i] << ' ';
        }
    }

    file << "\n\t\t\t\t</DataArray>\n";

    file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

    uint offset = 0;

    for (auto it = cells.begin(); it != cells.end(); ++it) {
        switch ((*it)[0]) {
            case VTKElementTypes::straight_triangle:
                offset += 3;
                file << offset << ' ';
                break;
            default:
                printf("\n");
                printf("MESH InitializeVTK - Fatal error!\n");
                printf("Undefined cell type = %d\n", (*it)[0]);
                exit(1);
        }
    }

    file << "\n\t\t\t\t</DataArray>\n";

    file << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";

    for (auto it = cells.begin(); it != cells.end(); ++it) {
        file << (*it)[0] << ' ';
    }

    file << "\n\t\t\t\t</DataArray>\n";

    file << "\t\t\t</Cells>\n";

    file << "\t\t</Piece>\n";
    file << "\t</UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}

#endif