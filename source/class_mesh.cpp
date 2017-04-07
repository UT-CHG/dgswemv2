#include "class_mesh.h"

MESH::MESH(int p, int p_geom) {
    this->p = p;
    this->p_geom = p_geom;

    this->InitializeElements();
    this->InitializeInterfaces();
    this->InitializeVTK();
}

MESH::~MESH() {
    for (auto it = this->interfaces.begin(); it != this->interfaces.end(); it++) {
        for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
            delete *itt;
        }
        it->second.clear();
    }
    this->interfaces.clear();

    for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
        delete it->second;
    }
    this->elements.clear();

    for (auto it = this->bases_2D.begin(); it != this->bases_2D.end(); it++) {
        delete it->second;
    }
    this->bases_2D.clear();

    for (auto it = this->geometric_bases_2D.begin(); it != this->geometric_bases_2D.end(); it++) {
        delete it->second;
    }
    this->geometric_bases_2D.clear();

    for (auto it = this->line_rules.begin(); it != this->line_rules.end(); it++) {
        delete it->second;
    }
    this->line_rules.clear();

    for (auto it = this->area_rules.begin(); it != this->area_rules.end(); it++) {
        delete it->second;
    }
    this->area_rules.clear();
}

void MESH::InitializeElements() {
    this->line_rules[GAUSS_LEGENDRE] = new INTEGRATION_1D(GAUSS_LEGENDRE, 2 * this->p);
    this->area_rules[DUNAVANT] = new INTEGRATION_2D(DUNAVANT, 2 * this->p);

    INTEGRATION_1D* line_rule;

    if (this->line_rules.find(GAUSS_LEGENDRE)->first == GAUSS_LEGENDRE) {
        line_rule = this->line_rules.find(GAUSS_LEGENDRE)->second;
    }
    else {
        this->line_rules[GAUSS_LEGENDRE] = new INTEGRATION_1D(GAUSS_LEGENDRE, 2 * p);
        line_rule = this->line_rules.find(GAUSS_LEGENDRE)->second;
    }

    INTEGRATION_2D* area_rule;

    if (this->area_rules.find(DUNAVANT)->first == DUNAVANT) {
        area_rule = this->area_rules.find(DUNAVANT)->second;
    }
    else {
        this->area_rules[DUNAVANT] = new INTEGRATION_2D(DUNAVANT, 2 * p);
        area_rule = this->area_rules.find(DUNAVANT)->second;
    }

    this->bases_2D[DUBINER] = new BASIS_2D(DUBINER, p, line_rule, area_rule);

    BASIS_2D* basis;

    if (this->bases_2D.find(DUBINER)->first == DUBINER) {
        basis = this->bases_2D.find(DUBINER)->second;
    }
    else {
        this->bases_2D[DUBINER] = new BASIS_2D(DUBINER, p, line_rule, area_rule);
        basis = this->bases_2D.find(DUBINER)->second;
    }

    unsigned int ID = 0;

    double x[3] = { -1,  1, -1 };
    double y[3] = { -1, -1,  1 };

    unsigned int neighbors[3] = { 1, DEFAULT_ID, DEFAULT_ID };
    unsigned char boundaries[3] = { INTERNAL, LAND, LAND };

    this->elements[ID] = new ELEMENT_2D(ID, neighbors, boundaries, x, y, basis);

    ID = 1;
    x[0] = 1;
    y[0] = 1;
    neighbors[0] = 0;

    this->elements[ID] = new ELEMENT_2D(ID, neighbors, boundaries, x, y, basis);
}

void MESH::InitializeInterfaces() {
    for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
        std::map<unsigned int, INTERFACE*> internal_interfaces(it->second->CreateInterfaces());

        for (auto itt = internal_interfaces.begin(); itt != internal_interfaces.end(); itt++) {
            this->elements.find(itt->first)->second->AppendInterface(it->first, itt->second);
        }

        std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces(it->second->GetOwnInterfaces());

        for (auto itt = own_interfaces.begin(); itt != own_interfaces.end(); itt++) {
            if (this->interfaces.find(itt->first) != this->interfaces.end()) {
                this->interfaces.find(itt->first)->second.push_back(itt->second);
            }
            else {
                this->interfaces[itt->first] = std::vector<INTERFACE*>{ itt->second };
            }
        }
    }
}

void MESH::InitializeVTK() {
    std::vector<double*> points;
    std::vector<unsigned int*> cells;

    for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
        it->second->InitializeVTK(points, cells);
    }

    std::string file_name = "geometry.vtk";
    std::ofstream file(file_name);

    file << "# vtk DataFile Version 3.0\n";
    file << "OUTPUT DATA\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << points.size() << " float\n";

    for (auto it = points.begin(); it != points.end(); it++) {
        file << (*it)[0] << '\t' << (*it)[1] << '\t' << (*it)[2] << '\n';
    }

    int n_cell_entries = 0;
    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
        case 5: n_cell_entries += 4; break;
        default:
            printf("\n");
            printf("MESH VTK - Fatal error!\n");
            printf("Undefined cell type = %d\n", (*it)[0]);
            exit(1);
        }
    }

    file << "CELLS " << cells.size() << ' ' << n_cell_entries << '\n';
    
    int n_nodes;

    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
        case 5: file << 3 << '\t'; n_nodes = 3; break;
        default:
            printf("\n");
            printf("MESH VTK - Fatal error!\n");
            printf("Undefined cell type = %d\n", (*it)[0]);
            exit(1);
        }

        for (int i = 1; i <= n_nodes; i++) {
            file << (*it)[i] << '\t';
        }
        file << '\n';
    }

    file << "CELL_TYPES " << cells.size() << '\n';

    for (auto it = cells.begin(); it != cells.end(); it++) {
        file << (*it)[0] << '\n';
    }
    
    for (auto it = cells.begin(); it != cells.end(); it++) delete[] *it;
    for (auto it = points.begin(); it != points.end(); it++) delete[] *it;

    points.clear();
    cells.clear();
}