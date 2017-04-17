#include "class_mesh.h"

MESH::MESH(int p, int p_geom) {
	this->p = p;
	this->p_geom = p_geom;

	INTEGRATION_1D* line_rule = new INTEGRATION_1D(GAUSS_LEGENDRE, 2 * this->p);
	this->line_rules[TRIANGLE] = line_rule;

	INTEGRATION_2D* area_rule = new INTEGRATION_2D(DUNAVANT, 2 * this->p);
	this->area_rules[TRIANGLE] = area_rule;
	
	this->bases_2D[TRIANGLE] = new BASIS_2D(DUBINER, p, line_rule, area_rule);
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

void MESH::RectangularDomainTest(double L, double W, int m, int n, int element_type) {
	double dx = L / m;
	double dy = W / n;

	unsigned int ID;
	
	switch (element_type) {
	case TRIANGLE:
		unsigned char boundaries[3];
		unsigned int neighbors[3];
		double x[3];
		double y[3];
		BASIS_2D* basis;

		if (!(this->bases_2D.empty()) && (this->bases_2D.find(TRIANGLE)->first == TRIANGLE)) {
			basis = this->bases_2D.find(TRIANGLE)->second;
		}
		else {
			printf("\n");
			printf("MESH RectangularDomainTest - Fatal error!\n");
			printf("Triangular basis not defined");
			exit(1);
		}

		// SIMPLE PATTERN
		//for (int i = 0; i < n; i++) {
		//	for (int j = 0; j < m; j++) {
		//		ID = 2 * j + 2 * m * i;

		//		neighbors[0] = ID + 1;
		//		neighbors[1] = ID + 1 + 2 * m;
		//		neighbors[2] = ID - 1;

		//		boundaries[0] = INTERNAL;
		//		boundaries[1] = INTERNAL;
		//		boundaries[2] = INTERNAL;

		//		if (i == n - 1) {
		//			neighbors[1] = DEFAULT_ID;
		//			boundaries[1] = LAND;
		//		}
		//		if (j == 0) {
		//			neighbors[2] = DEFAULT_ID;
		//			boundaries[2] = LAND;
		//		}

		//		x[0] = j*dx;
		//		x[1] = x[0];
		//		x[2] = x[0] + dx;

		//		y[0] = (i + 1)*dy;
		//		y[1] = y[0] - dy;
		//		y[2] = y[0];

		//		this->elements[ID] = new ELEMENT_2D(TRIANGLE, ID, neighbors, boundaries, x, y, basis);

		//		ID = ID + 1;

		//		neighbors[0] = ID - 1;
		//		neighbors[1] = ID + 1;
		//		neighbors[2] = ID - 1 - 2 * m;

		//		boundaries[0] = INTERNAL;
		//		boundaries[1] = INTERNAL;
		//		boundaries[2] = INTERNAL;
		//		
		//		if (i == 0) {
		//			neighbors[2] = DEFAULT_ID;
		//			boundaries[2] = LAND;
		//		}
		//		if (j == m - 1) {
		//			neighbors[1] = DEFAULT_ID;
		//			boundaries[1] = OCEAN;
		//		}

		//		x[0] = x[0] + dx;
		//		y[0] = y[0] - dy;

		//		this->elements[ID] = new ELEMENT_2D(TRIANGLE, ID, neighbors, boundaries, x, y, basis);
		//	}
		//}

		// CHECKERS PATTERN
		for (int i = 0; i < n; i++) {
			for (int j = i % 2; j < m; j += 2) {
				ID = 2 * j + 2 * m * i;

				neighbors[0] = ID + 1;
				neighbors[1] = ID + 2 * m;
				neighbors[2] = ID - 1;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == n - 1) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}
				if (j == 0) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}

				x[0] = j*dx;
				x[1] = x[0];
				x[2] = x[0] + dx;

				y[0] = (i + 1)*dy;
				y[1] = y[0] - dy;
				y[2] = y[0];

				this->elements[ID] = new ELEMENT_2D(TRIANGLE, ID, neighbors, boundaries, x, y, basis);

				ID = ID + 1;

				neighbors[0] = ID - 1;
				neighbors[1] = ID + 1;
				neighbors[2] = ID - 2 * m;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == 0) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}
				if (j == m - 1) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}

				x[0] = x[0] + dx;
				y[0] = y[0] - dy;

				this->elements[ID] = new ELEMENT_2D(TRIANGLE, ID, neighbors, boundaries, x, y, basis);
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = (i + 1) % 2; j < m; j += 2) {
				ID = 2 * j + 2 * m * i;

				neighbors[0] = ID + 1;
				neighbors[1] = ID - 1;
				neighbors[2] = ID - 2 * m;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == 0) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}
				if (j == 0) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}

				x[0] = j*dx;
				x[1] = x[0] + dx;
				x[2] = x[0];

				y[0] = i*dy;
				y[1] = y[0];
				y[2] = y[0] + dy;

				this->elements[ID] = new ELEMENT_2D(TRIANGLE, ID, neighbors, boundaries, x, y, basis);

				ID = ID + 1;

				neighbors[0] = ID - 1;
				neighbors[1] = ID + 2 * m;
				neighbors[2] = ID + 1;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == n - 1) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}
				if (j == m - 1) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}

				x[0] = x[0] + dx;
				y[0] = y[0] + dy;

				this->elements[ID] = new ELEMENT_2D(TRIANGLE, ID, neighbors, boundaries, x, y, basis);
			}
		}

		break;
	default:
		printf("\n");
		printf("MESH RectangularDomainTest - Fatal error!\n");
		printf("Undefined element type = %d\n", element_type);
		exit(1);
	}

	this->InitializeInterfaces();
	
	this->InitializeVTK();
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
        case TRIANGLE: n_cell_entries += 4; break;
        default:
            printf("\n");
            printf("MESH InitializeVTK - Fatal error!\n");
            printf("Undefined cell type = %d\n", (*it)[0]);
            exit(1);
        }
    }

    file << "CELLS " << cells.size() << ' ' << n_cell_entries << '\n';
    
    int n_nodes;

    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
        case TRIANGLE: file << 3 << '\t'; n_nodes = 3; break;
        default:
            printf("\n");
            printf("MESH InitializeVTK - Fatal error!\n");
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