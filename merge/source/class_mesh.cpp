#include "class_mesh.h"

MESH::MESH(int p, int p_geom) {
	this->p = p;
	this->p_geom = p_geom;
}

MESH::~MESH() {
    for (auto it = this->interfaces.begin(); it != this->interfaces.end(); it++) {
        for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
            delete *itt;
        }
        it->second.clear();
    }
}

void MESH::RectangularDomainTest
(double L, double W, int m, int n, int element_type) 
{
	double dx = L / m;
	double dy = W / n;

	switch (element_type) {
	case TRIANGLE: {
		unsigned int ID;

		Array2D<double> nodal_coordinates(2);
		nodal_coordinates[0].resize(3);
		nodal_coordinates[1].resize(3);

		std::vector<unsigned char> boundaries(3);
		std::vector<unsigned int> neighbors(3);

		//SIMPLE PATTERN
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

		//// CHECKERS PATTERN
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

				nodal_coordinates[0][0] = j*dx;
				nodal_coordinates[0][1] = nodal_coordinates[0][0];
				nodal_coordinates[0][2] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[1][0] = (i + 1)*dy;
				nodal_coordinates[1][1] = nodal_coordinates[1][0] - dy;
				nodal_coordinates[1][2] = nodal_coordinates[1][0];

				this->elements[ID] = new ELEMENT(triangle, ID, neighbors, boundaries, nodal_coordinates);

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
					boundaries[1] = OCEAN;
				}

				nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[1][0] = nodal_coordinates[1][0] - dy;

				this->elements[ID] = new ELEMENT(triangle, ID, neighbors, boundaries, nodal_coordinates);
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

				nodal_coordinates[0][0] = j*dx;
				nodal_coordinates[0][1] = nodal_coordinates[0][0] + dx;
				nodal_coordinates[0][2] = nodal_coordinates[0][0];

				nodal_coordinates[1][0] = i*dy;
				nodal_coordinates[1][1] = nodal_coordinates[1][0];
				nodal_coordinates[1][2] = nodal_coordinates[1][0] + dy;

				this->elements[ID] = new ELEMENT(triangle, ID, neighbors, boundaries, nodal_coordinates);

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
					boundaries[2] = OCEAN;
				}

				nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[1][0] = nodal_coordinates[1][0] + dy;

				this->elements[ID] = new ELEMENT(triangle, ID, neighbors, boundaries, nodal_coordinates);
			}
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
    std::vector<Point<3>> points;
    Array2D<unsigned int> cells;

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
}