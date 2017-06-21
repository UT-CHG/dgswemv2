#include "class_mesh.h"

MESH::~MESH() {
	delete shape; 
	delete triangle;

	for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
		delete it->second;
	}
	this->elements.clear();
}

void MESH::RectangularDomainTest
(double L, double W, int m, int n, int element_type) 
{
	double dx = L / m;
	double dy = W / n;

	switch (element_type) {
	case TRIANGLE: {
		unsigned int ID;

		std::vector<Point<2>> nodal_coordinates(3);

		std::vector<unsigned char> boundaries(3);
		std::vector<unsigned int> neighbors(3);
		
		triangle = new Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>(this->p); 

		shape = new Shape::StraightTriangle();

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
				nodal_coordinates[1][0] = nodal_coordinates[0][0];
				nodal_coordinates[2][0] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[0][1] = (i + 1)*dy;
				nodal_coordinates[1][1] = nodal_coordinates[0][1] - dy;
				nodal_coordinates[2][1] = nodal_coordinates[0][1];

				this->elements[ID] = new Element<>(*triangle, *shape, ID, neighbors, boundaries, nodal_coordinates);

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

				nodal_coordinates[0][1] = nodal_coordinates[0][1] - dy;

				this->elements[ID] = new Element<>(*triangle, *shape, ID, neighbors, boundaries, nodal_coordinates);
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
				nodal_coordinates[1][0] = nodal_coordinates[0][0] + dx;
				nodal_coordinates[2][0] = nodal_coordinates[0][0];

				nodal_coordinates[0][1] = i*dy;
				nodal_coordinates[1][1] = nodal_coordinates[0][1];
				nodal_coordinates[2][1] = nodal_coordinates[0][1] + dy;

				this->elements[ID] = new Element<>(*triangle, *shape, ID, neighbors, boundaries, nodal_coordinates);

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

				nodal_coordinates[0][1] = nodal_coordinates[0][1] + dy;

				this->elements[ID] = new Element<>(*triangle, *shape, ID, neighbors, boundaries, nodal_coordinates);
			}
		}

		//delete triangle;
	}
	break;
	default:
		printf("\n");
		printf("MESH RectangularDomainTest - Fatal error!\n");
		printf("Undefined element type = %d\n", element_type);
		exit(1);
	}

	this->InitializeBoundariesInterfaces();

	//printf("%d\n", this->boundaries[OCEAN].size());
	//printf("%d\n", this->interfaces_.size());

	this->InitializeVTK();
}

void MESH::InitializeBoundariesInterfaces() {
	std::map<unsigned char, std::vector<RawBoundary<>*>> pre_boundaries;
	std::map<unsigned int, std::map<unsigned int, RawBoundary<>*>> pre_interfaces;

	for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
		std::vector<RawBoundary<>*> raw_boundaries = it->second->CreateBoundaries();

		for (auto itt = raw_boundaries.begin(); itt != raw_boundaries.end(); itt++) {
			if ((*itt)->type != INTERNAL) {
				if (pre_boundaries.find((*itt)->type) != pre_boundaries.end()) {
					pre_boundaries.find((*itt)->type)->second.push_back(*itt);
				}
				else {
					pre_boundaries[(*itt)->type] = std::vector<RawBoundary<>*>{ *itt };
				}
			}
			else if ((*itt)->type == INTERNAL) {
				pre_interfaces[it->first][(*itt)->neighbor_ID] = *itt;
			}
		}
	}

	for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
		this->boundaries[it->first] = std::vector<Boundary<>*>();

		for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
			this->boundaries[it->first].push_back(new Boundary<>(*(*itt)));
			delete *itt;
		}
	}

	RawBoundary<>* raw_in;
	RawBoundary<>* raw_ex;

	for (auto it = pre_interfaces.begin(); it != pre_interfaces.end(); it++) {
		for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
			raw_in = itt->second;
			raw_ex = pre_interfaces[itt->second->neighbor_ID][it->first];

			this->interfaces.push_back(new Interface<>(*raw_in, *raw_ex));

			pre_interfaces[itt->second->neighbor_ID].erase(it->first);
			it->second.erase(itt);

			delete raw_in;
			delete raw_ex;
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

void MESH::Solve(){
	Stepper stepper(2, 2, 10.);

	int stages = stepper.get_num_stages();

	for (auto it = elements.begin(); it != elements.end(); it++){
		it->second->data.state = std::vector<SWE::State>(stages, 
			*it->second->data.state.begin());
	}

	for (auto it = elements.begin(); it != elements.end(); it++){
		SWE::volume_kernel(stepper, it->second);
	}

	for (auto it = boundaries.begin(); it != boundaries.end(); it++){
		for (auto itt = it->second.begin(); itt != it->second.end(); itt++){
			SWE::boundary_kernel(stepper, *itt);
		}
	}

	for (auto it = interfaces.begin(); it != interfaces.end(); it++){
		SWE::interface_kernel(stepper, *it);
	}

	for (auto it = elements.begin(); it != elements.end(); it++){
		SWE::update_kernel(stepper, it->second);
	}

	for (auto it = elements.begin(); it != elements.end(); it++){
		//SWE::swap_states(stepper, it->second);
	}
}