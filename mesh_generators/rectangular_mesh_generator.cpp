#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iomanip>

template<int dim>
using Point = std::array<double, dim>;

using uint = unsigned int;

struct node {
	uint ID;
	std::array<double, 3> coord;
};

struct element {
	uint ID;
	uint type;
	std::vector<uint> nodes;
};

int main() {
	double L = 100;
	double W = 100;
	uint m = 6;
	uint n = 8;

	double dx = L / m;
	double dy = W / n;

	std::vector<node> nodes((m+1)*(n+1));

	for (uint i=0;i<=m;i++){
		for(uint j=0;j<=n;j++){
			nodes[j*(m+1)+i].ID = j*(m+1)+i;

			nodes[j*(m+1)+i].coord[0] = dx*i;
			nodes[j*(m+1)+i].coord[1] = dy*j;
			nodes[j*(m+1)+i].coord[2] = 2.0; 
		}
	}

	std::vector<element> elements(2*m*n);

	//mesh with triangular elements checker pattern
	for (uint j = 0; j < n; j++) {
		for (uint i = j % 2; i < m; i += 2) {
			elements[2 * i + 2 * m * j].ID = 2 * i + 2 * m * j;
			elements[2 * i + 2 * m * j].type = 3;
			elements[2 * i + 2 * m * j].nodes.resize(3);
			elements[2 * i + 2 * m * j].nodes = std::vector<uint>{i+(m+1)*j, i+m+1+(m+1)*j, i+m+2+(m+1)*j};
			
			elements[2 * i + 2 * m * j + 1].ID = 2 * i + 2 * m * j + 1;
			elements[2 * i + 2 * m * j + 1].type = 3;
			elements[2 * i + 2 * m * j + 1].nodes = std::vector<uint>{i+(m+1)*j, i+1+(m+1)*j, i+m+2+(m+1)*j};
		}
	}

	for (uint j = 0; j < n; j++) {
		for (uint i = (j + 1) % 2; i < m; i += 2) {
			elements[2 * i + 2 * m * j].ID = 2 * i + 2 * m * j;
			elements[2 * i + 2 * m * j].type = 3;
			elements[2 * i + 2 * m * j].nodes = std::vector<uint>{i+(m+1)*j, i+m+1+(m+1)*j, i+1+(m+1)*j};

			elements[2 * i + 2 * m * j + 1].ID = 2 * i + 2 * m * j + 1;
			elements[2 * i + 2 * m * j + 1].type = 3;
			elements[2 * i + 2 * m * j + 1].nodes = std::vector<uint>{i+m+1+(m+1)*j, i+m+2+(m+1)*j, i+1+(m+1)*j};
		}
	}

	std::string file_name = "rectangular_mesh.14";
	std::ofstream file(file_name);

	file << std::fixed << std::setprecision(6);
	
	for(uint node=0; node<nodes.size();node++){
		file << nodes[node].ID << "\t\t";
		file << nodes[node].coord[0] << "\t\t";
		file << nodes[node].coord[1] << "\t\t";
		file << nodes[node].coord[2] << '\n';
	}

	for(uint element=0; element<elements.size();element++){
		file << elements[element].ID << "\t\t";
		file << elements[element].type << "\t\t";
		file << elements[element].nodes[0] << "\t\t";
		file << elements[element].nodes[1] << "\t\t";
		file << elements[element].nodes[2] << '\n';
	}

	//SIMPLE PATTERN
	//for (int i = 0; i < n; i++) {
	//	for (int j = 0; j < m; j++) {
	//		ID = 2 * j + 2 * m * i;

	//		x[0] = j*dx;
	//		x[1] = x[0];
	//		x[2] = x[0] + dx;

	//		y[0] = (i + 1)*dy;
	//		y[1] = y[0] - dy;
	//		y[2] = y[0];

	//		ID = ID + 1;

	//		x[0] = x[0] + dx;
	//		y[0] = y[0] - dy;
	//	}
	//}

	//// CHECKERS PATTERN
/*	for (int i = 0; i < n; i++) {
		for (int j = i % 2; j < m; j += 2) {
			ID = 2 * j + 2 * m * i;

			nodal_coordinates[0][0] = j*dx;
			nodal_coordinates[1][0] = nodal_coordinates[0][0];
			nodal_coordinates[2][0] = nodal_coordinates[0][0] + dx;

			nodal_coordinates[0][1] = (i + 1)*dy;
			nodal_coordinates[1][1] = nodal_coordinates[0][1] - dy;
			nodal_coordinates[2][1] = nodal_coordinates[0][1];

			ID = ID + 1;

			nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

			nodal_coordinates[0][1] = nodal_coordinates[0][1] - dy;
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = (i + 1) % 2; j < m; j += 2) {
			ID = 2 * j + 2 * m * i;

			nodal_coordinates[0][0] = j*dx;
			nodal_coordinates[1][0] = nodal_coordinates[0][0] + dx;
			nodal_coordinates[2][0] = nodal_coordinates[0][0];

			nodal_coordinates[0][1] = i*dy;
			nodal_coordinates[1][1] = nodal_coordinates[0][1];
			nodal_coordinates[2][1] = nodal_coordinates[0][1] + dy;

			ID = ID + 1;

			nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

			nodal_coordinates[0][1] = nodal_coordinates[0][1] + dy;
		}
	}*/
}