#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iomanip>

template<int dim>
using Point = std::array<double, dim>;

using uint = unsigned int;
using uchar = unsigned char;

struct node {
	uint ID;
	std::array<double, 3> coord;
};

struct element {
	uint ID;
	uint type;
	std::vector<uint> nodes;
};

struct boundary {
	uchar type;
	std::vector<uint> nodes;
};

int main() {
	double L = 50000;
	double W = 8000;
	uint m = 25;
	uint n = 4;

	std::vector<uchar> boundary_type{ 0,0,0,1 }; // 0 - land, 1 - tidal

	double dx = L / m;
	double dy = W / n;

	std::vector<node> nodes((m + 1)*(n + 1));

	for (uint i = 0; i <= m; i++) {
		for (uint j = 0; j <= n; j++) {
			nodes[j*(m + 1) + i].ID = j*(m + 1) + i;

			nodes[j*(m + 1) + i].coord[0] = dx*i;
			nodes[j*(m + 1) + i].coord[1] = dy*j;
			nodes[j*(m + 1) + i].coord[2] = 2.0;
		}
	}

	std::vector<element> elements(2 * m*n);

	//mesh with triangular elements checker pattern
	for (uint j = 0; j < n; j++) {
		for (uint i = j % 2; i < m; i += 2) {
			elements[2 * i + 2 * m * j].ID = 2 * i + 2 * m * j;
			elements[2 * i + 2 * m * j].type = 3;
			elements[2 * i + 2 * m * j].nodes.resize(3);
			elements[2 * i + 2 * m * j].nodes =
				std::vector<uint>{ i + (m + 1)*j, i + m + 1 + (m + 1)*j, i + m + 2 + (m + 1)*j };

			elements[2 * i + 2 * m * j + 1].ID = 2 * i + 2 * m * j + 1;
			elements[2 * i + 2 * m * j + 1].type = 3;
			elements[2 * i + 2 * m * j + 1].nodes =
				std::vector<uint>{ i + (m + 1)*j, i + 1 + (m + 1)*j, i + m + 2 + (m + 1)*j };
		}
	}

	for (uint j = 0; j < n; j++) {
		for (uint i = (j + 1) % 2; i < m; i += 2) {
			elements[2 * i + 2 * m * j].ID = 2 * i + 2 * m * j;
			elements[2 * i + 2 * m * j].type = 3;
			elements[2 * i + 2 * m * j].nodes =
				std::vector<uint>{ i + (m + 1)*j, i + m + 1 + (m + 1)*j, i + 1 + (m + 1)*j };

			elements[2 * i + 2 * m * j + 1].ID = 2 * i + 2 * m * j + 1;
			elements[2 * i + 2 * m * j + 1].type = 3;
			elements[2 * i + 2 * m * j + 1].nodes =
				std::vector<uint>{ i + m + 1 + (m + 1)*j, i + m + 2 + (m + 1)*j, i + 1 + (m + 1)*j };
		}
	}

	std::vector<boundary> boundaries(4);

	boundaries[0].type = boundary_type[0];
	boundaries[1].type = boundary_type[1];
	boundaries[2].type = boundary_type[2];
	boundaries[3].type = boundary_type[3];

	boundaries[0].nodes.resize(m + 1);
	boundaries[1].nodes.resize(n + 1);
	boundaries[2].nodes.resize(m + 1);
	boundaries[3].nodes.resize(n + 1);

	for (uint i = 0; i <= m; i++) {
		boundaries[0].nodes[i] = i;
		boundaries[2].nodes[i] = i + (m + 1)*n;
	}

	for (uint i = 0; i <= n; i++) {
		boundaries[1].nodes[i] = m + i*(m + 1);
		boundaries[3].nodes[i] = i*(m + 1);
	}

	std::string file_name = "rectangular_mesh.14";
	std::ofstream file(file_name);

	file << std::fixed << std::setprecision(6);
	file << "ADCIRC Model\n";
	file << 2 * m*n << "    " << (m + 1)*(n + 1) << '\n';

	for (uint node = 0; node < nodes.size(); node++) {
		file << nodes[node].ID << ' ';
		file << nodes[node].coord[0] << ' ';
		file << nodes[node].coord[1] << ' ';
		file << nodes[node].coord[2] << ' ';
		file << '\n';
	}

	for (uint element = 0; element < elements.size(); element++) {
		file << elements[element].ID << ' ';
		file << elements[element].type << ' ';
		file << elements[element].nodes[0] << ' ';
		file << elements[element].nodes[1] << ' ';
		file << elements[element].nodes[2] << ' ';
		file << '\n';
	}

	uint n_land = 0;
	uint n_tidal = 0;

	uint n_land_node = 0;
	uint n_tidal_node = 0;

	for (uint n_bound = 0; n_bound < 4; n_bound++) {
		if (boundaries[n_bound].type == 0) {
			n_land++;
			n_land_node += boundaries[n_bound].nodes.size();
		}
		else if (boundaries[n_bound].type == 1) {
			n_tidal++;
			n_tidal_node += boundaries[n_bound].nodes.size();
		}
	}

	file << n_tidal << " = Number of open boundaries\n";
	file << n_tidal_node << " = Total number of open boundary nodes\n";

	uint i = 1;
	for (uint n_bound = 0; n_bound < 4; n_bound++) {
		if (boundaries[n_bound].type == 1) {
			file << boundaries[n_bound].nodes.size() << " = Number of nodes for open boundary " << i << '\n';

			std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(),
				[&file](uint val) {file << val << '\n'; }
			);

			i++;
		}
	}

	file << n_land << " = Number of land boundaries\n";
	file << n_land_node << " = Total number of land boundary nodes\n";

	i = 1;
	for (uint n_bound = 0; n_bound < 4; n_bound++) {
		if (boundaries[n_bound].type == 0) {
			file << boundaries[n_bound].nodes.size() << " 0 = Number of nodes for land boundary " << i << '\n';

			std::for_each(boundaries[n_bound].nodes.begin(), boundaries[n_bound].nodes.end(),
				[&file](uint val) { file << val << '\n'; }
			);

			i++;
		}
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