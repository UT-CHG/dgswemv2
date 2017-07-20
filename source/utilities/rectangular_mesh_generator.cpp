#include <vector>
#include <array>

template<int dim>
using Point = std::array<double, dim>;

int main() {
	double L = 0;
	double W = 0;
	int m;
	int n;

	double dx = L / m;
	double dy = W / n;

	unsigned int ID;

	std::vector<Point<2>> nodal_coordinates(3);

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
	for (int i = 0; i < n; i++) {
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
	}
}