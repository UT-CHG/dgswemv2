#ifndef SWE_PRE_PREPROCESS_MESH_HPP
#define SWE_PRE_PREPROCESS_MESH_HPP

namespace SWE {
void Problem::preprocess_mesh_data(InputParameters<ProblemInputType>& input) {
    // need to do projection to cartesian coordinate system by projecting on a cylinder
    if (input.mesh_input.mesh_coordinate_sys == CoordinateSystem::spherical) {
        if (input.problem_input.spherical_projection.type == SphericalProjectionType::None) {
            std::cerr
                << "Warning: Spherical projection input not provided! Using default values to perfrom projection!\n";
            std::cerr << "Consider using spherical projection to account for geometric distortion!\n";
        }

        double R           = input.problem_input.spherical_projection.R;
        double longitude_o = input.problem_input.spherical_projection.longitude_o;
        double lattitude_o = input.problem_input.spherical_projection.latitude_o;

        double R_o = R * cos(lattitude_o * PI / 180.0);

        for (auto& node : input.mesh_input.mesh_data.nodes) {
            node.second.coordinates[GlobalCoord::x] =
                R_o * (node.second.coordinates[GlobalCoord::x] - longitude_o) * PI / 180.0;
            node.second.coordinates[GlobalCoord::y] = R * node.second.coordinates[GlobalCoord::y] * PI / 180.0;
        }
    }
}
}

#endif