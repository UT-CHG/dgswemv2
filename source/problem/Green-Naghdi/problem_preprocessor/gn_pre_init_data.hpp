#ifndef GN_PRE_INIT_DATA_HPP
#define GN_PRE_INIT_DATA_HPP

namespace GN {
template <typename MeshType>
void initialize_data_serial(MeshType& mesh) {
#ifdef D_INTERPOLATION
    mesh.CallForEachElement([](auto& elt) {
        auto& derivative         = elt.data.derivative;
        derivative.area          = elt.GetShape().GetArea();
        derivative.baryctr_coord = elt.GetShape().GetBarycentricCoordinates();
        derivative.midpts_coord  = elt.GetShape().GetMidpointCoordinates();
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in                                    = intface.data_in.derivative;
        auto& derivative_ex                                    = intface.data_ex.derivative;
        derivative_in.baryctr_coord_neigh[intface.bound_id_in] = derivative_ex.baryctr_coord;
        derivative_ex.baryctr_coord_neigh[intface.bound_id_ex] = derivative_in.baryctr_coord;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative                               = bound.data.derivative;
        derivative.baryctr_coord_neigh[bound.bound_id] = derivative.midpts_coord[bound.bound_id];
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;

        StatMatrix<double, 2, 2> A;
        StatVector<double, 2> alpha;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            for (uint perm = 0; perm < 2; ++perm) {
                const uint element_1 = bound;
                const uint element_2 = (bound + perm + 1) % elt.data.get_nbound();
                if (!is_distributed(elt.GetBoundaryType()[element_1]) &&
                    !is_distributed(elt.GetBoundaryType()[element_2])) {
                    column(A, 0) = derivative.baryctr_coord_neigh[element_1] - derivative.baryctr_coord;
                    column(A, 1) = derivative.baryctr_coord_neigh[element_2] - derivative.baryctr_coord;
                    alpha        = derivative.midpts_coord[bound] - derivative.baryctr_coord;
                    solve_sle(A, alpha);
                    for (uint i = 0; i < 2; ++i)
                        if (Utilities::almost_equal(alpha[i], 0.0))
                            alpha[i] = 0.0;
                    if (alpha[0] >= 0.0 && alpha[1] >= 0.0) {
                        derivative.a[bound]              = alpha;
                        derivative.a_elem[2 * bound]     = element_1;
                        derivative.a_elem[2 * bound + 1] = element_2;
                        continue;
                    }
                }
            }
        }
#ifdef D_LEASTSQUARES
        HybMatrix<double, 2> D_transpose(2, elt.data.get_nbound());
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(D_transpose, bound) = derivative.baryctr_coord_neigh[bound] - derivative.baryctr_coord;
        }
        derivative.P = inverse(D_transpose * transpose(D_transpose)) * D_transpose;
#endif
    });
#endif
}

template <typename MeshType>
void initialize_data_parallel_pre_send(MeshType& mesh, uint comm_type) {
    initialize_data_serial(mesh);

#ifdef D_INTERPOLATION
    mesh.CallForEachDistributedBoundary([comm_type](auto& dbound) {
        auto& derivative = dbound.data.derivative;

        std::vector<double> message(GN::n_dimensions);
        for (uint dim = 0; dim < GN::n_dimensions; ++dim) {
            message[dim] = derivative.baryctr_coord[dim];
        }
        dbound.boundary_condition.exchanger.SetToSendBuffer(comm_type, message);
    });
#endif
}

template <typename MeshType>
void initialize_data_parallel_post_receive(MeshType& mesh, uint comm_type) {
#ifdef D_INTERPOLATION
    mesh.CallForEachDistributedBoundary([comm_type](auto& dbound) {
        auto& derivative = dbound.data.derivative;

        std::vector<double> message(GN::n_dimensions);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(comm_type, message);
        for (uint dim = 0; dim < GN::n_dimensions; ++dim) {
            derivative.baryctr_coord_neigh[dbound.bound_id][dim] = message[dim];
        }
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;

        StatMatrix<double, 2, 2> A;
        StatVector<double, 2> alpha;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            for (uint perm = 0; perm < 2; ++perm) {
                const uint element_1 = bound;
                const uint element_2 = (bound + perm + 1) % elt.data.get_nbound();
                if (is_distributed(elt.GetBoundaryType()[element_1]) ||
                    is_distributed(elt.GetBoundaryType()[element_2])) {
                    column(A, 0) = derivative.baryctr_coord_neigh[element_1] - derivative.baryctr_coord;
                    column(A, 1) = derivative.baryctr_coord_neigh[element_2] - derivative.baryctr_coord;
                    alpha        = derivative.midpts_coord[bound] - derivative.baryctr_coord;
                    solve_sle(A, alpha);
                    for (uint i = 0; i < 2; ++i)
                        if (Utilities::almost_equal(alpha[i], 0.0))
                            alpha[i] = 0.0;
                    if (alpha[0] >= 0.0 && alpha[1] >= 0.0) {
                        derivative.a[bound]              = alpha;
                        derivative.a_elem[2 * bound]     = element_1;
                        derivative.a_elem[2 * bound + 1] = element_2;
                        continue;
                    }
                }
            }
        }
#ifdef D_LEASTSQUARES
        HybMatrix<double, 2> D_transpose(2, elt.data.get_nbound());
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(D_transpose, bound) = derivative.baryctr_coord_neigh[bound] - derivative.baryctr_coord;
        }
        derivative.P = inverse(D_transpose * transpose(D_transpose)) * D_transpose;
#endif
    });
#endif
}
}

#endif
