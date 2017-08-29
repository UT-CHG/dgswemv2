#ifndef SWE_KERNELS_HPP
#define SWE_KERNELS_HPP

#include "swe_LLF_flux.hpp"
#include "swe_true_src_functions.hpp"

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_boundaries_kernel(mesh_type& mesh, std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries) {
    uint n_bound_old_land = 0;
    uint n_bound_old_tidal = 0;

    std::ofstream log_file("output/" + mesh.GetMeshName() + "_log", std::ofstream::app);

    using BoundaryTypes = Geometry::BoundaryTypeTuple<SWE::Data, SWE::Land, SWE::Tidal>;

    for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
        switch (it->first) {
            case SWE::land:
                using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;

                n_bound_old_land = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeLand>(*itt);
                }

                log_file << "Number of land boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_land << std::endl;

                break;
            case SWE::tidal:
                using BoundaryTypeTidal = typename std::tuple_element<1, BoundaryTypes>::type;

                n_bound_old_tidal = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeTidal>(*itt);
                }

                log_file << "Number of tidal boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_tidal
                         << std::endl;

                break;
        }
    }

    log_file << std::endl;
}

template <typename RawBoundaryType, typename Communicator>
void Problem::create_distributed_interfaces_kernel(mesh_type& mesh,
                                                   Communicator& communicator,
                                                   std::map<uint, std::map<uint, RawBoundaryType>>& pre_boundaries) {

    using DistributedInterfaceType =
        std::tuple_element<0, Geometry::DistributedInterface<SWE::Data, SWE::Distributed>>::type;

    using Integration = DistributedInterfaceType::BoundaryIntegrationType;

    Integration integ;

    communicator.ResizeBuffers(integ, 3);

    for (uint n = 0; n < communicator.GetNumNeighbors(); ++n) {
        uint curr_indx{0};

        std::vector<double>& send_buff_ref = communicator.GetSendBufferReference(n);
        std::vector<double>& recv_buff_ref = communicator.GetReceiveBufferReference(n);

        uint buff_size = send_buff_ref.size();

        for (uint e = 0; e < communicator.GetNumEdges(n); ++e) {
            uint elt;
            uint fid;
            uint p;
            {
                std::tuple<uint, uint, uint> tup = communicator.GetElt_FaceID_PolynomialOrder(n, e);
                elt = std::get<0>(tup);
                fid = std::get<1>(tup);
                p = std::get<2>(tup);
            }

            uint num_gp = integ.GetNumGP(p);

            uint ze_in_indx = curr_indx;
            uint qx_in_indx = num_gp + curr_indx;
            uint qy_in_indx = 2 * num_gp + curr_indx;

            uint ze_ex_rindx = buff_size - (num_gp + curr_indx);
            uint qx_ex_rindx = buff_size - (2 * num_gp + curr_indx);
            uint qy_ex_rindx = buff_size - (3 * num_gp + curr_indx);

            SWE::Distributed buffs(send_buff_ref,
                                   recv_buff_ref,
                                   ze_in_indx,
                                   qx_in_indx,
                                   qy_in_indx,
                                   ze_ex_rindx,
                                   qx_ex_rindx,
                                   qy_ex_rindx);

            auto tmp_it = pre_boundaries.find(elt);
            assert(tmp_it != pre_boundaries.end());
            auto raw_bdry_iter = tmp_it->second.find(fid);
            assert(raw_bdry_iter != tmp_it->second.end());

            mesh.template CreateDistributedInterface<DistributedInterfaceType>(raw_bdry_iter->second, buffs);

            curr_indx += 3 * num_gp;
        }
    }
    std::ofstream log_file("output/" + mesh.GetMeshName() + "_log", std::ofstream::app);

    log_file << "Number of distributed boundaries: " << mesh.GetNumberDistributedInterfaces() << std::endl;
}

void Problem::initialize_data_kernel(mesh_type& mesh, const MeshMetaData& mesh_data) {
    mesh.CallForEachElement([](auto& elt) { elt.data.initialize(); });

    std::unordered_map<uint, std::vector<double>> bathymetry;

    for (const auto& elt : mesh_data._elements) {
        bathymetry.insert({elt.first, mesh_data.GetBathymetry(elt.first)});
    }

    mesh.CallForEachElement([&bathymetry](auto& elt) {
        uint id = elt.GetID();

        if (!bathymetry.count(id)) {
            throw std::logic_error("Error: could not find bathymetry for element with id: " + id);
        }

        elt.data.state[0].bath = elt.L2Projection(bathymetry[id]);

        elt.ComputeUgp(elt.data.state[0].bath, elt.data.internal.bath_at_gp);

        elt.ComputeDUgp(GlobalCoord::x, elt.data.state[0].bath, elt.data.internal.bath_deriv_wrt_x_at_gp);
        elt.ComputeDUgp(GlobalCoord::y, elt.data.state[0].bath, elt.data.internal.bath_deriv_wrt_y_at_gp);

        auto ze_init = [](Point<2>& pt) { return SWE::true_ze(0, pt); };

        elt.data.state[0].ze = elt.L2Projection(ze_init);

        auto qx_init = [](Point<2>& pt) { return SWE::true_qx(0, pt); };

        elt.data.state[0].qx = elt.L2Projection(qx_init);

        auto qy_init = [](Point<2>& pt) { return SWE::true_qy(0, pt); };

        elt.data.state[0].qy = elt.L2Projection(qy_init);
    });

    mesh.CallForEachInterface([](auto& intface) {
        intface.ComputeUgpIN(intface.data_in.state[0].bath, intface.data_in.boundary[intface.nbound_in].bath_at_gp);
        intface.ComputeUgpEX(intface.data_ex.state[0].bath, intface.data_ex.boundary[intface.nbound_ex].bath_at_gp);
    });

    mesh.CallForEachBoundary([](auto& bound) {
        bound.ComputeUgp(bound.data.state[0].bath, bound.data.boundary[bound.nbound].bath_at_gp);
    });
}

template <typename ElementType>
void Problem::volume_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];
    auto& internal = elt.data.internal;

    // get state at Gauss points
    elt.ComputeUgp(state.ze, internal.ze_at_gp);
    elt.ComputeUgp(state.qx, internal.qx_at_gp);
    elt.ComputeUgp(state.qy, internal.qy_at_gp);

    // assemble flux
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.water_column_hgt_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];

        internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
        internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];

        internal.qx_flux_at_gp[GlobalCoord::x][gp] =
            std::pow(internal.qx_at_gp[gp], 2) / internal.water_column_hgt_at_gp[gp] +
            Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
        internal.qx_flux_at_gp[GlobalCoord::y][gp] =
            internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

        internal.qy_flux_at_gp[GlobalCoord::x][gp] =
            internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
        internal.qy_flux_at_gp[GlobalCoord::y][gp] =
            std::pow(internal.qy_at_gp[gp], 2) / internal.water_column_hgt_at_gp[gp] +
            Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
    }

    // skip dof = 0, which is a constant and thus trivially 0 NOT ALWAYS!
    for (uint dof = 1; dof < elt.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.ze_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.ze_flux_at_gp[GlobalCoord::y]);

        state.rhs_qx[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qx_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qx_flux_at_gp[GlobalCoord::y]);

        state.rhs_qy[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qy_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qy_flux_at_gp[GlobalCoord::y]);
    }
}

template <typename ElementType>
void Problem::source_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];
    auto& internal = elt.data.internal;

    double t = stepper.get_t_at_curr_stage();

    auto source_ze = [t](Point<2>& pt) { return SWE::source_ze(t, pt); };

    auto source_qx = [t](Point<2>& pt) { return SWE::source_qx(t, pt); };

    auto source_qy = [t](Point<2>& pt) { return SWE::source_qy(t, pt); };

    elt.ComputeFgp(source_ze, internal.ze_source_term_at_gp);
    elt.ComputeFgp(source_qx, internal.qx_source_term_at_gp);
    elt.ComputeFgp(source_qy, internal.qy_source_term_at_gp);

    // note we assume that the values at gauss points have already been computed
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        // compute contribution of hydrostatic pressure
        internal.qx_source_term_at_gp[gp] += Global::g * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
        internal.qy_source_term_at_gp[gp] += Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];

        double u_at_gp = internal.qx_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
        double v_at_gp = internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

        // compute bottom friction contribution
        double bottom_friction_stress = Global::Cf * std::hypot(u_at_gp, v_at_gp) / internal.water_column_hgt_at_gp[gp];

        internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
        internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
    }

    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
        state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
        state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
    }
}

template <typename InterfaceType>
void Problem::interface_kernel(const Stepper& stepper, InterfaceType& intface) {
    const uint stage = stepper.get_stage();

    auto& state_in = intface.data_in.state[stage];
    auto& boundary_in = intface.data_in.boundary[intface.nbound_in];

    auto& state_ex = intface.data_ex.state[stage];
    auto& boundary_ex = intface.data_ex.boundary[intface.nbound_ex];

    intface.ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
    intface.ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
    intface.ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

    intface.ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
    intface.ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
    intface.ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);

    // assemble numerical fluxes
    for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(intface.nbound_in); ++gp) {
        uint gp_ex = intface.data_in.get_ngp_boundary(intface.nbound_in) - gp - 1;

        LLF_flux(boundary_in.ze_at_gp[gp],
                 boundary_ex.ze_at_gp[gp_ex],
                 boundary_in.qx_at_gp[gp],
                 boundary_ex.qx_at_gp[gp_ex],
                 boundary_in.qy_at_gp[gp],
                 boundary_ex.qy_at_gp[gp_ex],
                 boundary_in.bath_at_gp[gp],
                 intface.surface_normal[gp],
                 boundary_in.ze_numerical_flux_at_gp[gp],
                 boundary_in.qx_numerical_flux_at_gp[gp],
                 boundary_in.qy_numerical_flux_at_gp[gp]);

        boundary_ex.ze_numerical_flux_at_gp[gp_ex] = -boundary_in.ze_numerical_flux_at_gp[gp];
        boundary_ex.qx_numerical_flux_at_gp[gp_ex] = -boundary_in.qx_numerical_flux_at_gp[gp];
        boundary_ex.qy_numerical_flux_at_gp[gp_ex] = -boundary_in.qy_numerical_flux_at_gp[gp];
    }

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
        state_in.rhs_ze[dof] -= intface.IntegrationPhiIN(dof, boundary_in.ze_numerical_flux_at_gp);
        state_in.rhs_qx[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qx_numerical_flux_at_gp);
        state_in.rhs_qy[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qy_numerical_flux_at_gp);
    }

    for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
        state_ex.rhs_ze[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.ze_numerical_flux_at_gp);
        state_ex.rhs_qx[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qx_numerical_flux_at_gp);
        state_ex.rhs_qy[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qy_numerical_flux_at_gp);
    }
}

template <typename BoundaryType>
void Problem::boundary_kernel(const Stepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.get_stage();

    auto& state = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.nbound];

    bound.ComputeUgp(state.ze, boundary.ze_at_gp);
    bound.ComputeUgp(state.qx, boundary.qx_at_gp);
    bound.ComputeUgp(state.qy, boundary.qy_at_gp);

    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.nbound); ++gp) {
        bound.boundary_condition.set_ex(stepper,
                                        bound.surface_normal[gp],
                                        boundary.ze_at_gp[gp],
                                        boundary.qx_at_gp[gp],
                                        boundary.qy_at_gp[gp],
                                        ze_ex,
                                        qx_ex,
                                        qy_ex);

        LLF_flux(boundary.ze_at_gp[gp],
                 ze_ex,
                 boundary.qx_at_gp[gp],
                 qx_ex,
                 boundary.qy_at_gp[gp],
                 qy_ex,
                 boundary.bath_at_gp[gp],
                 bound.surface_normal[gp],
                 boundary.ze_numerical_flux_at_gp[gp],
                 boundary.qx_numerical_flux_at_gp[gp],
                 boundary.qy_numerical_flux_at_gp[gp]);
    }

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] -= bound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
        state.rhs_qx[dof] -= bound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
        state.rhs_qy[dof] -= bound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
    }
}

template <typename ElementType>
void Problem::update_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();
    double dt = stepper.get_dt();

    auto& state = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    curr_state.rhs_ze = elt.SolveLSE(curr_state.rhs_ze);
    curr_state.rhs_qx = elt.SolveLSE(curr_state.rhs_qx);
    curr_state.rhs_qy = elt.SolveLSE(curr_state.rhs_qy);

    std::fill(next_state.ze.begin(), next_state.ze.end(), 0);
    std::fill(next_state.qx.begin(), next_state.qx.end(), 0);
    std::fill(next_state.qy.begin(), next_state.qy.end(), 0);

    for (uint s = 0; s <= stage; ++s) {
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            next_state.ze[dof] +=
                stepper.ark[stage][s] * state[s].ze[dof] + dt * stepper.brk[stage][s] * state[s].rhs_ze[dof];

            next_state.qx[dof] +=
                stepper.ark[stage][s] * state[s].qx[dof] + dt * stepper.brk[stage][s] * state[s].rhs_qx[dof];

            next_state.qy[dof] +=
                stepper.ark[stage][s] * state[s].qy[dof] + dt * stepper.brk[stage][s] * state[s].rhs_qy[dof];
        }
    }
}

template <typename ElementType>
void Problem::swap_states_kernel(const Stepper& stepper, ElementType& elt) {
    uint n_stages = stepper.get_num_stages();
    auto& state = elt.data.state;

    std::swap(state[0].ze, state[n_stages].ze);
    std::swap(state[0].qx, state[n_stages].qx);
    std::swap(state[0].qy, state[n_stages].qy);
}

template <typename ElementType>
void Problem::scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt) {
    uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];

    for (auto& ze_mode : state.ze) {
        if (isnan(ze_mode)) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qx_mode : state.qx) {
        if (isnan(qx_mode)) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& qy_mode : state.qy) {
        if (isnan(qy_mode)) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_ze_mode : state.rhs_ze) {
        if (isnan(rhs_ze_mode)) {
            std::cerr << "Error: found isnan rhs_ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qx_mode : state.rhs_qx) {
        if (isnan(rhs_qx_mode)) {
            std::cerr << "Error: found isnan rhs_qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }

    for (auto& rhs_qy_mode : state.rhs_qy) {
        if (isnan(rhs_qy_mode)) {
            std::cerr << "Error: found isnan rhs_qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";
        }
    }
}

template <typename ElementType>
void Problem::extract_VTK_data_kernel(ElementType& elt, Array2D<double>& cell_data, Array2D<double>& point_data) {
    elt.WriteCellDataVTK(elt.data.state[0].ze, cell_data[0]);
    elt.WriteCellDataVTK(elt.data.state[0].qx, cell_data[1]);
    elt.WriteCellDataVTK(elt.data.state[0].qy, cell_data[2]);
    elt.WriteCellDataVTK(elt.data.state[0].bath, cell_data[3]);

    elt.WritePointDataVTK(elt.data.state[0].ze, point_data[0]);
    elt.WritePointDataVTK(elt.data.state[0].qx, point_data[1]);
    elt.WritePointDataVTK(elt.data.state[0].qy, point_data[2]);
    elt.WritePointDataVTK(elt.data.state[0].bath, point_data[3]);
}

template <typename MeshType>
void Problem::write_VTK_data_kernel(const Stepper& stepper, MeshType& mesh) {
    Array2D<double> cell_data;
    Array2D<double> point_data;

    cell_data.resize(4);
    point_data.resize(4);

    auto extract_VTK_data_kernel = [&cell_data, &point_data](auto& elt) {
        Problem::extract_VTK_data_kernel(elt, cell_data, point_data);
    };

    mesh.CallForEachElement(extract_VTK_data_kernel);

    std::string file_name = "output/" + mesh.GetMeshName() + "_raw_data.vtk";
    std::ofstream file(file_name);

    file << "CELL_DATA " << (*cell_data.begin()).size() << '\n';
    file << "SCALARS ze_cell float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = cell_data[0].begin(); it != cell_data[0].end(); it++)
        file << *it << '\n';

    file << "SCALARS qx_cell float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = cell_data[1].begin(); it != cell_data[1].end(); it++)
        file << *it << '\n';

    file << "SCALARS qy_cell float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = cell_data[2].begin(); it != cell_data[2].end(); it++)
        file << *it << '\n';

    file << "SCALARS bath_cell float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = cell_data[3].begin(); it != cell_data[3].end(); it++)
        file << *it << '\n';

    file << "POINT_DATA " << (*point_data.begin()).size() << '\n';
    file << "SCALARS ze_point float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = point_data[0].begin(); it != point_data[0].end(); it++)
        file << *it << '\n';

    file << "SCALARS qx_point float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = point_data[1].begin(); it != point_data[1].end(); it++)
        file << *it << '\n';

    file << "SCALARS qy_point float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = point_data[2].begin(); it != point_data[2].end(); it++)
        file << *it << '\n';

    file << "SCALARS bath_point float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto it = point_data[3].begin(); it != point_data[3].end(); it++)
        file << *it << '\n';

    file.close();

    std::string file_name_geom = "output/" + mesh.GetMeshName() + "_geometry.vtk";
    std::string file_name_data = "output/" + mesh.GetMeshName() + "_raw_data.vtk";

    std::ifstream file_geom(file_name_geom, std::ios_base::binary);
    std::ifstream file_data(file_name_data, std::ios_base::binary);

    uint n_step = (uint)(stepper.get_t_at_curr_stage() / stepper.get_dt());

    std::string file_name_merge = "output/" + mesh.GetMeshName() + "_data_" + std::to_string(n_step) + ".vtk";
    std::ofstream file_merge(file_name_merge, std::ios_base::binary);

    file_merge << file_geom.rdbuf() << file_data.rdbuf();
    file_merge.close();
}

template <typename ElementType>
void Problem::extract_modal_data_kernel(ElementType& elt, std::vector<std::pair<uint, Array2D<double>>>& modal_data) {
    modal_data.push_back(std::make_pair(elt.GetID(), Array2D<double>{elt.data.state[0].ze, elt.data.state[0].qx,
                                                                     elt.data.state[0].qy, elt.data.state[0].bath}));
}

template <typename MeshType>
void Problem::write_modal_data_kernel(const Stepper& stepper, MeshType& mesh) {
    std::vector<std::pair<uint, Array2D<double>>> modal_data;

    auto extract_modal_data_kernel = [&modal_data](auto& elt) { Problem::extract_modal_data_kernel(elt, modal_data); };

    mesh.CallForEachElement(extract_modal_data_kernel);

    std::ofstream file;

    std::string file_name = "output/" + mesh.GetMeshName() + "_modal_ze.txt";
    if (stepper.get_t_at_curr_stage() == 0.0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[0].begin(); itt != (*it).second[0].end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << '\n';
        }
    }

    file.close();

    file_name = "output/" + mesh.GetMeshName() + "_modal_qx.txt";
    if (stepper.get_t_at_curr_stage() == 0.0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[1].begin(); itt != (*it).second[1].end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << '\n';
        }
    }

    file.close();

    file_name = "output/" + mesh.GetMeshName() + "_modal_qy.txt";
    if (stepper.get_t_at_curr_stage() == 0.0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[2].begin(); itt != (*it).second[2].end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << '\n';
        }
    }

    file.close();

    file_name = "output/" + mesh.GetMeshName() + "_modal_bath.txt";
    if (stepper.get_t_at_curr_stage() == 0.0) {
        file = std::ofstream(file_name);
    } else {
        file = std::ofstream(file_name, std::ios::app);
    }

    file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
    for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
        for (auto itt = (*it).second[3].begin(); itt != (*it).second[3].end(); itt++) {
            file << (*it).first << ' ' << std::scientific << (*itt) << '\n';
        }
    }

    file.close();
}
}

#endif