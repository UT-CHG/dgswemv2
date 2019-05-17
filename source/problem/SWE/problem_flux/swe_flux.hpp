#ifndef SWE_FLUX_HPP
#define SWE_FLUX_HPP

namespace SWE {
void get_F(const std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables>& q,
           const std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries>& aux,
           DynMatrix<double>& F_x,
           DynMatrix<double>& F_y) {
    auto u = vec_cw_div(q[SWE::Variables::qx], aux[SWE::Auxiliaries::h]);
    auto v = vec_cw_div(q[SWE::Variables::qy], aux[SWE::Auxiliaries::h]);

    auto uuh = vec_cw_mult(u, q[SWE::Variables::qx]);
    auto vvh = vec_cw_mult(v, q[SWE::Variables::qy]);
    auto uvh = vec_cw_mult(u, q[SWE::Variables::qy]);
    auto pe  = Global::g * (0.5 * vec_cw_mult(q[SWE::Variables::ze], q[SWE::Variables::ze]) +
                           vec_cw_mult(q[SWE::Variables::ze], aux[SWE::Auxiliaries::bath]));

    row(F_x, SWE::Variables::ze) = q[SWE::Variables::qx];
    row(F_x, SWE::Variables::qx) = uuh + pe;
    row(F_x, SWE::Variables::qy) = uvh;

    row(F_y, SWE::Variables::ze) = q[SWE::Variables::qy];
    row(F_y, SWE::Variables::qx) = uvh;
    row(F_y, SWE::Variables::qy) = vvh + pe;
}

void get_Fn(const HybMatrix<double, SWE::n_variables>& q,
            const HybMatrix<double, SWE::n_auxiliaries>& aux,
            const std::array<DynRowVector<double>, SWE::n_dimensions>& surface_normal,
            std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables>& Fn) {
    auto u = vec_cw_div(row(q, SWE::Variables::qx), row(aux, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(q, SWE::Variables::qy), row(aux, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(q, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(q, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(q, SWE::Variables::qy));
    auto pe  = Global::g * (0.5 * vec_cw_mult(row(q, SWE::Variables::ze), row(q, SWE::Variables::ze)) +
                           vec_cw_mult(row(q, SWE::Variables::ze), row(aux, SWE::Auxiliaries::bath)));

    Fn[SWE::Variables::ze] = vec_cw_mult(row(q, SWE::Variables::qx), surface_normal[GlobalCoord::x]) +
                                  vec_cw_mult(row(q, SWE::Variables::qy), surface_normal[GlobalCoord::y]);
    Fn[SWE::Variables::qx] = vec_cw_mult(uuh + pe, surface_normal[GlobalCoord::x]) +
                                  vec_cw_mult(uvh, surface_normal[GlobalCoord::y]);
    Fn[SWE::Variables::qy] = vec_cw_mult(uvh, surface_normal[GlobalCoord::x]) +
                                  vec_cw_mult(vvh + pe, surface_normal[GlobalCoord::y]);
}

void get_Fn(double gravity,
            const Column<HybMatrix<double, SWE::n_variables>>& q,
            const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux,
            const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
            Column<HybMatrix<double, SWE::n_variables>>&& Fn) {
    double u = q[SWE::Variables::qx] / aux[SWE::Auxiliaries::h];
    double v = q[SWE::Variables::qy] / aux[SWE::Auxiliaries::h];

    double uuh = u * q[SWE::Variables::qx];
    double vvh = v * q[SWE::Variables::qy];
    double uvh = u * q[SWE::Variables::qy];
    double pe =
        gravity * (std::pow(q[SWE::Variables::ze], 2) / 2 + q[SWE::Variables::ze] * aux[SWE::Auxiliaries::bath]);

    double nx = surface_normal[GlobalCoord::x];
    double ny = surface_normal[GlobalCoord::y];

    Fn[SWE::Variables::ze] = q[SWE::Variables::qx] * nx + q[SWE::Variables::qy] * ny;
    Fn[SWE::Variables::qx] = (uuh + pe) * nx + uvh * ny;
    Fn[SWE::Variables::qy] = uvh * nx + (vvh + pe) * ny;
}


void get_dF_dq(const std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables>& q,
               const std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries>& aux,
               HybMatrix<double, SWE::n_variables * SWE::n_variables>& dFx_dq,
               HybMatrix<double, SWE::n_variables * SWE::n_variables>& dFy_dq) {
    auto u = vec_cw_div(q[SWE::Variables::qx], aux[SWE::Auxiliaries::h]);
    auto v = vec_cw_div(q[SWE::Variables::qy], aux[SWE::Auxiliaries::h]);

    // dFx/dq terms
    set_constant(row(dFx_dq, JacobianVariables::ze_ze), 0.0);
    set_constant(row(dFx_dq, JacobianVariables::ze_qx), 1.0);
    set_constant(row(dFx_dq, JacobianVariables::ze_qy), 0.0);

    row(dFx_dq, JacobianVariables::qx_ze) = -vec_cw_mult(u, u) + Global::g * aux[SWE::Auxiliaries::h];
    row(dFx_dq, JacobianVariables::qx_qx) = 2.0 * u;
    set_constant(row(dFx_dq, JacobianVariables::qx_qy), 0.0);

    row(dFx_dq, JacobianVariables::qy_ze) = -vec_cw_mult(u, v);
    row(dFx_dq, JacobianVariables::qy_qx) = v;
    row(dFx_dq, JacobianVariables::qy_qy) = u;

    // dFy/dq terms
    set_constant(row(dFy_dq, JacobianVariables::ze_ze), 0.0);
    set_constant(row(dFy_dq, JacobianVariables::ze_qx), 0.0);
    set_constant(row(dFy_dq, JacobianVariables::ze_qy), 1.0);

    row(dFy_dq, JacobianVariables::qx_ze) = -vec_cw_mult(u, v);
    row(dFy_dq, JacobianVariables::qx_qx) = v;
    row(dFy_dq, JacobianVariables::qx_qy) = u;

    row(dFy_dq, JacobianVariables::qy_ze) = -vec_cw_mult(v, v) + Global::g * aux[SWE::Auxiliaries::h];
    set_constant(row(dFy_dq, JacobianVariables::qy_qx), 0.0);
    row(dFy_dq, JacobianVariables::qy_qy) = 2.0 * v;
}

void get_dFn_dq(const HybMatrix<double, SWE::n_variables>& q,
                const HybMatrix<double, SWE::n_auxiliaries>& aux,
                const std::array<DynRowVector<double>, SWE::n_dimensions>& surface_normal,
                HybMatrix<double, SWE::n_variables * SWE::n_variables>& dFn_dq) {
    auto u = vec_cw_div(row(q, SWE::Variables::qx), row(aux, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(q, SWE::Variables::qy), row(aux, SWE::Auxiliaries::h));

    const DynRowVector<double>& nx = surface_normal[0];
    const DynRowVector<double>& ny = surface_normal[1];

    set_constant(row(dFn_dq, JacobianVariables::ze_ze), 0.0);
    row(dFn_dq, JacobianVariables::ze_qx) = nx;
    row(dFn_dq, JacobianVariables::ze_qy) = ny;

    row(dFn_dq, JacobianVariables::qx_ze) = vec_cw_mult(-vec_cw_mult(u, u) + Global::g * row(aux, SWE::Auxiliaries::h),
                                                        nx) -
                                            vec_cw_mult(vec_cw_mult(u, v), ny);
    row(dFn_dq, JacobianVariables::qx_qx) =
        2.0 * vec_cw_mult(u, nx) + vec_cw_mult(v, ny);
    row(dFn_dq, JacobianVariables::qx_qy) = vec_cw_mult(u, ny);

    row(dFn_dq, JacobianVariables::qy_ze) = -vec_cw_mult(vec_cw_mult(u, v), nx) +
                                            vec_cw_mult(-vec_cw_mult(v, v) + Global::g * row(aux, SWE::Auxiliaries::h),
                                                        ny);
    row(dFn_dq, JacobianVariables::qy_qx) = vec_cw_mult(v, nx);
    row(dFn_dq, JacobianVariables::qy_qy) =
        vec_cw_mult(u, nx) + 2.0 * vec_cw_mult(v, ny);
}
}

#endif