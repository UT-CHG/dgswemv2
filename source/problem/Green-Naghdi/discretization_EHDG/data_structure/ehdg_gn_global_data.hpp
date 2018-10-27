#ifndef EHDG_GN_GLOBAL_DATA_HPP
#define EHDG_GN_GLOBAL_DATA_HPP

namespace GN {
namespace EHDG {
struct GlobalData : SWE::EHDG::GlobalData {
#ifndef HAS_PETSC
    SparseMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;
#endif

#ifdef HAS_PETSC
    Mat w1_hat_w1_hat;
    Vec w1_hat_rhs;
    KSP dc_ksp;
    PC dc_pc;

    IS dc_from, dc_to;
    VecScatter dc_scatter;
    Vec dc_sol;

    DynVector<double> dc_solution;

    void destroy() {
        SWE::EHDG::GlobalData::destroy();

        MatDestroy(&w1_hat_w1_hat);
        VecDestroy(&w1_hat_rhs);
        KSPDestroy(&dc_ksp);

        ISDestroy(&dc_from);
        ISDestroy(&dc_to);
        VecScatterDestroy(&dc_scatter);
        VecDestroy(&dc_sol);
    }
#endif
};
}
}

#endif
