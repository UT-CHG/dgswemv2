#ifndef GN_GLOBAL_DATA_HPP
#define GN_GLOBAL_DATA_HPP

namespace GN {
struct GlobalData : SWE::GlobalData {
#ifndef HAS_PETSC
    SparseMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;
#endif

#ifdef HAS_PETSC
    Mat w1_hat_w1_hat;
    Vec w1_hat_rhs;
    KSP dc_ksp;
    PC dc_pc;

    std::vector<double> w1_hat_w1_hat_flat;

    PetscLogStage con_stage;
    PetscLogStage sol_stage;
    PetscLogStage prop_stage;
    PetscLogStage swe_stage;
    PetscLogStage d_stage;

    IS dc_from, dc_to;
    VecScatter dc_scatter;
    Vec dc_sol;

    DynVector<double> dc_solution;

    void destroy() {
        MatDestroy(&w1_hat_w1_hat);
        VecDestroy(&w1_hat_rhs);
        KSPDestroy(&dc_ksp);

        ISDestroy(&dc_from);
        ISDestroy(&dc_to);
        VecScatterDestroy(&dc_scatter);
        VecDestroy(&dc_sol);

#ifdef IHDG_SWE
        SWE::GlobalData::destroy();
#endif
    }
#endif
};
}

#endif
