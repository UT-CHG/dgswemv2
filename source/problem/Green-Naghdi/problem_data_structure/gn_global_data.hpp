#ifndef GN_GLOBAL_DATA_HPP
#define GN_GLOBAL_DATA_HPP

namespace GN {
struct GlobalData : SWE::GlobalData {
#ifndef HAS_PETSC
    SparseMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;
    DynVector<double> derivatives_at_node;
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

    Vec global_derivatives_at_node;
    Vec global_node_mult;
    DynVector<double> derivatives_at_node;
    std::vector<int> local_nodeIDs;

    IS derivatives_from, derivatives_to;
    VecScatter derivatives_scatter;
    Vec local_derivatives_at_node;

    IS node_mult_from, node_mult_to;
    VecScatter node_mult_scatter;
    Vec local_node_mult;

    PetscLogStage con_stage;
    PetscLogStage sol_stage;
    PetscLogStage prop_stage;
    PetscLogStage swe_stage;
    PetscLogStage d_stage;

    void destroy() {
        MatDestroy(&w1_hat_w1_hat);
        VecDestroy(&w1_hat_rhs);
        KSPDestroy(&dc_ksp);

        ISDestroy(&dc_from);
        ISDestroy(&dc_to);
        VecScatterDestroy(&dc_scatter);
        VecDestroy(&dc_sol);

#ifdef D_RECON_AVG
        VecDestroy(&global_derivatives_at_node);
        VecDestroy(&global_node_mult);

        ISDestroy(&derivatives_from);
        ISDestroy(&derivatives_to);
        VecScatterDestroy(&derivatives_scatter);
        VecDestroy(&local_derivatives_at_node);

        ISDestroy(&node_mult_from);
        ISDestroy(&node_mult_to);
        VecScatterDestroy(&node_mult_scatter);
        VecDestroy(&local_node_mult);
#endif

#ifdef IHDG_SWE
        SWE::GlobalData::destroy();
#endif
    }
#endif
};
}

#endif
