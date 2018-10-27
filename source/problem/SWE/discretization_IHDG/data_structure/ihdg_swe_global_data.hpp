#ifndef IHDG_SWE_GLOBAL_DATA_HPP
#define IHDG_SWE_GLOBAL_DATA_HPP

namespace SWE {
namespace IHDG {
struct GlobalData {
#ifndef HAS_PETSC
    SparseMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;
#endif

#ifdef HAS_PETSC
    Mat delta_hat_global;
    Vec rhs_global;
    KSP ksp;
    PC pc;

    IS from, to;
    VecScatter scatter;
    Vec sol;

    DynVector<double> solution;
    bool converged = false;

    void destroy() {
        MatDestroy(&delta_hat_global);
        VecDestroy(&rhs_global);
        KSPDestroy(&ksp);

        ISDestroy(&from);
        ISDestroy(&to);
        VecScatterDestroy(&scatter);
        VecDestroy(&sol);
    }
#endif
};
}
}

#endif
