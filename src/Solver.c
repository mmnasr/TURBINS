#include "definitions.h"
#include "DataTypes.h"
#include "Solver.h"

/* Function implementations */


/* This function creates a Petsc type matrix */
Mat Solver_create_LHS_matrix(int size) {

    Mat A;
    int ierr;

    /* Sparse parallel  matrix creation */
    ierr = MatCreate(PETSC_COMM_WORLD, &A); PETScErrAct(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, size, size); PETScErrAct(ierr);
    ierr = MatSetType(A, MATMPIAIJ);
    ierr = MatSetFromOptions(A);

    return A;
}
/*******************************************************************************************************/

/* This function creates a Petsc type vector */
Vec Solver_create_vector(int size) {

    Vec v;
    int ierr;

    ierr = VecCreate(PETSC_COMM_WORLD, &v); PETScErrAct(ierr);
    ierr = VecSetSizes(v, PETSC_DECIDE, size); PETScErrAct(ierr);
    ierr = VecSetFromOptions(v); PETScErrAct(ierr);

    return v;
}
/*******************************************************************************************************/

/* This function returns the suitable solver for the conc. linear system */
KSP Solver_get_concentration_solver(Mat lhs, Parameters *params) {

    KSP solver;
    PC pc;
    double rtol;
    int ierr;

    rtol = params->resmax;

    ierr = KSPCreate(PETSC_COMM_WORLD, &solver); PETScErrAct(ierr);

    ierr = KSPSetOperators(solver, lhs, lhs, SAME_NONZERO_PATTERN); PETScErrAct(ierr);

    //ierr = KSPSetType(solver, KSPGMRES); PETScErrAct(ierr);
    ierr = KSPSetType(solver, KSPIBCGS); PETScErrAct(ierr);

    ierr = KSPSetTolerances(solver, rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); PETScErrAct(ierr);

    ierr = KSPGetPC(solver, &pc); PETScErrAct(ierr);
    ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCILU); PETScErrAct(ierr);

    //ierr = PCSetType(pc, PCBJACOBI); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCLU); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCEISENSTAT); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCSOR); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCJACOBI); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);
    /*
 ierr = KSPGetPC(solver, &pc); PETScErrAct(ierr);
 ierr = PCSetType(pc, PCILU); PETScErrAct(ierr);
 ierr = PCFactorSetLevels(pc, 3); PETScErrAct(ierr);

 //ierr = PCFactorSetUseDropTolerance(pc, 1e-3, .1, 50); PETScErrAct(ierr);
 //ierr = PCFactorSetFill(pc, 30.7648); PETScErrAct(ierr);

 ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); PETScErrAct(ierr);
 ierr = PCFactorSetReuseFill(pc, PETSC_TRUE); PETScErrAct(ierr);
 ierr = PCFactorSetAllowDiagonalFill(pc); PETScErrAct(ierr);

 //ierr = PCFactorSetUseInPlace(pc); PETScErrAct(ierr);
*/
    ierr = KSPSetInitialGuessNonzero(solver, PETSC_TRUE); PETScErrAct(ierr);
    ierr = KSPSetFromOptions(solver); PETScErrAct(ierr);

    return solver;
}
/*******************************************************************************************************/

/* This function returns the suitable solver for pressure. linear system */
KSP Solver_get_pressure_solver(Mat lhs, Parameters *params) {

    KSP solver;
    PC pc;
    double rtol;
    int ierr;
    PetscLogDouble T1, T2;

    rtol = params->resmax;
    ierr = KSPCreate(PETSC_COMM_WORLD, &solver); PETScErrAct(ierr);

    short int hasnullspace = NO;
    if (hasnullspace) {

        MatNullSpace nullsp;
        MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
        KSPSetNullSpace(solver, nullsp);
        //MatNullSpaceDestroy(nullsp);
    }

    ierr = KSPSetOperators(solver, lhs, lhs, SAME_PRECONDITIONER); PETScErrAct(ierr);
    ierr = KSPSetType(solver, KSPGMRES); PETScErrAct(ierr);
    ierr = KSPGMRESSetRestart(solver, 20); PETScErrAct(ierr);
    ierr = KSPSetTolerances(solver, rtol, PETSC_DEFAULT, PETSC_DEFAULT, 300); PETScErrAct(ierr);

    ierr = KSPGetPC(solver, &pc); PETScErrAct(ierr);
    //PCSetType(pc, PCNONE);
    //PCSetType(pc, PCASM);

    PCSetType(pc, PCHYPRE);
    PCHYPRESetType(pc,"boomeramg");

    ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels", "25"); PETScErrAct(ierr);
    ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0"); PETScErrAct(ierr);
    ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all", "SOR/Jacobi"); PETScErrAct(ierr);

    //ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type", ""); PETScErrAct(ierr);
    //ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type", "V"); PETScErrAct(ierr);
    //ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type", "PMIS"); PETScErrAct(ierr);
    //ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_truncfactor", "0.9"); PETScErrAct(ierr);

    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /* Hypre-Petsc Interface.
The most important parameters to be set are
1- Strong Threshold
2- Truncation Factor
3- Coarsennig Type
*/ 

    /* Between 0 to 1 */
    /* "0 "gives better convergence rate (in 3D). */
    /* Suggested values (By Hypre manual): 0.25 for 2D, 0.5 for 3D

 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0"); PETScErrAct(ierr);
*/
    /*******************************************************************************************************/

    /* Available Options:
 "CLJP","Ruge-Stueben","modifiedRuge-Stueben","Falgout", "PMIS", "HMIS"
 Falgout is usually the best.
 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type", "Falgout"); PETScErrAct(ierr);
*/
    /*******************************************************************************************************/

    /* Availble options:
 "local", "global"

 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_measure_type", "local"); PETScErrAct(ierr);
*/
    /*******************************************************************************************************/

    /* Availble options:
 Jacobi,sequential-Gauss-Seidel,
 SOR/Jacobi,backward-SOR/Jacobi,symmetric-SOR/Jacobi,Gaussian-elimination

 Important: If you are using a symmetric KSP solver (like CG), you should use a symmetric smoother
 here.

 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all", "symmetric-SOR/Jacobi"); PETScErrAct(ierr);
*/
    /*******************************************************************************************************/

    /* Available options:
 "V", "W"

 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type", "V"); PETScErrAct(ierr);
*/
    /*******************************************************************************************************/

    /* Availble options:
 "classical", "", "", "direct", "multipass", "multipass-wts", "ext+i", "ext+i-cc", "standard", "standard-wts", "", "", "FF", "FF1"

 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type", ""); PETScErrAct(ierr);
*/
    /*******************************************************************************************************/
    /* Available options:
Greater than zero.
Use zero for the best convergence. However, if you have memory problems, use greate than zero to save some
memory.

 ierr = PetscOptionsSetValue("-pc_hypre_boomeramg_truncfactor", "0.0"); PETScErrAct(ierr);
 */




    /*
Preconditioner Generation Options
PCSetType(pc,PCHYPRE) or -pc_type hypre
-pc_hypre_boomeramg_max_levels nmax
-pc_hypre_boomeramg_truncfactor
-pc_hypre_boomeramg_strong_threshold
-pc_hypre_boomeramg_max_row_sum
-pc_hypre_boomeramg_no_CF
-pc_hypre_boomeramg_coarsen_type CLJP,Ruge-Stueben,modifiedRuge-Stueben,
-pc_hypre_boomeramg_measure_type local,global


Preconditioner Iteration Options
-pc_hypre_boomeramg_relax_type_all Jacobi,sequential-Gauss-Seidel,
   SOR/Jacobi,backward-SOR/Jacobi,symmetric-SOR/Jacobi,Gaussian-eliminat
-pc_hypre_boomeramg_relax_type_fine
-pc_hypre_boomeramg_relax_type_down
-pc_hypre_boomeramg_relax_type_up
-pc_hypre_boomeramg_relax_weight_all r
-pc_hypre_boomeramg_outer_relax_weight_all r
-pc_hypre_boomeramg_grid_sweeps_down n
-pc_hypre_boomeramg_grid_sweeps_up n
-pc_hypre_boomeramg_grid_sweeps_coarse n
-pc_hypre_boomeramg_tol tol
-pc_hypre_boomeramg_max_iter it

*/



    /*
 //ierr = PCSetType(pc, PCASM); PETScErrAct(ierr);
 ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);

 //ierr = PCSetType(pc, PCILU); PETScErrAct(ierr);

 //ierr = PCSetType(pc, PCBJACOBI); PETScErrAct(ierr);
 //ierr = PCSetType(pc, PCLU); PETScErrAct(ierr);
 //ierr = PCSetType(pc, PCEISENSTAT); PETScErrAct(ierr);
 //ierr = PCSetType(pc, PCSOR); PETScErrAct(ierr);
 //ierr = PCSetType(pc, PCJACOBI); PETScErrAct(ierr);
 //ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);

*/
    /*
 ierr = KSPGetPC(solver, &pc); PETScErrAct(ierr);
 ierr = PCSetType(pc, PCILU); PETScErrAct(ierr);
 ierr = PCFactorSetLevels(pc, 3); PETScErrAct(ierr);

 //ierr = PCFactorSetUseDropTolerance(pc, 1e-3, .1, 50); PETScErrAct(ierr);
 //ierr = PCFactorSetFill(pc, 30.7648); PETScErrAct(ierr);

 ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); PETScErrAct(ierr);
 ierr = PCFactorSetReuseFill(pc, PETSC_TRUE); PETScErrAct(ierr);
 ierr = PCFactorSetAllowDiagonalFill(pc); PETScErrAct(ierr);

 //ierr = PCFactorSetUseInPlace(pc); PETScErrAct(ierr);
*/
    ierr = KSPSetInitialGuessNonzero(solver, PETSC_TRUE); PETScErrAct(ierr);
    ierr = PetscGetTime(&T1);PETScErrAct(ierr);
    ierr = KSPSetFromOptions(solver); PETScErrAct(ierr);

    ierr = PetscGetTime(&T2);PETScErrAct(ierr);
    PetscPrintf(PCW, "Setup time for the Pressure solver was:%f\n", (T2 - T1));
    ierr = KSPView(solver, PETSC_VIEWER_STDOUT_WORLD); PETScErrAct(ierr);
    return solver;

}
/*******************************************************************************************************/

/* This function returns the suitable solver for the velocity linear system */
KSP Solver_get_velocity_solver(Mat lhs, Parameters *params) {

    KSP solver;
    PC pc;
    double rtol;
    int ierr;

    rtol = params->resmax;

    ierr = KSPCreate(PETSC_COMM_WORLD, &solver); PETScErrAct(ierr);
    ierr = KSPSetOperators(solver, lhs, lhs, SAME_NONZERO_PATTERN); PETScErrAct(ierr);

    //ierr = KSPSetType(solver, KSPGMRES); PETScErrAct(ierr);
    ierr = KSPSetType(solver, KSPIBCGS); PETScErrAct(ierr);

    ierr = KSPSetTolerances(solver, rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); PETScErrAct(ierr);

    ierr = KSPGetPC(solver, &pc); PETScErrAct(ierr);
    ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCJACOBI); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);

    //ierr = PCSetType(pc, PCASM); PETScErrAct(ierr);

    //ierr = PCSetType(pc, PCBJACOBI); PETScErrAct(ierr);

    //ierr = PCSetType(pc, PCBJACOBI); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCLU); PETScErrAct(ierr);

    //ierr = PCSetType(pc, PCEISENSTAT); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCSOR); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCJACOBI); PETScErrAct(ierr);
    //ierr = PCSetType(pc, PCNONE); PETScErrAct(ierr);

    /*
 ierr = KSPGetPC(solver, &pc); PETScErrAct(ierr);
 ierr = PCSetType(pc, PCILU); PETScErrAct(ierr);
 ierr = PCFactorSetLevels(pc, 3); PETScErrAct(ierr);

 //ierr = PCFactorSetUseDropTolerance(pc, 1e-3, .1, 50); PETScErrAct(ierr);
 //ierr = PCFactorSetFill(pc, 30.7648); PETScErrAct(ierr);

 ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); PETScErrAct(ierr);
 ierr = PCFactorSetReuseFill(pc, PETSC_TRUE); PETScErrAct(ierr);
 ierr = PCFactorSetAllowDiagonalFill(pc); PETScErrAct(ierr);

 //ierr = PCFactorSetUseInPlace(pc); PETScErrAct(ierr);
*/
    ierr = KSPSetInitialGuessNonzero(solver, PETSC_TRUE); PETScErrAct(ierr);
    ierr = KSPSetFromOptions(solver); PETScErrAct(ierr);


    return solver;
}
/*******************************************************************************************************/

/* This function checks if the linear system did converge */
short int Solver_is_converged(KSP solver, char *lsys_name) {

    KSPConvergedReason reason;
    int ierr = KSPGetConvergedReason(solver, &reason); PETScErrAct(ierr);
    if (reason < 0) {

        PetscPrintf(PCW, "Solver.c/ Warning! %s-linear system did not converged! Reason:%d\n", lsys_name, reason);
        return NO;
    } /* if */
    /* diverged */
/*
    KSP_DIVERGED_NULL                = -2,
    KSP_DIVERGED_ITS                 = -3,
    KSP_DIVERGED_DTOL                = -4,
    KSP_DIVERGED_BREAKDOWN           = -5,
    KSP_DIVERGED_BREAKDOWN_BICG      = -6,
    KSP_DIVERGED_NONSYMMETRIC        = -7,
    KSP_DIVERGED_INDEFINITE_PC       = -8,
    KSP_DIVERGED_NAN                 = -9,
    KSP_DIVERGED_INDEFINITE_MAT      = -10,
*/
    return YES;
}
