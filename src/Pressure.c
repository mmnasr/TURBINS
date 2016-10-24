#include "definitions.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Solver.h"
#include "Grid.h"
#include "Pressure.h"
#include "Communication.h"
#include "gvg.h"
#include "Output.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* This function allocates memory for pressure structure */
Pressure *Pressure_create( MAC_grid *grid, Parameters *params) {

    Pressure *new_p;
    int NX, NY, NZ;
    int NT;
    int NX_cell, NY_cell, NZ_cell;
    int ierr;

    new_p = (Pressure *)malloc(sizeof(Pressure));
    Memory_check_allocation(new_p);

    /* Number of real physical cells */
    NX_cell = params->NX;
    NY_cell = params->NY;
    NZ_cell = params->NZ;

    /* Add one cell in each direction to the far most end */
    NX = NX_cell + 1; /* Number of grid points in x-direction (half cell added in the far most direction) */
    NY = NY_cell + 1; /* Number of grid points in y-direction (half cell added in the far most direction)*/
    NZ = NZ_cell + 1; /* Number of grid points in z-direction (half cell added in the far most direction)*/
    NT = NX*NY*NZ;

    /***************************************/
    /* Now, create the vectors that hold the local and global data for all the 3D properties */
    ierr = DACreateGlobalVector(grid->DA_3D, &new_p->G_data); PETScErrAct(ierr);
    ierr = DACreateLocalVector(grid->DA_3D, &new_p->L_data); PETScErrAct(ierr);

    /* Pressure gradients */
    ierr = VecDuplicate(new_p->G_data, &new_p->G_dp_dx); PETScErrAct(ierr);
    ierr = VecDuplicate(new_p->G_data, &new_p->G_dp_dy); PETScErrAct(ierr);
    ierr = VecDuplicate(new_p->G_data, &new_p->G_dp_dz); PETScErrAct(ierr);

    ierr = VecDuplicate(new_p->L_data, &new_p->L_dp_dx); PETScErrAct(ierr);
    ierr = VecDuplicate(new_p->L_data, &new_p->L_dp_dy); PETScErrAct(ierr);
    ierr = VecDuplicate(new_p->L_data, &new_p->L_dp_dz); PETScErrAct(ierr);

    ierr = VecDuplicate(new_p->G_data, &new_p->G_divergence); PETScErrAct(ierr);

    /* Local and global {rhs}={b} vectors */
    ierr = VecDuplicate(new_p->G_data, &new_p->G_b); PETScErrAct(ierr);
    ierr = VecDuplicate(new_p->L_data, &new_p->L_b); PETScErrAct(ierr);

    /*new_p->A = Solver_create_LHS_matrix(NT);*/
    /* Create parallel sparse matrix for the LHS of the matrix */
    ierr = DAGetMatrix(grid->DA_3D_Lsys, MATMPIAIJ, &new_p->A); PETScErrAct(ierr);

    /* Index of left-bottom-back corner on current processor for the DA data layout */
    new_p->G_Is = grid->G_Is;
    new_p->G_Js = grid->G_Js;
    new_p->G_Ks = grid->G_Ks;

    /* Index of right-top-front corner on current processor for the DA data layout */
    new_p->G_Ie = grid->G_Ie;
    new_p->G_Je = grid->G_Je;
    new_p->G_Ke = grid->G_Ke;

    /* Index of left-bottom-back corner on current processor for the DA data layout including ghost nodes	*/
    new_p->L_Is = grid->L_Is;
    new_p->L_Js = grid->L_Js;
    new_p->L_Ks = grid->L_Ks;

    /* Index of right-top-front corner on current processor for the DA data layout including ghost nodes*/
    new_p->L_Ie = grid->L_Ie;
    new_p->L_Je = grid->L_Je;
    new_p->L_Ke = grid->L_Ke;

    /* Global number of cells */
    new_p->NX = NX;
    new_p->NY = NY;
    new_p->NZ = NZ;
    new_p->NT = NT;

    /* The last interior index of cell */
    new_p->NI = NX-1;
    new_p->NJ = NY-1;
    new_p->NK = NZ-1;


    return new_p;
}
/**************************************************************************************************************/

/* This function releases the allocated memory for pressure structure. */
void Pressure_destroy(Pressure *p) {

    int ierr;

    ierr = VecDestroy(p->G_data); PETScErrAct(ierr);
    ierr = VecDestroy(p->L_data); PETScErrAct(ierr);

    ierr = VecDestroy(p->G_dp_dx); PETScErrAct(ierr);
    ierr = VecDestroy(p->G_dp_dy); PETScErrAct(ierr);
    ierr = VecDestroy(p->G_dp_dz); PETScErrAct(ierr);

    ierr = VecDestroy(p->L_dp_dx); PETScErrAct(ierr);
    ierr = VecDestroy(p->L_dp_dy); PETScErrAct(ierr);
    ierr = VecDestroy(p->L_dp_dz); PETScErrAct(ierr);

    ierr = VecDestroy(p->G_divergence); PETScErrAct(ierr);

    /* Free linear system variables */
    ierr = VecDestroy(p->G_b); PETScErrAct(ierr);
    ierr = VecDestroy(p->L_b); PETScErrAct(ierr);

    ierr = MatDestroy(p->A); PETScErrAct(ierr);
    ierr = KSPDestroy(p->solver); PETScErrAct(ierr);

    free(p);

}
/***************************************************************************************************/

/* This function creates the vectors and matrix [A] for the pressure linear system (poisson equation), also it sets up the KSP
solver for the this linear system. */
void Pressure_setup_lsys_accounting_geometry(Pressure *p, MAC_grid *grid, Parameters *params) {

    int ierr;

    if (0) {
    //PetscViewer reader;
    //ierr = PetscViewerBinaryOpen(PCW, "../Resume_pressure_matrix.bin.gz", FILE_MODE_READ, &reader); PETScErrAct(ierr);
    //MatLoad(reader, PETSC_DECIDE, &(p->A));
    //ierr = PetscViewerDestroy(reader); PETScErrAct(ierr);
    }

    Pressure_form_LHS_matrix(p, grid);


#ifdef MEMORY_PROFILING
    PetscMallocGetCurrentUsage(&mem1);
#endif
    p->solver = Solver_get_pressure_solver(p->A, params);
    PetscPrintf(PCW, "Pressure.c/ Solver has been returned successfully\n");
    /* Since LHS matrix does not change for Pressure, we just create the Preconditioner once */
    ierr = KSPSetOperators(p->solver, p->A, p->A, SAME_PRECONDITIONER); PETScErrAct(ierr);
    ierr = KSPSetUp(p->solver); PETScErrAct(ierr);
#ifdef MEMORY_PROFILING
    PetscMallocGetCurrentUsage(&mem2);
    p_solver_mem += (int)(mem2 - mem1);
    p_mem += (int)(mem2 - mem1);
    printf("Pressure.c/ memory:%d\n", (int)(mem2-mem1)); getchar();
#endif

    if (params->output_lsys_matrix) {

        char filename[FILENAME_MAX];
        sprintf(filename, "Resume_pressure_matrix.bin.gz");
        Output_petsc_mat(&(p->A), filename, BINARY);

        PetscPrintf(PCW, "Pressure.c/ p-matrix has been written to %s successfully.\n", filename);
    } /* if */

}
/**************************************************************************************************************/

/* This function computes the RHS vector {b} for the pressure linear set of equations */
void Pressure_set_RHS(Pressure *p, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, double dt) {

    int i, j, k;
    int NI, NJ, NK;
    int ierr;
    double dudx, dvdy, dwdz;
    double ***u_data, ***v_data, ***w_data;
    double ***rhs_vec;
    double **u_outflow, **u_inflow;
    double a_dt;
    double *metric_xc, *metric_yc, *metric_zc;
    double a_dEta, a_dXi, a_dPhi;
    double C_x, C_y, C_z;
    double rhs_value;
    int Is, Js, Ks;
    int Ie, Je, Ke;


    /* Get regular data array for velocities data */
    /* Get the local velocity data on each processor including the ghost nodes */
    ierr = DAVecGetArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

    /* Now, got the RHS vector on current processor */
    ierr = DAVecGetArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    NI = p->NI;
    NJ = p->NJ;
    NK = p->NK;

    u_outflow = u->G_outflow;
    u_inflow  = u->G_inflow;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_zc = grid->metric_zc;

    a_dEta = 1.0/grid->dEta;
    a_dXi  = 1.0/grid->dXi;
    a_dPhi = 1.0/grid->dPhi;

    a_dt   = 1.0/dt;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* go over all the nodes and find div.Vstar */
    /* Exclude the last half cell added to the very end */
    for (k=Ks; k<Ke; k++) {

        C_z = metric_zc[k] * a_dPhi;
        for (j=Js; j<Je; j++) {

            C_y = metric_yc[j] * a_dXi;
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) != BOUNDARY)  {

                    C_x = metric_xc[i] * a_dEta;
                    /* Div_V_star at cell center p(i,j,k) */
                    /* Regular node */
                    dudx = (u_data[k][j][i+1] - u_data[k][j][i]) * C_x; /*use central differences*/

                    if ( (params->inflow) && (i == 0) ) {

                        dudx = (u_data[k][j][i+1] - u_inflow[k][j]) * C_x; /*use central differences*/
                    } /* if inflow */
                    if ( (params->outflow) && (i == NI-1) ) {

                        dudx = (u_outflow[k][j] - u_data[k][j][i]) * C_x; /*use central differences*/ /* get rid of */
                    } /* if outflow */

                    dvdy = (v_data[k][j+1][i] - v_data[k][j][i]) * C_y; /*use central differences*/
                    dwdz = (w_data[k+1][j][i] - w_data[k][j][i]) * C_z; /*use central differences*/

                    /*multiply by -1 because LHS is multiplied by -1 for Positive-Definite lsys */
                    rhs_value = -(dudx + dvdy + dwdz) * a_dt;

                } else {

                    rhs_value = 0.0;
                } /* else */

                /* insert the rhs_value into the rhs_vector */
                rhs_vec[k][j][i] = rhs_value;

            } /* for i*/
        } /* for j*/
    } /* for k*/

    /* Always restore the array after using the vector */
    ierr = DAVecRestoreArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

}
/**************************************************************************************************************/

/* This function updates the velocity field using projection method to get a divergence free velocity field */
void Pressure_project_velocity(Pressure *p, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, double dt) {

    int i, j, k;
    double ***phi;
    double ***u_data, ***v_data, ***w_data;
    double ***dp_dx, ***dp_dy, ***dp_dz;
    double ***div_V_star;
    double a_dEta, a_dXi, a_dPhi;
    double *metric_xu, *metric_yv, *metric_zw;
    double metric_ex;
    double C_y, C_z;
    double dphi_dx_ugrid;
    double dphi_dy_vgrid;
    double dphi_dz_wgrid;
    double coef, Re;
    double laplacian_term;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Get velocity data. No ghost node is required. So, Global data is retrieved */
    ierr = DAVecGetArray(grid->DA_3D, u->G_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->G_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->G_data, (void ***)&w_data); PETScErrAct(ierr);

    /* Phi. Get the ghost nodes from the neighboring processors */
    Communication_update_ghost_nodes(&grid->DA_3D, &p->G_data, &p->L_data, 'I');

    /* Now, get "phi" (Poisson equation variable) local data since ghost nodes are required */
    ierr = DAVecGetArray(grid->DA_3D, p->L_data, (void ***)&phi); PETScErrAct(ierr);

    /* get p-gradients */
    ierr = DAVecGetArray(grid->DA_3D, p->G_dp_dx, (void ***)&dp_dx); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, p->G_dp_dy, (void ***)&dp_dy); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, p->G_dp_dz, (void ***)&dp_dz); PETScErrAct(ierr);

    /* Update local ghost nodes for rhs of poisson equation: div_v_star. We will need this to update the v_* velocity field */
    Communication_update_ghost_nodes(&grid->DA_3D, &p->G_b, &p->L_b, 'I');
    ierr = DAVecGetArray(grid->DA_3D, p->L_b, (void ***)&div_V_star); PETScErrAct(ierr);

    Re = params->Re;

    /* To save time */
    coef = dt/Re;

    a_dEta = 1.0/grid->dEta;
    a_dXi  = 1.0/grid->dXi;
    a_dPhi = 1.0/grid->dPhi;

    metric_xu = grid->metric_xu;
    metric_yv = grid->metric_yv;
    metric_zw = grid->metric_zw;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* Update u_star to u_new (divergence free velocity field) */
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                short int status = Grid_get_u_status(grid, i, j, k);
                if ( (status != BOUNDARY) ) {//&& (status != SOLID) ){

                    /* (dx/dEta)^-1 at u-grid (i+1/2,j,k) */
                    metric_ex = metric_xu[i];
                    dphi_dx_ugrid = (phi[k][j][i] - phi[k][j][i-1]) * a_dEta * metric_ex;

                    /* Update u-vel to get divergence free velocity field */
                    u_data[k][j][i] -= dt * dphi_dx_ugrid;


                    /* Now, update the pressure gradient */
                    /* Corresponds to d/dx( div.V_star )/dt */
                    laplacian_term   = -(div_V_star[k][j][i] - div_V_star[k][j][i-1]) * metric_ex * a_dEta;

                    /* Update dp_dx)new = dp_dx_old + dphi_dx -dt/Re d/dx (Div_V_star) */
                    dp_dx[k][j][i]  += dphi_dx_ugrid - coef * laplacian_term;

                } /* if */
            } /* for i */
        } /* for j*/
    } /* for k*/

    /* Update v_star to v_new (divergence free velocity field) */
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {

            /* (dy/dXi)^-1 at v-grid (i,j+1/2,k) */
            C_y = a_dXi * metric_yv[j];
            for (i=Is; i<Ie; i++) {

                short int status = Grid_get_v_status(grid, i, j, k);
                if ( (status != BOUNDARY) ) {//&& (status != SOLID) ){

                    dphi_dy_vgrid = (phi[k][j][i] - phi[k][j-1][i]) * C_y;

                    /* Update u-vel to get divergence free velocity field */
                    v_data[k][j][i] -= dt * dphi_dy_vgrid;

                    /* Corresponds to d/dy( div.V_star )/dt */
                    laplacian_term   = -(div_V_star[k][j][i] - div_V_star[k][j-1][i]) * C_y;

                    /* Update dp_dy)new = dp_dy)old + dphi_dy -dt/Re d/dy (Div_V_star) */
                    dp_dy[k][j][i]  += dphi_dy_vgrid - coef * laplacian_term;

                } /* if */
            } /* for i */
        } /* for j*/
    } /* for k*/


    /* Update w_star to w_new (divergence free velocity field) */
    for (k=Ks; k<Ke; k++) {

        /* (dz/dPhi)^-1 at w-grid (i,j,k+1/2) */
        C_z = metric_zw[k] * a_dPhi;
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                short int status = Grid_get_w_status(grid, i, j, k);
                if ( (status != BOUNDARY) ) {//&& (status != SOLID) ){

                    dphi_dz_wgrid = (phi[k][j][i] - phi[k-1][j][i]) * C_z;

                    /* Update w-vel to get divergence free velocity field */
                    w_data[k][j][i] -= dt * dphi_dz_wgrid;

                    /* Corresponds to d/dz( div.V_star )/dt */
                    laplacian_term   = -(div_V_star[k][j][i] - div_V_star[k-1][j][i]) * C_z;

                    /* Update dp_dy)new = dp_dy)old + dphi_dy -dt/Re d/dy (Div_V_star) */
                    dp_dz[k][j][i]  += dphi_dz_wgrid - coef * laplacian_term;

                } /* if */
            } /* for i */
        } /* for j*/
    } /* for k*/

    /* restore arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, p->L_b, (void ***)&div_V_star); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_data, (void ***)&w_data); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, p->G_dp_dx, (void ***)&dp_dx); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->G_dp_dy, (void ***)&dp_dy); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->G_dp_dz, (void ***)&dp_dz); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, p->L_data, (void ***)&phi); PETScErrAct(ierr);

}
/**************************************************************************************************************/

/* This function solves the linear system to find pressure (phi) at each cell center.*/
int Pressure_solve(Pressure *p) {

    int iters;
    int ierr;

    /* Since LHS matrix is constant, no update is needed for the preconditioner */
    /* The solution is stored in data array */
    ierr = KSPSolve(p->solver, p->G_b, p->G_data); PETScErrAct(ierr);

    (void)Solver_is_converged(p->solver, (char *)"pressue");

    ierr = KSPGetIterationNumber(p->solver, &iters); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, " ** Pressure soln took %d iters.\n", iters); PETScErrAct(ierr);

    return iters;
}
/**************************************************************************************************************/

/* This function finds the maximum divergence of the velocity field to verify the solution, i.e. div_max --> 0.0 */
double Pressure_compute_velocity_divergence(Pressure *p, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) {

    int i, j, k;
    int NI, NJ, NK;
    double dudx, dvdy, dwdz;
    double ***divergence;
    double ***u_data, ***v_data, ***w_data;
    double **u_outflow, **u_inflow;
    double a_dEta, a_dXi, a_dPhi;
    double *metric_xc, *metric_yc, *metric_zc;
    double metric_px, metric_py, metric_pz;
    double div_cell;
    double div_max   = 0.0; /* local */
    double W_div_max = 0.0; /* (World) global */
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Generate the local u, v, w data (including the ghost nodes from neighboring processors) */
    /* u-vel */
    Communication_update_ghost_nodes(&grid->DA_3D, &u->G_data, &u->L_data, 'I');

    /* v-vel */
    Communication_update_ghost_nodes(&grid->DA_3D, &v->G_data, &v->L_data, 'I');

    /* w-vel */
    Communication_update_ghost_nodes(&grid->DA_3D, &w->G_data, &w->L_data, 'I');

    /* Get velocity data. Including ghost nodes */
    ierr = DAVecGetArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

    /* Get divergence data. No ghost node is required. So, Global data is retrieved */
    ierr = DAVecGetArray(grid->DA_3D, p->G_divergence, (void ***)&divergence); PETScErrAct(ierr);

    u_outflow  = u->G_outflow;
    u_inflow   = u->G_inflow;

    NI = p->NI;
    NJ = p->NJ;
    NK = p->NK;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_zc = grid->metric_zc;

    a_dEta = 1.0/grid->dEta;
    a_dXi  = 1.0/grid->dXi;
    a_dPhi = 1.0/grid->dPhi;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    int i_max=-1;
    int j_max=-1;
    int k_max=-1;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) != BOUNDARY) {

                    /* (dx/dEta)^-1 at (i,j,k) */
                    metric_px = metric_xc[i];
                    /* (dy/dXi)^-1 at (i,j,k) */
                    metric_py = metric_yc[j];
                    /* (dx/dPhi)^-1 at (i,j,k) */
                    metric_pz = metric_zc[k];

                    dudx = (u_data[k][j][i+1] - u_data[k][j][i]) * a_dEta * metric_px; /*use central differences*/
                    if ( (params->inflow) && (i == 0) ) {

                        dudx = (u_data[k][j][i+1] - u_inflow[k][j]) * a_dEta * metric_px; /*use central differences*/
                    }
                    if ( (params->outflow) && (i == NI-1) ) {

                        dudx = (u_outflow[k][j] - u_data[k][j][i]) * a_dEta * metric_px; /*use central differences*/
                    }
                    dvdy = ( v_data[k][j+1][i] - v_data[k][j][i]) * a_dXi  * metric_py;
                    dwdz = ( w_data[k+1][j][i] - w_data[k][j][i]) * a_dPhi * metric_pz;

                    div_cell = fabs(dudx + dvdy + dwdz);

                } else {

                    div_cell = 0.0;

                }
                divergence[k][j][i] = div_cell;

                /* Now check the maximum divergence in the whole domain */
                if ( div_cell >= div_max) {

                    div_max = div_cell;
                    i_max = i;
                    j_max = j;
                    k_max = k;
                }

            } /* for i*/
        } /* for j*/
    } /* for k*/

    /* restore arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, p->G_divergence, (void ***)&divergence); PETScErrAct(ierr);

    /* Now, find the maximum global divergence on all processors */
    (void)MPI_Allreduce (&div_max, &W_div_max, 1, MPI_DOUBLE, MPI_MAX, PCW);

    PetscPrintf(PCW, "Pressure.c/ div:%2.14f at (i,j,k)=(%d,%d,%d)\n", W_div_max, i_max, j_max, k_max);
    return (W_div_max);
}
/****************************************************************************************************************/

/* This function updates the real pressure (for output purposes) */
void Pressure_compute_real_pressure(Pressure *p, MAC_grid *grid) {

    int i, j, k;
    int NI, NJ, NK;
    int ierr;
    double ***rhs_vec;
    double ***dp_dx, ***dp_dy, ***dp_dz;
    double *metric_xc, *metric_yc, *metric_zc;
    double a_dEta, a_dXi, a_dPhi;
    double C_x, C_y, C_z;
    double rhs_value;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    NI = p->NI;
    NJ = p->NJ;
    NK = p->NK;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_zc = grid->metric_zc;

    a_dEta = 1.0/grid->dEta;
    a_dXi  = 1.0/grid->dXi;
    a_dPhi = 1.0/grid->dPhi;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* update the local version to get the ghost nodes values from the neighboring processors */
    Communication_update_ghost_nodes(&grid->DA_3D, &p->G_dp_dx, &p->L_dp_dx, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D, &p->G_dp_dy, &p->L_dp_dy, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D, &p->G_dp_dz, &p->L_dp_dz, 'I');

    /* rhs vector: div. Grad p */
    ierr = DAVecGetArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    /* get p gradients */
    ierr = DAVecGetArray(grid->DA_3D, p->L_dp_dx, (void ***)&dp_dx); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, p->L_dp_dy, (void ***)&dp_dy); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, p->L_dp_dz, (void ***)&dp_dz); PETScErrAct(ierr);

    /* go over all the nodes and find div.(grap_p) */
    /* Exclude the last half cell added to the very end */
    double ddx=0.0;
    double ddy=0.0;
    double ddz=0.0;

    for (k=Ks; k<Ke; k++) {

        C_z = metric_zc[k] * a_dPhi;
        for (j=Js; j<Je; j++) {

            C_y = metric_yc[j] * a_dXi;
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) != BOUNDARY)  {

                    C_x = metric_xc[i] * a_dEta;
                    /* Div_grad_p at cell center p(i,j,k) */
                    /* Regular node */
                    ddx = (dp_dx[k][j][i+1] - dp_dx[k][j][i]) * C_x; /*use central differences*/
                    ddy = (dp_dy[k][j+1][i] - dp_dy[k][j][i]) * C_y; /*use central differences*/
                    ddz = (dp_dz[k+1][j][i] - dp_dz[k][j][i]) * C_z; /*use central differences*/

                    /*multiply by -1 because LHS is multiplied by -1 for Positive-Definite lsys */
                    rhs_value = -(ddx + ddy + ddz);

                } else {

                    rhs_value = 0.0;
                } /* else */

                /* insert the rhs_value into the rhs_vector */
                rhs_vec[k][j][i] = rhs_value;

            } /* for i*/
        } /* for j*/
    } /* for k*/

    /* restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, p->L_dp_dx, (void ***)&dp_dx); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->L_dp_dy, (void ***)&dp_dy); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->L_dp_dz, (void ***)&dp_dz); PETScErrAct(ierr);

    Pressure_solve(p);
    PetscPrintf(PCW, "Pressure.c/ Solving Poisson equation to compute real pressure. (For output purposes). \n");

}
/****************************************************************************************************************/

/* This functions setups up the matrix for the solution of pressure Poisson equation */
void Pressure_form_LHS_matrix(Pressure *p, MAC_grid *grid) {

    PetscReal ae, aw, ap, as, an, af, ab;
    PetscReal ap_solid;
    int i, j, k;
    int NI, NJ, NK;
    int ierr;
    double dXi, dEta, dPhi;
    double Cx, Cy, Cz;
    double metric_px, metric_py, metric_pz;
    double metric_ex, metric_wx;
    double metric_ny, metric_sy;
    double metric_bz, metric_fz;
    double *metric_xc, *metric_yc, *metric_zc;
    double *metric_xu, *metric_yv, *metric_zw;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    MatStencil cols[STENCIL], row;
    double values[STENCIL];
    int nnz;

    NI  = p->NI;
    NJ  = p->NJ;
    NK  = p->NK;

    /*Assemble matrix*/
    /*
  a_? =[A]
  -p: diagonal
  -s: south
  -w: west
  -n: north
  -e: east
  -b: back
  -f: front
*/

    dEta = grid->dEta;
    dXi  = grid->dXi;
    dPhi = grid->dPhi;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_zc = grid->metric_zc;
    metric_xu = grid->metric_xu;
    metric_yv = grid->metric_yv;
    metric_zw = grid->metric_zw;

    /*Assemble matrix*/
    Cx  = -1.0/(dEta * dEta);
    Cy  = -1.0/(dXi * dXi);
    Cz  = -1.0/(dPhi * dPhi);

    /* To be inserted into the diagonal part of matrix for a solid node */
    /* The value is not important. It is chosen in a way to be the same order of the other diagonal elements */
    ap_solid = -Cx*metric_xc[0]-Cy*metric_yc[0]-Cz*metric_zc[0];

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                /* The 3D index of current row in the Matrix */
                row.i = i;
                row.j = j;
                row.k = k;
                nnz   = 0; /* Number of nonzero elements in the row */

                if (Grid_get_c_status(grid, i, j, k) != BOUNDARY) { /*setup the linear system for all the points except for the half cell added to the end */

                    /* get the metric coefficients for non-uniform grid formulation */
                    /* For more detail, check the manual */
                    metric_px = metric_xc[i];     /* (dx/deta)^-1 at (i,j,k) */
                    metric_py = metric_yc[j];     /* (dy/dxi)^-1 at (i,j,k)      */
                    metric_pz = metric_zc[k];     /* (dz/dphi)^-1 at (i,j,k)     */

                    metric_ex = metric_xu[i+1];   /* (dx/deta)^-1 at (i+1/2,j,k)   */
                    metric_wx = metric_xu[i];     /* (dx/deta)^-1 at (i-1/2,j,k) */

                    metric_ny = metric_yv[j+1];   /* (dy/dxi)^-1 at (i,j+1/2,k)  */
                    metric_sy = metric_yv[j];     /* (dy/dxi)^-1 at (i,j-1/2,k)  */

                    metric_fz = metric_zw[k+1];   /* (dz/dphi)^-1 at (i,j,k+1/2)  */
                    metric_bz = metric_zw[k];     /* (dz/dphi)^-1 at (i,j,k-1/2)  */

                    ap = -Cx * metric_px * ( metric_ex + metric_wx )
                            -Cy * metric_py * ( metric_ny + metric_sy )
                            -Cz * metric_pz * ( metric_bz + metric_fz );

                    aw = Cx * metric_px * metric_wx ;
                    ae = Cx * metric_px * metric_ex ;
                    an = Cy * metric_py * metric_ny ;
                    as = Cy * metric_py * metric_sy ;
                    ab = Cz * metric_pz * metric_bz ;
                    af = Cz * metric_pz * metric_fz ;

                    /* Now, we are setting the values for the current row of the LHS matrix [A] */

                    /* Check for Back neighbor */
                    if (k > 0) { /*Back neighbour is fluid, insert ab into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ab;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k-1;
                        nnz++;

                    } else { /* apply Neumann B.C. on the wall */

                        ap += ab;
                    }

                    /* Check for south neighbor */
                    if (j > 0) { /*south neighbour is fluid, insert as into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = as;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j-1;
                        cols[nnz].k = k;
                        nnz++;

                    } else { /* apply Neumann B.C. on the wall */

                        ap += as;
                    }


                    /* Check for west neighbor */
                    if (i > 0) { /*west neighbour is fluid, insert aw into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = aw;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i-1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                    } else { /* apply Neumann B.C. on the wall */

                        ap += aw;
                    }

                    /* East neighbor */
                    if (i < NI-1) { /*East neighbour is fluid, insert ae into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ae;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i+1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                    } else { /* apply Neumann B.C. on the wall */

                        /* Compatibility condition, Assume the east-node to the top-right-front corner has p=0 */
                        if (! ( (i == NI-1) && (j == NJ-1) && (k == NK-1) ) ) {

                            ap += ae;
                        } /* if */
                    }

                    /* Check for north neighbor */
                    if (j < NJ-1) { /*north neighbour is fluid, insert an into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = an;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j+1;
                        cols[nnz].k = k;
                        nnz++;
                    } else { /* apply Neumann B.C. on the wall */

                        ap += an;
                    }


                    if (k < NK-1) { /*Front neighbour is fluid, insert af into matrix*/
                        /* Check for front neighbor */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = af;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k+1;
                        nnz++;

                    } else { /* apply Neumann B.C. on the wall */

                        ap += af;
                    }


                    /* current node */
                    /* the value of nonzero element in the matrix */
                    values[nnz] = ap;

                    /* 3D index of nonzero element in current row */
                    cols[nnz].i = i;
                    cols[nnz].j = j;
                    cols[nnz].k = k;
                    nnz++;

                } else { /* current node is a solid node */


                    values[nnz] = ap_solid; /* the value is not important */
                    cols[nnz].i = i;
                    cols[nnz].j = j;
                    cols[nnz].k = k;
                    nnz++;

                } /* else */

                /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                ierr = MatSetValuesStencil(p->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);

            } /* for: i*/
        } /*for: j */
    }/* for :k*/

    ierr = MatAssemblyBegin(p->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
    ierr = MatAssemblyEnd(p->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
}
/**************************************************************************************************************/

