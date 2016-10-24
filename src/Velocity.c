#include "definitions.h"
#include "Boundary.h"
#include "DataTypes.h"
#include "Solver.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Conc.h"
#include "gvg.h"
#include "Immersed.h"
#include "Output.h"
#include "Velocity.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* This function allocates memory for u-velocity structure */
Velocity *Velocity_create(MAC_grid *grid, Parameters *params, char which_velocity) {

    Velocity *new_vel;
    int NX, NY, NZ, NT;
    int NX_cell, NY_cell, NZ_cell;
    int ierr;

    /* Number of real physical cells */
    NX_cell = params->NX;
    NY_cell = params->NY;
    NZ_cell = params->NZ;

    /* Add one more in each direction. This would help uniforming the for-loops and data distribution layout */
    NX = NX_cell + 1;
    NY = NY_cell + 1;
    NZ = NZ_cell + 1;
    NT = NX*NY*NZ;

    new_vel = (Velocity *)malloc(sizeof(Velocity));
    Memory_check_allocation(new_vel);

    new_vel->component  = which_velocity; /* 'u', 'v' or 'w' */
    /* Inlfux on current processor */
    new_vel->G_influx   = 0.0;

    /* Total influx (same on all processors) */
    new_vel->W_influx   = 0.0;

    /*
        - data    : value of the velocity
        - data_old: value of the velocity of the old time step
        - data_bc : value of the velocity at the cell_center
*/
    /* Create 3d Distributed array for the convective terms. This is due to the fact that we need more than
 one ghost nodes for calculating convective terms */
    /***************************************/

    /* Now, create the vectors that hold the local and global data for all the 3D properties */
    /* First, data at u or v or w grid */
    ierr = DACreateGlobalVector(grid->DA_3D, &new_vel->G_data); PETScErrAct(ierr);
    ierr = DACreateLocalVector(grid->DA_3D, &new_vel->L_data); PETScErrAct(ierr);

    /* Duplicate other vectors from existing Global and Local vector */
    ierr = VecDuplicate(new_vel->G_data, &new_vel->G_data_old); PETScErrAct(ierr);

    /* Cell centered velocity is mainly used for concentration convective terms VdcdX */
    ierr = VecDuplicate(new_vel->G_data, &new_vel->G_data_bc); PETScErrAct(ierr);
    ierr = VecDuplicate(new_vel->L_data, &new_vel->L_data_bc); PETScErrAct(ierr);

    /* Convective terms: uddx+vddy+wddz */
    ierr = VecDuplicate(new_vel->G_data, &new_vel->G_conv); PETScErrAct(ierr);


    /* Create a local vector based on a box stencil. This array might be useful in computing shear stress */
    ierr = DACreateLocalVector(grid->DA_3D_BOX, &new_vel->L_data_box); PETScErrAct(ierr);

    switch (which_velocity) {

    case 'u':

        /* Index of the last node excluding the half-cell added to the end */
        new_vel->NI = NX;
        new_vel->NJ = NY-1;
        new_vel->NK = NZ-1;


        /* Pass the pointers to the allocated memory for the transposed velocities at the current grid. Note that only
   two components are required. u -> Transposed: v, w  */

        new_vel->G_v_transposed = &grid->G_dummy1;
        new_vel->G_w_transposed = &grid->G_dummy2;
        new_vel->L_v_transposed = &grid->L_dummy1;
        new_vel->L_w_transposed = &grid->L_dummy2;
        break;

    case 'v':

        /* Index of the last node excluding the half-cell added to the end */
        new_vel->NI = NX-1;
        new_vel->NJ = NY;
        new_vel->NK = NZ-1;

        /* Pass the pointers to the allocated memory for the transposed velocities at the current grid. Note that only
   two components are required. v -> Transposed: u, w  */
        new_vel->G_u_transposed = &grid->G_dummy1;
        new_vel->G_w_transposed = &grid->G_dummy2;
        new_vel->L_u_transposed = &grid->L_dummy1;
        new_vel->L_w_transposed = &grid->L_dummy2;

        break;

    case 'w':

        /* Index of the last node excluding the half-cell added to the end */
        new_vel->NI = NX-1;
        new_vel->NJ = NY-1;
        new_vel->NK = NZ;

        /* Pass the pointers to the allocated memory for the transposed velocities at the current grid. Note that only
   two components are required. w -> Transposed: u, v  */
        new_vel->G_u_transposed = &grid->G_dummy1;
        new_vel->G_v_transposed = &grid->G_dummy2;
        new_vel->L_u_transposed = &grid->L_dummy1;
        new_vel->L_v_transposed = &grid->L_dummy2;

        break;

    } /* switch */

    new_vel->G_outflow     = NULL;
    new_vel->G_outflow_old = NULL;
    new_vel->G_inflow      = NULL;

    if (params->outflow) {

        /* create the array on the processors which contain the domain at the end of the regular box */
        if (grid->G_Ie == grid->NX) {

            /* Note that we create the whole array on the processor but the processor operates
   only on the part which contains the G_data, i.e. (Js to Je)*(Ks to Ke). This is done
   since a 2D array does not need much memory (Also, this way makes life a lot easier :) ) */
            new_vel->G_outflow     = Memory_allocate_2D_double(NY, NZ);
            new_vel->G_outflow_old = Memory_allocate_2D_double(NY, NZ);
        } /* if */
    } /* if */

    /* create inflow only for u-velocity */
    if ( (params->inflow) && (which_velocity == 'u') ){

        /* create the array on the processors which contain the domain at the beginning of the regular box */
        if (grid->G_Is == 0) {

            /* Note that we create the whole array on the processor but the processor operates
   only on the part which contains the G_data, i.e. (Js to Je)*(Ks to Ke). This is done
   since a 2D array does not need much memory (Also, this way makes life a lot easier :) ) */

            new_vel->G_inflow = Memory_allocate_2D_double(NY, NZ);
        } /* if */
    } /* if */


    /* Allocate memory for the shear stress on the bottom. tau = 1/Re sqrt(dudy^2 + dwdy^2) */
    if ( (which_velocity == 'u') && (params->shear_stress_output) ){

        new_vel->G_shear_stress_bottom  = Memory_allocate_2D_double(NX, NZ);
        new_vel->W_shear_stress_bottom  = Memory_allocate_2D_double(NX, NZ);

        new_vel->G_u_shear = Memory_allocate_2D_double(NX, NZ);
        new_vel->G_v_shear = Memory_allocate_2D_double(NX, NZ);
        new_vel->G_w_shear = Memory_allocate_2D_double(NX, NZ);

        new_vel->W_u_shear = Memory_allocate_2D_double(NX, NZ);
        new_vel->W_v_shear = Memory_allocate_2D_double(NX, NZ);
        new_vel->W_w_shear = Memory_allocate_2D_double(NX, NZ);

    } else {

        new_vel->G_shear_stress_bottom = NULL;
        new_vel->W_shear_stress_bottom = NULL;

        new_vel->G_u_shear = NULL;
        new_vel->G_v_shear = NULL;
        new_vel->G_w_shear = NULL;

        new_vel->W_u_shear = NULL;
        new_vel->W_v_shear = NULL;
        new_vel->W_w_shear = NULL;

    } /* else */

    /* Index of left-bottom-back corner on current processor for the DA data layout */
    new_vel->G_Is = grid->G_Is;
    new_vel->G_Js = grid->G_Js;
    new_vel->G_Ks = grid->G_Ks;

    /* Index of right-top-front corner on current processor for the DA data layout */
    new_vel->G_Ie = grid->G_Ie;
    new_vel->G_Je = grid->G_Je;
    new_vel->G_Ke = grid->G_Ke;

    /* Index of left-bottom-back corner on current processor for the DA data layout including ghost nodes	*/
    new_vel->L_Is = grid->L_Is;
    new_vel->L_Js = grid->L_Js;
    new_vel->L_Ks = grid->L_Ks;

    /* Index of right-top-front corner on current processor for the DA data layout including ghost nodes*/
    new_vel->L_Ie = grid->L_Ie;
    new_vel->L_Je = grid->L_Je;
    new_vel->L_Ke = grid->L_Ke;

    /* Global number of vel cells */
    new_vel->NX   = NX;
    new_vel->NY   = NY;
    new_vel->NZ   = NZ;
    new_vel->NT   = NT;

    /* Now, for the linear system */
    /* Create a 1D vector based on the Distributed Array setup */

    /* {rhs}={b} vector, */
    ierr = DACreateGlobalVector(grid->DA_3D_Lsys, &new_vel->G_b); PETScErrAct(ierr);

    /* Create the sparse lhs matrix (lsys setup) */
    ierr = DAGetMatrix(grid->DA_3D_Lsys, MATMPIAIJ, &new_vel->A); PETScErrAct(ierr);

    /* vector used to update the diagonal part of matrix A */
    /* Get the pointer to it so we can save memory */
    new_vel->lsys_diagonal = &grid->lsys_diagonal;

    //Velocity_nonzero_initialize(new_vel, grid, params);


    /* Pointer to store anything temporarily */
    new_vel->G_anything = &grid->G_dummy1;
    return new_vel;
}
/***************************************************************************************************/


/* This function releases the allocated memory for velocity structure. */
void Velocity_destroy(Velocity *vel) {

    int ierr;
    int NZ = vel->NZ;

    /* Velocity data */
    ierr = VecDestroy(vel->G_data); PETScErrAct(ierr);
    ierr = VecDestroy(vel->L_data); PETScErrAct(ierr);

    ierr = VecDestroy(vel->G_data_old); PETScErrAct(ierr);

    ierr = VecDestroy(vel->G_data_bc); PETScErrAct(ierr);
    ierr = VecDestroy(vel->L_data_bc); PETScErrAct(ierr);

    /* Convective terms */
    ierr = VecDestroy(vel->G_conv); PETScErrAct(ierr);

    /* Local data using a BOX_STENCIL DA */
    ierr = VecDestroy(vel->L_data_box); PETScErrAct(ierr);

    /* release the memory for outflow and inflow (if allocated) */
    if (vel->G_outflow != NULL) {

        Memory_free_2D_array(NZ, (void **)vel->G_outflow);
        Memory_free_2D_array(NZ, (void **)vel->G_outflow_old);
    }
    if (vel->G_inflow != NULL) {

        Memory_free_2D_array(NZ, (void **)vel->G_inflow);
    }

    if (vel->G_u_shear != NULL) {

        Memory_free_2D_array(NZ, (void **)vel->G_u_shear);
        Memory_free_2D_array(NZ, (void **)vel->G_v_shear);
        Memory_free_2D_array(NZ, (void **)vel->G_w_shear);

        Memory_free_2D_array(NZ, (void **)vel->W_u_shear);
        Memory_free_2D_array(NZ, (void **)vel->W_v_shear);
        Memory_free_2D_array(NZ, (void **)vel->W_w_shear);
    }
    if (vel->G_shear_stress_bottom != NULL) {

        Memory_free_2D_array(NZ, (void **)vel->G_shear_stress_bottom);
        Memory_free_2D_array(NZ, (void **)vel->W_shear_stress_bottom);
    }
    /* Free linear system variables */
    ierr = VecDestroy(vel->G_b); PETScErrAct(ierr);

    ierr = MatDestroy(vel->A); PETScErrAct(ierr);
    ierr = KSPDestroy(vel->solver); PETScErrAct(ierr);

    free(vel);

}
/***************************************************************************************************/

/* This function stores the old value for velocity, i.e. copying "data" into "data_old". */
void Velocity_store_old_data(Velocity *vel) {

    /* copy data into data_old */
    int ierr = VecCopy(vel->G_data, vel->G_data_old); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function solves the linear set of equations for velocity using the Petsc KSP solver */
int Velocity_solve(Velocity *vel) {

    int iters;
    int ierr;

    /* Call Petsc KSP-solver routine to solve linear system. */
    ierr = KSPSetOperators(vel->solver, vel->A, vel->A, SAME_NONZERO_PATTERN); PETScErrAct(ierr);

    /* The solution is stored in data array */
    ierr = KSPSolve(vel->solver, vel->G_b, vel->G_data); PETScErrAct(ierr);

    char *name = &vel->component;
    (void)Solver_is_converged(vel->solver, name);

    ierr = KSPGetIterationNumber(vel->solver, &iters); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, " ** %c-Velocity soln took %d iters.\n", vel->component, iters); PETScErrAct(ierr);

    return iters;
}
/***************************************************************************************************/

/* This function creates the vectors and matrix [A] for the velocity linear system, also it sets up the KSP
solver for the first time. */
void Velocity_setup_lsys_accounting_geometry(Velocity *vel, MAC_grid *grid, Parameters *params) {

    int ierr;

    /* Now, generate the LHS matrix for the (u, v, w) momentum equation using the 2nd order finite difference method */
    switch (vel->component) {

    case 'u': Velocity_u_form_LHS_matrix(vel, grid, params); break;

    case 'v': Velocity_v_form_LHS_matrix(vel, grid, params); break;

    case 'w': Velocity_w_form_LHS_matrix(vel, grid, params); break;

    }

    /* Now, using the given LHS matrix {and some known parameters like grid}, appropriate KSP solver (and
preconditioner) is generated using Petsc package */
#ifdef MEMORY_PROFILING
    PetscMallocGetCurrentUsage(&mem1);
#endif
    vel->solver = Solver_get_velocity_solver(vel->A, params);
    ierr = KSPSetOperators(vel->solver, vel->A, vel->A, SAME_PRECONDITIONER); PETScErrAct(ierr);
    ierr = KSPSetUp(vel->solver); PETScErrAct(ierr);
#ifdef MEMORY_PROFILING
    PetscMallocGetCurrentUsage(&mem2);
    switch (vel->component) {
    case 'u':  u_solver_mem += (int)(mem2 - mem1); u_mem += (int)(mem2 - mem1); break;
    case 'v':  v_solver_mem += (int)(mem2 - mem1); v_mem += (int)(mem2 - mem1); break;
    case 'w':  w_solver_mem += (int)(mem2 - mem1); w_mem += (int)(mem2 - mem1); break;
    } /* switch */
#endif

    if (params->output_lsys_matrix) {

        char filename[FILENAME_MAX];
        sprintf(filename, "Resume_%c_matrix.bin.gz", vel->component);
        Output_petsc_mat(&(vel->A), filename, BINARY);
        PetscPrintf(PCW, "Velocity.c/ %c-matrix has been written to %s successfully.\n", vel->component, filename);

    } /* if */
}
/***************************************************************************************************/

/* This function updates the diagonal part of matrix [A] for velocities by the new value of "dt" */
/* It will only update the nodes that are NOT immersed nodes */
void Velocity_modify_diagonal(Velocity *vel, MAC_grid *grid, Parameters *params, double dt, double dt_old) {

    int ierr;
    double mod;
    double mod_insert;
    double ***diagonal;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;


    if ( (params->viscosity_type == CONSTANT_VISCOSITY) && (fabs(dt_old) > 1.0e-7) ){

        mod = 1.0/dt - 1.0/dt_old;
    }
    else {

        mod = 1.0/dt;
    }

    /* To save time, add "mod" to the diagonal only if dt != dt_old */
    if (fabs(mod) > 1.0e-10) {

        /* Adds constant 'mod' to the diagonal part of matrix A*/
        /* Start index of bottom-left-back corner on current processor */
        Is = grid->G_Is;
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* End index of top-right-front corner on current processor */
        Ie = grid->G_Ie;
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        switch (vel->component) {

        case 'u': Grid_get_q_status = &Grid_get_u_status; break;
        case 'v': Grid_get_q_status = &Grid_get_v_status; break;
        case 'w': Grid_get_q_status = &Grid_get_w_status; break;
        } /* switch */

        /* Get the diagonal vector array */
        ierr = DAVecGetArray(grid->DA_3D, (*vel->lsys_diagonal), (void ***)&diagonal); PETScErrAct(ierr);

        /* Go through the entire domain and only insert the delta_t term for the nodes which are not immersed nodes */
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {

                    if (Grid_get_q_status(grid, i, j, k) == IMMERSED) {

                        mod_insert = 0.0;
                    } else {

                        mod_insert = mod;
                    } /* else */
                    diagonal[k][j][i] = mod_insert;
                } /* for i */
            } /* for j */
        } /* for k */

        ierr = DAVecRestoreArray(grid->DA_3D, (*vel->lsys_diagonal), (void ***)&diagonal); PETScErrAct(ierr);
        /* Now, add the diagonal part of the lsys matrix to update 1/dt */
        /* Maybe this could be done in a mroe efficietn way which requires less memory */
        ierr = MatDiagonalSet(vel->A, (*vel->lsys_diagonal), ADD_VALUES); PETScErrAct(ierr);
    } /* if */


}
/***************************************************************************************************/


/* This function creates the LHS matrix for the linear system of u-mom equations in the general case, i.e.
non-uniform mesh {and variable viscosity} */
void Velocity_u_form_LHS_matrix(Velocity *u, MAC_grid *grid, Parameters *params) {

    PetscReal ae, aw, ap, as, an, af, ab;
    PetscReal ap_solid;
    short int status;
    int ierr;
    int i, j, k;
    int NI, NJ, NK;
    double Re;
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
    MatStencil cols[STENCIL], row; /* 7 point stencil */
    double values[STENCIL];
    int nnz;
    ImmersedNode *ib_node=NULL;
    int ib_index;
    /* IMMERSED node setup */
    int g;
    int ii, jj, kk;
    int local_index, row_global_index;
    int im_global_indices[9];
    double im_values[9];
    int Is_g, Js_g, Ks_g;
    int nx_g, ny_g, nz_g, nt;
    int *global_indices;


    NI  = u->NI;
    NJ  = u->NJ;
    NK  = u->NK;

    Re  = params->Re;
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
    Cx  = -1.0/(Re * dEta * dEta);
    Cy  = -1.0/(Re * dXi * dXi);
    Cz  = -1.0/(Re * dPhi * dPhi);

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


    /* For the IMMERSED nodes since they do not use a 7-point finite difference stencil, to insert the nonzeron values in the matrix,
 we need to pay careful attention.
The function MatSetValueStencil() only inserts the values in the matrix if the current row of the matrix is owned within the local part (including the ghost nodes) on the current processor. 
So, we are NOT using this function to insert that values fot the IMMERSED nodes since the nonzero might be lying on the neighboring processor. 
MatSetValues() instead can be used to insert values on any processor. But we need to know the global_index of the node which is not easy to obtain using DA. 
To overcome that, we create a DA_BOX_STENCIL with 3 (or even more) layers of ghost nodes to make sure the immersed nodes' are bounded within the ghost node range. Now, we use DAGetGlobalIndices() to obtain the global indices correspodning to all the nodes (including the ghost nodes) on current processor. 
This way, we make sure that we can always get the global index of any nodes on the current processor. 
Note that, creating the matrix using 3-layers of ghost nodes and DA_BOX_STENCIL is also an easier solution, but that would make the memory requirement enormous, e.g. 400-500% more!
*/	
    ierr = DAGetGhostCorners(grid->DA_3D_BOX, &Is_g, &Js_g, &Ks_g, &nx_g, &ny_g, &nz_g); PETScErrAct(ierr);
    ierr = DAGetGlobalIndices(grid->DA_3D_BOX, &nt, &global_indices); PETScErrAct(ierr);

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                /* The 3D index of current row in the Matrix */
                row.i = i;
                row.j = j;
                row.k = k;
                nnz   = 0; /* Number of nonzero elements in the row */

                /* status of the current node */
                status = Grid_get_u_status(grid, i, j, k);

#ifdef SOLVE_INSIDE_SOLID
                if( ( status == FLUID) || (status == BUFFER_SOLID)) {
#else
                if( ( status == FLUID) ) {
#endif

                    /* get the metric coefficients for non-uniform grid formulation */
                    /* For more detail, check the manual */
                    metric_px = metric_xu[i];     /* (dx/deta)^-1 at (i+1/2,j,k) */
                    metric_py = metric_yc[j];     /* (dy/dxi)^-1 at (i,j,k)      */
                    metric_pz = metric_zc[k];     /* (dz/dphi)^-1 at (i,j,k)     */
                    metric_ex = metric_xc[i];     /* (dx/deta)^-1 at (i+1,j,k)   */
                    metric_wx = metric_xc[i-1];   /* (dx/deta)^-1 at (i-1/2,j,k) */
                    metric_ny = metric_yv[j+1];   /* (dy/dxi)^-1 at (i,j+1/2,k)  */
                    metric_sy = metric_yv[j];     /* (dy/dxi)^-1 at (i,j-1/2,k)  */
                    metric_bz = metric_zw[k];     /* (dz/dphi)^-1 at (i,j,k-1/2)  */
                    metric_fz = metric_zw[k+1];   /* (dz/dphi)^-1 at (i,j,k+1/2)  */

                    ap = -Cx * metric_px * ( metric_ex + metric_wx )
                            -Cy * metric_py * ( metric_ny + metric_sy )
                            -Cz * metric_pz * ( metric_fz + metric_bz );

                    aw = Cx * metric_px * metric_wx ;
                    ae = Cx * metric_px * metric_ex ;
                    an = Cy * metric_py * metric_ny ;
                    as = Cy * metric_py * metric_sy ;
                    ab = Cz * metric_pz * metric_bz ;
                    af = Cz * metric_pz * metric_fz ;

                    /* Now, we are setting the values for the current row of the LHS matrix [A] */

                    /* Check for Back neighbor */
                    if (k == 0) {
                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef BACK_WALL_VELOCITY_NOSLIP
                        ap -= ab;
#endif
#ifdef BACK_WALL_VELOCITY_FREESLIP
                        ap += ab;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ab;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k-1;
                        nnz++;
                    }

                    /* Check for South neighbor */
                    if (j == 0) {
                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
                        ap -= as;
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
                        ap += as;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = as;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j-1;
                        cols[nnz].k = k;
                        nnz++;
                    } /* else */

                    /* Check for West neighbor */
                    if (Grid_get_u_status(grid, i-1, j, k) != BOUNDARY) { /* Do not include the box boundaries node at i=0 */
                        /* the value of nonzero element in the matrix */
                        values[nnz] = aw;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i-1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;
                    } /* if */

                    /* East neighbor */
                    if (Grid_get_u_status(grid, i+1, j, k) != BOUNDARY) { /* Do not include the box boundary node at i=NI */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ae;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i+1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                    }  /* if east */

                    /* Check for north neighbor */
                    if (j == NJ-1 ) {
                        /* Define the boundary conditions for the top-wall in 'Boundary.h' file */
#ifndef DRIVEN_CAVITY
#ifdef TOP_WALL_VELOCITY_NOSLIP
                        ap -= an;
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
                        ap += an;
#endif
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = an;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j+1;
                        cols[nnz].k = k;
                        nnz++;
                    } /* else */


                    /* Check for Front neighbor */
                    if (k == NK-1 ) {
                        /* Define the boundary conditions for the Front-wall in 'Boundary.h' file */
#ifdef FRONT_WALL_VELOCITY_NOSLIP
                        ap -= af;
#endif
#ifdef FRONT_WALL_VELOCITY_FREESLIP
                        ap += af;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = af;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k+1;
                        nnz++;
                    } /* else */

                    /* current node */
                    /* the value of nonzero element in the matrix */
                    values[nnz] = ap;

                    /* 3D index of nonzero element in current row */
                    cols[nnz].i = i;
                    cols[nnz].j = j;
                    cols[nnz].k = k;
                    nnz++;

                    /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                    ierr = MatSetValuesStencil(u->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);

                } else { /* current node is Immersed node or the extra-half cell or Boundary node */

                    if (status == IMMERSED) {

                        local_index = (k - Ks_g) * ny_g * nx_g + (j - Js_g) * nx_g + (i - Is_g);
                        /* Get the global index (row number in the Matrix) for current node at (i,j,k) */
                        row_global_index = global_indices[local_index];
                        im_global_indices[nnz] = global_indices[local_index];
                        im_values[nnz] = 1.0;
                        nnz++;

                        /* Now, go through the fluid nodes, find the global index in the matrix, copy the value and use MatSetValues() instead of MatSetValuesStencil. This is done b/c some of the fluid node might be on the neighboring processor and out of the range of the ghost layers. */
                        ib_index = Immersed_get_ib_global_index(grid->u_immersed, i, j, k);
                        ib_node  = Immersed_get_ib_node(grid->u_immersed, ib_index);

                        Immersed_is_valid(ib_node) ;
                        for (g=0; g<ib_node->n_fluid; g++) {

                            ii  = ib_node->fluid_index[g].x_index;
                            jj  = ib_node->fluid_index[g].y_index;
                            kk  = ib_node->fluid_index[g].z_index;
                            if ( (ii == i) && (jj == j) && (kk == k) ) { /* current immersed node is a function of itself value */

                                im_values[0] += ib_node->fluid_coef[g];
                            } else {


                                /* Only add the interpolation point if the coefficient is greater than zero. This will help with reducing the nonzero stencil */
                                if (fabs(ib_node->fluid_coef[g]) >= LSYS_DROP_TOLERANCE) {

                                    local_index = (kk - Ks_g) * ny_g * nx_g + (jj - Js_g) * nx_g + (ii - Is_g);
                                    im_values[nnz] = ib_node->fluid_coef[g];
                                    im_global_indices[nnz] = global_indices[local_index];
                                    nnz++;

                                }
                            }

                        } /* for */

                        /* All interpolation coefficients are scaled with a constant factor. Better convergence */
                        int s;
                        for (s=0; s<nnz; s++) im_values[s] *= ap_solid;

                        ierr = MatSetValues(u->A, 1, &row_global_index, nnz, im_global_indices, im_values, INSERT_VALUES); PETScErrAct(ierr);
                    } else {

                        values[nnz] = ap_solid; /* the value is not important */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                        /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                        ierr = MatSetValuesStencil(u->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);
                    } /* else BOUNDARY node */

                } /* else */

            } /* for: i*/
        } /*for: j */
    }/* for :k*/

    ierr = MatAssemblyBegin(u->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
    ierr = MatAssemblyEnd(u->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function creates the LHS matrix for the linear system of v-mom equations in the general case. */
void Velocity_v_form_LHS_matrix(Velocity *v, MAC_grid *grid, Parameters *params) {

    PetscReal ae, aw, ap, as, an, af, ab;
    PetscReal ap_solid;
    short int status;
    int ierr;
    int i, j, k;
    int NI, NJ, NK;
    double Re;
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
    MatStencil cols[2*STENCIL], row;
    double values[2*STENCIL];
    int nnz;
    ImmersedNode *ib_node;
    int ib_index;
    /* IMMERSED node setup */
    int g;
    int ii, jj, kk;
    int local_index, row_global_index;
    int im_global_indices[9];
    double im_values[9];
    int Is_g, Js_g, Ks_g;
    int nx_g, ny_g, nz_g, nt;
    int *global_indices;

    NI  = v->NI;
    NJ  = v->NJ;

    NK  = v->NK;

    Re  = params->Re;

    /*Assemble matrix*/
    /*
  a? =[A]
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
    Cx  = -1.0/(Re * dEta * dEta);
    Cy  = -1.0/(Re * dXi * dXi);
    Cz  = -1.0/(Re * dPhi * dPhi);

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

    /* For the IMMERSED nodes since they do not use a 7-point finite difference stencil, to insert the nonzeron values in the matrix,
 we need to pay careful attention.
The function MatSetValueStencil() only inserts the values in the matrix if the current row of the matrix is owned within the local part (including the ghost nodes) on the current processor. 
So, we are NOT using this function to insert that values fot the IMMERSED nodes since the nonzero might be lying on the neighboring processor. 
MatSetValues() instead can be used to insert values on any processor. But we need to know the global_index of the node which is not easy to obtain using DA. 
To overcome that, we create a DA_BOX_STENCIL with 3 (or even more) layers of ghost nodes to make sure the immersed nodes' are bounded within the ghost node range. Now, we use DAGetGlobalIndices() to obtain the global indices correspodning to all the nodes (including the ghost nodes) on current processor. 
This way, we make sure that we can always get the global index of any nodes on the current processor. 
Note that, creating the matrix using 3-layers of ghost nodes and DA_BOX_STENCIL is also an easier solution, but that would make the memory requirement enormous, e.g. 400-500% more!
*/	
    ierr = DAGetGhostCorners(grid->DA_3D_BOX, &Is_g, &Js_g, &Ks_g, &nx_g, &ny_g, &nz_g); PETScErrAct(ierr);
    ierr = DAGetGlobalIndices(grid->DA_3D_BOX, &nt, &global_indices); PETScErrAct(ierr);

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                /* The 3D index of current row in the Matrix */
                row.i = i;
                row.j = j;
                row.k = k;
                nnz   = 0; /* Number of nonzero elements in the row */

                /* status of the current node */
                status = Grid_get_v_status(grid, i, j, k);
#ifdef SOLVE_INSIDE_SOLID
                if( ( status == FLUID) || (status == BUFFER_SOLID)) {
#else
                if( ( status == FLUID) ) {
#endif

                    /* get the metric coefficients for non-uniform grid formulation */
                    /* For more detail, check the manual */
                    metric_px = metric_xc[i];     /* (dx/deta)^-1 at (i,j,k) */
                    metric_py = metric_yv[j];     /* (dy/dxi)^-1 at (i,j+1/2,k)      */
                    metric_pz = metric_zc[k];     /* (dz/dphi)^-1 at (i,j,k)     */
                    metric_ex = metric_xu[i+1];   /* (dx/deta)^-1 at (i+1/2,j,k)   */
                    metric_wx = metric_xu[i];     /* (dx/deta)^-1 at (i-1/2,j,k) */
                    metric_ny = metric_yc[j];     /* (dy/dxi)^-1 at (i,j+1,k)  */
                    metric_sy = metric_yc[j-1];   /* (dy/dxi)^-1 at (i,j,k)  */
                    metric_bz = metric_zw[k];     /* (dz/dphi)^-1 at (i,j,k-1/2)  */
                    metric_fz = metric_zw[k+1];   /* (dz/dphi)^-1 at (i,j,k+1/2)  */

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
                    if (k == 0) {
                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef BACK_WALL_VELOCITY_NOSLIP
                        ap -= ab;
#endif



#ifdef BACK_WALL_VELOCITY_FREESLIP
                        ap += ab;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ab;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k-1;
                        nnz++;
                    }

                    /* Check for South neighbor */
                    if (Grid_get_v_status(grid, i, j-1, k) != BOUNDARY) { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = as;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j-1;
                        cols[nnz].k = k;
                        nnz++;
                    }

                    /* Check for West neighbor */
                    if (i == 0)  { /* Do not include the box boundaries node at i=0 */

                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef LEFT_WALL_VELOCITY_NOSLIP
                        ap -= aw;
#endif
#ifdef LEFT_WALL_VELOCITY_FREESLIP
                        ap += aw;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = aw;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i-1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;
                    }

                    /* Check for East neighbor */
                    if ( (i == NI-1) && (!params->outflow) )  { /* Do not include the box boundaries node at i=NI */

                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef RIGHT_WALL_VELOCITY_NOSLIP
                        ap -= ae;
#endif
#ifdef RIGHT_WALL_VELOCITY_FREESLIP
                        ap += ae;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        if ( i < (NI - 1) ) {

                            /* the value of nonzero element in the matrix */
                            values[nnz] = ae;

                            /* 3D index of nonzero element in current row */
                            cols[nnz].i = i+1;
                            cols[nnz].j = j;
                            cols[nnz].k = k;
                            nnz++;
                        } /* if i */
                    }

                    /* Check for north neighbor */
                    if (Grid_get_v_status(grid, i, j+1, k) != BOUNDARY ) {

                        /* the value of nonzero element in the matrix */
                        values[nnz] = an;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j+1;
                        cols[nnz].k = k;
                        nnz++;

                    } /* if */


                    /* Check for Front neighbor */
                    if (k == NK-1 ) {
                        /* Define the boundary conditions for the Front-wall in 'Boundary.h' file */
#ifdef FRONT_WALL_VELOCITY_NOSLIP
                        ap -= af;
#endif
#ifdef FRONT_WALL_VELOCITY_FREESLIP
                        ap += af;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = af;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k+1;
                        nnz++;
                    } /* else */

                    /* current node */
                    /* the value of nonzero element in the matrix */
                    values[nnz] = ap;

                    /* 3D index of nonzero element in current row */
                    cols[nnz].i = i;
                    cols[nnz].j = j;
                    cols[nnz].k = k;
                    nnz++;

                    /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                    ierr = MatSetValuesStencil(v->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);

                    //printf("Velocity.c/ (i,j,k)=(%d,%d,%d) ap:%f aw:%f ae:%f as:%f an:%f ab:%f af:%f\n", i, j, k, ap, aw, ae, as, an, ab, af);

                } else { /* current node is Immersed node or the extra-half cell or Boundary node */

                    if (status == IMMERSED) {

                        local_index = (k - Ks_g) * ny_g * nx_g + (j - Js_g) * nx_g + (i - Is_g);
                        /* Get the global index (row number in the Matrix) for current node at (i,j,k) */
                        row_global_index = global_indices[local_index];
                        im_global_indices[nnz] = global_indices[local_index];
                        im_values[nnz] = 1.0;
                        nnz++;

                        /* Now, go through the fluid nodes, find the global index in the matrix, copy the value and use MatSetValues() instead of MatSetValuesStencil. This is done b/c some of the fluid node might be on the neighboring processor and out of the range of the ghost layers. */
                        ib_index = Immersed_get_ib_global_index(grid->v_immersed, i, j, k);
                        ib_node  = Immersed_get_ib_node(grid->v_immersed, ib_index);
                        for (g=0; g<ib_node->n_fluid; g++) {

                            ii  = ib_node->fluid_index[g].x_index;
                            jj  = ib_node->fluid_index[g].y_index;
                            kk  = ib_node->fluid_index[g].z_index;

                            if ( (ii == i) && (jj == j) && (kk == k) ) { /* current immersed node is a function of itself value */

                                im_values[0] += ib_node->fluid_coef[g];
                            } else {


                                /* Only add the interpolation point if the coefficient is greater than zero. This will help with reducing the nonzero stencil */
                                if (fabs(ib_node->fluid_coef[g]) >= LSYS_DROP_TOLERANCE) {

                                    local_index = (kk - Ks_g) * ny_g * nx_g + (jj - Js_g) * nx_g + (ii - Is_g);
                                    im_values[nnz] = ib_node->fluid_coef[g];
                                    im_global_indices[nnz] = global_indices[local_index];
                                    nnz++;
                                } /* if */
                            }
                        } /* for */

                        /* All interpolation coefficients are scaled with a constant factor. Better convergence */
                        int s;
                        for (s=0; s<nnz; s++) im_values[s] *= ap_solid;

                        ierr = MatSetValues(v->A, 1, &row_global_index, nnz, im_global_indices, im_values, INSERT_VALUES); PETScErrAct(ierr);
                    } else {

                        values[nnz] = ap_solid; /* the value is not important */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                        /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                        ierr = MatSetValuesStencil(v->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);
                    } /* else BOUNDARY node */

                } /* else */


            } /* for: i*/
        } /*for: j */
    }/* for :k*/

    /* TBC */
    ierr = MatAssemblyBegin(v->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
    ierr = MatAssemblyEnd(v->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
}
/***************************************************************************************************/



/* This function creates the LHS matrix for the linear system of u-mom equations in the general case, i.e.
non-uniform mesha and variable viscosity */
void Velocity_w_form_LHS_matrix(Velocity *w, MAC_grid *grid, Parameters *params) {

    PetscReal ae, aw, ap, as, an, af, ab;
    PetscReal ap_solid;
    short int status;
    int ierr;
    int i, j, k;
    int NI, NJ, NK;
    double Re;
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
    MatStencil cols[2*STENCIL], row;
    double values[2*STENCIL];
    int nnz;
    ImmersedNode *ib_node;
    int ib_index;
    /* IMMERSED node setup */
    int g;
    int ii, jj, kk;
    int local_index, row_global_index;
    int im_global_indices[9];
    double im_values[9];
    int Is_g, Js_g, Ks_g;
    int nx_g, ny_g, nz_g, nt;
    int *global_indices;

    NI  = w->NI;
    NJ  = w->NJ;
    NK  = w->NK;

    Re  = params->Re;

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
    Cx  = -1.0/(Re * dEta * dEta);
    Cy  = -1.0/(Re * dXi * dXi);
    Cz  = -1.0/(Re * dPhi * dPhi);

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

    /* For the IMMERSED nodes since they do not use a 7-point finite difference stencil, to insert the nonzeron values in the matrix,
 we need to pay careful attention.
The function MatSetValueStencil() only inserts the values in the matrix if the current row of the matrix is owned within the local part (including the ghost nodes) on the current processor. 
So, we are NOT using this function to insert that values fot the IMMERSED nodes since the nonzero might be lying on the neighboring processor. 
MatSetValues() instead can be used to insert values on any processor. But we need to know the global_index of the node which is not easy to obtain using DA. 
To overcome that, we create a DA_BOX_STENCIL with 3 (or even more) layers of ghost nodes to make sure the immersed nodes' are bounded within the ghost node range. Now, we use DAGetGlobalIndices() to obtain the global indices correspodning to all the nodes (including the ghost nodes) on current processor. 
This way, we make sure that we can always get the global index of any nodes on the current processor. 
Note that, creating the matrix using 3-layers of ghost nodes and DA_BOX_STENCIL is also an easier solution, but that would make the memory requirement enormous, e.g. 400-500% more!
*/	
    ierr = DAGetGhostCorners(grid->DA_3D_BOX, &Is_g, &Js_g, &Ks_g, &nx_g, &ny_g, &nz_g); PETScErrAct(ierr);
    ierr = DAGetGlobalIndices(grid->DA_3D_BOX, &nt, &global_indices); PETScErrAct(ierr);


    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                /* The 3D index of current row in the Matrix */
                row.i = i;
                row.j = j;
                row.k = k;
                nnz   = 0; /* Number of nonzero elements in the row */

                status = Grid_get_w_status(grid, i, j, k);
#ifdef SOLVE_INSIDE_SOLID
                if( ( status == FLUID) || (status == BUFFER_SOLID)) {
#else
                if( ( status == FLUID) ) {
#endif

                    /* get the metric coefficients for non-uniform grid formulation */
                    /* For more detail, check the manual */
                    metric_px = metric_xc[i];     /* (dx/deta)^-1 at (i,j,k) */
                    metric_py = metric_yc[j];     /* (dy/dxi)^-1 at (i,j,k)      */
                    metric_pz = metric_zw[k];     /* (dz/dphi)^-1 at (i,j,k+1/2)     */
                    metric_ex = metric_xu[i+1];   /* (dx/deta)^-1 at (i+1/2,j,k)   */
                    metric_wx = metric_xu[i];     /* (dx/deta)^-1 at (i-1/2,j,k) */
                    metric_ny = metric_yv[j+1];   /* (dy/dxi)^-1 at (i,j+1/2,k)  */
                    metric_sy = metric_yv[j];     /* (dy/dxi)^-1 at (i,j+1/2,k)  */
                    metric_bz = metric_zc[k-1];   /* (dz/dphi)^-1 at (i,j,k)  */
                    metric_fz = metric_zc[k];     /* (dz/dphi)^-1 at (i,j,k+1)  */

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
                    if (Grid_get_w_status(grid, i, j, k-1) != BOUNDARY) { /*Back neighbour is fluid, insert ab into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ab;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k-1;
                        nnz++;
                    }

                    /* south neighbor */
                    if (j == 0) { /*South neighbour is fluid, insert as into matrix*/

                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
                        ap -= as;
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
                        ap += as;
#endif

                    } else { /* Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = as;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j-1;
                        cols[nnz].k = k;
                        nnz++;
                    }

                    /* Check for west neighbor */
                    if (i == 0) {

                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef LEFT_WALL_VELOCITY_NOSLIP
                        ap -= aw;
#endif
#ifdef LEFT_WALL_VELOCITY_FREESLIP
                        ap += aw;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = aw;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i-1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;
                    }

                    /* East neighbor */
                    if ( (i == NI-1) && (!params->outflow) ) {
                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef RIGHT_WALL_VELOCITY_NOSLIP
                        ap -= ae;
#endif
#ifdef RIGHT_WALL_VELOCITY_FREESLIP
                        ap += ae;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        if ( i < (NI-1) ) {

                            /* the value of nonzero element in the matrix */
                            values[nnz] = ae;

                            /* 3D index of nonzero element in current row */
                            cols[nnz].i = i+1;
                            cols[nnz].j = j;
                            cols[nnz].k = k;
                            nnz++;
                        } /* if i */
                    }

                    /* North neighbor */
                    if (j == NJ-1) {

                        /* Define the boundary conditions for the side-wall in 'Boundary.h' file */
#ifdef TOP_WALL_VELOCITY_NOSLIP
                        ap -= an;
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
                        ap += an;
#endif

                    } else { /* The neighbor is either an Immersed or a Fluid node */

                        /* the value of nonzero element in the matrix */
                        values[nnz] = an;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j+1;
                        cols[nnz].k = k;
                        nnz++;
                    }

                    /* Front neighbor */
                    if (Grid_get_w_status(grid, i, j, k+1) != BOUNDARY) { /*Front neighbour is fluid, insert af into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = af;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k+1;
                        nnz++;

                    }

                    /* current node */
                    /* the value of nonzero element in the matrix */
                    values[nnz] = ap;

                    /* 3D index of nonzero element in current row */
                    cols[nnz].i = i;
                    cols[nnz].j = j;
                    cols[nnz].k = k;
                    nnz++;

                    /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                    ierr = MatSetValuesStencil(w->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);

                } else { /* current node is Immersed node or the extra-half cell or Boundary node */

                    if (status == IMMERSED) {

                        local_index = (k - Ks_g) * ny_g * nx_g + (j - Js_g) * nx_g + (i - Is_g);
                        /* Get the global index (row number in the Matrix) for current node at (i,j,k) */
                        row_global_index = global_indices[local_index];
                        im_global_indices[nnz] = global_indices[local_index];
                        im_values[nnz] = 1.0;
                        nnz++;

                        /* Now, go through the fluid nodes, find the global index in the matrix, copy the value and use MatSetValues() instead of MatSetValuesStencil. This is done b/c some of the fluid node might be on the neighboring processor and out of the range of the ghost layers. */
                        ib_index = Immersed_get_ib_global_index(grid->w_immersed, i, j, k);
                        ib_node  = Immersed_get_ib_node(grid->w_immersed, ib_index);
                        for (g=0; g<ib_node->n_fluid; g++) {

                            ii  = ib_node->fluid_index[g].x_index;
                            jj  = ib_node->fluid_index[g].y_index;
                            kk  = ib_node->fluid_index[g].z_index;

                            if ( (ii == i) && (jj == j) && (kk == k) ) { /* current immersed node is a function of itself value */

                                im_values[0] += ib_node->fluid_coef[g];
                            } else {

                                /* Only add the interpolation point if the coefficient is greater than zero. This will help with reducing the nonzero stencil */
                                if (fabs(ib_node->fluid_coef[g]) >= LSYS_DROP_TOLERANCE) {

                                    local_index = (kk - Ks_g) * ny_g * nx_g + (jj - Js_g) * nx_g + (ii - Is_g);
                                    im_values[nnz] = ib_node->fluid_coef[g];
                                    im_global_indices[nnz] = global_indices[local_index];
                                    nnz++;
                                }
                            }
                        } /* for */

                        /* All interpolation coefficients are scaled with a constant factor. Better convergence */
                        int s;
                        for (s=0; s<nnz; s++) im_values[s] *= ap_solid;

                        ierr = MatSetValues(w->A, 1, &row_global_index, nnz, im_global_indices, im_values, INSERT_VALUES); PETScErrAct(ierr);
                    } else {

                        values[nnz] = ap_solid; /* the value is not important */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                        /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                        ierr = MatSetValuesStencil(w->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);
                    } /* else BOUNDARY node */

                } /* else */
            } /* for: i*/
        } /*for: j */
    }/* for :k*/

    /* TBC */
    ierr = MatAssemblyBegin(w->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
    ierr = MatAssemblyEnd(w->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
}
/***************************************************************************************************/

/* Have to check this for inflow and outflow stuff*/
/* This function calculates velocities at the cell center using linear averaging of the 4 neighbors. */
void Velocity_cell_center( Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) {

    int i, j, k;
    int NI, NJ, NK;
    double ***u_data, ***v_data, ***w_data;
    double ***u_data_bc, ***v_data_bc, ***w_data_bc;
    double **u_outflow;
    double **u_inflow;
    double *xc, *yc, *zc;
    double *xu, *yv, *zw;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    Scalar velocities[2];
    Scalar velocity_bc;

    NI = grid->NI;
    NJ = grid->NJ;
    NK = grid->NK;

    /* First, update the local part of velocity u,v, w data generating ghost nodes */
    Communication_update_ghost_nodes(&grid->DA_3D, &u->G_data, &u->L_data, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D, &v->G_data, &v->L_data, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D, &w->G_data, &w->L_data, 'I');

    /* Now, get the array using a local pointer to a 3D array */
    ierr = DAVecGetArray(grid->DA_3D, u->L_data, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->L_data, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->L_data, (void ***)&w_data);PETScErrAct(ierr);

    /* Now get the global array for Velocity at cell center */
    ierr = DAVecGetArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc);PETScErrAct(ierr);

    u_inflow  = u->G_inflow;
    u_outflow = u->G_outflow;

    xc = grid->xc;
    yc = grid->yc;
    zc = grid->zc;

    xu = grid->xu;
    yv = grid->yv;
    zw = grid->zw;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* exclude the half cell added */
    Ie = min(Ie, NI);
    Je = min(Je, NJ);
    Ke = min(Ke, NK);

    /* Go through entire domain and interpolate the value of the velocity at the cell center */
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                //printf("C(i,j,k)=(%d,%d,%d) = %d\n", i, j, k , Grid_c_is_solid(grid, i, j, k));
                /* First, find u at cell center */
                if ( (params->inflow) && (i == 0) ) {

                    velocities[0].Value = u_inflow[k][j];

                } else {

                    velocities[0].Value = u_data[k][j][i];
                } /* else */
                if ( (params->outflow) && (i == NI-1) ) {

                    velocities[1].Value = u_outflow[k][j];

                } else {

                    velocities[1].Value = u_data[k][j][i+1];
                } /* else */

                velocities[0].X1    = xu[i];
                velocities[1].X1    = xu[i+1];
                velocity_bc.X1      = xc[i];
                u_data_bc[k][j][i]  = MyMath_interpolate_quantity(velocities, velocity_bc, 2);

                /* Now, find v at cell center */
                velocities[0].Value = v_data[k][j][i];
                velocities[1].Value = v_data[k][j+1][i];
                /* Note that for two point interpolation, always pass the position (no matter x or y)to the X1 */

                velocities[0].X1     = yv[j];
                velocities[1].X1     = yv[j+1];
                velocity_bc.X1       = yc[j];

                v_data_bc[k][j][i]   = MyMath_interpolate_quantity(velocities, velocity_bc, 2);

                /* Now, find w at cell center */
                velocities[0].Value = w_data[k][j][i];
                velocities[1].Value = w_data[k+1][j][i];

                /* Note that for two point interpolation, always pass the position (no matter x or y)to the X1 */
                velocities[0].X1     = zw[k];
                velocities[1].X1     = zw[k+1];
                velocity_bc.X1       = zc[k];

                w_data_bc[k][j][i]   = MyMath_interpolate_quantity(velocities, velocity_bc, 2);
            } /* for i*/
        } /* for j*/
    } /* for k*/

    /* Now free allocated memory. Always one should call RestoreArray() after done with the Vector
 read/write */

    ierr = DAVecRestoreArray(grid->DA_3D, u->L_data, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->L_data, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->L_data, (void ***)&w_data);PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc);PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc);PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc);PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function computes the RHS of the u-momentum linear system. Convective terms and inflow, outflow are
treated explicitly */
void Velocity_u_set_RHS( Velocity *u, Pressure *p, MAC_grid *grid, Parameters *params, double dt) {

    int i,j, k;
    int NI, NJ, NK;
    int ierr;
    double ***conv;
    double ***data;
    double ***dp_dx;
    double ***rhs_vec;
    double **u_outflow, **u_inflow;
    double *metric_xc, *metric_xu, *metric_yv, *metric_yc;
    double Re, ae, aw;
    double dEta, dXi, dPhi;
    double Cx, Cy, Cz;
    double a_dt;
    double rhs_value;
    double metric_px, metric_ex, metric_wx;
    double value;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int status;
    /* Pointer to a function */


    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D, u->G_data, (void ***)&data); PETScErrAct(ierr);

    /* Get convective terms global data  */
    ierr = DAVecGetArray(grid->DA_3D, u->G_conv, (void ***)&conv); PETScErrAct(ierr);

    /* Get pressure gradient from old time step */
    ierr = DAVecGetArray(grid->DA_3D, p->G_dp_dx, (void ***)&dp_dx); PETScErrAct(ierr);

    /* Now, got the RHS vector on current processor */
    ierr = DAVecGetArray(grid->DA_3D_Lsys, u->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    /***********/
    a_dt = 1.0/dt;

    u_outflow = u->G_outflow;
    u_inflow  = u->G_inflow;

    NI  = u->NI;
    NJ  = u->NJ;
    NK  = u->NK;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_xu = grid->metric_xu;
    metric_yv = grid->metric_yv;

    dEta = grid->dEta;
    dXi  = grid->dXi;
    dPhi = grid->dPhi;

    Re = params->Re;

    Cx = -1.0/(Re*dEta*dEta);
    Cy = -1.0/(Re*dXi*dXi);
    Cz = -1.0/(Re*dPhi*dPhi);

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

                status = Grid_get_u_status(grid, i, j, k);
#ifdef SOLVE_INSIDE_SOLID
                if( ( status == FLUID) || (status == BUFFER_SOLID)) {
#else
                if( ( status == FLUID) ) {
#endif
                    rhs_value = data[k][j][i]*a_dt
                            -conv[k][j][i] - dp_dx[k][j][i];

                    if (params->inflow) {

                        if (i == 1) {

                            /* get the metric coeeffients for non-uniform grid formulation */
                            metric_px = metric_xu[i];     /* (dx/deta)^-1 at (i+1/2,j)   */
                            metric_wx = metric_xc[i-1];   /* (dx/deta)^-1 at (i,j)     */

                            aw = -Cx * metric_px * metric_wx;

                            value = aw*u_inflow[k][j];

                            rhs_value +=value;
                        }
                    } /* if inflow */

                    if (params->outflow) {

                        if (i == NI-2) {

                            /* get the metric coeeffients for non-uniform grid formulation */
                            metric_px = metric_xu[i];   /* (dx/deta)^-1 at (i+1/2,j)   */
                            metric_ex = metric_xc[i];   /* (dx/deta)^-1 at (i+1,j)     */

                            ae = -Cx * metric_px * metric_ex;

                            value = ae*u_outflow[k][j];
                            rhs_value +=value;
                            //PetscPrintf(PCW, "(i,j,k)=(%d,%d,%d) ae:%f u_outflow:%f value:%f\n", i, j, k, ae, u_outflow[k][j], value);
                        }
                    } /* if outflow */

#ifdef DRIVEN_CAVITY
                    if (j == NJ-1) {

                        an = -Cy * metric_yv[j+1] * metric_yc[j];

                        value = an*1.0;

                        rhs_value +=value;
                    }
#endif

                } else { /* An Immersed node, set rhs equal to zero. This does not affect the real solution */

                    rhs_value = 0.0;
                }/* else */

                rhs_vec[k][j][i] = rhs_value;
            } /* for i*/
        } /* for j*/
    }/* for k*/

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data, (void ***)&data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, u->G_conv, (void ***)&conv); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->G_dp_dx, (void ***)&dp_dx); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_Lsys, u->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function computes the RHS of the v-momentum linear system. Convective terms and inflow, outflow are
treated explicitly */
void Velocity_v_set_RHS(Velocity *v, Pressure *p, Concentration **c, MAC_grid *grid, Parameters *params, double dt)  {

    int i, j, k;
    int NI, NJ, NK;
    int NConc;
    int ierr;
    double *yc, *yv;
    double ***conv;
    double ***data;
    double ***conc_total;
    double ***rhs_vec;
    double c_term;
    double *metric_xc, *metric_xu;
    double ***dp_dy, **v_outflow;
    double Re;
    double a_dt;
    double rhs_value;
    double dEta, dXi, dPhi;
    double Cx, Cy, Cz;
    double ae, value;
    double metric_px, metric_ex;
    double c_total, c_north_total;
    Scalar concentrations[2], concentration_at_v;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int status;


    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D, v->G_data, (void ***)&data); PETScErrAct(ierr);

    /* Get convective terms global data  */
    ierr = DAVecGetArray(grid->DA_3D, v->G_conv, (void ***)&conv); PETScErrAct(ierr);

    /* Get pressure gradient from old time step */
    ierr = DAVecGetArray(grid->DA_3D, p->G_dp_dy, (void ***)&dp_dy); PETScErrAct(ierr);

    /* Now, got the RHS vector on current processor */
    ierr = DAVecGetArray(grid->DA_3D_Lsys, v->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    Re    = params->Re;
    NConc = params->NConc;

    /* If we have more than one concentration field, c_term in y-momentum is store in c_total (always in
 c[0] */
    if (NConc > 1) {

        /* Generate local ghost nodes on current processor */
        Communication_update_ghost_nodes(&grid->DA_3D, &c[0]->G_c_total, &c[0]->L_c_total, 'I');


        /* Get total_conc array store in c[0] */
        ierr = DAVecGetArray(grid->DA_3D, c[0]->L_c_total, (void ***)&conc_total); PETScErrAct(ierr);

    } else { /* only one concentration field. Use c[0]->conc in v-momentum equation */

        /* Generate local ghost nodes on current processor */
        Communication_update_ghost_nodes(&grid->DA_3D, &c[0]->G_data, &c[0]->L_data, 'I');

        /* Get regular data array store in c[0] */
        ierr = DAVecGetArray(grid->DA_3D, c[0]->L_data, (void ***)&conc_total); PETScErrAct(ierr);
    } /* else */

    yc = grid->yc;
    yv = grid->yv;

    a_dt = 1.0/dt;

    dEta = grid->dEta;
    dXi  = grid->dXi;
    dPhi = grid->dPhi;

    metric_xc = grid->metric_xc;
    metric_xu = grid->metric_xu;

    Cx = -1.0/(Re*dEta*dEta);
    Cy = -1.0/(Re*dXi*dXi);
    Cz = -1.0/(Re*dPhi*dPhi);

    v_outflow = v->G_outflow;

    NI = v->NI;
    NJ = v->NJ;
    NK = v->NK;

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

                status = Grid_get_v_status(grid, i, j, k);
#ifdef SOLVE_INSIDE_SOLID
                if( ( status == FLUID) || (status == BUFFER_SOLID)) {
#else
                if( ( status == FLUID) ) {
#endif

                    c_total       = conc_total[k][j-1][i];;
                    c_north_total = conc_total[k][j][i];;

                    /* Find the total concentration term based on each individual concentration field */
                    /* C_t = SUM_i ( Conc_alpha(i) * C(i) ) i=0,1, ..., N-1*/

                    /* Linear interpolation of c_north and c_south to find c at v_grid*/
                    /* for the two point interpolation, always pass the position (no matter x or y or z) to X1 */
                    concentrations[0].Value  = c_total;
                    concentrations[1].Value  = c_north_total;
                    concentrations[0].X1     = yc[j-1];
                    concentrations[1].X1     = yc[j];
                    concentration_at_v.X1    = yv[j];

                    /* Just linear interpolation for the north and south conc to find c_term at the v_grid node */
                    c_term = 0.0;
#ifdef SOLVE_CONC
                    c_term = MyMath_interpolate_quantity(concentrations, concentration_at_v, 2);
#endif
                    //printf("c_term(i:%d,j:%d,k:%d)=%f\n", i, j, k, c_term);
                    //getchar();
                    rhs_value  = data[k][j][i]*a_dt
                            - conv[k][j][i] - dp_dy[k][j][i] - 0.0*c_term;

                    if (params->outflow) {

                        if (i==NI-1) {

                            metric_px = metric_xc[i];
                            metric_ex = metric_xu[i+1];

                            ae    = -Cx * metric_px * metric_ex;
                            value = ae*v_outflow[k][j];

                            rhs_value+=value;
                        }
                    }/* if outflow */

                    /*if ( (i==1) && (k==2) ){

      PetscPrintf(PCW, "v-q(i,j,k)=(%d,%d,%d) uddx:%f vddy:%f wddz:%f dpdy:%f rhs:%f v_old:%f\n", i, j, k, uddx[k][j][i], vddy[k][j][i], wddz[k][j][i], dp_dy[k][j][i], rhs_value, data[k][j][i]);

     }*/
                } else { /* An Immersed or Boundary node */

                    rhs_value = 0.0;
                }/* if solid */

                rhs_vec[k][j][i] = rhs_value;
            } /* for i*/
        } /* for j*/
    } /* for k */
    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data, (void ***)&data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_conv, (void ***)&conv); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->G_dp_dy, (void ***)&dp_dy); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_Lsys, v->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    /* Now, restore the allocated arrays */
    if (NConc > 1) {

        ierr = DAVecRestoreArray(grid->DA_3D, c[0]->L_c_total, (void ***)&conc_total); PETScErrAct(ierr);
    } else {

        ierr = DAVecRestoreArray(grid->DA_3D, c[0]->L_data, (void ***)&conc_total); PETScErrAct(ierr);
    }
}
/***************************************************************************************************/

/* This function computes the RHS of the w-momentum linear system. Convective terms and inflow, outflow are
treated explicitly */
void Velocity_w_set_RHS( Velocity *w, Pressure *p, MAC_grid *grid, Parameters *params, double dt) {

    int i, j, k;
    int NI, NJ, NK;
    int ierr;
    double ***conv;
    double ***data;
    double ***dp_dz;
    double ***rhs_vec;
    double **w_outflow;
    double *metric_xc, *metric_xu, *metric_yv, *metric_yc;
    double Re, ae;
    double dEta, dXi, dPhi;
    double Cx, Cy, Cz;
    double a_dt;
    double rhs_value, value;
    double metric_px, metric_ex;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int status;

    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D, w->G_data, (void ***)&data); PETScErrAct(ierr);

    /* Get convective terms global data  */
    ierr = DAVecGetArray(grid->DA_3D, w->G_conv, (void ***)&conv); PETScErrAct(ierr);

    /* Get pressure gradient from old time step */
    ierr = DAVecGetArray(grid->DA_3D, p->G_dp_dz, (void ***)&dp_dz); PETScErrAct(ierr);

    /* Now, got the RHS vector on current processor */
    ierr = DAVecGetArray(grid->DA_3D_Lsys, w->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    a_dt = 1.0/dt;

    w_outflow = w->G_outflow;

    NI  = w->NI;
    NJ  = w->NJ;
    NK  = w->NK;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_xu = grid->metric_xu;
    metric_yv = grid->metric_yv;

    dEta = grid->dEta;
    dXi  = grid->dXi;
    dPhi = grid->dPhi;

    Re = params->Re;

    Cx = -1.0/(Re*dEta*dEta);
    Cy = -1.0/(Re*dXi*dXi);
    Cz = -1.0/(Re*dPhi*dPhi);

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

                status = Grid_get_w_status(grid, i, j, k);
#ifdef SOLVE_INSIDE_SOLID
                if( ( status == FLUID) || (status == BUFFER_SOLID)) {
#else
                if( ( status == FLUID) ) {
#endif
                    rhs_value = data[k][j][i]*a_dt
                            - conv[k][j][i] - dp_dz[k][j][i];

                    if (params->outflow) {

                        if (i==NI-1) {

                            metric_px = metric_xc[i];
                            metric_ex = metric_xu[i+1];

                            ae = -Cx * metric_px * metric_ex;
                            value = ae*w_outflow[k][j];

                            rhs_value +=value;
                        }
                    } /* if outflow */


                } else { /* immersed node */

                    rhs_value = 0.0;
                }/* if !solid */

                rhs_vec[k][j][i] = rhs_value;

            } /* for i*/
        } /* for j*/
    }/* for k*/

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_data, (void ***)&data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_conv, (void ***)&conv); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, p->G_dp_dz, (void ***)&dp_dz); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_Lsys, w->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function does the velocity RK average */
void Velocity_rk_average(Velocity *vel, double old_frac, double new_frac) {

    int ierr;
    /* use Petsc vec product and summation */
    /*  data[k][j][i] = old_frac * data_old[k][j][i] + new_frac * data[k][j][i] */

    ierr = VecAXPBY(vel->G_data, old_frac, new_frac, vel->G_data_old); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function stores the old data for the outflow velocity. This is done for time integration purposes */
void Velocity_store_old_outflow(Velocity *vel) {

    int j, k;
    int Js, Ks;
    int Je, Ke;
    double **outflow, **outflow_old;

    /* Only apply for the processors including the yz plane at the end */
    if (vel->G_Ie == vel->NX) {

        outflow     = vel->G_outflow;
        outflow_old = vel->G_outflow_old;

        /* start index on current processor */
        Js = vel->G_Js;
        Ks = vel->G_Ks;

        /* end index on current processor */
        Je = vel->G_Je;
        Ke = vel->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                outflow_old[k][j] = outflow[k][j];
            } /* for j */
        } /* for k */
    } /* if */
}
/********************************************************************************/

/* This function takes the rk-average for velocity outflow */
void Velocity_rk_average_outflow(Velocity *vel, double old_frac, double new_frac) {

    int j, k;
    int Js, Ks;
    int Je, Ke;
    double **outflow, **outflow_old;

    /* Only apply for the processors including the yz plane at the end */
    if (vel->G_Ie == vel->NX) {

        outflow     = vel->G_outflow;
        outflow_old = vel->G_outflow_old;

        /* start index on current processor */
        Js = vel->G_Js;
        Ks = vel->G_Ks;

        /* end index on current processor */
        Je = vel->G_Je;
        Ke = vel->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                outflow[k][j] = old_frac * outflow_old[k][j] + new_frac * outflow[k][j];
            } /* for j */
        } /* for k */
    } /* if */
}
/***************************************************************************************************/


/* This function generates the ghost nodes for calculation of the u-convective terms using ENO scheme */
void Velocity_u_ENO_ghost_cells(ENO_Scheme *ENO, Velocity *u, MAC_grid *grid, Parameters *params, Indices
                                start_cell, Indices end_cell, char which_direction) {

    double *D0, *Position, *Velocity;
    double *xu, *yv, *zw;
    double **inflow, **outflow;
    double x_start, x_end, z_start;
    double y_start, y_end, z_end;
    int NI, NJ, NK;
    int i_start, j_start, k_start;
    int i_end, j_end, k_end;
    int ghost_index, interior_index, ex;
    int ENO_order;
    int Nmax_local;

    D0         = ENO->D0;
    Position   = ENO->Position;
    Velocity   = ENO->Velocity;
    ENO_order  = ENO->ENO_order;

    NI = u->NI;
    NJ = u->NJ;
    NK = u->NK;

    inflow  = u->G_inflow ;
    outflow = u->G_outflow;

    xu = grid->xu;
    yv = grid->yv;
    zw = grid->zw;

    /* Get the indices of the first fluid node */
    i_start = start_cell.x_index;
    j_start = start_cell.y_index;
    k_start = start_cell.z_index;

    /* Get the indices of the last fluid node */
    i_end   = end_cell.x_index;
    j_end   = end_cell.y_index;
    k_end   = end_cell.z_index;


    switch (which_direction) {

    case 'x':

        Nmax_local = i_end - i_start + 2*ENO_order + 1;
        /* generate the ghost nodes if the current processor is in the vicinity of the boundaries */
        if (i_start == 1) {

            /* First, for the left side of the array D0, corresponds to the left side of the domain */
            /* One sided zero u-velocity on the boundaries except for the inflow and outflow case */
            x_start = xu[i_start-1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 2;

                if (ex == ENO_order-1) { /* The last ghost node is exaclty on the wall. D0 is exacly zero.  */

                    Position[ghost_index] = x_start;
                    D0[ghost_index]       = 0.0;
                }
                else { /* Mirror the postion */

                    Position[ghost_index] = 2.0*x_start - Position[interior_index];
                    /* central */
                    D0[ghost_index]       = -D0[ENO_order];//-D0[interior_index];
                }
                /* Node right next to the inflow region */
                if (params->inflow) {

                    D0[ghost_index]       = inflow[k_start][j_start];
                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = inflow[k_start][j_start];
#endif
                }
            } /* for ex */
        } /* if i_start */

        /* generate the ghost nodes if the current processor is in the vicinity of the boundaries */
        if (i_end == NI-2) {

            /* Now, for the right side of the array D0*/
            x_end = xu[i_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex + 1;

                if (ex == ENO_order-1) { /* The last ghost node is exaclty on the wall. D0 is exaclty zero */

                    Position[ghost_index] = x_end;
                    D0[ghost_index]       = 0.0;
                }
                else { /* Mirror the postion */

                    Position[ghost_index] = 2.0*x_end - Position[interior_index];
                    D0[ghost_index]       = 0.0;//-D0[Nmax_local - ENO_order -1];//-D0[interior_index];

                }
                /* Node right before the outlfow region */
                if (params->outflow) {

                    D0[ghost_index]       = outflow[k_end][j_end];
                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = outflow[k_end][j_end];
#endif
                } /* if */
            } /* for ex */
        } /* if *i_end */
        break; /* case x */
        /**********************************/

        /* for vdudy */
    case 'y':

        Nmax_local = j_end - j_start + 2*ENO_order + 1;
        if (j_start == 0) {

            /* Apply Dirichlet B.C.. Zero u-velocity on the bottom boundary */
            /* Left side of D0, corresponds to the bottom boundary */

            y_start = yv[j_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*y_start - Position[interior_index];
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
                D0[ghost_index]       = -D0[ENO_order];//-D0[interior_index];
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
                D0[ghost_index]       = +D0[ENO_order];//-D0[interior_index];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[ENO_order];
#endif

            } /* for ex */
        } /* if j_start */

        if (j_end == NJ-1) {

            /* Now, for the right side of the array D0, corresponds to the top boundary. No variable geomtery on top. */
            y_end = yv[j_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*y_end - Position[interior_index];
#ifdef TOP_WALL_VELOCITY_NOSLIP
                D0[ghost_index]       = -D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
                D0[ghost_index]       = +D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
#ifdef DRIVEN_CAVITY
                D0[ghost_index]       = 1.0;
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[Nmax_local - ENO_order -1];
#endif

            } /* for ex */
        } /* if */
        break; /* case y */
        /**********************************/

        /* for wdudz */
    case 'z':

        Nmax_local = k_end - k_start + 2*ENO_order + 1;
        if (k_start == 0) {

            /* Left side of D0, corresponds to the back boundary */
            z_start = zw[k_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*z_start - Position[interior_index];

#ifdef BACK_WALL_VELOCITY_NOSLIP
                D0[ghost_index]       = -D0[ENO_order];
#endif
#ifdef BACK_WALL_VELOCITY_FREESLIP
                D0[ghost_index]       = +D0[ENO_order];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[ENO_order];
#endif
            } /* for ex */
        } /* if k_start */

        if (k_end == NK-1) {

            /* Now, for the right side of the array D0, corresponds to the side boundary */
            z_end = zw[k_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*z_end - Position[interior_index];

#ifdef FRONT_WALL_VELOCITY_NOSLIP
                D0[ghost_index]       = -D0[Nmax_local - ENO_order -1];
#endif
#ifdef FRONT_WALL_VELOCITY_FREESLIP
                D0[ghost_index]       = +D0[Nmax_local - ENO_order -1];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[Nmax_local - ENO_order -1];
#endif
            } /* for ex */
        } /* if k_end */
        break; /* case z */
    } /* switch */
}
/***************************************************************************************************/

/* This function generates the ghost nodes for calculation of the v-convective terms using ENO scheme */
void Velocity_v_ENO_ghost_cells(ENO_Scheme *ENO, Velocity *v, MAC_grid *grid, Parameters *params, Indices
                                start_cell, Indices end_cell, char which_direction) {

    double *D0, *Position, *Velocity;
    double *xu, *yv, *zw;
    double **outflow;
    double x_start, x_end;
    double y_start, y_end;
    double z_start, z_end;
    int NI, NJ, NK;
    int ghost_index, interior_index, ex;
    int ENO_order;
    int Nmax_local;
    int i_start, j_start, k_start;
    int i_end, j_end, k_end;

    D0         = ENO->D0;
    Position   = ENO->Position;
    Velocity   = ENO->Velocity;
    ENO_order  = ENO->ENO_order;

    NI = v->NI;
    NJ = v->NJ;
    NK = v->NK;

    outflow  = v->G_outflow;

    xu = grid->xu;
    yv = grid->yv;
    zw = grid->zw;

    /* Get the indices of the first fluid node */
    i_start = start_cell.x_index;
    j_start = start_cell.y_index;
    k_start = start_cell.z_index;

    /* Get the indices of the last fluid node */
    i_end   = end_cell.x_index;
    j_end   = end_cell.y_index;
    k_end   = end_cell.z_index;

    switch (which_direction) {

    /* udvdx */
    case 'x':

        Nmax_local = i_end - i_start + 2*ENO_order + 1;
        if (i_start == 0) {

            /* First, for the left side of the array D0, corresponds to the left side of the domain*/
            x_start = xu[i_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*x_start - Position[interior_index];
#ifdef LEFT_WALL_VELOCITY_NOSLIP
                D0[ghost_index]       = -D0[ENO_order];
#endif
#ifdef LEFT_WALL_VELOCITY_FREESLIP
                D0[ghost_index]       = +D0[ENO_order];
#endif
                //My debug: central
                //Velocity[ghost_index] = -Velocity[ENO_order];
                if (params->inflow) {

                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = Velocity[ENO_order];
#endif
                }
            } /* for ex */
        }/* if i_start */

        if (i_end == NI-1) {

            /* Now, for the right side of the array D0*/
            x_end = xu[i_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*x_end - Position[interior_index];

                if (params->outflow) {

                    D0[ghost_index]       = outflow[k_end][j_end];
                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = Velocity[Nmax_local - ENO_order -1];
#endif
                } else {

#ifdef RIGHT_WALL_VELOCITY_NOSLIP
                    D0[ghost_index]       = -D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
#ifdef RIGHT_WALL_VELOCITY_FREESLIP

                    D0[ghost_index]       = +D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = -Velocity[Nmax_local - ENO_order -1];
#endif
                }
            } /* for ex */
        } /* if i_end */
        break; /* case x */
        /*************************/

        /* vdvdy */
    case 'y':

        Nmax_local = j_end - j_start + 2*ENO_order + 1;

        if (j_start == 1) {

            /* Apply Dirichlet B.C.. Zero v-velocity on the bottom boundary */
            y_start = yv[j_start-1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 2;

                if (ex == ENO_order-1) { /* last ghost node is exactly on the wall. v is zero */

                    Position[ghost_index] = y_start;
                    D0[ghost_index]       = 0.0;

                } else {

                    Position[ghost_index] = 2.0*y_start - Position[interior_index];
                    D0[ghost_index]       = 0.0;//-D0[ENO_order];//-D0[interior_index];
                } /* else */
            } /* for ex */
        } /* if j */

        /* Now, for the right side of the array D0, corresponds to the top boundary. No variable geomtery. */
        if (j_end == NJ-2) {

            y_end = yv[j_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex + 1;

                if (ex == ENO_order-1) { /* last ghost node is exactly on the wall. v is zero */

                    Position[ghost_index] = y_end;
                    D0[ghost_index]       = 0.0;

                } else {

                    Position[ghost_index] = 2.0*y_end - Position[interior_index];
                    D0[ghost_index]       = 0.0;//-D0[Nmax_local - ENO_order -1];//-D0[interior_index];
                } /* else */
            } /* for ex */
        } /* if j */
        break; /* case y */
        /*************************/

        /* for wdvdz */
    case 'z':

        Nmax_local = k_end - k_start + 2*ENO_order + 1;
        if ( k_start == 0 ) {

            z_start = zw[k_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*z_start - Position[interior_index];

#ifdef BACK_WALL_VELOCITY_NOSLIP

                D0[ghost_index]       = -D0[ENO_order];
#endif
#ifdef BACK_WALL_VELOCITY_FREESLIP

                D0[ghost_index]       = +D0[ENO_order];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[ENO_order];
#endif
            } /* for ex */
        } /* if k */

        /* Now, for the right side of the array D0, corresponds to the side boundary */
        if (k_end == NK-1) {

            z_end = zw[k_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*z_end - Position[interior_index];

#ifdef FRONT_WALL_VELOCITY_NOSLIP

                D0[ghost_index]       = -D0[Nmax_local - ENO_order -1];
#endif
#ifdef FRONT_WALL_VELOCITY_FREESLIP

                D0[ghost_index]       = +D0[Nmax_local - ENO_order -1];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[Nmax_local - ENO_order -1];
#endif
            } /* for ex */
        } /* if */
        break;
    } /* switch */
}
/***************************************************************************************************/

/* This function generates the ghost nodes for calculation of the w-convective terms using ENO scheme */
void Velocity_w_ENO_ghost_cells(ENO_Scheme *ENO, Velocity *w, MAC_grid *grid, Parameters *params, Indices
                                start_cell, Indices end_cell, char which_direction) {

    double *D0, *Position, *Velocity;
    double *xu, *yv, *zw;
    double **outflow;
    double x_start, x_end;
    double y_start, y_end;
    double z_start, z_end;
    int NI, NJ, NK;
    int ghost_index, interior_index, ex;
    int ENO_order;
    int Nmax_local;
    int i_start, j_start, k_start;
    int i_end, j_end, k_end;


    D0         = ENO->D0;
    Position   = ENO->Position;
    Velocity   = ENO->Velocity;
    ENO_order  = ENO->ENO_order;

    NI = w->NI;
    NJ = w->NJ;
    NK = w->NK;

    outflow = w->G_outflow;

    xu = grid->xu;
    yv = grid->yv;
    zw = grid->zw;

    /* Get the indices of the first fluid node */
    i_start = start_cell.x_index;
    j_start = start_cell.y_index;
    k_start = start_cell.z_index;

    /* Get the indices of the last fluid node */
    i_end   = end_cell.x_index;
    j_end   = end_cell.y_index;
    k_end   = end_cell.z_index;

    switch (which_direction) {

    /* udwdx */
    case 'x':

        Nmax_local = i_end - i_start + 2*ENO_order + 1;
        if (i_start == 0) {

            x_start = xu[i_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*x_start - Position[interior_index];
#ifdef LEFT_WALL_VELOCITY_NOSLIP

                D0[ghost_index]       = -D0[ENO_order];//-D0[interior_index];
#endif
#ifdef LEFT_WALL_VELOCITY_FREESLIP

                D0[ghost_index]       = +D0[ENO_order];//-D0[interior_index];
#endif

#ifdef VELOCITY_CONVECTIVE_CENTRAL
                if (params->inflow) {
                    //My debug: central
                    Velocity[ghost_index] = Velocity[ENO_order];
                } else {

                    //My debug: central
                    Velocity[ghost_index] = -Velocity[ENO_order];
                }
#endif
            } /* for ex */
        } /* if i_start */

        /* Now, for the right side of the array D0*/
        if (i_end == NI-1) {

            x_end = xu[i_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*x_end - Position[interior_index];

                if (params->outflow) {

                    D0[ghost_index]       = outflow[k_end][j_end];
                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = Velocity[Nmax_local - ENO_order -1];
#endif
                } else {

#ifdef RIGHT_WALL_VELOCITY_NOSLIP

                    D0[ghost_index]       = -D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
#ifdef RIGHT_WALL_VELOCITY_FREESLIP

                    D0[ghost_index]       = +D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
                    //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                    Velocity[ghost_index] = -Velocity[Nmax_local - ENO_order -1];
#endif
                } /* else */
            } /* for ex */
        } /* for i_end */
        break; /* case x*/

        /* vdwdy */
    case 'y':

        Nmax_local = j_end - j_start + 2*ENO_order + 1;
        /* Apply Dirichlet B.C.. Zero v-velocity on the bottom boundary */

        if (j_start == 0 ) {

            y_start = yv[j_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*y_start - Position[interior_index];
#ifdef BOTTOM_WALL_VELOCITY_NOSLIP

                D0[ghost_index]       = -D0[ENO_order];//-D0[interior_index];
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP

                D0[ghost_index]       = +D0[ENO_order];//-D0[interior_index];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[ENO_order];
#endif
            } /* for ex */
        } /* if j_start */

        /* Now, for the right side of the array D0, corresponds to the top boundary. No variable geomtery. */
        if (j_end == NJ-1) {

            y_end = yv[j_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*y_end - Position[interior_index];
#ifdef TOP_WALL_VELOCITY_NOSLIP

                D0[ghost_index]       = -D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP

                D0[ghost_index]       = +D0[Nmax_local - ENO_order -1];//-D0[interior_index];
#endif
                //My debug: central
#ifdef VELOCITY_CONVECTIVE_CENTRAL
                Velocity[ghost_index] = -Velocity[Nmax_local - ENO_order -1];
#endif
            } /* for ex */
        } /* if j_end */
        break; /* case y*/

        /* for wdwdz */
    case 'z':

        Nmax_local = k_end - k_start + 2*ENO_order + 1;
        /* Apply Dirichlet B.C.. Zero w-velocity on the bottom boundary */

        if (k_start == 1) {

            z_start = zw[k_start-1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 2;

                if (ex == ENO_order-1) { /* last ghost node is exactly on the wall. v is zero */

                    Position[ghost_index] = z_start;
                    D0[ghost_index]       = 0.0;

                } else {

                    Position[ghost_index] = 2.0*z_start - Position[interior_index];
                    D0[ghost_index]       = 0.0;//-D0[ENO_order];//-D0[interior_index];
                } /* else */
            } /* for ex */
        } /* if k_start */

        /* Now, for the right side of the array D0, corresponds to the top boundary. No variable geomtery. */
        if (k_end == NK-2) {

            z_end = zw[k_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex + 1;

                if (ex == ENO_order-1) { /* last ghost node is exactly on the wall. v is zero */

                    Position[ghost_index] = z_end;
                    D0[ghost_index]       = 0.0;

                } else {

                    Position[ghost_index] = 2.0*z_end - Position[interior_index];
                    D0[ghost_index]       = 0.0;//-D0[Nmax_local - ENO_order -1];//-D0[interior_index];
                } /* else */
            } /* for ex */
        } /* if */
        break; /* case z*/
    } /* switch */
}
/***************************************************************************************************/

/* This function computes the total vicous dissipation due to the fluid motion */
/* e= 2.0/Re S_ij S_ij */
void Velocity_compute_viscous_dissipation_rate(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) {

    double ***u_data=NULL;
    double ***v_data=NULL;
    double ***w_data=NULL;

    /* It is defined as
  e = 2/Re S_ij S_ij
 */
    double G_dissipation_rate = 0.0;

    int ierr;
    double ***strain_rate = NULL;
    ierr = DAVecGetArray(grid->DA_3D, (*u->G_anything), (void ***)&strain_rate); PETScErrAct(ierr);

    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D, u->L_data_bc, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->L_data_bc, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->L_data_bc, (void ***)&w_data);PETScErrAct(ierr);


    /* Grid coordinates */
    double *xc = grid->xc;
    double *yc = grid->yc;
    double *zc = grid->zc;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    int NI = grid->NI;
    int NJ = grid->NJ;
    int NK = grid->NK;

    /* Coefficients used to find the values of the ghost nodes across the box boundaries */
    double k_bottom = 0.0;
    double k_top = 0.0;
    double k_left = 0.0;
    double k_right = 0.0;
    double k_back = 0.0;
    double k_front = 0.0;


#ifdef BOTTOM_WALL_VELOCITY_NOSLIP
    k_bottom = -1.0;
#endif
#ifdef BOTTOM_WALL_VELOCITY_FREESLIP
    k_bottom = +1.0;
#endif
#ifdef TOP_WALL_VELOCITY_NOSLIP
    k_top = -1.0;
#endif
#ifdef TOP_WALL_VELOCITY_FREESLIP
    k_top = +1.0;
#endif
#ifdef LEFT_WALL_VELOCITY_NOSLIP
    k_left = -1.0;
#endif
#ifdef LEFT_WALL_VELOCITY_FREESLIP
    k_left = +1.0;
#endif
#ifdef RIGHT_WALL_VELOCITY_NOSLIP
    k_right = -1.0;
#endif
#ifdef RIGHT_WALL_VELOCITY_FREESLIP
    k_right = +1.0;
#endif
#ifdef BACK_WALL_VELOCITY_NOSLIP
    k_back = -1.0;
#endif
#ifdef BACK_WALL_VELOCITY_FREESLIP
    k_back = +1.0;
#endif
#ifdef FRONT_WALL_VELOCITY_NOSLIP
    k_front = -1.0;
#endif
#ifdef FRONT_WALL_VELOCITY_FREESLIP
    k_front = +1.0;
#endif

    double coef = -2.0/params->Re;
    double **u_inflow=u->G_inflow;
    double **u_outflow=u->G_outflow;
    double **v_outflow=v->G_outflow;
    double **w_outflow=w->G_outflow;
    double U_w=0.0; double U_e=0.0;
    double V_w=0.0; double V_e=0.0;
    double W_w=0.0; double W_e=0.0;
    double U_s=0.0; double U_n=0.0;
    double V_s=0.0; double V_n=0.0;
    double W_s=0.0; double W_n=0.0;
    double U_f=0.0; double U_b=0.0;
    double V_f=0.0; double V_b=0.0;
    double W_f=0.0; double W_b=0.0;
    int i, j, k;
    for (k=Ks; k<Ke; k++) {

        int kf = min(k+1, NK-1); int kb = max(k-1, 0);
        double a_dz = 1.0/(zc[kf] - zc[kb]);
        if ( (k == 0) || (k == NK-1) ) a_dz *= 0.5;

        for (j=Js; j<Je; j++) {

            int jn = min(j+1, NJ-1); int js = max(j-1, 0);
            double a_dy = 1.0/(yc[jn] - yc[js]);
            if ( (j == 0) || (j == NJ-1) ) a_dy *= 0.5;

            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /*only include if point is fluid*/

                    int ie = min(i+1, NI-1); int iw = max(i-1, 0);
                    double a_dx = 1.0/(xc[ie] - xc[iw]);
                    if ( (i == 0) || (i == NI-1) ) a_dx *= 0.5;

                    if (i == 0) {

                        if (params->inflow) {

                            U_w = u_inflow[k][j];
                            V_w = k_left*v_data[k][j][i];
                            W_w = k_left*w_data[k][j][i];
                        } else {

                            U_w = k_left*u_data[k][j][i];
                            V_w = k_left*v_data[k][j][i];
                            W_w = k_left*w_data[k][j][i];
                        }

                    } else {

                        U_w = u_data[k][j][i-1];
                        V_w = v_data[k][j][i-1];
                        W_w = w_data[k][j][i-1];

                    } /* else i == 0 */
                    if (j == 0) {

                        U_s = k_bottom*u_data[k][j][i];
                        V_s = k_bottom*v_data[k][j][i];
                        W_s = k_bottom*w_data[k][j][i];
                    } else {

                        U_s = u_data[k][j-1][i];
                        V_s = v_data[k][j-1][i];
                        W_s = w_data[k][j-1][i];
                    } /* else j == 0 */

                    if (k == 0) {

                        U_b = k_back*u_data[k][j][i];
                        V_b = k_back*v_data[k][j][i];
                        W_b = k_back*w_data[k][j][i];
                    } else {

                        U_b = u_data[k-1][j][i];
                        V_b = v_data[k-1][j][i];
                        W_b = w_data[k-1][j][i];
                    } /* else k == 0 */

                    if (i == NI-1) {

                        if (params->outflow) {

                            U_e = u_outflow[k][j];
                            V_e = v_outflow[k][j];
                            W_e = w_outflow[k][j];
                        } else {

                            U_e = k_right*u_data[k][j][i];
                            V_e = k_right*v_data[k][j][i];
                            W_e = k_right*w_data[k][j][i];
                        }
                    } else {

                        U_e = u_data[k][j][i+1];
                        V_e = v_data[k][j][i+1];
                        W_e = w_data[k][j][i+1];
                    } /* else i == NI-1 */

                    if (j == NJ-1) {

                        U_n = k_top*u_data[k][j][i];
                        V_n = k_top*v_data[k][j][i];
                        W_n = k_top*w_data[k][j][i];
                    } else {

                        U_n = u_data[k][j+1][i];
                        V_n = v_data[k][j+1][i];
                        W_n = w_data[k][j+1][i];
                    } /* else j == NJ-1 */

                    if (k == NK-1) {

                        U_f = k_front*u_data[k][j][i];
                        V_f = k_front*v_data[k][j][i];
                        W_f = k_front*w_data[k][j][i];
                    } else {

                        U_f = u_data[k+1][j][i];
                        V_f = v_data[k+1][j][i];
                        W_f = w_data[k+1][j][i];
                    } /* else k == NK-1 */

                    double dudx = (U_e -U_w) * a_dx;
                    double dudy = (U_n -U_s) * a_dy;
                    double dudz = (U_f -U_b) * a_dz;

                    double dvdx = (V_e -V_w) * a_dx;
                    double dvdy = (V_n -V_s) * a_dy;
                    double dvdz = (V_f -V_b) * a_dz;

                    double dwdx = (W_e -W_w) * a_dx;
                    double dwdy = (W_n -W_s) * a_dy;
                    double dwdz = (W_f -W_b) * a_dz;

                    /* Strain rate tensor: Symmetric */
                    /* S_ij = 0.5 * (du_i/dx_j + du_j/dx_i) */
                    double S_xx = dudx;
                    double S_yy = dvdy;
                    double S_zz = dwdz;
                    double S_xy = 0.5*(dudy + dvdx);
                    double S_xz = 0.5*(dudz + dwdx);
                    double S_yz = 0.5*(dvdz + dwdy);

                    G_dissipation_rate = S_xx*S_xx + S_yy*S_yy + S_zz*S_zz + 2.0*(S_xy*S_xy + S_xz*S_xz + S_yz*S_yz) ;

                    strain_rate[k][j][i] = G_dissipation_rate * coef;
                } /* if */
            } /* for i*/
        } /* for j*/
    }/* for k*/

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, (*u->G_anything), (void ***)&strain_rate); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, u->L_data_bc, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->L_data_bc, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->L_data_bc, (void ***)&w_data); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function initializes the velocity field to a nonzero value */
/* Make sure the initial velocity field is divergence free */
void Velocity_nonzero_initialize(Velocity *vel, MAC_grid *grid, Parameters *params){ 

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, j, k;
    double ***data;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = min(grid->G_Ie, grid->NI);
    Je = min(grid->G_Je, grid->NJ);
    Ke = min(grid->G_Ke, grid->NK);

    ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&data); PETScErrAct(ierr);


    double *xq = NULL;
    double *yq = NULL;
    double *zq = NULL;

    switch (vel->component) {

    case 'u':
        xq = grid->xu; yq = grid->yu; zq = grid->zu; break;
    case 'v':
        xq = grid->xv; yq = grid->yv; zq = grid->zv; break;
    case 'w':
        xq = grid->xw; yq = grid->yw; zq = grid->zw; break;
    }

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                double value=0.0;
                switch (vel->component) {

                case 'u':
                    value = sin(2.0*PI*xq[i]/params->Lx)*cos(2.0*PI*yq[j]/params->Ly)*cos(2.0*PI*zq[k]/params->Lz); break;

                case 'v':
                    value = cos(2.0*PI*xq[i]/params->Lx)*sin(2.0*PI*yq[j]/params->Ly)*cos(2.0*PI*zq[k]/params->Lz); break;

                case 'w':
                    value = cos(2.0*PI*xq[i]/params->Lx)*cos(2.0*PI*yq[j]/params->Ly)*sin(2.0*PI*zq[k]/params->Lz); break;
                }

                data[k][j][i] = value;
            }
        }
    }

    ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&data); PETScErrAct(ierr);
}
/***************************************************************************************************/

/* just to debug */
void Velocity_debug_flow(Velocity *vel, MAC_grid *grid, Parameters *params) {

    double yy;
    double an_sol, num_sol;
    double error;
    double norm2=0.0;
    double norm1=0.0;
    double norm_inf=0.0;
    int Nt=0;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je-1;
    int Ke = grid->G_Ke-1;

    double ***data;
    int ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&data); PETScErrAct(ierr);

    double *yu = grid->yu;

    double u_max = 1.0;

    int i, j, k;
    int i_inf=-1;
    int j_inf=-1;
    int k_inf=-1;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                num_sol = data[k][j][i];
                if (vel->component == 'u') {

                    yy = yu[j]/params->Ly;

                    an_sol = 4.0*u_max*(yy - yy*yy);
                    if (i == 0) {

                        num_sol = vel->G_inflow[k][j];
                    }
                } /* if */ else {

                    an_sol = 0.0;
                }
                if (i == grid->NI) {

                    num_sol = vel->G_outflow[k][j];
                }

                error = fabs(an_sol - num_sol);
                norm2 += error*error;
                norm1 += error;
                if (error > norm_inf) {

                    norm_inf = error;
                    i_inf = i;
                    j_inf = j;
                    k_inf = k;
                } /* if */
                Nt++;
            } /* for i */
        } /* for j */
    } /* for k */

    norm2  = sqrt(norm2/Nt);
    norm1 /= Nt;

    PetscPrintf(PCW, "Velocity.c/ Error analysis for %c-velocity norm1:%f norm2:%f norm_inf:%f at (i,j,k)=(%d,%d,%d)\n", vel->component, norm1, norm2, norm_inf, i_inf, j_inf, k_inf);
    //getchar();
    ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&data); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function computes one vorticity component */
/* w = w_xI + w_yJ + w_zK */
/* which_component:
  'x': w_x
  'y': w_y
  'z': w_z
*/
void Velocity_compute_vorticity(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, char which_component, Vec *G_vort) {

    int ierr;
    double ***u_data=NULL;
    double ***v_data=NULL;
    double ***w_data=NULL;
    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D, u->L_data_bc, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->L_data_bc, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->L_data_bc, (void ***)&w_data); PETScErrAct(ierr);

    double ***vort=NULL;
    ierr = DAVecGetArray(grid->DA_3D, (*G_vort), (void ***)&vort); PETScErrAct(ierr);

    int NI  = grid->NI;
    int NJ  = grid->NJ;
    int NK  = grid->NK;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    double *xc = grid->xc;
    double *yc = grid->yc;
    double *zc = grid->zc;

    int i, j, k;
    switch (which_component) {
    case 'x':

        for (k=Ks; k<Ke; k++) {

            int k_f = min(k+1, NK-1);
            int k_b = max(k-1, 0);
            double a_dz = 1.0/(zc[k_f] - zc[k_b]);

            for (j=Js; j<Je; j++) {

                int j_n = min(j+1, NJ-1);
                int j_s = max(j-1, 0);
                double a_dy = 1.0/(yc[j_n] - yc[j_s]);
                for (i=Is; i<Ie; i++) {

                    if ( Grid_get_c_status(grid, i, j, k) == FLUID ) {

                        double dwdy = (w_data[k][j_n][i] - w_data[k][j_s][i]) * a_dy;
                        double dvdz = (v_data[k_f][j][i] - v_data[k_b][j][i]) * a_dz;

                        vort[k][j][i] = dwdy - dvdz;

                    } /* if */

                } /* for i */
            } /* for j */
        } /* for k */
        break;
        /*************************************************/
        /*************************************************/
    case 'y':
        for (k=Ks; k<Ke; k++) {

            int k_f = min(k+1, NK-1);
            int k_b = max(k-1, 0);
            double a_dz = 1.0/(zc[k_f] - zc[k_b]);

            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {

                    if ( Grid_get_c_status(grid, i, j, k) == FLUID ) {

                        int i_e = min(i+1, NI-1);
                        int i_w = max(i-1, 0);
                        double a_dx = 1.0/(xc[i_e] - xc[i_w]);

                        double dwdx = (w_data[k][j][i_e] - w_data[k][j][i_w]) * a_dx;
                        double dudz = (u_data[k_f][j][i] - u_data[k_b][j][i]) * a_dz;

                        vort[k][j][i] = dudz - dwdx;

                    } /* if */

                } /* for i */
            } /* for j */
        } /* for k */
        break;
        /*************************************************/
        /*************************************************/
    case 'z':

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                int j_n = min(j+1, NJ-1);
                int j_s = max(j-1, 0);
                double a_dy = 1.0/(yc[j_n] - yc[j_s]);
                for (i=Is; i<Ie; i++) {

                    if ( Grid_get_c_status(grid, i, j, k) == FLUID ) {

                        int i_e = min(i+1, NI-1);
                        int i_w = max(i-1, 0);
                        double a_dx = 1.0/(xc[i_e] - xc[i_w]);

                        double dvdx = (v_data[k][j][i_e] - v_data[k][j][i_w]) * a_dx;
                        double dudy = (u_data[k][j_n][i] - u_data[k][j_s][i]) * a_dy;

                        vort[k][j][i] = dvdx - dudy;

                    } /* if */
                } /* for i */
            } /* for j */
        } /* for k */
        break;
    default:
        PetscPrintf(PCW, "Velocity.c/ No such vorticity component %c\n", which_component);

    } /* switch */

    /* Restore arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, (*G_vort), (void ***)&vort); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, u->L_data_bc, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->L_data_bc, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->L_data_bc, (void ***)&w_data); PETScErrAct(ierr);
}
/***************************************************************************************************/

/*  This function finds the velocity components at each others locations...This is done for computing
convective terms */
void Velocity_transpose_velocities(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, char which_quantity) {

    int i, j ,k;
    int uNI, vNI, wNI;
    double ***v_at_ugrid, ***w_at_ugrid;
    double ***u_at_vgrid, ***w_at_vgrid;
    double ***u_at_wgrid, ***v_at_wgrid;
    double ***u_data, ***v_data, ***w_data;
    double *xu, *xv, *xw;
    double *yu, *yv, *yw;
    double *zu, *zv, *zw;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    Scalar velocities[4];
    Scalar vel_transpose;

    /*get sizes of velocity fields */
    uNI = u->NI;
    vNI = v->NI;
    wNI = w->NI;

    /* Coordinates of MAC-grid velocities */
    xu = grid->xu;
    xv = grid->xv;
    xw = grid->xw;

    yu = grid->yu;
    yv = grid->yv;
    yw = grid->yw;


    zu = grid->zu;
    zv = grid->zv;
    zw = grid->zw;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;


    switch (which_quantity) {

    case 'u':

        ierr = DAVecGetArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
        ierr = DAVecGetArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

        ierr = DAVecGetArray(grid->DA_3D, (*u->G_v_transposed), (void ***)&v_at_ugrid); PETScErrAct(ierr);
        ierr = DAVecGetArray(grid->DA_3D, (*u->G_w_transposed), (void ***)&w_at_ugrid); PETScErrAct(ierr);
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {

                    /********************************************************************/
                    /* Find v at u_grid */
                    /*
  v2------v3
   |  u   |
  v0------v1

     */
                    if (Grid_get_u_status(grid, i, j, k) != BOUNDARY) {

                        velocities[0].Value = v_data[k][j][i-1];
                        velocities[1].Value = v_data[k][j][i];
                        velocities[2].Value = v_data[k][j+1][i-1];
                        velocities[3].Value = v_data[k][j+1][i];

                        velocities[0].X1 = xv[i-1];
                        velocities[1].X1 = xv[ i ];
                        velocities[2].X1 = xv[i-1];
                        velocities[3].X1 = xv[ i ];

                        velocities[0].X2 = yv[ j ];
                        velocities[1].X2 = yv[ j ];
                        velocities[2].X2 = yv[j+1];
                        velocities[3].X2 = yv[j+1];

                        vel_transpose.X1 = xu[i];
                        vel_transpose.X2 = yu[j];

                        //printf("i:%d j:%d k:%d\n", i, j, k);
                        v_at_ugrid[k][j][i] = MyMath_interpolate_quantity(velocities, vel_transpose, 4);

                        /* Find w at u_grid */

                        /*
  w2------w3
   |  u   |
  w0------w1

     */
                        velocities[0].Value = w_data[k][j][i-1];
                        velocities[1].Value = w_data[k][j][i];
                        velocities[2].Value = w_data[k+1][j][i-1];
                        velocities[3].Value = w_data[k+1][j][i];

                        velocities[0].X1 = xw[i-1];
                        velocities[1].X1 = xw[ i ];
                        velocities[2].X1 = xw[i-1];
                        velocities[3].X1 = xw[ i ];

                        velocities[0].X2 = zw[ k ];
                        velocities[1].X2 = zw[ k ];
                        velocities[2].X2 = zw[k+1];
                        velocities[3].X2 = zw[k+1];

                        vel_transpose.X1 = xu[i];
                        vel_transpose.X2 = zu[k];

                        w_at_ugrid[k][j][i] = MyMath_interpolate_quantity(velocities, vel_transpose, 4);
                    } /* if u_status */
                } /* for i */
            } /* for j */
        } /* for k */

        ierr = DAVecRestoreArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
        ierr = DAVecRestoreArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

        ierr = DAVecRestoreArray(grid->DA_3D, (*u->G_v_transposed), (void ***)&v_at_ugrid); PETScErrAct(ierr);
        ierr = DAVecRestoreArray(grid->DA_3D, (*u->G_w_transposed), (void ***)&w_at_ugrid); PETScErrAct(ierr);
        break;

    case 'v':

        ierr = DAVecGetArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
        ierr = DAVecGetArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

        ierr = DAVecGetArray(grid->DA_3D, (*v->G_u_transposed), (void ***)&u_at_vgrid); PETScErrAct(ierr);
        ierr = DAVecGetArray(grid->DA_3D, (*v->G_w_transposed), (void ***)&w_at_vgrid); PETScErrAct(ierr);
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {
                    /********************************************************************/
                    /* Find u at v_grid */
                    /*
  u2------u3
   |  v   |
  u0------u1

     */
                    if (Grid_get_v_status(grid, i, j, k) != BOUNDARY) {

                        if ( (params->inflow) && (i == 0) ){


                            velocities[0].Value = u->G_inflow[k][j-1];
                            velocities[2].Value = u->G_inflow[k][ j ];
                        } else {

                            velocities[0].Value = u_data[k][j-1][i];
                            velocities[2].Value = u_data[k][j][i];
                        }

                        if ( (params->outflow) && (i == vNI-1) ){

                            velocities[1].Value = u->G_outflow[k][j-1];
                            velocities[3].Value = u->G_outflow[k][ j ];
                        } else {

                            velocities[1].Value = u_data[k][j-1][i+1];
                            velocities[3].Value = u_data[k][j][i+1];
                        }

                        velocities[0].X1 = xu[ i ];
                        velocities[1].X1 = xu[i+1];
                        velocities[2].X1 = xu[ i ];
                        velocities[3].X1 = xu[i+1];

                        velocities[0].X2 = yu[j-1];
                        velocities[1].X2 = yu[j-1];
                        velocities[2].X2 = yu[ j ];
                        velocities[3].X2 = yu[ j ];

                        vel_transpose.X1 = xv[i];
                        vel_transpose.X2 = yv[j];

                        u_at_vgrid[k][j][i] = MyMath_interpolate_quantity(velocities, vel_transpose, 4);

                        /* Find w at v_grid */

                        /*
  w2------w3
   |  v   |
  w0------w1


     */
                        velocities[0].Value = w_data[k][j-1][i];
                        velocities[1].Value = w_data[k+1][j-1][i];
                        velocities[2].Value = w_data[k][j][i];
                        velocities[3].Value = w_data[k+1][j][i];

                        velocities[0].X1 = zw[ k ];
                        velocities[1].X1 = zw[k+1];
                        velocities[2].X1 = zw[ k ];
                        velocities[3].X1 = zw[k+1];

                        velocities[0].X2 = yw[j-1];
                        velocities[1].X2 = yw[j-1];
                        velocities[2].X2 = yw[ j ];
                        velocities[3].X2 = yw[ j ];

                        vel_transpose.X1 = zv[k];
                        vel_transpose.X2 = yv[j];

                        w_at_vgrid[k][j][i] = MyMath_interpolate_quantity(velocities, vel_transpose, 4);
                    } /* if v-status */

                } /* for i */
            } /* for j */
        } /* for k */
        ierr = DAVecRestoreArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
        ierr = DAVecRestoreArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

        ierr = DAVecRestoreArray(grid->DA_3D, (*v->G_u_transposed), (void ***)&u_at_vgrid); PETScErrAct(ierr);
        ierr = DAVecRestoreArray(grid->DA_3D, (*v->G_w_transposed), (void ***)&w_at_vgrid); PETScErrAct(ierr);
        break;

    case 'w':

        ierr = DAVecGetArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
        ierr = DAVecGetArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);

        ierr = DAVecGetArray(grid->DA_3D, (*w->G_u_transposed), (void ***)&u_at_wgrid); PETScErrAct(ierr);
        ierr = DAVecGetArray(grid->DA_3D, (*w->G_v_transposed), (void ***)&v_at_wgrid); PETScErrAct(ierr);
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {

                    /********************************************************************/
                    /* Find u at w_grid */
                    /*
  u2------u3
   |  w   |
  u0------u1

     */

                    if (Grid_get_w_status(grid, i, j, k) != BOUNDARY) {

                        if ( (params->inflow) && (i == 0) ){

                            velocities[0].Value = u->G_inflow[k-1][j];
                            velocities[2].Value = u->G_inflow[k  ][j];
                        } else {

                            velocities[0].Value = u_data[k-1][j][i];
                            velocities[2].Value = u_data[k][j][i];
                        }

                        if ( (params->outflow) && (i == wNI-1) ){

                            velocities[1].Value = u->G_outflow[k-1][j];
                            velocities[3].Value = u->G_outflow[k][ j ];
                        } else {

                            velocities[1].Value = u_data[k-1][j][i+1];
                            velocities[3].Value = u_data[k][j][i+1];
                        }

                        velocities[0].X1 = xu[ i ];
                        velocities[1].X1 = xu[i+1];
                        velocities[2].X1 = xu[ i ];
                        velocities[3].X1 = xu[i+1];

                        velocities[0].X2 = zu[k-1];
                        velocities[1].X2 = zu[k-1];
                        velocities[2].X2 = zu[ k ];
                        velocities[3].X2 = zu[ k ];

                        vel_transpose.X1 = xw[i];
                        vel_transpose.X2 = zw[k];

                        u_at_wgrid[k][j][i] = MyMath_interpolate_quantity(velocities, vel_transpose, 4);

                        /* Find v at w_grid */

                        /*
  v2------v3
   |  w   |
  v0------v1

     */
                        velocities[0].Value = v_data[k-1][j][i];
                        velocities[1].Value = v_data[k][j][i];
                        velocities[2].Value = v_data[k-1][j+1][i];
                        velocities[3].Value = v_data[k][j+1][i];

                        velocities[0].X1 = zv[k-1];
                        velocities[1].X1 = zv[ k ];
                        velocities[2].X1 = zv[k-1];
                        velocities[3].X1 = zv[ k ];

                        velocities[0].X2 = yv[ j ];
                        velocities[1].X2 = yv[ j ];
                        velocities[2].X2 = yv[j+1];
                        velocities[3].X2 = yv[j+1];

                        vel_transpose.X1 = zw[k];
                        vel_transpose.X2 = yw[j];

                        v_at_wgrid[k][j][i] = MyMath_interpolate_quantity(velocities, vel_transpose, 4);
                    } /* if w_status */

                } /* for i*/
            } /*for j*/
        } /* for k*/
        ierr = DAVecRestoreArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
        ierr = DAVecRestoreArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);

        ierr = DAVecRestoreArray(grid->DA_3D, (*w->G_u_transposed), (void ***)&u_at_wgrid); PETScErrAct(ierr);
        ierr = DAVecRestoreArray(grid->DA_3D, (*w->G_v_transposed), (void ***)&v_at_wgrid); PETScErrAct(ierr);
        break;
    } /* switch */
}
/***************************************************************************************************/
