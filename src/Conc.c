#include "definitions.h"
#include "DataTypes.h"
#include "Solver.h"
#include "Conc.h"
#include "Velocity.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Immersed.h"
#include "Display.h"
#include "Output.h"
#include "gvg.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Function implementations */
/**********************************************************************************************/
/* This function allocates enough memory for the concentration structure based on parameters defined in "*params" */
Concentration *Conc_create(MAC_grid *grid, Parameters *params, int conc_index) {

    Concentration *new_conc;
    int NX, NY, NZ, NT;
    int NX_cell, NY_cell, NZ_cell;
    int ierr;

    /* Number of real physical cells */
    NX_cell = params->NX;
    NY_cell = params->NY;
    NZ_cell = params->NZ;

    /* Add one cell in each direction to the far most end */
    NX = NX_cell + 1; /* Number of grid points in x-direction (half cell added in the far most direction) */
    NY = NY_cell + 1; /* Number of grid points in y-direction (half cell added in the far most direction)*/
    NZ = NZ_cell + 1; /* Number of grid points in z-direction (half cell added in the far most direction)*/
    NT = NX*NY*NZ;

    /* We extend the Field region by one column to the left and right and also one row below and above the regular region. Reason: ASK BRENDON!!*/
    new_conc = (Concentration *)malloc(sizeof(Concentration));
    Memory_check_allocation(new_conc);

    /* index of the associated concentraion field */
    new_conc->conc_index = conc_index;

    /***************************************/
    /* Now, create the vectors that hold the local and global data for all the 3D properties */
    ierr = DACreateGlobalVector(grid->DA_3D, &new_conc->G_data); PETScErrAct(ierr);
    ierr = DACreateLocalVector(grid->DA_3D, &new_conc->L_data); PETScErrAct(ierr);

    ierr = VecDuplicate(new_conc->G_data, &new_conc->G_data_old); PETScErrAct(ierr);

    /* Convective terms: uddx+vddy+wddz */
    ierr = VecDuplicate(new_conc->G_data, &new_conc->G_conv); PETScErrAct(ierr);


    /* Total concentration: To be used in v-momentum equation. Always store the total concentration in c[0]->G_c_total */
    if ( (params->NConc > 1) && (conc_index == 0) ){

        ierr = VecDuplicate(new_conc->G_data, &new_conc->G_c_total); PETScErrAct(ierr);
        ierr = VecDuplicate(new_conc->L_data, &new_conc->L_c_total); PETScErrAct(ierr);

    } /* else, only one concentration field c[0] */


    new_conc->G_outflow     = NULL;
    new_conc->G_outflow_old = NULL;
    new_conc->G_inflow      = NULL;

    if (params->outflow) {

        /* create the array on the processors which contain the domain at the end of the regular box */
        if (grid->G_Ie == grid->NX) {

            /* Note that we create the whole array on the processor but the processor operates
   only on the part which contains the G_data, i.e. (Js to Je)*(Ks to Ke). This is done
   since a 2D array does not need much memory (Also, this way makes life a lot easier :) ) */

            new_conc->G_outflow     = Memory_allocate_2D_double(NY, NZ);
            new_conc->G_outflow_old = Memory_allocate_2D_double(NY, NZ);
        } /* if */
    } /* if */

    /* create inflow */
    if (params->inflow) {

        /* create the array on the processors which contain the domain at the beginning of the regular box */
        if (grid->G_Is == 0) {

            /* Note that we create the whole array on the processor but the processor operates
   only on the part which contains the G_data, i.e. (Js to Je)*(Ks to Ke). This is done
   since a 2D array does not need much memory (Also, this way makes life a lot easier :) ) */

            new_conc->G_inflow = Memory_allocate_2D_double(NY, NZ);
        } /* if */
    } /* if */


    /* concentration averaged height (only a function of x) */
    if (params->ave_height_output) {

        new_conc->G_ave_height_x = Memory_allocate_1D_double(NX);
        new_conc->W_ave_height_x = Memory_allocate_1D_double(NX);

    }/* if */

    if (params->sedim_rate_output) {

        new_conc->G_sed_rate_bottom = Memory_allocate_2D_double(NX, NZ);
        new_conc->W_sed_rate_bottom = Memory_allocate_2D_double(NX, NZ);
    } else {

        new_conc->G_sed_rate_bottom = NULL;
        new_conc->W_sed_rate_bottom = NULL;
    }

    /* Index of left-bottom-back corner on current processor for the DA data layout */
    new_conc->G_Is = grid->G_Is;
    new_conc->G_Js = grid->G_Js;
    new_conc->G_Ks = grid->G_Ks;

    /* Index of right-top-front corner on current processor for the DA data layout */
    new_conc->G_Ie = grid->G_Ie;
    new_conc->G_Je = grid->G_Je;
    new_conc->G_Ke = grid->G_Ke;

    /* Index of left-bottom-back corner on current processor for the DA data layout including ghost nodes	*/
    new_conc->L_Is = grid->L_Is;
    new_conc->L_Js = grid->L_Js;
    new_conc->L_Ks = grid->L_Ks;

    /* Index of right-top-front corner on current processor for the DA data layout including ghost nodes*/
    new_conc->L_Ie = grid->L_Ie;
    new_conc->L_Je = grid->L_Je;
    new_conc->L_Ke = grid->L_Ke;


    /* Global number of cells */
    new_conc->NX = NX;
    new_conc->NY = NY;
    new_conc->NZ = NZ;
    new_conc->NT = NT;

    /* The last interior index of conc cell excluding the half cell */
    new_conc->NI = NX-1;
    new_conc->NJ = NY-1;
    new_conc->NK = NZ-1;

    if (params->hindered_settling) {

    }

    if (params->resuspension == YES) {

    }

    /* Type of conc field: SALINITY, PARTICLE, TEMPERATURE */
    new_conc->Type   = params->conc_type[conc_index];

    /* Volume Fraction: phi = (rho_bulk - rho_H2O)/(rho_i - rho_H2O) * 1.0    */
    /* Max volume fraction of Conc conc_index */
    new_conc->Phi    = params->Phi[conc_index];

    /* Corresponing Peclet Number */
    new_conc->Pe     = params->Pe[conc_index];

    /* PARTICLE: Settling Speed */
    if (new_conc->Type == PARTICLE) {

        new_conc->v_settl0       = params->V_s0[conc_index];

        /* V_p)y-direction = V_s + V_fluid */
        ierr = VecDuplicate(new_conc->G_data, &new_conc->G_v_particle); PETScErrAct(ierr);
        ierr = VecDuplicate(new_conc->L_data, &new_conc->L_v_particle); PETScErrAct(ierr);

        if (params->hindered_settling) {

            ierr = VecDuplicate(new_conc->G_data, &new_conc->G_v_settl); PETScErrAct(ierr);
        }

    } else { /* To make sure, */

        new_conc->v_settl0   = 0.0;
    }

    /* particle deposit height */
    if ( (new_conc->Type == PARTICLE) && (params->deposit_height_output) ) {

        new_conc->G_deposit_height = Memory_allocate_2D_double(NX, NZ);
        new_conc->W_deposit_height = Memory_allocate_2D_double(NX, NZ);

        new_conc->G_deposit_height_dumped = Memory_allocate_2D_double(NX, NZ);
        new_conc->W_deposit_height_dumped = Memory_allocate_2D_double(NX, NZ);

    } else {

        new_conc->G_deposit_height = NULL;
        new_conc->W_deposit_height = NULL;

        new_conc->G_deposit_height_dumped = NULL;
        new_conc->W_deposit_height_dumped = NULL;
    } /* else */


    /* {rhs}={b} vector, */
    ierr = DACreateGlobalVector(grid->DA_3D_Lsys, &new_conc->G_b); PETScErrAct(ierr);

    /* Create parallel sparse matrix for the LHS of the matrix */
    ierr = DAGetMatrix(grid->DA_3D_Lsys, MATMPIAIJ, &new_conc->A); PETScErrAct(ierr);

    /* vector used to update the diagonal part of matrix A */
    /* Get the pointer to it so we can save memory */
    new_conc->lsys_diagonal = &grid->lsys_diagonal;

    return new_conc;
}
/**********************************************************************************************/

/* This function releases the allocated memory for velocity structure. */
void Conc_destroy(Concentration *c, Parameters *params, int conc_index) {

    int ierr;
    int NZ = c->NZ;

    ierr = VecDestroy(c->G_data); PETScErrAct(ierr);
    ierr = VecDestroy(c->L_data); PETScErrAct(ierr);
    ierr = VecDestroy(c->G_data_old); PETScErrAct(ierr);
    ierr = VecDestroy(c->G_conv); PETScErrAct(ierr);

    ierr = VecDestroy(c->G_b); PETScErrAct(ierr);

    if ( (params->NConc > 1) && (conc_index == 0) ) {

        ierr = VecDestroy(c->G_c_total); PETScErrAct(ierr);
        ierr = VecDestroy(c->L_c_total); PETScErrAct(ierr);
    } /* if */
    if (c->Type == PARTICLE) {

        ierr = VecDestroy(c->G_v_particle); PETScErrAct(ierr);
        ierr = VecDestroy(c->L_v_particle); PETScErrAct(ierr);
        if (params->hindered_settling) {

            ierr = VecDestroy(c->G_v_settl); PETScErrAct(ierr);
        }
    } /* if */

    /* release the memory for outflow and inflow (if allocated) */
    if (c->G_outflow != NULL) {

        Memory_free_2D_array(NZ, (void **)c->G_outflow);
        Memory_free_2D_array(NZ, (void **)c->G_outflow_old);
    } /* if */
    if (c->G_inflow != NULL) {

        Memory_free_2D_array(NZ, (void **)c->G_inflow);
    }/* if */

    if (params->ave_height_output) {

        free(c->G_ave_height_x);
        free(c->W_ave_height_x);
    }/* if */

    if (c->G_deposit_height != NULL) {

        Memory_free_2D_array(NZ, (void **)c->G_deposit_height);
        Memory_free_2D_array(NZ, (void **)c->W_deposit_height);
    } /* if */

    if (c->G_deposit_height_dumped != NULL) {

        Memory_free_2D_array(NZ, (void **)c->G_deposit_height_dumped);
        Memory_free_2D_array(NZ, (void **)c->W_deposit_height_dumped);
    } /* if */

    if (c->G_sed_rate_bottom != NULL) {

        Memory_free_2D_array(NZ, (void **)c->G_sed_rate_bottom);
        Memory_free_2D_array(NZ, (void **)c->W_sed_rate_bottom);
    }
    ierr = MatDestroy(c->A); PETScErrAct(ierr);
    ierr = KSPDestroy(c->solver); PETScErrAct(ierr);

    free(c);

}
/***************************************************************************************************/


/* This function updates the diagonal part of matrix [A] for velocities by the new value of "dt" */
/* It will only update the nodes that are NOT immersed nodes */
void Conc_modify_diagonal(Concentration *c, MAC_grid *grid, double dt, double dt_old) {

    int ierr;
    double mod;
    double mod_insert;
    double ***diagonal;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;


    if (fabs(dt_old) > 1.0e-7) {

        mod = 1.0/dt - 1.0/dt_old;
    }
    else {

        mod = 1.0/dt;
    }

    /* To save time, add "mod" to the diagonal only if dt != dt_old */
    if (fabs(mod) > 1.0e-7) {

        /* Adds constant 'mod' to the diagonal part of matrix A*/
        /* Start index of bottom-left-back corner on current processor */
        Is = grid->G_Is;
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* End index of top-right-front corner on current processor */
        Ie = grid->G_Ie;
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        /* Get the diagonal vector array */
        ierr = DAVecGetArray(grid->DA_3D, (*c->lsys_diagonal), (void ***)&diagonal); PETScErrAct(ierr);

        /* Go through the entire domain and only insert the delta_t term for the nodes which are not immersed nodes */
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {

                    if (Grid_get_c_status(grid, i, j, k) == IMMERSED) {

                        mod_insert = 0.0;
                    } else {

                        mod_insert = mod;
                    } /* else */
                    diagonal[k][j][i] = mod_insert;
                } /* for i */
            } /* for j */
        } /* for k */

        ierr = DAVecRestoreArray(grid->DA_3D, (*c->lsys_diagonal), (void ***)&diagonal); PETScErrAct(ierr);
        /* Now, add the diagonal part of the lsys matrix to update 1/dt */
        /* Maybe this could be done in a mroe efficietn way which requires less memory */
        ierr = MatDiagonalSet(c->A, (*c->lsys_diagonal), ADD_VALUES); PETScErrAct(ierr);
    } /* if */


}
/***************************************************************************************************/

/* This function creates the matrices and solver for conc linear set of equations */
void Conc_setup_lsys_accounting_geometry(Concentration *c, MAC_grid *grid, Parameters *params) {

    /* Now, generate the LHS matrix for the Transport equation using the 2nd order finite difference method */
    Conc_form_LHS_matrix(c, grid, params);

    /* Now, using the given LHS matrix {and some known parameters like grid}, appropriate KSP solver (and
preconditioner) is generated using Petsc package */
#ifdef MEMORY_PROFILING
    PetscMallocGetCurrentUsage(&mem1);
#endif
    c->solver = Solver_get_concentration_solver(c->A, params);
    int ierr = KSPSetOperators(c->solver, c->A, c->A, SAME_PRECONDITIONER); PETScErrAct(ierr);
    ierr = KSPSetUp(c->solver); PETScErrAct(ierr);
#ifdef MEMORY_PROFILING
    PetscMallocGetCurrentUsage(&mem2);
    c_solver_mem += (int)(mem2 - mem1);
    c_mem += (int)(mem2 - mem1);
#endif

    if (params->output_lsys_matrix) {

        char filename[FILENAME_MAX];
        sprintf(filename, "Resume_conc%d_matrix.bin.gz", c->conc_index);
        Output_petsc_mat(&(c->A), filename, BINARY);

        PetscPrintf(PCW, "Conc.c/ conc%d-matrix has been written to %s successfully.\n", c->conc_index, filename);
    } /* if */

}
/**********************************************************************************************/

/* This function generates LHS matrix of the linear system ignoring the temporal diagonal term */
void Conc_form_LHS_matrix( Concentration *c, MAC_grid *grid,  Parameters *params) {

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
    double Pe;
    double V_s0;
    double coef, b;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    MatStencil cols[2*STENCIL], row;
    double values[2*STENCIL]; /* STENCIL = 7 */
    int nnz; /* number of nonzeros on each row of the Matrix */
    int status;
    Immersed *c_immersed;
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

    c_immersed = grid->c_immersed;


    NI  = c->NI;
    NJ  = c->NJ;
    NK  = c->NK;

    /* Constant settling speed of particle. Zero for temperature and Salinity */
    V_s0 = c->v_settl0;

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

    Pe   = c->Pe;
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
    Cx  = -1.0/(Pe * dEta * dEta);
    Cy  = -1.0/(Pe * dXi * dXi);
    Cz  = -1.0/(Pe * dPhi * dPhi);

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
                status = Grid_get_c_status(grid, i, j, k);
                //PetscPrintf(PCW, "Conc.c/ (i,j,k)=(%d,%d,%d) status:%d\n", i, j, k, status);
                if (status == FLUID) {

                    /* get the metric coefficients for non-uniform grid formulation */
                    /* For more detail, check the manual */
                    metric_px = metric_xc[i];     /* (dx/deta)^-1 at (i,j,k) */
                    metric_py = metric_yc[j];     /* (dy/dxi)^-1 at (i,j,k)      */
                    metric_pz = metric_zc[k];     /* (dz/dphi)^-1 at (i,j,k)     */

                    metric_ex = metric_xu[i+1];   /* (dx/deta)^-1 at (i+1/2,j,k)   */
                    metric_wx = metric_xu[i];     /* (dx/deta)^-1 at (i-1/2,j,k) */

                    metric_ny = metric_yv[j+1];   /* (dy/dxi)^-1 at (i,j+1/2,k)  */
                    metric_sy = metric_yv[j];     /* (dy/dxi)^-1 at (i,j-1/2,k)  */

                    metric_bz = metric_zw[k];   /* (dz/dphi)^-1 at (i,j,k+1/2)  */
                    metric_fz = metric_zw[k+1];     /* (dz/dphi)^-1 at (i,j,k-1/2)  */

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
                    if (k > 0) {

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

                    } /* else */

                    if (i > 0) { /*west neighbour is fluid, insert aw into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = aw;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i-1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;
                    } else { /* apply Neumann B.C. on the wall */

                        /* If we have inflow, do not impose Neumann B.C. */
                        if (! params->inflow ){

                            ap += aw;
                        }
                    } /* else */

                    /* East neighbor */
                    if (i < NI-1) {/*East neighbour is fluid, insert as into matrix*/

                        /* the value of nonzero element in the matrix */
                        values[nnz] = ae;

                        /* 3D index of nonzero element in current row */
                        cols[nnz].i = i+1;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                    } else { /* apply Neumann B.C. on the wall */

                        /* If outflow, do not impose Neumann B.C. */
                        if ( ! ( (params->outflow) && (i==NI-1) ) ) {

                            ap += ae;
                        } /* if */
                    } /* else */

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

                        /* the no flux boundary condition on the top wall, i.e. c.U_s + 1/Pe.dc/dy = 0 */
                        if ( ( j == NJ-1 ) && (c->Type == PARTICLE) ) {

                            b    = 2.0/(Pe * dXi) * metric_ny;
                            coef = an * (b - fabs(V_s0)) / (b + fabs(V_s0) );
                        } else { /* Neumann B.C. */

                            coef = an;
                        }

                        ap += coef;
                    } /* else */

                    /* Check for front neighbor */
                    if (k < NK-1) { /*Front neighbour is fluid, insert af into matrix*/

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

                    /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                    ierr = MatSetValuesStencil(c->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);

                } else { /* current node is Immersed node or the extra-half cell or Boundary node */

                    if (status == IMMERSED) {

                        local_index = (k - Ks_g) * ny_g * nx_g + (j - Js_g) * nx_g + (i - Is_g);
                        /* Get the global index (row number in the Matrix) for current node at (i,j,k) */
                        row_global_index = global_indices[local_index];
                        im_global_indices[nnz] = global_indices[local_index];
                        im_values[nnz] = 1.0;
                        nnz++;

                        /* Now, go through the fluid nodes, find the global index in the matrix, copy the value and use MatSetValues() instead of MatSetValuesStencil. This is done b/c some of the fluid node might be on the neighboring processor and out of the range of the ghost layers. */
                        ib_index = Immersed_get_ib_global_index(c_immersed, i, j, k);
                        ib_node  = Immersed_get_ib_node(c_immersed, ib_index);
                        for (g=0; g<ib_node->n_fluid; g++) {

                            ii  = ib_node->fluid_index[g].x_index;
                            jj  = ib_node->fluid_index[g].y_index;
                            kk  = ib_node->fluid_index[g].z_index;

                            local_index = (kk - Ks_g) * ny_g * nx_g + (jj - Js_g) * nx_g + (ii - Is_g);

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
                            } /* else */
                        } /* for */

                        /* All interpolation coefficients are scaled with a constant factor. Better convergence */
                        int s;
                        for (s=0; s<nnz; s++) im_values[s] *= ap_solid;

                        ierr = MatSetValues(c->A, 1, &row_global_index, nnz, im_global_indices, im_values, INSERT_VALUES); PETScErrAct(ierr);
                    } else {

                        values[nnz] = ap_solid; /* the value is not important */
                        cols[nnz].i = i;
                        cols[nnz].j = j;
                        cols[nnz].k = k;
                        nnz++;

                        /* Insert the values for the current node, maximum 7 point stencil. Minimum 1 (for solid nodes) */
                        ierr = MatSetValuesStencil(c->A, 1, &row, nnz, cols, values, INSERT_VALUES); PETScErrAct(ierr);
                    } /* else BOUNDARY node */

                } /* SOLID or BOUNDARY node */

            } /* for: i*/
        } /*for: j */
    }/* for :k*/

    ierr = MatAssemblyBegin(c->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);
    ierr = MatAssemblyEnd(c->A, MAT_FINAL_ASSEMBLY); PETScErrAct(ierr);

}
/**************************************************************************************************************/

/* This function adds the settling velocity to "v" velocity of the fluid calculated at the center of each cell. */
void Conc_set_particle_v_velocity(Concentration *c, Velocity *v, Parameters *params) {

    int ierr;

    /* Only add the settling speed to the "particle" concentration fields */
    if (c->Type == PARTICLE) {

        if (params->hindered_settling) {

            /* use Petsc vec summation */

            /* v_particle = v_data_bc + v_settl */
            ierr = VecWAXPY(c->G_v_particle, 1.0, v->G_data_bc, c->G_v_settl); PETScErrAct(ierr);

        } /* if hindered */
        else { /* Just add constant settling speed */

            /* First, copy
     v_particle <-- v_data_bc */
            ierr = VecCopy(v->G_data_bc, c->G_v_particle); PETScErrAct(ierr);

            /* Now, add v_settl0 to the the v_bc to find particle's velocity */
            ierr = VecShift(c->G_v_particle, c->v_settl0); PETScErrAct(ierr);

        } /* else */
    } /* if particle */
    /* else, concentraion velocity is equal to fluid velocity field */
}
/************************************************************************************************/
/* This Function, calculates the settling speed due to hindered settling effect for the whole domain.*/
void Conc_compute_particle_settling_speed(Concentration *c, MAC_grid *grid) {

    double ***v_settl, ***conc;
    int i, j, k;
    double phi, v_settl0;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Partcile Volume fraction */
    phi  = c->Phi;

    /* Particle Settling Speed in absence of hindered settling */
    v_settl0 = c->v_settl0;

    ierr = DAVecGetArray(grid->DA_3D, c->G_v_settl, (void ***)&v_settl); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

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

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {

                    v_settl[k][j][i] = Conc_settling_speed_function(conc[k][j][i], phi, v_settl0);

                }
                else { /* Check it for bottom boundary!!!!! */

                    v_settl[k][j][i] = 0.0; /* Is this true */
                }
            } /* for i*/
        } /* for j*/
    } /* for k*/

    ierr = DAVecRestoreArray(grid->DA_3D, c->G_v_settl, (void ***)&v_settl); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

}
/************************************************************************************************/

/* This Function, calculates the settling speed due to hindered settling effect. */
double Conc_settling_speed_function(double conc, double phi_max, double V_s0) {

    /* Various correlations could be used:

   1- Batchelers' model: 1977 (Theoretical equation)
        V_s / V_s0 = (1.0 - 6.55 * phi ) :: Valid for phi <= 0.02

   2- Ham, Homsy: 1988 (experimental relation)
        V_s / V_s0 = (1.0 - 4.0 * phi + 8 * phi^2) :: Valid for phi <= 0.10

   3- Zaki, Richardson type: 1954 (experimental relation)
        V_s / V_s0 = (1.0 - phi)^n :: Valid for phi <= 0.60 (To be verified) ::    2.39 <= n <= 5.
        Typical value n = 4.
*/	       
    double phi = conc * phi_max; /* particle volume fraction */
    double n   = 4.0;

    /*return ( V_s0 * ( 1.0 - 6.55 * phi) );  */
    /*return ( V_s0 * ( 1.0 - 4.0 * phi + 8.0 * phi * phi ) ); */
    return ( V_s0 * pow( 1.0 - phi , n) );
}
/************************************************************************************************/

/* This function solves the linear set of equations for concentration using the Petsc KSP solver */
int Conc_solve(Concentration *c) {

    int iters;
    int ierr;

    /* Call Petsc KSP-solver routine to solve linear system. */
    ierr = KSPSetOperators(c->solver, c->A, c->A, SAME_NONZERO_PATTERN); PETScErrAct(ierr);

    /* The solution is stored in data array */
    ierr = KSPSolve(c->solver, c->G_b, c->G_data); PETScErrAct(ierr);

    char name[5];
    sprintf(name, "c%d", c->conc_index);
    (void)Solver_is_converged(c->solver, name);

    ierr = KSPGetIterationNumber(c->solver, &iters); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, " ** Conc soln took %d iters.\n", iters); PETScErrAct(ierr);

    return iters;
}
/***************************************************************************************************/


/* This function computes the RHS vector {b} for the conc linear set of equations */
void Conc_set_RHS(Concentration *c, MAC_grid *grid, Parameters *params, double dt) {

    int i, j, k;
    int NI, NJ, NK;
    int ierr;
    double ***conv;
    double ***data;
    double ***rhs_vec;
    double **c_outflow, **c_inflow;
    double *metric_xc, *metric_xu, *metric_yv, *metric_yc;
    double ae, aw;
    double dEta, dXi, dPhi;
    double Cx, Cy, Cz;
    double a_dt;
    double metric_px, metric_ex, metric_wx;
    double rhs_value, value;
    double Pe;
    int Is, Js, Ks;
    int Ie, Je, Ke;


    /* Get regular data array for c data */
    ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&data); PETScErrAct(ierr);

    /* Get convective terms global data  */
    ierr = DAVecGetArray(grid->DA_3D, c->G_conv, (void ***)&conv); PETScErrAct(ierr);

    /* Now, got the RHS vector on current processor */
    ierr = DAVecGetArray(grid->DA_3D_Lsys, c->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    /*********************/
    a_dt = 1.0/dt;

    c_outflow = c->G_outflow;
    c_inflow  = c->G_inflow;

    NI  = c->NI;
    NJ  = c->NJ;
    NK  = c->NK;

    metric_xc = grid->metric_xc;
    metric_yc = grid->metric_yc;
    metric_xu = grid->metric_xu;
    metric_yv = grid->metric_yv;

    dEta = grid->dEta;
    dXi  = grid->dXi;
    dPhi = grid->dPhi;

    Pe = c->Pe;

    Cx = -1.0/(Pe*dEta*dEta);
    Cy = -1.0/(Pe*dXi*dXi);
    Cz = -1.0/(Pe*dPhi*dPhi);

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

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /*only include if point is fluid*/


                    rhs_value = data[k][j][i]*a_dt - conv[k][j][i];

                    if (params->inflow) {

                        if (i == 0) {

                            metric_px = metric_xc[i];   /* (dx/deta)^-1 at (i,j)     */
                            metric_wx = metric_xu[i];   /* (dx/deta)^-1 at (i-1/2,j) */
                            aw        = -Cx * metric_px * metric_wx;

                            value      = aw * c_inflow[k][j];
                            rhs_value +=value;
                        } /* if */
                    } /* if inflow */

                    if (params->outflow) {

                        /* only for the last node before to the outflow boundary */
                        if (i == NI-1) {

                            /* get the metric coeeffients for non-uniform grid formulation */
                            metric_px = metric_xc[i];   /* (dx/deta)^-1 at (i,j)     */
                            metric_ex = metric_xu[i];   /* (dx/deta)^-1 at (i+1/2,j) */
                            ae        = -Cx * metric_px * metric_ex;

                            value      = ae * c_outflow[k][j];
                            rhs_value += value;
                        }
                    }

                    if (params->hindered_settling) {

                    }

                    if (params->resuspension) {

                    }


                } else { /* Just a solid node */

                    rhs_value = 0.0;
                }/* else */

                rhs_vec[k][j][i] = rhs_value;

            } /* for i*/
        } /* for j*/
    }/* for k*/

    //getchar();
    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_conv, (void ***)&conv); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_Lsys, c->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function copiies the value of c->data into c->data_old. This is done for the time integration
purposes. */
void Conc_store_old_data(Concentration *c) {

    int ierr;

    /* copy data into data_old */
    ierr = VecCopy(c->G_data, c->G_data_old); PETScErrAct(ierr);
}
/************************************************************************************************/

/* This function finds the conc data by applying the linear RK- Averaging (see TVD rk3 method) */
void Conc_rk_average(Concentration *c, double old_frac, double new_frac) {


    int ierr;

    /* use Petsc vec product and summation */
    /*  data[k][j][i] = old_frac * data_old[k][j][i] + new_frac * data[k][j][i] */

    ierr = VecAXPBY(c->G_data, old_frac, new_frac, c->G_data_old); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function generates the ghost nodes for calculation of the w-convective terms using ENO scheme */
void Conc_ENO_ghost_cells(ENO_Scheme *ENO, Concentration *c, MAC_grid *grid, Parameters *params, Indices
                          start_cell, Indices end_cell, char which_direction) {

    double *D0, *Position, *Velocity;
    double *xu, *xv, *yv, *zw;
    double *metric_yv;
    int **interface_y_index;
    double **inflow;
    double **outflow;
    double **E_s;
    double x_start, x_end;
    double y_start, y_end;
    double z_start, z_end;
    double dXi;
    int NI, NJ, NK;
    int ghost_index, interior_index, ex;
    int ENO_order;
    int Nmax_local;
    int i_start, j_start, k_start;
    int i_end, j_end, k_end;
    double Pe;
    double conc_value, coef, Vs;

    outflow   = c->G_outflow;
    inflow    = c->G_inflow;

    D0         = ENO->D0;
    Position   = ENO->Position;
    Velocity   = ENO->Velocity;
    ENO_order  = ENO->ENO_order;

    NI = c->NI;
    NJ = c->NJ;
    NK = c->NK;

    xv = grid->xv;
    xu = grid->xu;
    yv = grid->yv;
    zw = grid->zw;

    metric_yv = grid->metric_yv;

    dXi = grid->dXi;

    Pe = c->Pe;

    /* (Resuspension) Erosion flux */
    E_s = c->E_s;

    /* Bottom geometry first j-index */
    interface_y_index = grid->interface_y_index;

    /* Get the indices of the first fluid node */
    i_start = start_cell.x_index;
    j_start = start_cell.y_index;
    k_start = start_cell.z_index;

    /* Get the indices of the last fluid node */
    i_end   = end_cell.x_index;
    j_end   = end_cell.y_index;
    k_end   = end_cell.z_index;


    /* udcdx */
    /* Apply Neumann B.C. except for the nodes neighboring inflow and outflow */
    if (which_direction == 'x') {

        if (i_start == 0) {

            /* First, for the left side of the array D0, corresponds to the left side of the domain*/
            x_start = xu[i_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*x_start - Position[interior_index];

                if (params->inflow) {

                    D0[ghost_index]       = inflow[k_start][j_start];
                } else {

                    D0[ghost_index]       = D0[ENO_order];//D0[interior_index];
                }
            } /* for ex */
        } /* if i_start */

        /* Now, for the right side of the array D0*/
        Nmax_local = i_end - i_start + 2*ENO_order + 1;
        if (i_end == NI-1) {

            x_end = xu[i_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*x_end - Position[interior_index];
                if ( params->outflow ){

                    D0[ghost_index]       = outflow[k_end][j_end];
                } else {

                    D0[ghost_index]       = D0[Nmax_local - ENO_order - 1];//D0[interior_index];
                }
            } /* for */
        } /* if i_end */
    } /* if x */

    /* vdcdy */
    /* Apply Neumann B.C. except for the nodes neighboring bottom and top boundaries (in case of
 particle-sedimentation */
    if (which_direction == 'y') {

        /* First, for the left side of the array D0, corresponds to the left side of the domain*/
        if (j_start == 0) {

            y_start = yv[j_start];

            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*y_start - Position[interior_index];

                /* Fix this for the case of resuspension */
                D0[ghost_index]       = D0[ENO_order];//D0[interior_index];
            } /* for ex */

        } /* j_start */

        /* Now, for the right side of the array D0, corresponds to the top boundary. No variable geomtery. */
        Nmax_local = j_end - j_start + 2*ENO_order + 1;
        if (j_end == NJ-1) {
            y_end = yv[j_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*y_end - Position[interior_index];
                if (c->Type == PARTICLE) {

                    if ( (params->sedimentation) ){

                        /* Have to check (fix) this (for hindered settling) */
                        coef        = 2.0/(Pe * dXi) * metric_yv[j_end+1];
                        Vs          = fabs(c->v_settl0);
                        conc_value  = (coef - Vs) / (coef + Vs) * D0[Nmax_local - ENO_order - 1];

                        D0[ghost_index]       = conc_value;

                    } else {

                        D0[ghost_index]       = D0[Nmax_local - ENO_order - 1];//D0[interior_index];
                    } /* else */
                } /* if Particle */
            } /* for ex */
        } /* if j_end */
    } /* if y*/

    /* for wdcdz */

    /* Apply Neumann B.C. */
    if (which_direction == 'z') {

        if (k_start == 0) {

            z_start = zw[k_start];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = ex;
                interior_index = 2*ENO_order - ex - 1;

                Position[ghost_index] = 2.0*z_start - Position[interior_index];
                D0[ghost_index]       = D0[ENO_order];//D0[interior_index];

            } /* for ex */
        } /* if k_start */

        /* Now, for the right side of the array D0, corresponds to the top boundary. No variable geomtery. */
        Nmax_local = k_end - k_start + 2*ENO_order + 1;

        if (k_end == NK-1) {

            z_end = zw[k_end+1];
            for (ex=0; ex<ENO_order; ex++) {

                ghost_index    = Nmax_local - ex - 1;
                interior_index = Nmax_local - 2*ENO_order + ex;

                Position[ghost_index] = 2.0*z_end - Position[interior_index];
                D0[ghost_index]       = D0[Nmax_local - ENO_order - 1];//D0[interior_index];

            } /* for */
        } /* if */
    } /* if z */
}
/***************************************************************************************************/

/* This function stores the old data for the outflow concentration. This is done for time integration purposes */
void Conc_store_old_outflow(Concentration *c) {

    int j, k;
    int Js, Ks;
    int Je, Ke;
    double **outflow, **outflow_old;

    /* Only apply for the processors including the yz plane at the end */
    if (c->G_Ie == c->NX) {

        outflow     = c->G_outflow;
        outflow_old = c->G_outflow_old;

        /* start index on current processor */
        Js = c->G_Js;
        Ks = c->G_Ks;

        /* end index on current processor */
        Je = c->G_Je;
        Ke = c->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                outflow_old[k][j] = outflow[k][j];
            } /* for j */
        } /* for k */
    } /* if */

}
/**********************************************************************************************/

/* This function takes the rk-average for concentration outflow */
void Conc_rk_average_outflow(Concentration *c, double old_frac, double new_frac) {

    int j, k;
    int Js, Ks;
    int Je, Ke;
    double **outflow, **outflow_old;

    /* Only apply for the processors including the yz plane at the end */
    if (c->G_Ie == c->NX) {

        outflow     = c->G_outflow;
        outflow_old = c->G_outflow_old;

        /* start index on current processor */
        Js = c->G_Js;
        Ks = c->G_Ks;

        /* end index on current processor */
        Je = c->G_Je;
        Ke = c->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                outflow[k][j] = old_frac * outflow_old[k][j] + new_frac * outflow[k][j];
            } /* for j */
        } /* for k */
    } /* if */
}
/***************************************************************************************************/

/* This function computes the total concentration for more than 1 concentration field */
/* This is to be used in y-momentum equation */
/* Note that conc_total is always stored in c[0] */
void Conc_compute_total_concentration(Concentration **c, Parameters *params) {

    int NConc, iconc;
    int ierr;

    NConc = params->NConc;

    /* Copy first concentration field data into "conc_toal" */
    /* Store the total concentration in c[0] */

    /* c_total = c[0]->data; */
    ierr = VecCopy(c[0]->G_data, c[0]->G_c_total); PETScErrAct(ierr);

    /* c_total = c_total * Conc_alpha[0] */
    ierr = VecScale(c[0]->G_c_total, params->Conc_alpha[0]); PETScErrAct(ierr);

    for (iconc=1; iconc<NConc; iconc++) {

        /* compute total concentration (to be used in y-momentum equation */
        /* c_total = c_total + alpha_i * c[i]->data */
        ierr = VecAXPY(c[0]->G_c_total, params->Conc_alpha[iconc], c[iconc]->G_data); PETScErrAct(ierr);
    } /* for */

}
/***************************************************************************************************/

/*This function sets the initial lock profile */
void Conc_initialize ( Concentration *c, MAC_grid *grid, Parameters *params ) {

    int i, j, k;
    int lock_smooth;
    double x_fr, y_fr, z_fr;
    double *xc, *yc, *zc;
    double ***conc;
    double c_test;
    double Re;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    short int lock_cube;
    double xb_s, xb_e, yb_s, yb_e, zb_s, zb_e;


    /* Get concentration data  */
    ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

    Re = params->Re;

    x_fr = params->x_fr;
    y_fr = params->y_fr;
    z_fr = params->z_fr;

    xc   = grid->xc;
    yc   = grid->yc;
    zc   = grid->zc;

    lock_smooth = params->lock_smooth;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;


    /* Assume the lock is a block. It could be located anywhere within the domain. */
    /* One has to define the coordinates of the bottom-left-back and top-right-front corner */
    lock_cube   = NO;
    xb_s = 0.0;
    xb_e = 0.2*params->Lx;
    yb_s = 0.0;
    yb_e = 0.5*params->Ly;
    zb_s = 0.3*params->Lz;
    zb_e = 0.7*params->Lz;

    /* To perturb the interface */
    short int AddPerturbation = YES;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) { /*only set initial profile if cell is fluid*/

                    if (params->sim_type == LOCK_EXCHANGE) {

                        if (lock_smooth) {

                            c_test = 1.0/4.0*
                                    (1.0-erf((xc[i]-x_fr)*sqrt(Re))) *
                                    (1.0-erf((yc[j]-y_fr)*sqrt(Re))) ;

                            c_test = min(c_test, 1.0);
                            c_test = max(c_test, 0.0);

/* Perturb the interface */
                            if (AddPerturbation) {
                                if ( fabs(x_fr - xc[i]) < 0.05) {
                                    if (c_test >= 0.9){

                                        double c_perturb = MyMath_random(0.1);
                                        int sign        = (int)MyMath_random(10);
                                        if (sign < 5) {

                                            c_perturb *= -1.0;
                                        }
                                        c_test += c_perturb;
                                    }
                                    c_test = min(c_test, 1.0);
                                }
                            } /* if perturbed */

                            conc[k][j][i] = c_test;

                        } else { /*otherwise just do step-function concentration profile*/

                            if ( (xc[i] <= x_fr) && (yc[j] <= y_fr) && (zc[k] <=z_fr)) {

                                conc[k][j][i] = 1.0;
                            }
                            else {

                                conc[k][j][i] = 0.0;
                            }
                        } /* else */
                    }
                }
            }
        }
    }

    if (lock_cube) {

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {
                for (i=Is; i<Ie; i++) {

                    if (!Grid_get_c_status(grid, i, j, k) == FLUID) { /*only set initial profile if cell is fluid*/

                        if ( (xc[i] >= xb_s) && (xc[i] <= xb_e) &&
                             (yc[j] >= yb_s) && (yc[j] <= yb_e) &&
                             (zc[k] >= zb_s) && (zc[k] <= zb_e) ) {

                            conc[k][j][i] = 1.0;
                        } else {

                            conc[k][j][i] = 0.0;
                        }
                    } /* if grid */
                } /* for i*/
            } /* for j*/
        } /* for k */
    } /* if lock cube */



    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

}
/**********************************************************************************************/

/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/* Post processing functions */
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
/* This function computes total suspended mass within the domain */
void Conc_compute_total_suspended_mass(Concentration *c, MAC_grid *grid) {

    /* This would only work if we have uniform grid */
    /* int ierr = VecSum(c->G_data, &(c->W_total_susp_mass)); PETScErrAct(ierr); */

    double ***conc = NULL;
    int ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    double c_tot = 0.0;

    int i, j, k;
    for (k=Ks; k<Ke; k++) {
        double dz = grid->dz_c[k];
        for (j=Js; j<Je; j++) {

            double dy = grid->dy_c[j];
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {

                    double dV = dz*dy*grid->dx_c[i];;
                    c_tot += conc[k][j][i]*dV;
                } /* if */
            } /* for i */
        } /* for j */
    } /* for k */

    /* restore array now */
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

    double W_c_tot = 0.0;
    /* now, communicate to find the total amount of suspended mass in the entire domain */
    (void)MPI_Allreduce (&c_tot, &W_c_tot, 1, MPI_DOUBLE, MPI_SUM, PCW);

    c->W_total_susp_mass = W_c_tot;


}
/**********************************************************************************************/

/* This function computes the (concentration based) average height of the flow. It is only a function of x
since it has been averaged in y and z directions */
/* Note that, it is stored only in processor zero */
void Conc_compute_ave_height_x(Concentration **c, MAC_grid *grid, Parameters *params) {


    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, j, k;
    int NConc, iconc;
    int NX;
    double h_temp, coef;
    double dy, dz;
    double *yv, *zw;
    double *dy_c, *dz_c;
    double ***conc;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    yv = grid->yv;
    zw = grid->zw;

    /* Grid dimensions */
    dy_c = grid->dy_c;
    dz_c = grid->dz_c;

    NX = grid->NX;

    coef = 1.0/(params->Ly * params->Lz);

    NConc = params->NConc;
    for (iconc=0; iconc<NConc; iconc++) {

        /* Get concentration data  */
        ierr = DAVecGetArray(grid->DA_3D, c[iconc]->G_data, (void ***)&conc); PETScErrAct(ierr);

        /* At each x location (i-index), find

    h_ave = 1/(Ly*Lz) sum(over j and k or (y,z) indices) {C(k,j) * (dy*dz)} */

        /* This will give us the current height only as a function of x*/
        for (i=Is; i<Ie; i++) {

            h_temp = 0.0;
            for (k=Ks; k<Ke; k++) {

                dz = dz_c[k];
                for (j=Js; j<Je; j++) {

                    if (Grid_get_c_status(grid, i, j, k) == FLUID) {

                        dy = dy_c[j];
                        h_temp += conc[k][j][i]*dy*dz;
                    } /* if */

                } /* for j*/
            } /* for k*/

            c[iconc]->G_ave_height_x[i] = h_temp * coef;
        } /* for i*/

        /* Now, reduce all G_ave_height on processor zero (or all processors) and add them all to find the correct ave_height .*/
        /* Note that the trick is that on all processors, we create the full array of
  ave_height_x[0->NX-1] and set all the parts not located on current processor equal to zero
  (all the time) and just compute the part located on current processor */

        /* Store the final result on processor zero W_ave_height */
        MPI_Reduce( (void *)c[iconc]->G_ave_height_x, (void *)c[iconc]->W_ave_height_x, NX, MPI_DOUBLE, MPI_SUM, MASTER, PCW);

        /* This will send to all processors after computing ave_height */
        //MPI_Allreduce( (void *)c[iconc]->G_ave_height_x, (void *)c[iconc]->W_ave_height_x, NX, MPI_DOUBLE, MPI_SUM, PCW);
        /* Restore array for iconc */
        ierr = DAVecRestoreArray(grid->DA_3D, c[iconc]->G_data, (void ***)&conc); PETScErrAct(ierr);

    } /* for iconc */
}
/**********************************************************************************************/

/* This function finds the front location using spanwise averaged height */
/* It is done only on processor zero */
double Conc_find_front_location(Concentration **c, MAC_grid *grid, Parameters *params, double front_limit) {

    Scalar heights[2];
    Scalar height_front;
    int iconc, NConc;
    int i;
    int NX;
    double *xc;
    double front_location = -1.0;
    double ave_height_t_0;
    double ave_height_t_1;


    if (params->rank == MASTER) {

        xc    = grid->xc;

        NConc = params->NConc;
        NX    = c[0]->NX;

        for (i=NX-2; i>0; i--) {


            ave_height_t_0 = 0.0;
            ave_height_t_1 = 0.0;
            for(iconc=0; iconc<NConc; iconc++) {

                ave_height_t_0 += c[iconc]->W_ave_height_x[i-1] * params->Conc_alpha[iconc] ;
                ave_height_t_1 += c[iconc]->W_ave_height_x[i] * params->Conc_alpha[iconc] ;

            } /*for iconc */
            /* Found front location */
            if ( (ave_height_t_1<= front_limit) && (ave_height_t_0>= front_limit) ) {

                /* for the two point interpolation, Play a little trick since we have
    front limit, but we want to find the location X. */
                /* h0 is found at x0. h1 is found at x1.
       Now, we want to find xf corresponding to the front limit */

                heights[0].X1  = ave_height_t_0;
                heights[1].X1  = ave_height_t_1;

                heights[0].Value = xc[i-1];
                heights[1].Value = xc[i];

                height_front.X1 = front_limit;

                /* Just linear interpolation for the north and south conc to find c_term at the v_grid node */
                front_location = MyMath_interpolate_quantity(heights, height_front, 2);

                c[0]->W_front_location = front_location;
                break;

            } /* if */
        } /* for i*/

    } /* if rank */

    return front_location;

}
/**********************************************************************************************/

/* This function integrates (in time) the deposited sediment height for the particle concentration fields */
/* For now, just use it for cases without hindered settling */
void Conc_integrate_deposited_height(Concentration *c, MAC_grid *grid, double weight_factor, double dt) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, k;
    int j_index;
    int **interface_y_index;
    int NX, NZ;
    double delta_h;
    double dx, dz;
    double *xu, *zw;
    double *dx_c, *dz_c;
    double ***conc, ***conc_old;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    xu = grid->xu;
    zw = grid->zw;

    /* Grid dimensions */
    dx_c = grid->dx_c;
    dz_c = grid->dz_c;

    /* y-index of first fluid node from bottom, (y=0) */
    interface_y_index = grid->interface_y_index;

    /* Get concentration data  */
    ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, c->G_data_old, (void ***)&conc_old); PETScErrAct(ierr);

    NX = grid->NX;
    NZ = grid->NZ;

    /* Now go through 2D array and check if the first fluid node lies on current processor */
    for (k=Ks; k<Ke; k++) {

        dz = dz_c[k];
        for (i=Is; i<Ie; i++) {

            /* Index of first interior fluid node (on the bottom boundary */
            j_index = interface_y_index[k][i];

            dx = dx_c[i];
            /* Check if current prcossor includes the bottom geometry */
            if ( ( j_index >= Js) && (j_index < Je) ){

                /* delta_h added height at each time step: second order in time */
                delta_h = dt * fabs(c->v_settl0) * weight_factor * 0.5 * (conc[k][j_index][i]+conc_old[k][j_index][i]);
                //printf("conc:%f conc_old:%f\n", conc[k][j_index][i], conc_old[k][j_index][i]);
                //printf("dt:%f dx:%f dz:%f v_set:%f weigh:%f Delta_h: %2.12f\n", dt, dx, dz, fabs(c->v_settl0), weight_factor, delta_h);
                //getchar();
                /* Update the height of deposited sediment */
                c->G_deposit_height[k][i] += delta_h;

            } /* if */

        } /* for */
    } /* for k*/

    /* Restore array */
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data_old, (void ***)&conc_old); PETScErrAct(ierr);
}
/**********************************************************************************************/

/* This function communicates all the G_deposit_height amongst processors to find W_deposit height */
/* Result could be stored either on processor zero or all processors */
void Conc_update_world_deposited_height(Concentration **c, MAC_grid *grid, Parameters *params) {

    int iconc, NConc;
    Indices G_s, G_e, W_e;

    NConc = params->NConc;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Ke = grid->G_Ke;

    /*  Since the 2D arrays are assumed to have [Y][X] index order, pass always the min, max indices
 based on this rule */
    /* Start and End indices of the 2D array on the current processor */
    G_s.x_index = Is;
    G_s.y_index = Ks;

    G_e.x_index = Ie;
    G_e.y_index = Ke;

    /* Total number of grid point in the W_ array */
    W_e.x_index = grid->NX;
    W_e.y_index = grid->NZ;

    /* Go through each particle concentration field and Reduce all the 2D arrays on processor zero */
    for (iconc=0; iconc<NConc; iconc++) {

        if (c[iconc]->Type == PARTICLE) {

            /* Communicate to add all the G_deposit_height[] to processor zero for each concentration field */
            Communication_reduce_2D_arrays(c[iconc]->G_deposit_height, c[iconc]->W_deposit_height, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
        } /* if */
    } /* for iconc */

}
/**********************************************************************************************/

/* This function communicates all the G_deposit_height_dumped amongst processors to find W_deposit height_dumped */
/* Result could be stored either on processor zero or all processors */
void Conc_update_world_deposited_height_dumped(Concentration **c, MAC_grid *grid, Parameters *params) {

    int iconc, NConc;
    Indices G_s, G_e, W_e;

    NConc = params->NConc;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Ke = grid->G_Ke;

    /*  Since the 2D arrays are assumed to have [Y][X] index order, pass always the min, max indices
 based on this rule */
    /* Start and End indices of the 2D array on the current processor */
    G_s.x_index = Is;
    G_s.y_index = Ks;

    G_e.x_index = Ie;
    G_e.y_index = Ke;

    /* Total number of grid point in the W_ array */
    W_e.x_index = grid->NX;
    W_e.y_index = grid->NZ;

    /* Go through each particle concentration field and Reduce all the 2D arrays on processor zero */
    for (iconc=0; iconc<NConc; iconc++) {

        if (c[iconc]->Type == PARTICLE) {

            /* Communicate to add all the G_deposit_height[] to processor zero for each concentration field */
            Communication_reduce_2D_arrays(c[iconc]->G_deposit_height_dumped, c[iconc]->W_deposit_height_dumped, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
        } /* if */
    } /* for iconc */

}
/**********************************************************************************************/


/* This function adds all the available particles in the domain to the deposited height. 
   This is useful when the velocity field is small and the only important thing happening in the domain is
   the settling of particles. This would somehow predicts how the deposited height would be after a long run
   calculation */
void Conc_dump_particles(Concentration **c, MAC_grid *grid, Parameters *params) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    int ierr;
    int iconc, NConc;
    int NX, NZ;
    double ***conc;
    double *xu, *yv, *zw;
    double *dx_c, *dy_c, *dz_c;
    double dx, dy, dz;
    double h_temp;
    double weight_factor;

    xu = grid->xu;
    yv = grid->yv;
    zw = grid->zw;

    /* Grid dimensions */
    dx_c = grid->dx_c;
    dy_c = grid->dy_c;
    dz_c = grid->dz_c;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    NX = grid->NX;
    NZ = grid->NZ;

    NConc = params->NConc;

    /* Integrate all the particles in y-direction and add them to the sed_height for the particle fields */
    for (iconc=0; iconc<NConc; iconc++) {


        if (c[iconc]->Type == PARTICLE) {


            /* Get concentration data  */
            ierr = DAVecGetArray(grid->DA_3D, c[iconc]->G_data, (void ***)&conc); PETScErrAct(ierr);

            /* At each location (x,z) find c_t = integral (0 to Ly) of c_dx.dy.dz */

            weight_factor = params->Conc_alpha[iconc];
            /* This will give us the current height only as a function of x*/
            for (k=Ks; k<Ke; k++) {

                dz = dz_c[k];
                for (i=Is; i<Ie; i++) {

                    dx     = dx_c[i];
                    h_temp = 0.0;

                    c[iconc]->G_deposit_height_dumped[k][i] = c[iconc]->G_deposit_height[k][i];
                    for (j=Js; j<Je; j++) {

                        if (Grid_get_c_status(grid, i, j, k) == FLUID) {

                            dy = dy_c[j];
                            h_temp += conc[k][j][i]*dy;
                        } /* if */
                    } /* for j */

                    /* Now add the integrated value for the shrinked value of concentration to the deposited one and store the result into the two dimensional array locally on each processor */
                    c[iconc]->G_deposit_height_dumped[k][i] += h_temp * weight_factor;

                } /* for i*/
            } /* for k*/

            ierr = DAVecRestoreArray(grid->DA_3D, c[iconc]->G_data, (void ***)&conc); PETScErrAct(ierr);
        } /* if particle */

    } /* for iconc*/

    /* Now, update the W_dep... on processor zero */
    Conc_update_world_deposited_height_dumped(c, grid, params);
}
/**********************************************************************************************/

/* This function computes the potential energies for the concentration field:
This potential includes:
1- Active potential energy
2- Active potential eneryg grad
3- Passive potentital energy: particles that have settled out. This
potential energy could not be converted to kinetic energy and is computed just for the cases where we have
variable topography. */
void Conc_compute_potential_energies(Concentration *c, Velocity *v, MAC_grid *grid, Parameters *params) {

    double ***conc=NULL;
    /* Get concentration data  */
    int ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

    /* Grid coordinates */
    double *yc = grid->yc;

    /* Grid dimensions */
    double *dx_c = grid->dx_c;
    double *dz_c = grid->dz_c;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    /* Reference height to define potential energy */
    double y_ref = grid->yv[0];

    /* Add the potential energy due to the deposited particles * on the bottom geometry */
    double G_Ep_passive = 0.0;
    if (c->Type == PARTICLE) {

        double dx, dz;
        int i, k;
        int j_index;
        for (k=Ks; k<Ke; k++) {

            dz = dz_c[k];
            for (i=Is; i<Ie; i++) {

                dx = dx_c[i];

                /* y-index of first interior fluid node */
                j_index = grid->interface_y_index[k][i];

                /* Ony add the contribution if the bottom surface lies on the current processor */
                if ( (j_index >= Js) && (j_index < Je)) {

                    G_Ep_passive += c->G_deposit_height[k][i] * dx*dz * (yc[j_index] - y_ref);

                } /* if */

            } /* for i*/
        }/* for k*/
    } /* if */

    /* Get regular data array for vel data at cell center */
    double ***v_data_bc=NULL;
    ierr = DAVecGetArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc); PETScErrAct(ierr);

    /* Available potential energy in the domain */
    double G_Ep = 0.0;
    double G_Stokes_dissipation_rate = 0.0; /* Stokes dissipation rate associated with the current particle field, integral (c*u_s *dVol) */
    double G_kinetic_potential_conversion_rate = 0.0; /* Potential to kinetic energy conversion rate for the current concentration field */
    int i, j, k;
    for (k=Ks; k<Ke; k++) {

        /* to save computational time, compute 1/dz */
        double dz = grid->dz_c[k];
        for (j=Js; j<Je; j++) {

            double dy = grid->dy_c[j];
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /*only include if point is fluid*/

                    double dx = grid->dx_c[i];

                    /* differntial volume of the current grid */
                    double dV = dx*dy*dz;

                    G_Ep  += conc[k][j][i] * dV * (yc[j]-y_ref);
                    G_Stokes_dissipation_rate += conc[k][j][i] * dV;

                    double v_ = v_data_bc[k][j][i];

                    /* Add to the conversion rate on current processor */
                    G_kinetic_potential_conversion_rate += v_* conc[k][j][i] * (dx*dy*dz);

                } /* if */
            } /* for i*/
        } /* for j*/
    }/* for k*/

    /* restore arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc); PETScErrAct(ierr);

    G_Ep                  *= params->Conc_alpha[c->conc_index];
    G_Stokes_dissipation_rate *= fabs(c->v_settl0) * params->Conc_alpha[c->conc_index];
    G_kinetic_potential_conversion_rate *= params->Conc_alpha[c->conc_index];

    /* Now, communicate to obtain the reduced world values */
    double send[4];
    double recv[4];
    int n_comm = 4;
    send[0] = G_Ep;
    send[1] = G_Ep_passive;
    send[2] = G_Stokes_dissipation_rate;
    send[3] = G_kinetic_potential_conversion_rate;

    (void)MPI_Allreduce( (void *)send, (void *)recv, n_comm, MPI_DOUBLE, MPI_SUM, PCW);

    /* Store the total potential energies associated with the current concentration field */
    c->W_Ep              = recv[0];
    c->W_Ep_passive      = recv[1];
    c->W_Stokes_dissipation_rate = recv[2];
    c->W_kinetic_potential_conversion_rate = recv[3];

    c->W_Ep_passive_rate = c->W_Ep_passive * fabs(c->v_settl0);
    c->W_Ep_rate = c->W_kinetic_potential_conversion_rate + c->W_Ep_passive_rate - c->W_Stokes_dissipation_rate;

}
/**********************************************************************************************/

/* This function resets that nonlocal values of the G_deposit_height[][].
This is to ensure the correct MPI_Allreduce() when updating for the W_deposit_height is performed.
This task is used when we resume the simulations from an older run */
void Conc_reset_nonlocal_deposit_height(Concentration *c, MAC_grid *grid) {

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    int NX = grid->NX;
    int NZ = grid->NZ;

    /* Go through the entire 2D array and reset the values that do not lie within the current processor including the bottom boundary */
    /* This is essential for the global MPI_Allreduce() function when updating the W_deposit_height */
    int i, k;
    for (k=0; k<NZ; k++) {
        for (i=0; i<NX; i++) {

            /* Index of first interior fluid node (on the bottom boundary */
            int j_index = grid->interface_y_index[k][i];

            /* Check if current prcossor includes the bottom geometry */
            if (!  (
                      (j_index >= Js) && (j_index < Je) &&
                      (i >= Is) && (i < Ie) &&
                      (k >= Ks) && (k < Ke)
                   )
                ) {
                c->G_deposit_height[k][i] = 0.0;
            } /* if */
        } /* for i */
    } /* for k */
}
/**********************************************************************************************/

/* This function computes the sedimentation rate for the current particle field */
void Conc_compute_sedimentation_rate(Concentration *c, MAC_grid *grid, Parameters *params) {

    /* Only if it is particle */
    if (c->Type == PARTICLE) {

        /* Start index of bottom-left-back corner on current processor */
        int Is = grid->G_Is;
        int Js = grid->G_Js;
        int Ks = grid->G_Ks;

        /* End index of top-right-front corner on current processor */
        int Ie = grid->G_Ie;
        int Je = grid->G_Je;
        int Ke = grid->G_Ke;

        /* y-index of first fluid node from bottom, (y=0) */
        int **interface_y_index = grid->interface_y_index;

        /* Get concentration data  */
        double ***conc=NULL;
        int ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

        double u_s = fabs(c->v_settl0);
        double G_sed_rate = 0.0;
        /* Now go through 2D array and check if the first fluid node lies on current processor */
        int k, i;
        for (k=Ks; k<Ke; k++) {

            double dz = grid->dz_c[k];
            for (i=Is; i<Ie; i++) {

                /* Index of first interior fluid node (on the bottom boundary */
                int j_index = interface_y_index[k][i];

                double dx = grid->dx_c[i];
                /* Check if current prcossor includes the bottom geometry */
                if ( ( ( j_index >= Js) && (j_index < Je) ) ||
                   ( (j_index == 0) && (Grid_get_c_status(grid, i, j_index, k) == FLUID) ) ) {

                    /* add the sedimentation rate */
                    G_sed_rate += u_s * conc[k][j_index][i] * dx*dz;
                }

            } /* for i*/
        } /* for k*/

        /* Restore array */
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

        double W_sed_rate = 0.0;
        (void)MPI_Allreduce (&G_sed_rate, &W_sed_rate, 1, MPI_DOUBLE, MPI_SUM, PCW);

        c->W_sed_rate = W_sed_rate;
    } /* if particle */ else {

        c->W_sed_rate = 0.0;
    }
}

/* This function gets the sedimenation rate on the bottom boundary */
void Conc_bottom_sed_rate(Concentration *c, MAC_grid *grid, Parameters *params) {

    /* Only if it is particle */
    if (c->Type == PARTICLE) {

        /* Start index of bottom-left-back corner on current processor */
        int Is = grid->G_Is;
        int Js = grid->G_Js;
        int Ks = grid->G_Ks;

        /* End index of top-right-front corner on current processor */
        int Ie = grid->G_Ie;
        int Je = grid->G_Je;
        int Ke = grid->G_Ke;

        /* y-index of first fluid node from bottom, (y=0) */
        int **interface_y_index = grid->interface_y_index;

        /* Get concentration data  */
        double ***conc=NULL;
        int ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

        double u_s = fabs(c->v_settl0);

        /* Now go through 2D array and check if the first fluid node lies on current processor */
        int k, i;
        for (k=Ks; k<Ke; k++) {

            double dz = grid->dz_c[k];
            for (i=Is; i<Ie; i++) {

                /* Index of first interior fluid node (on the bottom boundary */
                int j_index = interface_y_index[k][i];

                double dx = grid->dx_c[i];
                /* Check if current prcossor includes the bottom geometry */
                if ( ( ( j_index >= Js) && (j_index < Je) ) ||
                   ( (j_index == 0) && (Grid_get_c_status(grid, i, j_index, k) == FLUID) ) ) {

                    /* add the sedimentation rate */
                    c->G_sed_rate_bottom[k][i] = u_s * conc[k][j_index][i] * dx*dz;
                } else {

                    c->G_sed_rate_bottom[k][i] = 0.0;
                }

            } /* for i*/
        } /* for k*/

        /* Restore array */
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

        /* Now, update the the global values */
        /*  Since the 2D arrays are assumed to have [Y][X] index order, pass always the min, max indices
     based on this rule */
        /* Start and End indices of the 2D array on the current processor */
        Indices G_s, G_e, W_e;
        G_s.x_index = Is;
        G_s.y_index = Ks;

        G_e.x_index = Ie;
        G_e.y_index = Ke;

        /* Total number of grid point in the W_ array */
        W_e.x_index = grid->NX;
        W_e.y_index = grid->NZ;

        /* Now, sum all the computed shear stresses from all the processors (zero for the ones which don't
     contain the bottom boundary and send the result to processor zero: Store on W_.... */
        Communication_reduce_2D_arrays(c->G_sed_rate_bottom, c->W_sed_rate_bottom, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    }
}

