/* Note: Function names are in alphabetical order */
#include "definitions.h"
#include "DataTypes.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Levelset.h"
#include "Display.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* This function allocates memory for the levelset function */
LevelsetType *Levelset_create(MAC_grid *grid) {

    int ierr;
    LevelsetType *new_levelset;

    new_levelset = (LevelsetType *)calloc(1, sizeof(LevelsetType));
    Memory_check_allocation(new_levelset);

    /* phi data */
    ierr = DACreateGlobalVector(grid->DA_3D, &new_levelset->G_data); PETScErrAct(ierr);
    ierr = DACreateLocalVector(grid->DA_3D, &new_levelset->L_data); PETScErrAct(ierr);

    /* Duplicate other vectors from existing Global and Local vector */
    /* phi data_old */
    ierr = VecDuplicate(new_levelset->G_data, &new_levelset->G_data_old); PETScErrAct(ierr);

    /* initial levelset function: phi_0(x,y,z). To remain unchanged */
    ierr = VecDuplicate(new_levelset->G_data, &new_levelset->G_data_0); PETScErrAct(ierr);
    ierr = VecDuplicate(new_levelset->L_data, &new_levelset->L_data_0); PETScErrAct(ierr);

    /* global vectors: First order upwind derivatives */
    ierr = VecDuplicate(new_levelset->G_data, &new_levelset->G_Dx); PETScErrAct(ierr);
    ierr = VecDuplicate(new_levelset->G_data, &new_levelset->G_Dy); PETScErrAct(ierr);
    ierr = VecDuplicate(new_levelset->G_data, &new_levelset->G_Dz); PETScErrAct(ierr);

    /* local vectors: First order upwind derivatives */
    ierr = VecDuplicate(new_levelset->L_data, &new_levelset->L_Dx); PETScErrAct(ierr);
    ierr = VecDuplicate(new_levelset->L_data, &new_levelset->L_Dy); PETScErrAct(ierr);
    ierr = VecDuplicate(new_levelset->L_data, &new_levelset->L_Dz); PETScErrAct(ierr);

    /* Total number of the grid points including the half cell added to the end */
    new_levelset->NX = grid->NX;
    new_levelset->NY = grid->NY;
    new_levelset->NZ = grid->NZ;

    /* Total number of grid points excluding the half cell */
    new_levelset->NI = grid->NI;
    new_levelset->NJ = grid->NJ;
    new_levelset->NK = grid->NK;

    return new_levelset;
}
/***************************************************************************************************/

/* This function releases allocated memory for the levelset type */
void Levelset_destroy(LevelsetType *phi) {

    int ierr;

    ierr = VecDestroy(phi->G_data); PETScErrAct(ierr);
    ierr = VecDestroy(phi->L_data); PETScErrAct(ierr);

    ierr = VecDestroy(phi->G_data_0); PETScErrAct(ierr);
    ierr = VecDestroy(phi->L_data_0); PETScErrAct(ierr);

    ierr = VecDestroy(phi->G_data_old); PETScErrAct(ierr);

    ierr = VecDestroy(phi->G_Dx); PETScErrAct(ierr);
    ierr = VecDestroy(phi->G_Dy); PETScErrAct(ierr);
    ierr = VecDestroy(phi->G_Dz); PETScErrAct(ierr);

    ierr = VecDestroy(phi->L_Dx); PETScErrAct(ierr);
    ierr = VecDestroy(phi->L_Dy); PETScErrAct(ierr);
    ierr = VecDestroy(phi->L_Dz); PETScErrAct(ierr);

    free(phi);
}
/***************************************************************************************************/

/* This function integrates the initial levelset function in fictitous time to get the sign distance function*/
void Levelset_reinitialize(LevelsetType *phi, MAC_grid *grid, int max_iterations) {

    int ierr;
    int converged = NO;
    int c=0;
    int flag;
    double dt;
    short int Need_iteration = YES;
    /* For the cases where the initial sdf is analytical solution, no iteration is needed */
    if (!Need_iteration) {
        converged = YES;
    }

    //PetscPrintf(PCW, "Levelset.c/ done-11\n");
    /* copy data_0 into data */
    ierr = VecCopy(phi->G_data_0, phi->G_data); PETScErrAct(ierr);
    //PetscPrintf(PCW, "Levelset.c/ done-12\n");

    /* update ghost nodes for phi_0. Since phi_0 is constant in time, we only do this procedure once */
    Communication_update_ghost_nodes(&grid->DA_3D, &phi->G_data_0, &phi->L_data_0, 'I');
    //PetscPrintf(PCW, "Levelset.c/ done-13\n");

    /* find "dt" fictitious time step based on the grid size */
    dt = Levelset_compute_fictitious_time(grid);
    //PetscPrintf(PCW, "Levelset.c/ done-14\n");

    /* Not integrate in time to reach steady solution or maximum number of iterations have reached */
    while (!converged) {

        /* update ghost nodes for phi */
        Communication_update_ghost_nodes(&grid->DA_3D, &phi->G_data, &phi->L_data, 'I');
        //PetscPrintf(PCW, "Levelset.c/ done2\n");

        /* first compute the derivatives Dx, Dy and Dz using simple upwind method */
        Levelset_compute_derivatives(phi, grid);
        //PetscPrintf(PCW, "Levelset.c/ done3\n");

        /* update the derivative vectors to get the ghost node values */
        Communication_update_ghost_nodes(&grid->DA_3D, &phi->G_Dx, &phi->L_Dx, 'I');
        Communication_update_ghost_nodes(&grid->DA_3D, &phi->G_Dy, &phi->L_Dy, 'I');
        Communication_update_ghost_nodes(&grid->DA_3D, &phi->G_Dz, &phi->L_Dz, 'I');


        /* now, integrate in time (fictious time) for one time step using Godunov method */
        Levelset_do_Godunov(phi, grid, dt);

        //Display_DA_3D_data(phi->G_data, grid, params, "phi after", 'c');
        //getchar();



        //VecView(phi->G_data, 0);
        //getchar();
        /* One could use NORM1 or NORM_INF for the other types of error */
        flag = Levelset_is_converged(phi, grid, NORM2);
        if ( (c > max_iterations+1) || (flag) ) {

            converged = YES;
        } /* if */

        /* copy phi into phi_old */
        Levelset_store_old_data(phi);

        PetscPrintf(PCW, "Levelset.c/ counter:%d\n", c);

        //getchar();
        c++;
    } /*  while */

    /* update the last planes */
    Levelset_update_last_planes(phi, grid);


}
/***************************************************************************************************/

/* This function computes the first divided differences, Dx, Dy and Dz */
/* For any given node (i,j,k) the indexing is such that; Dx+ = Dx(i,j,k), Dx- = Dx(i-1,j,k); */
void Levelset_compute_derivatives(LevelsetType *phi, MAC_grid *grid) {

    double ***phi_;
    double ***Dx, ***Dy, ***Dz;
    double a_dx, a_dy, a_dz;
    double a_dEta, a_dXi, a_dPhi;
    double *metric_xu, *metric_yv, *metric_zw;
    int NI, NJ, NK;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, j, k;

    /* For non-uniform grid formulation */
    a_dEta  = 1.0/grid->dEta;
    a_dXi   = 1.0/grid->dXi;
    a_dPhi  = 1.0/grid->dPhi;
    /* metric coefficients used for nonuniform grid */
    metric_xu = grid->metric_xu;
    metric_yv = grid->metric_yv;
    metric_zw = grid->metric_zw;

    /* grid points excluding half cell */
    NI = phi->NI;
    NJ = phi->NJ;
    NK = phi->NK;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* get the local version of phi_ */
    ierr = DAVecGetArray(grid->DA_3D, phi->L_data, (void ***)&phi_); PETScErrAct(ierr);

    /* get the global versions of Dx, Dy and Dz */
    ierr = DAVecGetArray(grid->DA_3D, phi->G_Dx, (void ***)&Dx); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, phi->G_Dy, (void ***)&Dy); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, phi->G_Dz, (void ***)&Dz); PETScErrAct(ierr);

    /* go through all the nodes and compute first order Dx, Dy and Dz */
    /* Note that for any phi(i,j,k), Dx(i,j,k) is Dx+, i.e. Dx(i+1/2,j,k) */
    /* Same rule applies in other directions */
    /* One can impose boundary conditions for the ones at (i==0) and (i==NI-1) */
    for (k=Ks; k<Ke; k++) {

        a_dz = 1.0/grid->dz_c[k];
        for (j=Js; j<Je; j++) {

            a_dy = 1.0/grid->dy_c[j];
            for (i=Is; i<Ie; i++) {

                a_dx = 1.0/grid->dx_c[i];
                if (i < NI-1) {

                    //Dx[k][j][i] = (phi_[k][j][i+1] - phi_[k][j][i]) * a_dx;
                    Dx[k][j][i] = (phi_[k][j][i+1] - phi_[k][j][i]) * metric_xu[i] * a_dEta;
                } else {

                    //Dx[k][j][i] = (phi_[k][j][i] - phi_[k][j][i-1]) * a_dx;
                    Dx[k][j][i] = (phi_[k][j][i] - phi_[k][j][i-1]) * metric_xu[i-1] * a_dEta;
                } /* else */

                if ( j < NJ-1) {

                    Dy[k][j][i] = (phi_[k][j+1][i] - phi_[k][j][i]) * metric_yv[j] * a_dXi;
                    //Dy[k][j][i] = (phi_[k][j+1][i] - phi_[k][j][i]) * a_dy;
                } else {

                    Dy[k][j][i] = (phi_[k][j][i] - phi_[k][j-1][i]) * metric_yv[j-1] * a_dXi;
                    //Dy[k][j][i] = (phi_[k][j][i] - phi_[k][j-1][i]) * a_dy;
                } /* else */

                if ( k < NK-1) {

                    Dz[k][j][i] = (phi_[k+1][j][i] - phi_[k][j][i]) * metric_zw[k] * a_dPhi;
                    //Dz[k][j][i] = (phi_[k+1][j][i] - phi_[k][j][i]) * a_dz;
                } else {

                    Dz[k][j][i] = (phi_[k][j][i] - phi_[k-1][j][i]) * metric_zw[k-1] * a_dPhi;
                    //Dz[k][j][i] = (phi_[k][j][i] - phi_[k-1][j][i]) * a_dz;
                } /* else */

            } /* for i */
        } /* for j */
    } /* for k */

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, phi->L_data, (void ***)&phi_); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_Dx, (void ***)&Dx); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_Dy, (void ***)&Dy); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_Dz, (void ***)&Dz); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function updates the levelset function from old time step phi_n to phi_{n+1} using Godunov method. */
/* It uses an improved method to preserve the exact location of the interface */
/* See paper, "A Remark on Computing Distance Functions" Authors: Russo, Giovanni; Smereka, Peter. JCP (2000) */
void Levelset_do_Godunov(LevelsetType *phi, MAC_grid *grid, double dt) {

    double ***p, ***p0;
    double ***Dx, ***Dy, ***Dz;
    double *xc, *yc, *zc;
    double dx, dy, dz;
    double Cx, Cy, Cz;
    double sign_0;
    double Ds;
    double ap, am, bp, bm, cp, cm, dp, dm, ep, em, fp, fm;
    double Dxm, Dxp, Dym, Dyp, Dzm, Dzp;
    double tx, ty, tz;
    double H;
    double p_p, p_e, p_w, p_n, p_s, p_b, p_f;
    double Dphi_dx, Dphi_dy, Dphi_dz;
    double Dt, tol = 0.01;
    double D;
    int index;
    int NI, NJ, NK;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, j, k;
    double DeltaX, DeltaY, DeltaZ;


    /* grid coordinates */
    xc = grid->xc;
    yc = grid->yc;
    zc = grid->zc;

    /* For non-uniform grid formulation */
    Cx  = 1.0/(grid->dEta);
    Cy  = 1.0/(grid->dXi);
    Cz  = 1.0/(grid->dPhi);


    /* grid points excluding the added half cell */
    NI = phi->NI;
    NJ = phi->NJ;
    NK = phi->NK;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    /* exclude the last half cell added to save computation time */
    Ie = min(grid->G_Ie, NI);
    Je = min(grid->G_Je, NJ);
    Ke = min(grid->G_Ke, NK);

    /* get the global version of phi */
    ierr = DAVecGetArray(grid->DA_3D, phi->G_data, (void ***)&p); PETScErrAct(ierr);

    /* get the local version of phi_0: constant through time */
    ierr = DAVecGetArray(grid->DA_3D, phi->L_data_0, (void ***)&p0); PETScErrAct(ierr);

    /* get the local versions of Dx, Dy and Dz */
    ierr = DAVecGetArray(grid->DA_3D, phi->L_Dx, (void ***)&Dx); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, phi->L_Dy, (void ***)&Dy); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, phi->L_Dz, (void ***)&Dz); PETScErrAct(ierr);

    /* go through all the nodes and compute first order Dx, Dy and Dz */
    /* Note that for any phi(i,j,k), Dx(i,j,k) is Dx+, i.e. Dx(i+1/2,j,k) */
    /* Same rule applies in other directions */
    /* One can impose boundary conditions for the ones at (i==0) and (i==NI-1) */
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                dx = grid->dx_c[i];

                sign_0 = p0[k][j][i] / sqrt (p0[k][j][i]*p0[k][j][i] + dx*dx);

                index = max(0,i-1);
                Dxp = Dx[k][j][i];
                Dxm = Dx[k][j][index];

                index = max(0,j-1);
                Dyp = Dy[k][j][i];
                Dym = Dy[k][index][i];

                index = max(0,k-1);
                Dzp = Dz[k][j][i];
                Dzm = Dz[index][j][i];

                ap = max(Dxm, 0.0);
                am = min(Dxm, 0.0);

                bp = max(Dxp, 0.0);
                bm = min(Dxp, 0.0);

                cp = max(Dym, 0.0);
                cm = min(Dym, 0.0);

                dp = max(Dyp, 0.0);
                dm = min(Dyp, 0.0);

                ep = max(Dzm, 0.0);
                em = min(Dzm, 0.0);

                fp = max(Dzp, 0.0);
                fm = min(Dzp, 0.0);

                if (sign_0 > 0.0) {

                    tx = max(ap*ap, bm*bm);
                    ty = max(cp*cp, dm*dm);
                    tz = max(ep*ep, fm*fm);

                } else {

                    tx = max(am*am, bp*bp);
                    ty = max(cm*cm, dp*dp);
                    tz = max(em*em, fp*fp);
                } /* else */
                H = sign_0 * (sqrt(tx + ty + tz) - 1.0);


                /* special treatment for the nodes next to the boundary */
                /* This is based on the paper, "A remark on Computing Distance Functions" by Russo, Giovanni; Smereka, Peter JCP (2000) */
                p_p = p0[k][j][i];

                /* Find the phi_0 of the neighboring nodes */
                index = min(i+1, NI-1);
                p_e = p0[k][j][index];
                DeltaX = xc[index];
                
                index = min(j+1, NJ-1);
                p_n = p0[k][index][i];
                DeltaY = yc[index];

                index = min(k+1, NK-1);
                p_f = p0[index][j][i];
                DeltaZ = zc[index];
                
                index = max(i-1, 0);
                p_w = p0[k][j][index];
                DeltaX -= xc[index];
                
                index = max(j-1, 0);
                p_s = p0[k][index][i];
                DeltaY -= yc[index];

                index = max(k-1, 0);
                p_b = p0[index][j][i];
                DeltaZ -= zc[index];
                
                /* if one of these cases happens, the node is right next to the interface */
                /* We keep the distance D of the current node to the interface using a Taylor series expansion */
                if (
                        ( (p_p*p_w) < 0.0) ||
                        ( (p_p*p_e) < 0.0) ||
                        ( (p_p*p_n) < 0.0) ||
                        ( (p_p*p_s) < 0.0) ||
                        ( (p_p*p_b) < 0.0) ||
                        ( (p_p*p_f) < 0.0) ) {

                    /* 2nd order centeral difference approximation of the Dx, Dy and Dz at node (i,j,k) */
                    if (fabs(DeltaX) > tol) {

                        Dphi_dx = (p_e - p_w) /DeltaX;
                    } else {

                        Dphi_dx = 0.0;
                    } /* else */
                    if (fabs(DeltaY) > tol) {

                        Dphi_dy = (p_n - p_s) /DeltaY;
                    } else {

                        Dphi_dy = 0.0;
                    } /* else */

                    if (fabs(DeltaZ) > tol) {

                        Dphi_dz = (p_f - p_b) /DeltaZ;
                    } else {

                        Dphi_dz = 0.0;
                    } /* else */

                    Dt = sqrt(Dphi_dx*Dphi_dx + Dphi_dy*Dphi_dy + Dphi_dz*Dphi_dz);
                    if (Dt < tol) {

                        Dt = tol;
                    } /* to avoid division by zero */
                    /* Distance of the current node from the interface */
                    D  = p_p / Dt;

                    /* Use the original definition of sign function */
                    /* For some reason the "smoothed" defition of the function does not work for the nodes right next to the boundary */
                    if (fabs(sign_0) < 1.0e-8) {

                        sign_0 = 0.0;
                    } else if (sign_0 > 0.0) {

                        sign_0 = 1.0;
                    } else {

                        sign_0 = -1.0;
                    }  /* else */

                    dx = grid->dx_c[i];
                    dy = grid->dy_c[j];
                    dz = grid->dz_c[k];

                    Ds = sqrt(dx*dx + dy*dy + dz*dz);

                    H = (sign_0 * fabs(p[k][j][i]) - D) / Ds;

                } /* if boundary node */

                /* Euler integration in time */
                p[k][j][i] += -H*dt;
            } /* for i */
        } /* for j */
    } /* for k */
    /* restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, phi->L_data_0, (void ***)&p0); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_data, (void ***)&p); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->L_Dx, (void ***)&Dx); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->L_Dy, (void ***)&Dy); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->L_Dz, (void ***)&Dz); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function computes the error in the entire domain. One should define which error to be computed */
/* NORM_1, NORM_2 or NORM_INF */
int Levelset_is_converged(LevelsetType *phi, MAC_grid *grid, int which_error) {

    double ***p, ***p_old;
    double dx, dy, dz;
    double W_E=0.0;
    double E1=0.0, E2=0.0, E_inf=0.0;
    double tol=1.0e-5; /* convergence tolerance */
    double error;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, j, k;
    int mpi_flag;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* get the global version of phi_ */
    ierr = DAVecGetArray(grid->DA_3D, phi->G_data, (void ***)&p); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, phi->G_data_old, (void ***)&p_old); PETScErrAct(ierr);

    for (k=Ks; k<Ke; k++) {

        dz = grid->dz_c[k];
        for (j=Js; j<Je; j++) {

            dy = grid->dy_c[j];
            for (i=Is; i<Ie; i++) {

                dx = grid->dx_c[i];
                error = fabs(p[k][j][i] - p_old[k][j][i]);
                //PetscPrintf(PCW, "Levelset.c/ error(i,j,k)=(%d,%d,%d) p:%f p_old:%f\n", i, j, k, p[k][j][i], p_old[k][j][i]);
                E1 += error * dx*dy*dz; /* L1 norm of the error */
                E2 += error * error * dx*dy*dz; /* L2 norm of the error */
                if (error > E_inf) { /* L_infinity norm of the error */

                    E_inf = error;
                } /* if */
            } /* for i */
        } /* for j */
    } /* for k */

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_data, (void ***)&p); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_data_old, (void ***)&p_old); PETScErrAct(ierr);

    /* Now, add all the errors on all the processors (for E1 and E2) or find the Maximum (for E_inf) */
    switch (which_error) {

    case NORM1:
        mpi_flag = MPI_Allreduce( &E1, &W_E, 1, MPI_DOUBLE, MPI_SUM, PCW);
        break;
    case NORM2:
        mpi_flag = MPI_Allreduce( &E2, &W_E, 1, MPI_DOUBLE, MPI_SUM, PCW);
        W_E = sqrt(W_E); /* first add E2's and then take the sqrt of the summation */
        break;
    case NORM_INF:
        mpi_flag = MPI_Allreduce( &E_inf, &W_E, 1, MPI_DOUBLE, MPI_MAX, PCW);
        break;
    default:
        PetscPrintf(PCW, "Levelset.c/ Error in convergence. Unknown type for the \"error\" \n");

    } /* switch */

    PetscPrintf(PCW, "Levelset.c/ Convergence Error:%2.14f\n", W_E);
    if (W_E < tol) {

        return YES;
    } /* converged */
    return NO;
}
/***************************************************************************************************/

/* This function stores the old value for phi, i.e. copying "data" into "data_old". */
void Levelset_store_old_data(LevelsetType *phi) {

    int ierr;

    /* copy data into data_old */
    ierr = VecCopy(phi->G_data, phi->G_data_old); PETScErrAct(ierr);
}
/***************************************************************************************************/

/* This function computes the fictitious time based on the minimum grid size */
double Levelset_compute_fictitious_time(MAC_grid *grid) {

    double dx, dy, dz;
    double dx_min = 1000.0;
    double dy_min = 1000.0;
    double dz_min = 1000.0;
    double d_min, dt;
    int i, j, k;

    /* min(dx) */
    for (i=0; i<grid->NI; i++) {

        dx = grid->dx_c[i];
        if (dx < dx_min) {

            dx_min = dx;
        } /* if */
    } /* for i*/

    /* min(dy) */
    for (j=0; j<grid->NJ; j++) {

        dy = grid->dy_c[j];
        if (dy < dy_min) {

            dy_min = dy;
        } /* if */
    } /* for j*/

    /* min(dz) */
    for (k=0; k<grid->NK; k++) {

        dz = grid->dz_c[k];
        if (dz < dz_min) {

            dz_min = dz;
        } /* if */
    } /* for k*/

    /* find the minimum dx (or dy or dz) in the whole domain */
    d_min = min(dx_min, dy_min);
    d_min = min(d_min, dz_min);

    /* fictitious time: for stability always less than d_min */
    dt = 0.5*d_min;

    return dt;
}
/***************************************************************************************************/

/* This function copies the values of the last interior phi to the last plane (added cells to the very end). 
This might not do anything useful */
void Levelset_update_last_planes(LevelsetType *phi, MAC_grid *grid) {


    double ***p, ***p_old;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int NX, NY, NZ;
    int i, j, k;
    int ierr;

    NX = grid->NX;
    NY = grid->NY;
    NZ = grid->NZ;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* get the global version of phi_ */
    ierr = DAVecGetArray(grid->DA_3D, phi->G_data, (void ***)&p); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, phi->G_data_old, (void ***)&p_old); PETScErrAct(ierr);

    if (Ie == NX) {

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                p[k][j][Ie-1]     = p[k][j][Ie-2];
                p_old[k][j][Ie-1] = p_old[k][j][Ie-2];
            } /* for j */
        } /* for k */
    } /* if Ie */


    if (Je == NY) {

        for (k=Ks; k<Ke; k++) {
            for (i=Is; i<Ie; i++) {

                p[k][Je-1][i]     = p[k][Je-2][i];
                p_old[k][Je-1][i] = p_old[k][Je-2][i];
            } /* for i */
        } /* for k */
    } /* if Je */

    if (Ke == NZ) {

        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                p[Ke-1][j][i]     = p[Ke-2][j][i];
                p_old[Ke-1][j][i] = p_old[Ke-2][j][i];
            } /* for i*/
        } /* for j */
    } /* if Ke */

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_data, (void ***)&p); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, phi->G_data_old, (void ***)&p_old); PETScErrAct(ierr);


}
