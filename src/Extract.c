#include "definitions.h"
#include "DataTypes.h"
#include "Boundary.h"
#include "Solver.h"
#include "Conc.h"
#include "Velocity.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Communication.h"
#include "Immersed.h"
#include "Display.h"
#include "Extract.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void Extract_wall_shear(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) {

    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    int ib_index;
    double ***u_data, ***v_data, ***w_data;
    double ***u_data_bc, ***v_data_bc, ***w_data_bc;
    double *xc, *yc, *zc;
    double *yv;
    double **shear_stress;
    double Re;
    Indices G_s, G_e, W_e;
    double V_n;
    double nx, ny, nz;
    double u_int, v_int, w_int;
    double **u_shear, **v_shear, **w_shear;
    PointType *p_int;
    ImmersedNode *ib_node;


    /* allocate memory for the interpolation point */
    p_int = (PointType *)calloc(1, sizeof(PointType));

    /* Update the values of the ghost nodes based on a BOX stencil */
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &u->G_data, &u->L_data_box, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &v->G_data, &v->L_data_box, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &w->G_data, &w->L_data_box, 'I');

    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D_BOX, v->L_data_box, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D_BOX, w->L_data_box, (void ***)&w_data);PETScErrAct(ierr);

    /* Get regular data array for vel data at cell center*/
    ierr = DAVecGetArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc);PETScErrAct(ierr);

    /* Grid coordinates */
    xc = grid->xc;
    yc = grid->yc;
    zc = grid->zc;
    yv = grid->yv;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    Re = params->Re;
    shear_stress = u->G_shear_stress_bottom;

    /* Tangential velocity on the bottom */
    u_shear = u->G_u_shear;
    v_shear = u->G_v_shear;
    w_shear = u->G_w_shear;

    /* Go through the x-z plane and find the shear stress based on the tangential velocity right next to the boundary */
    for (k=Ks; k<Ke; k++) {
        for (i=Is; i<Ie; i++) {

            /* index of first IMMERSED (or FLUID for no solid boundary on the bottom) node */
            j = grid->interface_y_index[k][i];

            double u_tan = 0.0;
            double v_tan = 0.0;
            double w_tan = 0.0;
            double tau   = 0.0;
            double d = 0.0;
            //PetscPrintf(PCW, "Extract.c/ Here1 (i,j,k)=(%d,%d,%d)\n", i, j, k);
            /* Calculate shear stress only if the current processor contains the bottom boundary */
            if ( (j >= Js) && (j < Je) && (j != -1) ) {

                if (Grid_get_c_status(grid, i, j, k) == IMMERSED) {

                    ib_index = Immersed_get_ib_global_index(grid->c_immersed, i, j, k);
                    ib_node  = Immersed_get_ib_node(grid->c_immersed, ib_index);

                    /* normal vector to the boundary */
                    nx = ib_node->n.vx;
                    ny = ib_node->n.vy;
                    nz = ib_node->n.vz;

                    //p_int->x = ib_node->control_point.x + nx*d;
                    //p_int->y = ib_node->control_point.y + ny*d;
                    //p_int->z = ib_node->control_point.z + nz*d;

                    /* coordinates of the node at which the tangential velcoity is computed. */
                    p_int->x = ib_node->image_point1.x;
                    p_int->y = ib_node->image_point1.y;
                    p_int->z = ib_node->image_point1.z;

                    /* distance between interpolation point and surface */
                    d = MyMath_get_point_point_distance(p_int, &ib_node->control_point);

                    /* Use a trilinear interpolation to find the velocity at the interpolation node */
                    u_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, u_data, grid, 'u');
                    v_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, v_data, grid, 'v');
                    w_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, w_data, grid, 'w');

                    /* V = V_t.e_t + V_n.e_n */
                    /* e_n = nx.i + ny.j + nz.k */
                    /* magnitude of the normal velocity vector */
                    V_n = u_int*nx + v_int*ny + w_int*nz;

                    /* Tangential component: V_t = V - V_n */
                    u_tan = u_int - V_n*nx;
                    v_tan = v_int - V_n*ny;
                    w_tan = w_int - V_n*nz;

                } /* if status */
                if ( (j == 0) && (Grid_get_c_status(grid, i, j, k) == FLUID) ){

                    u_tan = u_data_bc[k][j][i];
                    v_tan = 0.0;
                    w_tan = w_data_bc[k][j][i];
                    d = yc[j] - yv[j];

                } /* if */

                tau = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan) / (d * Re);
            } /* if */

            /* Store the shear stress on the current processor */
            shear_stress[k][i] = tau;

            /* Tangential velocity */
            u_shear[k][i] = u_tan;
            v_shear[k][i] = v_tan;
            w_shear[k][i] = w_tan;

        } /* for i*/
    } /* for k*/

    /* Always restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, v->L_data_box, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, w->L_data_box, (void ***)&w_data); PETScErrAct(ierr);

    /* Always restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc); PETScErrAct(ierr);


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

    /* Now, sum all the computed shear stresses from all the processors (zero for the ones which don't
 contain the bottom boundary and send the result to processor zero: Store on W_.... */
    Communication_reduce_2D_arrays(u->G_shear_stress_bottom, u->W_shear_stress_bottom, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);

    /* Get the bottom tangential velocity */
    Communication_reduce_2D_arrays(u->G_u_shear, u->W_u_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(u->G_v_shear, u->W_v_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(u->G_w_shear, u->W_w_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);

    free(p_int);
}
/***************************************************************************************************/

/* This function computes the bottom surfac shear stress using a quadratic interpolation */
void Extract_wall_shear_quadratic(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) {

    /* Update the values of the ghost nodes based on a BOX stencil */
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &u->G_data, &u->L_data_box, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &v->G_data, &v->L_data_box, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &w->G_data, &w->L_data_box, 'I');

    /* Get regular data array for vel data */
    int ierr;
    double ***u_data=NULL, ***v_data=NULL, ***w_data=NULL;
    ierr = DAVecGetArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D_BOX, v->L_data_box, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D_BOX, w->L_data_box, (void ***)&w_data);PETScErrAct(ierr);

    /* Get regular data array for vel data at cell center*/
    double ***u_data_bc=NULL, ***v_data_bc=NULL, ***w_data_bc=NULL;
    ierr = DAVecGetArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc);PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    /* Grid coordinates */
    double *yc = grid->yc;
    double *yv = grid->yv;

    /* Shear stress */
    double **shear_stress = u->G_shear_stress_bottom;

    /* Tangential velocity on the bottom */
    double **u_shear = u->G_u_shear;
    double **v_shear = u->G_v_shear;
    double **w_shear = u->G_w_shear;

    /* Allocate memory for the matrix. Used for quadratric interpolation */
    MatrixType *M = (MatrixType *)calloc(1, sizeof(MatrixType));
    M->size = 3;

    double dVt_dn=0.0;
    double a_Re = 1.0/params->Re;
    /* allocate memory for the interpolation point */
    PointType *p_int;
    p_int = (PointType *)calloc(1, sizeof(PointType));

    int k, i;
    /* Go through the x-z plane and find the shear stress based on the tangential velocity right next to the boundary */
    for (k=Ks; k<Ke; k++) {
        for (i=Is; i<Ie; i++) {

            /* index of first IMMERSED (or FLUID for no solid boundary on the bottom) node */
            int j = grid->interface_y_index[k][i];

            double u_tan = 0.0;
            double v_tan = 0.0;
            double w_tan = 0.0;
            double tau   = 0.0;
            //PetscPrintf(PCW, "Extract.c/ Here1 (i,j,k)=(%d,%d,%d)\n", i, j, k);
            /* Calculate shear stress only if the current processor contains the bottom boundary */
            if ( (j >= Js) && (j < Je) && (j != -1) ) {

                if (Grid_get_c_status(grid, i, j, k) == IMMERSED) {

                    int ib_index = Immersed_get_ib_global_index(grid->c_immersed, i, j, k);
                    ImmersedNode *ib_node  = Immersed_get_ib_node(grid->c_immersed, ib_index);

                    /* Control point: on the solid surface */
                    double d0 = 0.0;
                    double V_t0 = 0.0;

                    /******************************************************************************/
                    /* normal vector to the boundary */
                    double nx = ib_node->n.vx;
                    double ny = ib_node->n.vy;
                    double nz = ib_node->n.vz;

                    //double dtest = min(fabs(grid->dx_c[i]/nx), fabs(grid->dy_c[j]/ny));
                    //dtest = min(dtest, fabs(grid->dz_c[k]/nz));

                    /* coordinates of the node at which the tangential velcoity is computed. */
                    p_int->x = ib_node->image_point1.x;
                    p_int->y = ib_node->image_point1.y;
                    p_int->z = ib_node->image_point1.z;

                    /* Use a trilinear interpolation to find the value of the velocity */
                    double d1 = MyMath_get_point_point_distance(p_int, &ib_node->control_point);
                    double u_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, u_data, grid, 'u');
                    double v_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, v_data, grid, 'v');
                    double w_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, w_data, grid, 'w');

                    /* V = V_t.e_t + V_n.e_n */
                    /* e_n = nx.i + ny.j + nz.k */
                    /* magnitude of the normal velocity vector */
                    double V_n = fabs(u_int*nx + v_int*ny + w_int*nz);

                    /* Tangential component (vector relationship): V_t = V - V_n */
                    u_tan = u_int - V_n*nx;
                    v_tan = v_int - V_n*ny;
                    w_tan = w_int - V_n*nz;
                    double V_t1 = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan);

                    /******************************************************************************/
                    /* Use a trilinear interpolation to find the value of the velocity */
                    /* coordinates of the node at which the tangential velcoity is computed. */
                    p_int->x = ib_node->image_point2.x;
                    p_int->y = ib_node->image_point2.y;
                    p_int->z = ib_node->image_point2.z;

                    double d2 = MyMath_get_point_point_distance(p_int, &ib_node->control_point);
                    u_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, u_data, grid, 'u');
                    v_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, v_data, grid, 'v');
                    w_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, w_data, grid, 'w');

                    /* Check if the second point lies within the current processor range. If not, shift p2 to the midpoint between p1 and p2 */
                    if ( ( fabs(u_int-GVG_INF) < 1.0e-5) ||
                         ( fabs(v_int-GVG_INF) < 1.0e-5) ||
                         ( fabs(w_int-GVG_INF) < 1.0e-5) ) {

                        printf("Extract.c/ Warning! Interpolation point (%f,%f,%f) does not lie within the local range on the current processor.\n Shifting the point\n", p_int->x, p_int->y, p_int->z);
                        p_int->x = 0.5*(ib_node->image_point2.x + ib_node->image_point1.x);
                        p_int->y = 0.5*(ib_node->image_point2.y + ib_node->image_point1.y);
                        p_int->z = 0.5*(ib_node->image_point2.z + ib_node->image_point1.z);

                        d2 = MyMath_get_point_point_distance(p_int, &ib_node->control_point);
                        u_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, u_data, grid, 'u');
                        v_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, v_data, grid, 'v');
                        w_int = MyMath_do_trilinear_inter(p_int->x, p_int->y, p_int->z, w_data, grid, 'w');

                        /* Reduce the shear stress computation to first order */
                        if ( ( fabs(u_int-GVG_INF) < 1.0e-5) ||
                             ( fabs(v_int-GVG_INF) < 1.0e-5) ||
                             ( fabs(w_int-GVG_INF) < 1.0e-5) ) {

                            tau = (V_t1 /d1) * a_Re;
                            /* Store the shear stress on the current processor */
                            shear_stress[k][i] = tau;
                            /* Tangential velocity */
                            u_shear[k][i] = u_tan;
                            v_shear[k][i] = v_tan;
                            w_shear[k][i] = w_tan;

                            continue;
                        }
                    } /* if */

                    /* V = V_t.e_t + V_n.e_n */
                    /* e_n = nx.i + ny.j + nz.k */
                    /* magnitude of the normal velocity vector */
                    V_n = fabs(u_int*nx + v_int*ny + w_int*nz);

                    /* Tangential component (vector relationship): V_t = V - V_n */
                    u_tan = u_int - V_n*nx;
                    v_tan = v_int - V_n*ny;
                    w_tan = w_int - V_n*nz;
                    double V_t2 = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan);

                    /* Now, we employ a quadratic interpolation to find the velocity and ultimately the shear stress on the solid surface at the control point */
                    /* fluid point 1 */
                    M->A[0][0] = d1*d1;
                    M->A[0][1] = d1;
                    M->A[0][2] = 1.0;

                    /* fluid point 2 */
                    M->A[1][0] = d2*d2;
                    M->A[1][1] = d2;
                    M->A[1][2] = 1.0;

                    /* control point: d0=0.0 */
                    M->A[2][0] = d0*d0;
                    M->A[2][1] = d0;
                    M->A[2][2] = 1.0;

                    /* Invert the matrix to find the coefficients */
                    MyMath_GaussJordan_matrix_inv(M);

                    double a = V_t1*M->A_inv[0][0] + V_t2*M->A_inv[0][1] + V_t0*M->A_inv[0][2];
                    double b = V_t1*M->A_inv[1][0] + V_t2*M->A_inv[1][1] + V_t0*M->A_inv[1][2];

                    /* Derivative of the tangential velocity on the solid surface (control point) */
                    dVt_dn = 2.0*a*d0 + b;

                } /* if status */
                if ( (j == 0) && (Grid_get_c_status(grid, i, j, k) == FLUID) ){

                    u_tan = u_data_bc[k][j][i];
                    v_tan = 0.0;
                    w_tan = w_data_bc[k][j][i];
                    double a_d = 1.0/(yc[j] - yv[j]);

                    dVt_dn = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan) * a_d;
                } /* if */

                tau = dVt_dn * a_Re;
            } /* if */

            /* Store the shear stress on the current processor */
            shear_stress[k][i] = tau;

            /* Tangential velocity */
            u_shear[k][i] = u_tan;
            v_shear[k][i] = v_tan;
            w_shear[k][i] = w_tan;

        } /* for i*/
    } /* for k*/
    //getchar();
    /* Always restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, v->L_data_box, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, w->L_data_box, (void ***)&w_data); PETScErrAct(ierr);

    /* Always restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc); PETScErrAct(ierr);


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
    Communication_reduce_2D_arrays(u->G_shear_stress_bottom, u->W_shear_stress_bottom, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);

    /* Get the bottom tangential velocity */
    Communication_reduce_2D_arrays(u->G_u_shear, u->W_u_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(u->G_v_shear, u->W_v_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(u->G_w_shear, u->W_w_shear, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);

    free(p_int);
    free(M);
}
/***************************************************************************************************/


/* This function computes the total vicous dissipation rate due to the fluid motion */
/* e = 2.0/Re int( S_ij S_ij dV) (over entire domain) */
/* For now, only cases without in/outflow are working properly */
void Extract_viscous_dissipation_rate(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) {

    double ***u_data=NULL;
    double ***v_data=NULL;
    double ***w_data=NULL;

    int ierr;
    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D, u->L_data, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->L_data, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->L_data, (void ***)&w_data);PETScErrAct(ierr);

    /* Grid dimensions */
    double *dx_c = grid->dx_c;
    double *dy_c = grid->dy_c;
    double *dz_c = grid->dz_c;

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

    double **u_inflow=u->G_inflow;
    double **u_outflow=u->G_outflow;
    double **v_outflow=v->G_outflow;
    double **w_outflow=w->G_outflow;
    double U_w=0.0; double U_e=0.0;
    double U_n=0.0; double U_s=0.0;
    double U_f=0.0; double U_b=0.0;
    double V_e=0.0; double V_w=0.0;
    double V_s=0.0; double V_n=0.0;
    double V_f=0.0; double V_b=0.0;
    double W_e=0.0; double W_w=0.0;
    double W_s=0.0; double W_n=0.0;
    double W_f=0.0; double W_b=0.0;

    double Cx  = 1.0/(grid->dEta * grid->dEta);
    double Cy  = 1.0/(grid->dXi * grid->dXi);
    double Cz  = 1.0/(grid->dPhi * grid->dPhi);

    /* Viscous dissipation rate on current processor */
    double G_dissipation_rate = 0.0;

    int i, j, k;
    for (k=Ks; k<Ke; k++) {
       for (j=Js; j<Je; j++) {
           for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /* only include if the current point is fluid*/

                    /* Metric coefficients used for central difference discretization */
                    double U_ae   = Cx * grid->metric_xu[i] * grid->metric_xc[i];
                    double U_aw   = Cx * grid->metric_xu[i];
                    double U_an   = Cy * grid->metric_yc[j] * grid->metric_yv[j+1];
                    double U_as   = Cy * grid->metric_yc[j] * grid->metric_yv[j];;
                    double U_af   = Cz * grid->metric_zc[k] * grid->metric_zw[k+1];
                    double U_ab   = Cz * grid->metric_zc[k] * grid->metric_zw[k];;

                    double V_ae   = Cx * grid->metric_xc[i] * grid->metric_xu[i+1];
                    double V_aw   = Cx * grid->metric_xc[i] * grid->metric_xu[i];;
                    double V_an   = Cy * grid->metric_yv[j] * grid->metric_yc[j];
                    double V_as   = Cy * grid->metric_yv[j];
                    double V_af   = Cz * grid->metric_zc[k] * grid->metric_zw[k+1];
                    double V_ab   = Cz * grid->metric_zc[k] * grid->metric_zw[k];;

                    double W_ae   = Cx * grid->metric_xc[i] * grid->metric_xu[i+1];
                    double W_aw   = Cx * grid->metric_xc[i] * grid->metric_xu[i];;
                    double W_an   = Cy * grid->metric_yc[j] * grid->metric_yv[j+1];
                    double W_as   = Cy * grid->metric_yc[j] * grid->metric_yv[j];;
                    double W_af   = Cz * grid->metric_zw[k] * grid->metric_zc[k];
                    double W_ab   = Cz * grid->metric_zw[k];

                    if (i == 0) {

                        U_aw   *= grid->metric_xc[i];
                        if (params->inflow) {

                            U_w = u_inflow[k][j];
                            V_w = k_left*v_data[k][j][i];
                            W_w = k_left*w_data[k][j][i];
                        } else {

                            U_w = 0.0;
                            V_w = k_left*v_data[k][j][i];
                            W_w = k_left*w_data[k][j][i];
                        }

                    } else {

                        U_aw   *= grid->metric_xc[i-1];
                        U_w = u_data[k][j][i-1];
                        V_w = v_data[k][j][i-1];
                        W_w = w_data[k][j][i-1];

                    } /* else i == 0 */
                    if (j == 0) {

                        V_as   *= grid->metric_yc[j];

                        U_s = k_bottom*u_data[k][j][i];
                        V_s = 0.0;
                        W_s = k_bottom*w_data[k][j][i];
                    } else {

                        V_as   *= grid->metric_yc[j-1];

                        U_s = u_data[k][j-1][i];
                        V_s = v_data[k][j-1][i];
                        W_s = w_data[k][j-1][i];
                    } /* else j == 0 */

                    if (k == 0) {

                        W_ab   *= grid->metric_zc[k];
                        U_b = k_back*u_data[k][j][i];
                        V_b = k_back*v_data[k][j][i];
                        W_b = 0.0;
                    } else {

                        W_ab   *= grid->metric_zc[k-1];
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

                            U_e = 0.0;
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
                        V_n = 0.0;
                        W_n = k_top*w_data[k][j][i];
                    } else {

                        U_n = u_data[k][j+1][i];
                        V_n = v_data[k][j+1][i];
                        W_n = w_data[k][j+1][i];
                    } /* else j == NJ-1 */

                    if (k == NK-1) {

                        U_f = k_front*u_data[k][j][i];
                        V_f = k_front*v_data[k][j][i];
                        W_f = 0.0;
                    } else {

                        U_f = u_data[k+1][j][i];
                        V_f = v_data[k+1][j][i];
                        W_f = w_data[k+1][j][i];
                    } /* else k == NK-1 */

                    double u_ = u_data[k][j][i];
                    double v_ = v_data[k][j][i];
                    double w_ = w_data[k][j][i];

                    double U_ap = -(U_ae + U_aw + U_ab + U_af + U_as + U_an);
                    double V_ap = -(V_ae + V_aw + V_ab + V_af + V_as + V_an);
                    double W_ap = -(W_ae + W_aw + W_ab + W_af + W_as + W_an);

                    double laplace_u = (U_ae*U_e + U_aw*U_w) + (U_an*U_n + U_as*U_s) + (U_af*U_f + U_ab*U_b) + U_ap*u_;
                    double laplace_v = (V_ae*V_e + V_aw*V_w) + (V_an*V_n + V_as*V_s) + (V_af*V_f + V_ab*V_b) + V_ap*v_;
                    double laplace_w = (W_ae*W_e + W_aw*W_w) + (W_an*W_n + W_as*W_s) + (W_af*W_f + W_ab*W_b) + W_ap*w_;

                    double dV_u = (grid->xc[i+1] - grid->xc[i]) * dy_c[j] *dz_c[k];
                    double dV_v = (grid->yc[j+1] - grid->yc[j]) * dx_c[i] *dz_c[k];
                    double dV_w = (grid->zc[k+1] - grid->zc[k]) * dx_c[i] *dy_c[j];

                    /* We are using the trapezoidal rule integration. So, add the contribution of each term to both current cell with dV and the previous cells (west or south or back) */
                    G_dissipation_rate += u_*laplace_u*dV_u + v_*laplace_v*dV_v + w_*laplace_w*dV_w;

                    //printf("Extract.c/ (i,j,k)=(%d,%d,%d) (u,v,w)=(%f,%f,%f) laplace(u,v,w)=(%f,%f,%f)\n", i, j, k, u_, v_, w_, laplace_u, laplace_v, laplace_w);
                    //printf("Extract.c/ (i,j,k)=(%d,%d,%d) U(w,e,s,n,b,f,p)=(%f,%f,%f,%f,%f,%f,%f) U_a(w,e,s,n,b,f,p)=(%f,%f,%f,%f,%f,%f,%f)\n", i, j, k, U_w, U_e, U_s, U_n, U_b, U_f, u_, U_aw, U_ae, U_as, U_an, U_ab, U_af, U_ap);

                } /* if */
            } /* for i*/
        } /* for j*/
    }/* for k*/

    G_dissipation_rate *= 1.0/params->Re;

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, u->L_data, (void ***)&u_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->L_data, (void ***)&v_data); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->L_data, (void ***)&w_data); PETScErrAct(ierr);

    /* Get the global dissipation rate */
    double W_dissipation_rate = 0.0;
    (void)MPI_Allreduce( &G_dissipation_rate, &W_dissipation_rate, 1, MPI_DOUBLE, MPI_SUM, PCW);

    u->G_dissipation_rate = G_dissipation_rate;
    v->G_dissipation_rate = G_dissipation_rate;
    w->G_dissipation_rate = G_dissipation_rate;

    u->W_dissipation_rate = W_dissipation_rate;
    v->W_dissipation_rate = W_dissipation_rate;
    w->W_dissipation_rate = W_dissipation_rate;

}
/***************************************************************************************************/

/* This function computes the rate of change of kinetic to potential energy */
/* _rate = integral (v c dV) */
void Extract_kinetic_potential_conversion_rate(Velocity *v, Concentration **c, MAC_grid *grid, Parameters *params) {

    /* Get regular data array for vel data at cell center */
    double ***v_data_bc=NULL;
    int ierr;
    ierr = DAVecGetArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc);PETScErrAct(ierr);

    Vec *c_vec=NULL;
    /* Pointer to the total concentration field */
    if (params->NConc > 1) {

        c_vec = &c[0]->G_c_total;
    } else {

        c_vec = &c[0]->G_data;
    } /* else */

    double ***conc_total=NULL;
    ierr = DAVecGetArray(grid->DA_3D, *c_vec, (void ***)&conc_total);PETScErrAct(ierr);

    /* Grid dimensions */
    double *dx_c = grid->dx_c;
    double *dy_c = grid->dy_c;
    double *dz_c = grid->dz_c;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    double G_rate = 0.0;
    int i, j, k;
    for (k=Ks; k<Ke; k++) {

        double dz = dz_c[k];
        for (j=Js; j<Je; j++) {

            double dy = dy_c[j];
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /*only include if point is fluid*/

                    double dx = dx_c[i];

                    double v_ = v_data_bc[k][j][i];
                    double c_t = conc_total[k][j][i];

                    /* Add to the conversion rate on current processor */
                    G_rate += -v_* c_t * (dx*dy*dz);

                } /* if */

            } /* for i*/
        } /* for j*/
    }/* for k*/

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, *c_vec, (void ***)&conc_total);PETScErrAct(ierr);

    /* Get the global value */
    double W_rate = 0.0;
    (void)MPI_Allreduce( &G_rate, &W_rate, 1, MPI_DOUBLE, MPI_SUM, PCW);

    v->G_kinetic_potential_conversion_rate = G_rate;
    v->W_kinetic_potential_conversion_rate = W_rate;
}
/***************************************************************************************************/

/* This function computes the total kinetic energy within the interior of the domain */
/* K = 0.5 integral (u.u dV) */
void Extract_total_kinetic_energy(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid) {

    /* Get regular data array for vel data at cell center */
    double ***u_data_bc=NULL;
    double ***v_data_bc=NULL;
    double ***w_data_bc=NULL;
    int ierr;
    ierr = DAVecGetArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc);PETScErrAct(ierr);

    /* Grid dimensions */
    double *dx_c = grid->dx_c;
    double *dy_c = grid->dy_c;
    double *dz_c = grid->dz_c;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    /* Kinetic energy on the current processor */
    double G_kinetic_energy = 0.0;
    int i, j, k;
    for (k=Ks; k<Ke; k++) {

        double dz = dz_c[k];
        for (j=Js; j<Je; j++) {

            double dy = dy_c[j];
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /*only include if point is fluid*/

                    double  dx = dx_c[i];

                    double dV = dx * dy * dz;

                    double u_ = u_data_bc[k][j][i];
                    double v_ = v_data_bc[k][j][i];
                    double w_ = w_data_bc[k][j][i];

                    /* Add to the kinetic energy on current processor */
                    G_kinetic_energy += (u_*u_ + v_*v_ + w_*w_) * dV;

                } /* if */

            } /* for i*/
        } /* for j*/
    }/* for k*/

    G_kinetic_energy *= 0.5;

    /* Restore the arrays */
    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data_bc, (void ***)&u_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, v->G_data_bc, (void ***)&v_data_bc); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, w->G_data_bc, (void ***)&w_data_bc); PETScErrAct(ierr);

    /* Get the global value */
    double W_kinetic_energy = 0.0;
    (void)MPI_Allreduce( &G_kinetic_energy, &W_kinetic_energy, 1, MPI_DOUBLE, MPI_SUM, PCW);

    u->G_kinetic_energy = G_kinetic_energy;
    v->G_kinetic_energy = G_kinetic_energy;
    w->G_kinetic_energy = G_kinetic_energy;

    u->W_kinetic_energy = W_kinetic_energy;
    v->W_kinetic_energy = W_kinetic_energy;
    w->W_kinetic_energy = W_kinetic_energy;

}
/***************************************************************************************************/

/* This function computes the stokes dissipation for any particle concentration field.
   This dissipation is due to the fact that the Stokes flow around settling particles cause some microscopic
dissipation of tenergy
e_s = int( u_s * c dV) domain
*/
void Extract_Stokes_dissipation_rate(Concentration *c, MAC_grid *grid, Parameters *params) {

    /* Get concentration data  */
    int ierr;
    double ***conc=NULL;
    ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

    /* Grid dimensions */
    double *dx_c = grid->dx_c;
    double *dy_c = grid->dy_c;
    double *dz_c = grid->dz_c;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    double G_Stokes_dissipation_rate = 0.0;
    int i, j, k;
    /* Only compute if the concentration field is particle. Otherwise, it is zero */
    if (c->Type == PARTICLE) {

        for (k=Ks; k<Ke; k++) {

            double dz = dz_c[k];
            for (j=Js; j<Je; j++) {

                double dy = dy_c[j];
                for (i=Is; i<Ie; i++) {

                    if (Grid_get_c_status(grid, i, j, k) == FLUID) {  /*only include if point is fluid*/

                        double dx = dx_c[i];

                        /* differntial volume of the current grid */
                        double dV   = dx * dy * dz;

                        G_Stokes_dissipation_rate += conc[k][j][i] * dV;
                    } /* if */

                } /* for i*/
            } /* for j*/
        }/* for k*/
    } /* if */
    else {

        G_Stokes_dissipation_rate = 0.0;
    }

    /* Now, restore the array */
    ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc); PETScErrAct(ierr);

    /* Now, multiply by the constant settling speed and the coefficient contributing to the entire mixture */
    G_Stokes_dissipation_rate *= fabs(c->v_settl0) * params->Conc_alpha[c->conc_index];

    /* Store the rate for the current processor */
    c->G_Stokes_dissipation_rate = G_Stokes_dissipation_rate;

    /* Get the global value */
    double W_Stokes_dissipation_rate=0.0;
    (void)MPI_Allreduce( &G_Stokes_dissipation_rate, &W_Stokes_dissipation_rate, 1, MPI_DOUBLE, MPI_SUM, PCW);

    c->W_Stokes_dissipation_rate = W_Stokes_dissipation_rate;
}
/**********************************************************************************************/

void Extract_sphere_Drag_coefficient(Velocity *u, Velocity *v, Velocity *w, Pressure *p, MAC_grid *grid, Parameters *params, double time) {

	double f_D = 0.0; 
	double f_L = 0.0; 
	int ierr; 

    /* Update the values of the ghost nodes based on a BOX stencil */
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &u->G_data, &u->L_data_box, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &v->G_data, &v->L_data_box, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &w->G_data, &w->L_data_box, 'I');

    double ***u_data=NULL; 
	double ***v_data=NULL; 
	double ***w_data=NULL; 

    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D_BOX, v->L_data_box, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D_BOX, w->L_data_box, (void ***)&w_data);PETScErrAct(ierr);

	/* sphere properties */
	double x_cent = params->sphere_cent_x; 
	double y_cent = params->sphere_cent_y; 
	double z_cent = params->sphere_cent_z; 

	int n_segments = 200; 
	double dphi = 2.0*PI/(double)n_segments; 
	double dtheta = dphi;
	double dR = 1.0*0.015;//grid->dx_c[grid->NI/2]; 
	double R = params->sphere_radius; /* 0.5 */

	/* force due to shear */
	if (params->current_rank_has_sphere == YES) {

		double theta, phi; 
		for (theta=0.0; theta <= PI; theta+=dtheta) {
			for (phi=0.0; phi<= 2.0*PI; phi+=dphi) {

				double dA = R*R*sin(theta)*dtheta*dphi;
				PointType f_int; 

				f_int.x = x_cent + (R + dR) * sin(theta) * cos(phi); 
				f_int.y = y_cent + (R + dR) * sin(theta) * sin(phi); 
				f_int.z = z_cent + (R + dR) * cos(theta);

                /* Use a trilinear interpolation to find the velocity at the interpolation node */
				double u_int = MyMath_do_trilinear_inter(f_int.x, f_int.y, f_int.z, u_data, grid, 'u');
				double v_int = MyMath_do_trilinear_inter(f_int.x, f_int.y, f_int.z, v_data, grid, 'v');
				double w_int = MyMath_do_trilinear_inter(f_int.x, f_int.y, f_int.z, w_data, grid, 'w');

				printf("(x,y,z)=(%f,%f,%f) (u,v,w)_int=(%f,%f,%f)\n", f_int.x, f_int.y, f_int.z, u_int, v_int, w_int); 
				/* pressure */
				//double p_int = MyMath_do_trilinear_inter(f_int.x, f_int.y, f_int.z, p_data, grid, 'p');

				double nx = sin(theta) * cos(phi); 
				double ny = sin(theta) * sin(phi); 
				double nz = cos(theta); 

                /* V = V_t.e_t + V_n.e_n */
                /* e_n = nx.i + ny.j + nz.k */
                /* magnitude of the normal velocity vector */
				double V_n = u_int*nx + v_int*ny + w_int*nz;

				/* Tangential component: V_t = V - V_n */
				double u_tan = u_int - V_n*nx;
				double v_tan = v_int - V_n*ny;
				double w_tan = w_int - V_n*nz;

				/* tangential velocity component */
				double V_tan = sqrt(u_tan*u_tan + v_tan*v_tan + w_tan*w_tan);
                double tau = V_tan / (dR * params->Re);

				/* Pressure force */
				//f_P += dA * P_; 

				double F_x = 0.0; 
				double F_y = 0.0; 
				if ( (V_tan > 1.e-12) ){
					F_x = (tau * dA) * (u_tan / V_tan); 
					F_y = (tau * dA) * (v_tan / V_tan); 
				}

				/* Drag force */
				f_D += F_x; 
				f_L += F_y; 

				//printf("(x,y,z)=(%f,%f,%f) (u,v,w)_int=(%f,%f,%f) (fX,fY)=(%f,%f) \n", f_int.x, f_int.y, f_int.z, u_int, v_int, w_int, F_x, F_y); 
			}
		}
	}

    /* Get regular data array for vel data */
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&u_data);PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, v->L_data_box, (void ***)&v_data);PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, w->L_data_box, (void ***)&w_data);PETScErrAct(ierr);


    /* Update the values of the ghost nodes based on a BOX stencil */
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &p->G_data, &u->L_data_box, 'I');

    double ***p_data=NULL; 

    /* Get regular data array for vel data */
    ierr = DAVecGetArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&p_data);PETScErrAct(ierr);

	/* force due to pressure */
	if (params->current_rank_has_sphere == YES) {

		double theta, phi; 
		dR = 1.0*0.015; 
		for (theta=0.0; theta <= PI; theta+=dtheta) {
			for (phi=0.0; phi<= 2.0*PI; phi+=dphi) {

				double dA = R*R*sin(theta)*dtheta*dphi;
				PointType f_int; 

				f_int.x = x_cent + (R + dR) * sin(theta) * cos(phi); 
				f_int.y = y_cent + (R + dR) * sin(theta) * sin(phi); 
				f_int.z = z_cent + (R + dR) * cos(theta);

                /* Use a trilinear interpolation to find the pressure at the interpolation node */
				double p_int = MyMath_do_trilinear_inter(f_int.x, f_int.y, f_int.z, p_data, grid, 'c');

				double nx = sin(theta) * cos(phi); 
				double ny = sin(theta) * sin(phi); 
				double nz = cos(theta); 

				/* Pressure force */
				double f_P = dA * p_int; 

				double F_x = -f_P * nx; 
				double F_y = -f_P * ny;

				/* Drag force */
				f_D += F_x; 
				f_L += F_y; 
			}
		}

		f_D /= 0.5*PI*R*R; 
		f_L /= 0.5*PI*R*R; 

		FILE *f_out=NULL; 
		f_out = fopen("Cd_Cl.dat", "a"); 

		if (f_out == NULL) {

			printf("Error! Could not open Cd_Cl.dat file\n"); 
		}
		fprintf(f_out, "%f %2.14lf %2.14lf\n", time, f_D, f_L); 
		fclose(f_out); 
	}

    /* Get regular data array for vel data */
    ierr = DAVecRestoreArray(grid->DA_3D_BOX, u->L_data_box, (void ***)&p_data);PETScErrAct(ierr);
}


