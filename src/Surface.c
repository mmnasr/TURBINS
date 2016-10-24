#include "definitions.h"
#include "DataTypes.h"
#include "Communication.h"
#include "MyMath.h"
#include "Memory.h"
#include "Grid.h"
#include "Levelset.h"
#include "Surface.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* This function allocates memory for the surface of the solid object in the domain */
/* Note that surface is represented in an implicit way, i.e. using a levelset function with sdf(x,y,z) with which sdf(x,y,z) = 0 represents the surface of the solid boundaries */
SurfaceType *Surface_create(MAC_grid *grid) {

    int ierr;
    SurfaceType *new_surf;

    new_surf = (SurfaceType *)calloc(1, sizeof(SurfaceType));
    Memory_check_allocation(new_surf);

    /* total number of grid points including the half cell added to the right. */
    /* note that surface is represented implicitly, so we need a Signed Distance Function given at all (x,y,z) nodes to give the distance from the surface */
    new_surf->NX = grid->NX;
    new_surf->NY = grid->NY;
    new_surf->NZ = grid->NZ;

    /* Create levelset object associated with the representation of the surface. Note that the goal is to find the sign distance function for the whole */
    /* entire domain. That way, the value of sdf (Sign Distance Function) at each grid point will give us the distance to the interface of the solid object */
    /* after reinitializtion, we can destroy this part of the surface structure, however, if the surface is evolving in time, one can keep this part of the surface structure */
    new_surf->sdf_object = Levelset_create(grid);

    /* sdf at u_grid position */
    /* create the local and global versions of the arrays */
    ierr = DACreateGlobalVector(grid->DA_3D, &new_surf->G_u_sdf); PETScErrAct(ierr);
    ierr = DACreateLocalVector(grid->DA_3D, &new_surf->L_u_sdf); PETScErrAct(ierr);

    /* sdf at v_grid position */
    /* create the local and global versions of the arrays */
    ierr = VecDuplicate(new_surf->G_u_sdf, &new_surf->G_v_sdf); PETScErrAct(ierr);
    ierr = VecDuplicate(new_surf->L_u_sdf, &new_surf->L_v_sdf); PETScErrAct(ierr);

    /* sdf at w_grid position */
    /* create the local and global versions of the arrays */
    ierr = VecDuplicate(new_surf->G_u_sdf, &new_surf->G_w_sdf); PETScErrAct(ierr);
    ierr = VecDuplicate(new_surf->L_u_sdf, &new_surf->L_w_sdf); PETScErrAct(ierr);

    /* sdf at c_grid position */
    /* create the local and global versions of the arrays */
    ierr = VecDuplicate(new_surf->G_u_sdf, &new_surf->G_c_sdf); PETScErrAct(ierr);
    ierr = VecDuplicate(new_surf->L_u_sdf, &new_surf->L_c_sdf); PETScErrAct(ierr);

    /* get the pointers to the vec arrays. That will save a huge amount of time */
    ierr = DAVecGetArray(grid->DA_3D, new_surf->L_c_sdf, (void ***)&new_surf->c_sdf);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, new_surf->L_u_sdf, (void ***)&new_surf->u_sdf);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, new_surf->L_v_sdf, (void ***)&new_surf->v_sdf);PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, new_surf->L_w_sdf, (void ***)&new_surf->w_sdf);PETScErrAct(ierr);

    return new_surf;
}
/***************************************************************************************************/

/* This function releases the allocated memory for the surface */
void Surface_destroy(SurfaceType *surf, MAC_grid *grid)  {

    int ierr;

    /* First, restore the arrays from the local pointers */
    ierr = DAVecRestoreArray(grid->DA_3D, surf->L_c_sdf, (void ***)&surf->c_sdf); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, surf->L_u_sdf, (void ***)&surf->u_sdf); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, surf->L_v_sdf, (void ***)&surf->v_sdf); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, surf->L_w_sdf, (void ***)&surf->w_sdf); PETScErrAct(ierr);

    ierr = VecDestroy(surf->G_u_sdf); PETScErrAct(ierr);
    ierr = VecDestroy(surf->L_u_sdf); PETScErrAct(ierr);

    ierr = VecDestroy(surf->G_v_sdf); PETScErrAct(ierr);
    ierr = VecDestroy(surf->L_v_sdf); PETScErrAct(ierr);

    ierr = VecDestroy(surf->G_w_sdf); PETScErrAct(ierr);
    ierr = VecDestroy(surf->L_w_sdf); PETScErrAct(ierr);

    ierr = VecDestroy(surf->G_c_sdf); PETScErrAct(ierr);
    ierr = VecDestroy(surf->L_c_sdf); PETScErrAct(ierr);

    free(surf);
}
/***************************************************************************************************/

/* This function initializes the sdf based on the location of the interface and/or solid objects in the domain */
/* Note that this is NOT the correct signed distance function */
/* All you need to do is a "good" guess for that in a way that the location of inteface is correct */
/* Also, the fluid region has a positive sign and in the solid region, a negative sign */
void Surface_initialize(SurfaceType *surf, MAC_grid *grid, Parameters *params) {

    double ***sdf_guess;
    double x, y, z;
    double *xc, *yc, *zc;
    double Lx, Ly, Lz;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;
    int i, j, k;
    int Only_bottom_surface = YES;

    /* Length of the domain */
    Lx = params->Lx;
    Ly = params->Ly;
    Lz = params->Lz;

    /* cell centered coordinates of the grid points */
    xc = grid->xc;
    yc = grid->yc;
    zc = grid->zc;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* get the global version of first guess for the sdf */
    ierr = DAVecGetArray(grid->DA_3D, surf->sdf_object->G_data_0, (void ***)&sdf_guess); PETScErrAct(ierr);

    /* Go through all the node in the entire domain on this processor and set the first guess for the signed distance function which represents */
    /* the correct location of the inteface */
    /* We are doing this only for the cell center values */
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                x = xc[i];
                y = yc[j];
                z = zc[k];

                /* If we have only a solid surface on the bottom h=f(x,z) then use that grid->inteface_position */
                if (Only_bottom_surface) {

                    /* height of the current node minus the height of bottom topography at h=G(x,z) */
                    sdf_guess[k][j][i] = yc[j] - grid->interface_position[k][i];
                } else {

                    /* Else, define the inteface of the solid object(s) here as a function of (x,y,z) */
                    /* Insert box */
                    double d_x, d_y;
                    double XS = 0.5*Lx;
                    double YS = 0.5*Ly;
                    double ax = 1.0;
                    double ay = 1.0;
                    d_x = min( fabs(xc[i] - (XS -0.5*ax)), fabs(xc[i] + (XS -0.5*ax)));
                    d_y = min( fabs(yc[j] - (YS -0.5*ay)), fabs(yc[j] + (YS -0.5*ay)));
                    if ( (xc[i] >= (XS-0.5*ax) ) && (xc[i] <= (XS+0.5*ay) ) &&
                         (yc[j] >= (YS-0.5*ay) ) && (yc[j] <= (YS+0.5*ay) ) ) {

                        sdf_guess[k][j][i] = -min(d_x,d_y);
                    } else {

                        sdf_guess[k][j][i] = min(d_x,d_y)+2.0*grid->dx;
                    }


                    /* insert cyliner */
                    sdf_guess[k][j][i] = sqrt( (y - 0.5*Ly)*(y - 0.5*Ly) + (x - 0.5*Lx)*(x - 0.5*Lx) ) - 0.25;

                    //d_y = 0.5 - 0.2*x - 0.2*z;
                    //sdf_guess[k][j][i] = y - d_y;

                    /* For instance, for a sphere located at (Lx/2, Ly/2, Lz/2) and Radius of 0.5, we have */
                    sdf_guess[k][j][i] = sqrt( (x - 0.5*Lx)*(x - 0.5*Lx) + (y - 0.5*Ly)*(y - 0.5*Ly) + (z - 0.50*Lz)*(z - 0.50*Lz) ) - 0.5;

                    /* insert sphere with diameter 1.0 */
                    double x_cent = params->sphere_cent_x;
                    double y_cent = params->sphere_cent_y;
                    double z_cent = params->sphere_cent_z;
                    sdf_guess[k][j][i] = sqrt( (x - x_cent)*(x - x_cent) + (y - y_cent)*(y - y_cent) + (z - z_cent)*(z - z_cent) ) - 0.5;


                } /* else */
            }/* for i */
        } /* for j */
    } /* for k */

    /* Restore the array */
    ierr = DAVecRestoreArray(grid->DA_3D, surf->sdf_object->G_data_0, (void ***)&sdf_guess); PETScErrAct(ierr);

}
/***************************************************************************************************/

/* This function reinitializes the sdf to compute the correct sdf using a levelset approach */
/* Note that levelset function is computed at the cell centers */
void Surface_find_exact_sdf(SurfaceType *surf, MAC_grid *grid) {

    int max_iterations = 500;  /* maximum number of iterations for reinitialization scheme */
    int ierr;

    /* Reinitialize the inital signed distance function to compute the correct signed distance function */
    Levelset_reinitialize(surf->sdf_object, grid, max_iterations);

    /* Now, copy the levelset function into the data array. Note that levelset function was computed in the cell centers */
    ierr = VecCopy(surf->sdf_object->G_data, surf->G_c_sdf); PETScErrAct(ierr);

    /* Now, you can delete the "sdf_object" for a solid boundary which is not moving */
    Levelset_destroy(surf->sdf_object);
}
/***************************************************************************************************/

/* This function computes the value of the signed distance function for the other variables assuming that we have the signed distance function at grid cell center */
/* For now, we only use a 1-D interpolation */
/* q: u, v, and w */
void Surface_find_q_sdf(SurfaceType *surf, MAC_grid *grid) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    int ierr;
    double ***c_sdf;
    double ***u_sdf, ***v_sdf, ***w_sdf;
    double *xu, *yv, *zw;
    double *xc, *yc, *zc;
    Scalar q[2], q_int;

    /* grid coordinates */
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

    /* Get the local version of c_sdf since we need the ghost nodes for the interpolations */

    /* update the ghost nodes first on the current processor */
    Communication_update_ghost_nodes(&grid->DA_3D, &surf->G_c_sdf, &surf->L_c_sdf, 'I');
    ierr = DAVecGetArray(grid->DA_3D, surf->L_c_sdf, (void ***)&c_sdf); PETScErrAct(ierr);

    /* get the arrays for the u,v and w grid sdf functions */
    ierr = DAVecGetArray(grid->DA_3D, surf->G_u_sdf, (void ***)&u_sdf); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, surf->G_v_sdf, (void ***)&v_sdf); PETScErrAct(ierr);
    ierr = DAVecGetArray(grid->DA_3D, surf->G_w_sdf, (void ***)&w_sdf); PETScErrAct(ierr);

    /* Go throught the entire domain and use a two point interpolation to find the value of the sdf at the cell center onto the u,v,w grid nodes */

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                /* sdf at u_grid node */
                /* values of sdf to the west and east of the current u-node */
                q[0].Value = c_sdf[k][j][i];
                q[0].X1    = xc[i];
                if (i > 0) {

                    q[1].Value = c_sdf[k][j][i-1];
                    q[1].X1    = xc[i-1];
                } else { /* The boundary node. The value is not important */

                    q[1].Value = c_sdf[k][j][i];
                    q[1].X1    = xc[i];
                } /* else */

                q_int.X1   = xu[i];
                u_sdf[k][j][i] = MyMath_interpolate_quantity(q, q_int, 2);
                /****************************/

                /* sdf at v_grid node */
                /* values of sdf to the south and north of the current v-node */
                q[0].Value = c_sdf[k][j][i];
                q[0].X1    = yc[j];
                if (j > 0) {

                    q[1].Value = c_sdf[k][j-1][i];
                    q[1].X1    = yc[j-1];
                } else { /* The boundary node. The value is not important */

                    q[1].Value = c_sdf[k][j][i];
                    q[1].X1    = yc[j];
                } /* else */

                q_int.X1   = yv[j];
                v_sdf[k][j][i] = MyMath_interpolate_quantity(q, q_int, 2);
                /****************************/

                /* sdf at w_grid node */
                /* values of sdf to the bacl and front of the current w-node */
                q[0].Value = c_sdf[k][j][i];
                q[0].X1    = zc[k];
                if (k > 0) {

                    q[1].Value = c_sdf[k-1][j][i];
                    q[1].X1    = zc[k-1];
                } else { /* The boundary node. The value is not important */

                    q[1].Value = c_sdf[k][j][i];
                    q[1].X1    = zc[k];
                } /* else */

                q_int.X1   = zw[k];
                w_sdf[k][j][i] = MyMath_interpolate_quantity(q, q_int, 2);
                /****************************/


            } /* for i */
        } /* for j */
    } /* for k */

    ierr = DAVecRestoreArray(grid->DA_3D, surf->L_c_sdf, (void ***)&c_sdf); PETScErrAct(ierr);

    ierr = DAVecRestoreArray(grid->DA_3D, surf->G_u_sdf, (void ***)&u_sdf); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, surf->G_v_sdf, (void ***)&v_sdf); PETScErrAct(ierr);
    ierr = DAVecRestoreArray(grid->DA_3D, surf->G_w_sdf, (void ***)&w_sdf); PETScErrAct(ierr);

    /* update the arrays with the ghost nodes */
    Communication_update_ghost_nodes(&grid->DA_3D, &surf->G_u_sdf, &surf->L_u_sdf, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D, &surf->G_v_sdf, &surf->L_v_sdf, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D, &surf->G_w_sdf, &surf->L_w_sdf, 'I');

}
/***************************************************************************************************/

/* This function finds the distance (in x-direction and not the normal distance) between from the current node (i,j,k) and interface */
double Surface_find_Dx_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity) {

    double tol = 1.0e-2;
    double *xq=NULL;

    int NI = grid->NI;
    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        NI++;
        break;
    case 'v':
        xq = grid->xv;
        break;
    case 'w':
        xq = grid->xw;
        break;
    case 'c':
        xq = grid->xc;
        break;
    } /* switch */

    int index_e = max(i-1,0);
    int index_w = min(i+1,NI-1);
    double Delta_x = xq[index_e] - xq[index_w];

    /* D_Dx at (i,j,k) */
    /* One can use a higher order interpolation to compute Distance */
    double D_sdf_dx = fabs ( (sdf[k][j][index_e] - sdf[k][j][index_w]) / Delta_x ) ;

    if ( (fabs(D_sdf_dx) < tol) || (isnan(D_sdf_dx)) ){

        D_sdf_dx = tol;
    } /* if */

    /* The horizontal distance between current node and inteface */
    /* We have used a first order approximation of Taylor series to reach result */
    double Dx = sdf[k][j][i] / D_sdf_dx;

    return Dx;
}
/***************************************************************************************************/

/* This function finds the distance (in y-direction and not the normal distance) between from the current node (i,j,k) and interface */
double Surface_find_Dy_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity) {

    double tol = 1.0e-2;

    int NJ = grid->NJ;
    double *yq=NULL;

    switch (which_quantity) {

    case 'u':
        yq = grid->yu;
        break;
    case 'v':
        yq = grid->yv;
        NJ++;
        break;
    case 'w':
        yq = grid->yw;
        break;
    case 'c':
        yq = grid->yc;
        break;
    } /* switch */

    int index_s = max(j-1,0);
    int index_n = min(j+1,NJ-1);
    double Delta_y = yq[index_n] - yq[index_s];

    /* D_Dy at (i,j,k) */
    /* One can use a higher order interpolation to compute Distance */
    double D_sdf_dy = fabs ( (sdf[k][index_n][i] - sdf[k][index_s][i]) / Delta_y ) ;

    if ( (fabs(D_sdf_dy) < tol) || (isnan(D_sdf_dy) ) ){

        D_sdf_dy = tol;
    } /* if */

    /* The horizontal distance between current node and inteface */
    /* We have used a first order approximation of Taylor series to reach result */
    double Dy = sdf[k][j][i] / D_sdf_dy;

    return Dy;
}
/***************************************************************************************************/

/* This function finds the distance (in x-direction and not the normal distance) between from the current node (i,j,k) and interface */
double Surface_find_Dz_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity) {

    double tol = 1.0e-2;
    double *zq=NULL;

    int NK=grid->NK;
    switch (which_quantity) {

    case 'u':
        zq = grid->zu;
        break;
    case 'v':
        zq = grid->zv;
        break;
    case 'w':
        zq = grid->zw;
        NK++;
        break;
    case 'c':
        zq = grid->zc;
        break;
    } /* switch */

    int index_b = max(k-1,0);
    int index_f = min(k+1,NK-1);
    double Delta_z = zq[index_f] - zq[index_b];

    /* D_Dz at (i,j,k) */
    /* One can use a higher order interpolation to compute Distance */
    double D_sdf_dz = fabs ( (sdf[index_f][j][i] - sdf[index_b][j][i]) / Delta_z ) ;

    if ( (fabs(D_sdf_dz) < tol) || (isnan(D_sdf_dz) ) ){

        D_sdf_dz = tol;
    } /* if */

    /* The horizontal distance between current node and inteface */
    /* We have used a first order approximation of Taylor series to reach result */
    double Dz = sdf[k][j][i] / D_sdf_dz;

    return Dz;
}
/***************************************************************************************************/

/* This function finds the distance between the current node and inteface along with the given "n" vector direction */
/* "n" could be any arbitrary direction */
double Surface_find_Dn_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, VectorType *n, char which_quantity) {

    double tol = 1.0e-2;

    int NI=grid->NI;
    int NJ=grid->NJ;
    int NK=grid->NK;

    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;

    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        NI++;
        break;
    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        NJ++;
        break;
    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        NK++;
        break;
    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        break;
    } /* switch */

    int index_w = max(i-1,0);
    int index_e = min(i+1,NI-1);
    int index_s = max(j-1,0);
    int index_n = min(j+1,NJ-1);
    int index_b = max(k-1,0);
    int index_f = min(k+1,NK-1);

    double Delta_x = xq[index_e] - xq[index_w];
    double Delta_y = yq[index_n] - yq[index_s];
    double Delta_z = zq[index_f] - zq[index_b];

    /* One can use a higher order interpolation to compute Distance */
    double D_sdf_dx = (sdf[k][j][index_e] - sdf[k][j][index_w]) / Delta_x ;
    double D_sdf_dy = (sdf[k][index_n][i] - sdf[k][index_s][i]) / Delta_y ;
    double D_sdf_dz = (sdf[index_f][j][i] - sdf[index_b][j][i]) / Delta_z ;

    /* d_sdf_dn = grad(sdf) . n */
    double D_sdf_dn = n->vx*D_sdf_dx + n->vy*D_sdf_dy + n->vz*D_sdf_dz;
    if ( (fabs(D_sdf_dn) < tol) || (isnan(D_sdf_dn) ) ){

        D_sdf_dn = tol;
    } /* if */

    /* The distance between the current node and interface along with the "n" vector direction */
    /* "n" could be any arbitrary vector */
    double Dn = sdf[k][j][i] / D_sdf_dn;

    return Dn;
}
/***************************************************************************************************/

/* This function finds the normal vector on the surface of the boundary. 
We approximate the normal vector to be equal to the vector computed at the immersed node */
VectorType Surface_compute_normal(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity) {

    double Delta_x, Delta_y, Delta_z;
    double D_sdf_dx, D_sdf_dy, D_sdf_dz;
    double a_norm;
    int index_e, index_w;
    int index_n, index_s;
    int index_f, index_b;
    VectorType n;

    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;

    int NI = grid->NI;
    int NJ = grid->NJ;
    int NK = grid->NK;

    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        NI++;
        break;

    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        NJ++;
        break;

    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        NK++;
        break;

    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        break;
    } /* switch */

    index_w = max(i-1,0);
    index_e = min(i+1,NI-1);
    index_s = max(j-1,0);
    index_n = min(j+1,NJ-1);
    index_b = max(k-1,0);
    index_f = min(k+1,NK-1);

    Delta_x = xq[index_e] - xq[index_w];
    Delta_y = yq[index_n] - yq[index_s];
    Delta_z = zq[index_f] - zq[index_b];

    /* One can use a higher order interpolation to compute Distance */
    D_sdf_dx = (sdf[k][j][index_e] - sdf[k][j][index_w]) / Delta_x ;
    D_sdf_dy = (sdf[k][index_n][i] - sdf[k][index_s][i]) / Delta_y ;
    D_sdf_dz = (sdf[index_f][j][i] - sdf[index_b][j][i]) / Delta_z ;

    /* It should be equal to one, but for the nodes close to the boundary it might not be, so let us manually set it equal to the unit vector */
    a_norm = 1.0/sqrt(D_sdf_dx*D_sdf_dx + D_sdf_dy*D_sdf_dy + D_sdf_dz*D_sdf_dz);

    n.vx = D_sdf_dx * a_norm;
    n.vy = D_sdf_dy * a_norm;
    n.vz = D_sdf_dz * a_norm;

    return n;
}
