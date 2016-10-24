#include "definitions.h"
#include "DataTypes.h"
#include "Grid.h"
#include "Surface.h"
#include "Immersed.h"
#include "MyMath.h"
#include "Memory.h"
#include "Display.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This function allocates memory for the immersed structure based on the given number of the grids */
/* For now, we assume that only the bottom boundary is changing! No I have upgraded it to a cooler sexier version! :-) */
Immersed *Immersed_create(MAC_grid *grid, char which_quantity) {

    int N;
    Immersed *new_immersed;
    int NZ, NY, NX;
    int ierr;

    NX = grid->NX;
    NY = grid->NY;
    NZ = grid->NZ;

    new_immersed = (Immersed *)malloc(sizeof(Immersed));
    Memory_check_allocation(new_immersed);

    /* Total number of the immersed nodes. It should be more or equal to NX */
    /* TBC */
    N = Immersed_count_immersed_nodes(grid, which_quantity);
    PetscPrintf(PCW, "Found %d immeresed points in the domain on the current processor for the %c-quantity\n", N, which_quantity);

    /* Maximum number of immersed nodes in the whole domain */
    new_immersed->N = N;

    /* Allocate the memory for each individual immersed node */
    new_immersed->ib_nodes = (ImmersedNode *)calloc(N, sizeof(ImmersedNode));
    Memory_check_allocation(new_immersed->ib_nodes);

    /* Allocate memory for the global index of the current ib_node in the entire domain. If not an ib_node, the value is -1 */
    /* This will help finding the (i,j,k) location of the immersed node */
    ierr = DACreateGlobalVector(grid->DA_3D, &new_immersed->G_global_index); PETScErrAct(ierr);
    ierr = DACreateLocalVector(grid->DA_3D, &new_immersed->L_global_index); PETScErrAct(ierr);

    /* Pass the 3d int-pointer to the L_global_index. This will help save a huge amount of time */
    ierr = DAVecGetArray(grid->DA_3D, new_immersed->L_global_index, (void ***)&new_immersed->global_index); PETScErrAct(ierr);

    /* u, v, c (or p) */
    new_immersed->quantity = which_quantity;

    return new_immersed;
}
/**************************************************************************************************************/

/* This function releases the memory allocated for the immersed structure */
void Immersed_destroy(Immersed *q_immersed, MAC_grid *grid) {

    int ierr;

    free(q_immersed->ib_nodes);

    /* release the pointer to the Vector */
    ierr = DAVecRestoreArray(grid->DA_3D, q_immersed->L_global_index, (void ***)&q_immersed->global_index); PETScErrAct(ierr);

    ierr = VecDestroy(q_immersed->G_global_index); PETScErrAct(ierr);
    ierr = VecDestroy(q_immersed->L_global_index); PETScErrAct(ierr);

    free(q_immersed);
}
/**************************************************************************************************************/

/* This function counts the total number of immeresed nodes for each quantity */
int Immersed_count_immersed_nodes(MAC_grid *grid, char which_quantity) {

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status; break;
    case 'v':
        Grid_get_q_status = &Grid_get_v_status; break;
    case 'w':
        Grid_get_q_status = &Grid_get_w_status; break;
    case 'c':
        Grid_get_q_status = &Grid_get_c_status; break;
    default:
        PetscPrintf(PCW, "Immersed.c/ Cannot count immersed nodes. Invalid quantity\n");
    } /* switch */

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    /* Go through local part on each processor. Count number of immersed nodes on the current processor */
    int i, j, k;
    int N_immersed = 0;
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (Grid_get_q_status(grid, i, j, k) == IMMERSED ) {

                    N_immersed++;
                } /* if */
            } /* for i */
        } /* for j */
    } /* for k */

    return N_immersed;
}
/**************************************************************************************************************/

/* This function computes the interpolation coefficients of the IB node based on the poistion of the
neighboring fluid nodes and the position of the exact boundary point (control point) */
void Immersed_compute_interpolation_coef(ImmersedNode *ib_node, char which_quantity) {

    double ds, df;
    double coef;
    int g;
    int n_fluid;
    double weight_coef[IBM_MAX];
    double alpha, beta, gamma;

    /* Number of fluid nodes used for interpolation */
    n_fluid = ib_node->n_fluid;

    /* distance between the image point and the boundary (control point) */
    df = MyMath_get_point_point_distance(&ib_node->image_point, &ib_node->control_point);

    /* distance between the immersed point and the boundary (control point) */
    ds = MyMath_get_point_point_distance(&ib_node->im_point, &ib_node->control_point);

    /* should be multiplied to the coefficients of the interpolation */
    /* corresponding to a linear interpolation between the image node and the immersed node satisfying the boundary condition exaclty on the boundary (control point) */
    double tol = 1e-10;
    /* to avoid division by zero */
    if ( (fabs(ds) <= tol) && (fabs(df) <= tol) ) {
        coef = 1.0;
    } else {
        coef = ds/df;
    }


    /* Coefficient should be multiplied by the value at at the control point */
    /* For zero Dirichlet, no need to do anything */
    if ( (coef - 1.0) > 1e-5) {

        printf("Immersed.c/ Warning coef=%2.14f > 1. This might cause the system go unstable. ds:%f df:%f\n", coef, ds, df);
        //Display_immersed_node(ib_node);
        //getchar();
    } /* if coef */

    /* Perform a trilinear interpolation using the 8 nodes located at the corners of the box */
    /* Always follow this standard for the 8 neigboring fluid nodes. This standard is used to find the interpolation coefficients using a trilinear interpolation */

    /*
                 F4------------F6
                / |           /|
               /  |          / |
              /   |         /  |
            F5------------F7   |
            |   F0---------|---F2
            |  /	 I     |  /
            | /            | /
            |/             |/
            F1-------------F3

     z
     |
     |___y
     /
    /
   x
*/

    if (n_fluid == 8) {

        /* Coordinates of the top-right-front corner */
        double x7 = ib_node->fluid_point[7].x;
        double y7 = ib_node->fluid_point[7].y;
        double z7 = ib_node->fluid_point[7].z;

        double x0 = ib_node->fluid_point[0].x;
        double y0 = ib_node->fluid_point[0].y;
        double z0 = ib_node->fluid_point[0].z;


        /* Trilinear interpolation coefficients */
        alpha = (x7 - ib_node->image_point.x)/(x7 - x0);
        beta  = (y7 - ib_node->image_point.y)/(y7 - y0);
        gamma = (z7 - ib_node->image_point.z)/(z7 - z0);

        /* Get the interpolation coefficients to find the value of the image node based on the neigboring fluid nodes */
        weight_coef[0] = alpha * beta * gamma;
        weight_coef[1] = (1.0-alpha) * beta * gamma;
        weight_coef[2] = alpha * (1.0-beta) * gamma;
        weight_coef[3] = (1.0-alpha) * (1.0-beta) * gamma;
        weight_coef[4] = alpha * beta * (1.0-gamma);
        weight_coef[5] = (1.0-alpha) * beta * (1.0-gamma);
        weight_coef[6] = alpha * (1.0-beta) * (1.0-gamma);
        weight_coef[7] = (1.0-alpha)*(1.0-beta)*(1.0-gamma);

    } else if (n_fluid == 1){ /* Just copy the value of one fluid node into the image node */

        weight_coef[0] = 1.0;

    } else { /* Something is wrong */

        PetscPrintf(PCW, "Immersed.c/ Can not perform the interpolation for %d number of fluid nodes (should be either 8 or 1)\n", n_fluid);
    } /* else */

    /* Now, depending on the boundary condition, interpolate the Immersed node value based on the neighboring fluid nodes and the normal unit vector on the boundary (used for Neumann and Robin boundary conditon) */
    switch (ib_node->boundary_condition) {

    /* phi(xcont, ycont, zcont) = b */
    case DIRICHLET: {

        /* If non-zero Dirichlet B.C. is imposed (b!=0), this coefficient is used to account for the nonzero term in the RHS of the u, v, and w linear system */
        ib_node->boundary_coef = 1.0 + coef;
        for (g=0; g<n_fluid; g++) {

            ib_node->fluid_coef[g] = weight_coef[g] * coef * (+1.0);  /* q(xc,yc,zc) = b */
        } /* for g */
    }/* DIRICHLET */
        break;

        /* dphi_dn(xcont, ycont, zcont) = b */
    case NEUMANN: {


        /* If only one fluid node is found, set the coefficient and stop here */
        if (n_fluid == 1) {

            ib_node->boundary_coef = -(ds+df);
            ib_node->fluid_coef[0] = weight_coef[0]*(-1.0);
            break;
        } /* if 1 */


        /* If non-zero Neumann B.C. is imposed (b!=0), this coefficient is used to account for the nonzero term in the RHS of the c linear system */
        int which_method = 3;
        if (which_method == 1) {
            ib_node->boundary_coef = -(ds+df);
            for (g=0; g<n_fluid; g++) {

                ib_node->fluid_coef[g] = weight_coef[g] * (-1.0);  /* dq(xc,yc,zc)/dn = b */
            } /* for g */
        } /* if which_method = 1 */

        /* This method is more accurate since it includes the cases where the mirrored nodes is very close to the boundaries as well */
        if (which_method == 3) {

            /* Same as method 1, except for the case that mirrored point is surrounded by the immersed node itself. */
            /* In such case, we use the boundary condition in the interpolation itself */
            /* q = a*xyz + b*xy + c*xz + d*yz + e*x + f*y + g*z + h */

            MatrixType *M;
            M = (MatrixType *)calloc(1,sizeof(MatrixType));
            M->size = 8;
            short int self_included = NO;
            int r = 0;
            for (g=0; g<n_fluid; g++) {

                if ( (ib_node->fluid_index[g].x_index == ib_node->im_index.x_index)
                     && (ib_node->fluid_index[g].y_index == ib_node->im_index.y_index)
                     && (ib_node->fluid_index[g].z_index == ib_node->im_index.z_index) ) {

                    self_included = YES;
                } else {

                    double x = ib_node->fluid_point[g].x - ib_node->image_point.x;
                    double y = ib_node->fluid_point[g].y - ib_node->image_point.y;
                    double z = ib_node->fluid_point[g].z - ib_node->image_point.z;
                    M->A[r][0] = x*y*z;
                    M->A[r][1] = x*y;
                    M->A[r][2] = x*z;
                    M->A[r][3] = y*z;
                    M->A[r][4] = x;
                    M->A[r][5] = y;
                    M->A[r][6] = z;
                    M->A[r][7] = 1.0;
                    r++;
                }
            }
            /* If ib_node is also included, then just use the boundary condition as the 4th equation to solve for the bilinear coefficients */
            if (self_included) {

                double xc = ib_node->control_point.x - ib_node->image_point.x;
                double yc = ib_node->control_point.y - ib_node->image_point.y;
                double zc = ib_node->control_point.z - ib_node->image_point.z;
                double nx = ib_node->n.vx;
                double ny = ib_node->n.vy;
                double nz = ib_node->n.vz;

                M->A[r][0] = nx*yc*zc + ny*xc*zc + nz*xc*yc;
                M->A[r][1] = yc*nx + xc*ny;
                M->A[r][2] = zc*nx + xc*nz;
                M->A[r][3] = zc*ny + yc*nz;
                M->A[r][4] = nx;
                M->A[r][5] = ny;
                M->A[r][6] = nz;
                M->A[r][7] = 0.0;
            }
            MyMath_GaussJordan_matrix_inv(M);

            double rhs = 0.0;
            if (self_included) {

                rhs = M->A_inv[7][7];
            }
            ib_node->boundary_coef = -(ds+df) + rhs;
            r=0;
            for (g=0; g<n_fluid; g++) {

                if ( (ib_node->fluid_index[g].x_index == ib_node->im_index.x_index)
                     && (ib_node->fluid_index[g].y_index == ib_node->im_index.y_index)
                     && (ib_node->fluid_index[g].z_index == ib_node->im_index.z_index) ) {
                    ib_node->fluid_coef[g] = 0.0;
                } else {
                    ib_node->fluid_coef[g] = -M->A_inv[7][r];
                    r++;
                }
            } /* for g */
            free(M);
        }
    } /* NEUMANN */
        break;

        /* dphi/dn + K*phi = b; at (xcont,ycont,zcont) */
        /* where dphi_dn = dphi/dx*nx + dphi/dy*ny + dphi/dz*nz */
        /* n=(nx,ny,nz) unit normal vector pointing toward the fluid region */

        /* Set the coefficient used for the Robbin boundary condition */
        /* Seems like K should take a negative value for stability reasons */
    case ROBBIN: {

        /* For now, we assume a constant value for K */
        double K = -1.0;
        double delta = df + ds;
        double w = ds / delta;
        double temp = (1.0 + K*delta*w) / (K * delta * (1.0-w) - 1.0);

        /* If non-zero Robbin B.C. is imposed (b!=0), this coefficient is used to account for the nonzero term in the RHS of the linear system */
        ib_node->boundary_coef = delta / (K * delta * (1.0-w) - 1.0);
        for (g=0; g<n_fluid; g++) {

            ib_node->fluid_coef[g] = weight_coef[g] * temp; /*dphi/dn + K*phi = b at (xc,yc,zc) */
        } /* for g */
    } /* ROBBIN */
        break;

    default:
        PetscPrintf(PCW, "Immersed.c/ Unknown boundary condition type...\n");
    } /* switch */

    /* now check if all the coefficients were computed correctly. */
    for (g=0; g<n_fluid; g++) {
        if ( (isnan(ib_node->fluid_coef[g]) ) || ( isinf(ib_node->fluid_coef[g]) ) ) {

            printf("Immersed.c/ Warning... Invalid coefficient for the %c-immersed node. Setting it to zero.\n", which_quantity);

            ib_node->fluid_coef[g] = 0.0;
            ib_node->boundary_coef = 0.0;
            getchar();
        } /* if */
    } /* for g */
}
/**************************************************************************************************************/

/* This function sets both global and xyz indices of the immersed nodes in the domain on the current processor */
/* TBC */
void Immersed_set_q_immersed_indices(MAC_grid *grid, char which_quantity) {

    Immersed *q_immersed=NULL;
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed  = grid->u_immersed;
        break;

    case 'v':
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed  = grid->v_immersed;
        break;

    case 'w':
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed  = grid->w_immersed;
        break;

    case 'c':
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed  = grid->c_immersed;
        break;

    default:
        PetscPrintf(PCW, "Immersed.c/ Cannot count immersed nodes. Invalid quantity %c\n", which_quantity);
    } /* switch */

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    int i, j, k;
    int g_index = 0;
    /* Go through all the nodes in the physical domain on the current processor and set the physical indices for each node and also the */
    /* the global index of each immersed node. Global index is simply the number of the ib node */
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (Grid_get_q_status(grid, i, j, k) == IMMERSED) {

                    /* Set the global index of the current node */
                    Immersed_set_ib_global_index(q_immersed, i, j, k, g_index);

                    /* Set the x and y indices of the current ib node */
                    Immersed_set_ib_xyz_index(&q_immersed->ib_nodes[g_index], i, j, k);

                    g_index++;

                } else { /* Not an immersed node. Set the global index to -1 */

                    /* Set null for the current node which is not an Immersed node */
                    Immersed_set_ib_global_index(q_immersed, i, j, k, NOT_FOUND);
                } /* else */
            } /* for i */
        } /* for j*/
    } /* for k*/
}
/**************************************************************************************************************/

/* This function sets the global index for the current immersed node located at (x_index,y_index) in the physical domain */
/* Note that the Immersed nodes go from 0--->global_index-1 */
void Immersed_set_ib_global_index(Immersed *q_immersed, int x_index, int y_index, int z_index, int global_index) {

    q_immersed->global_index[z_index][y_index][x_index] = global_index;
}
/**************************************************************************************************************/


/* This function returns the global index for the current immersed node located at (x_index,y_index) in the physical domain */
/* Note that the Immersed nodes go from 0--->global_index-1 */
int Immersed_get_ib_global_index(Immersed *q_immersed, int x_index, int y_index, int z_index) {

    return q_immersed->global_index[z_index][y_index][x_index];
}
/*************************************************************************************************************/

/* This function sets the (x,y,z) index of the current IB node */
void Immersed_set_ib_xyz_index(ImmersedNode *ib_node, int x_index, int y_index, int z_index) {

    ib_node->im_index.x_index = x_index;
    ib_node->im_index.y_index = y_index;
    ib_node->im_index.z_index = z_index;
}
/**************************************************************************************************************/

/* This function sets the coordinates of the immersed nodes */
void Immersed_set_q_immersed_coordinates(MAC_grid *grid, char which_quantity) {

    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;

    Immersed *q_immersed=NULL;
    switch (which_quantity) {

    case 'u':

        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        q_immersed = grid->u_immersed;
        break;

    case 'v':

        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        q_immersed = grid->v_immersed;
        break;

    case 'w':

        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        q_immersed = grid->w_immersed;
        break;

    case 'c':

        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        q_immersed = grid->c_immersed;
        break;

    default:

        PetscPrintf(PCW, "Immersed.c/ Error setting coordinates. Unknown quantity\n");
    } /* switch */

    /* Total number of immersed nodes on current processor */
    int N = q_immersed->N;
    int g;
    ImmersedNode *ib_node;
    /* Go through all the immersed nodes and set the coordinates of the immersed nodes */
    for (g=0; g<N; g++) {

        /* Get the current immersed node */
        ib_node = Immersed_get_ib_node(q_immersed, g);

        /* Index of the current ib_node */
        int i_im = ib_node->im_index.x_index;
        int j_im = ib_node->im_index.y_index;
        int k_im = ib_node->im_index.z_index;

        /* Set the coordinates of the current ib_node based on the grid position of the given quantity */
        ib_node->im_point.x = xq[i_im];
        ib_node->im_point.y = yq[j_im];
        ib_node->im_point.z = zq[k_im];

    } /* for g */
} 
/**************************************************************************************************************/

/* This function sets up all the necessary functions for the immersed nodes */
void Immersed_setup_q_immersed_nodes(MAC_grid *grid, char which_quantity) {

    Immersed *q_immersed;

    /* First, allocate memory for the q_immersed nodes */
    q_immersed = Immersed_create(grid, which_quantity);

    /* Pass the pointer to the grid structure */
    switch (which_quantity) {

    case 'u':

        grid->u_immersed = q_immersed;
        break;

    case 'v':

        grid->v_immersed = q_immersed;
        break;

    case 'w':

        grid->w_immersed = q_immersed;
        break;

    case 'c':

        grid->c_immersed = q_immersed;
        break;

    default:
        PetscPrintf(PCW, "Immeresed.c /Can not setup immersed nodes. Unknown quantity %c.\n", which_quantity);
    } /* switch */

    /* Set the xyz and global indices of the immersed nodes in the physical domain */
    Immersed_set_q_immersed_indices(grid, which_quantity);

    /* Set the coordinates of the immersed nodes */
    Immersed_set_q_immersed_coordinates(grid, which_quantity);

    /* find the unit normal vector to the surface */
    Immersed_find_surface_normals(grid, which_quantity);

    /* Set the control points on the boundary associated with each ib nodes */
    Immersed_set_q_immersed_control_points(grid, which_quantity);
    //printf("Immersed.c/ setup -1\n");


    /* Find the image node for the immersed nodes in the interior region (fluid region) */
    Immersed_set_q_immersed_image_points(grid, which_quantity);
    //printf("Immersed.c/ setup 0\n");


    /* Set the appropriate boundary condition for the imersed nodes. u,v: Dirichlet. c: Neumann */
    Immersed_set_q_immersed_boundary_condition(grid, which_quantity);
   // printf("Immersed.c/ setup 1\n");

    /* Now, set the neighboring nodes for each ib node. Ideally, for each ib_node, we should be able to find two neighbors */
    Immersed_set_q_immersed_neighboring_fluid_nodes(grid, which_quantity);
    //printf("Immersed.c/ setup 2\n");

    /* Now, using all the data, use a linear interpolation to find the relationship between the */
    /* immersed nodes and the fluid nodes and the correct boundary condition */
    Immersed_set_q_interpolation_coef(grid, which_quantity);

    //printf("Immersed.c/ setup 3\n");
    if (which_quantity == 'c') {
        //printf("Immersed.c/ stage 5\n");
        //Display_immersed_node(&q_immersed->ib_nodes[0]);
        //getchar();
    }
}
/**************************************************************************************************************/

/* This function finds the fluid neighboring nodes for each immersed node. Stores the indices and the state of the node comparing to the current node */
void Immersed_set_q_immersed_neighboring_fluid_nodes(MAC_grid *grid, char which_quantity) {

    Immersed *q_immersed=NULL;
    switch (which_quantity) {

    case 'u':
        q_immersed = grid->u_immersed;
        break;

    case 'v':
        q_immersed = grid->v_immersed;
        break;

    case 'w':
        q_immersed = grid->w_immersed;
        break;

    case 'c':
        q_immersed = grid->c_immersed;
        break;

    default:

        PetscPrintf(PCW, "Immersed.c/ Could not set the neighboring immersed nodes. Unknown quantity\n");
    }  /* switch */

    /* Total number of the immersed nodes */
    int N = q_immersed->N;
    int g;
    ImmersedNode *ib_node;

    /* Go through all the ib nodes and find the neighboring fluid nodes. Store the indices and the coordinates */
    for (g=0; g<N; g++) {

        /* Get the current immersed node */
        ib_node = Immersed_get_ib_node(q_immersed, g);

        /* Set the fluid nodes neighboring the current ib node */
        Immersed_assign_q_neighboring_fluid_nodes(ib_node, grid, which_quantity);
    } /* for i */
} 
/**************************************************************************************************************/

/* This function computes the interpolation coefficients relating the value of the immersed node to the fluid neighboring nodes in a way that it satisfies the boundary conditions */
void Immersed_set_q_interpolation_coef(MAC_grid *grid, char which_quantity) {

    Immersed *q_immersed=NULL;
    switch (which_quantity) {

    case 'u':
        q_immersed = grid->u_immersed;
        break;

    case 'v':
        q_immersed = grid->v_immersed;
        break;

    case 'w':
        q_immersed = grid->w_immersed;
        break;

    case 'c':
        q_immersed = grid->c_immersed;
        break;

    default:

        PetscPrintf(PCW, "Immersed.c/ Could not find the inerpolation coefficient. Unknown quantity\n");
    }  /* switch */


    /* Total number of the immersed nodes */
    int N = q_immersed->N;
    int g;
    ImmersedNode *ib_node;

    /* Go through all the ib nodes and find the neighboring fluid nodes. Store the indices and the coordinates */
    for (g=0; g<N; g++) {

        /* Get the current immersed node */
        ib_node = Immersed_get_ib_node(q_immersed, g);

        /* Find the interpolation coefficient for the current in node */
        Immersed_compute_interpolation_coef(ib_node, which_quantity);

    } /* for g */
} 
/**************************************************************************************************************/

/* This function assigns the neighboring fluid nodes and positions of them for the current immersed boundary node */
void Immersed_assign_q_neighboring_fluid_nodes(ImmersedNode *ib_node, MAC_grid *grid, char which_quantity) {

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status;
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        break;

    case 'v':
        Grid_get_q_status = &Grid_get_v_status;
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        break;

    case 'w':
        Grid_get_q_status = &Grid_get_w_status;
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        break;

    case 'c':
        Grid_get_q_status = &Grid_get_c_status;
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        break;

    default:

        PetscPrintf(PCW, "Immersed.c/ Could not set the neighboring immersed nodes. Unknown quantity %c\n", which_quantity);
    }  /* switch */

    static const short int dir_x[IBM_MAX] = {0, 1, 0, 1, 0, 1, 0, 1};
    static const short int dir_y[IBM_MAX] = {0, 0, 1, 1, 0, 0, 1, 1};
    static const short int dir_z[IBM_MAX] = {0, 0, 0, 0, 1, 1, 1, 1};

    /* Always follow this standard for the 8 neigboring fluid nodes. This standard is used to find the interpolation coefficients using a trilinear interpolation */

    /*
        F4------------F6
       / |           /|
      /  |          / |
     /   |         /  |
   F5------------F7   |
   |   F0---------|---F2
   |  /	 I     |  /
            | /            | /
            |/	           |/
   F1-------------F3

     z
     |
  |___y
 /
   /
  x
*/
    /* Indices of grid nearest to the image point */
    /* Return the index of the F0 node (left-bottom-back) index which represnts the box surroudning the image node */
    Indices Fluid_index = Immersed_find_box(ib_node, &ib_node->image_point, grid, which_quantity);

    int i_ig = Fluid_index.x_index;
    int j_ig = Fluid_index.y_index;
    int k_ig = Fluid_index.z_index;

    //printf("Immersed.c/ Done box (i,j,k)=(%d,%d,%d) n_fluid:%d\n", i_ig, j_ig, k_ig, ib_node->n_fluid);
    if (i_ig != NOT_FOUND) {

        PointType Fluid_point;

        /* A full box is found which surrounds the image node associated with the current immersed node */
        if (ib_node->n_fluid == IBM_MAX) { /* IBM_MAX = 8 */

            int m;
            for (m=0; m<IBM_MAX; m++) {

                int i_f = i_ig + dir_x[m];
                int j_f = j_ig + dir_y[m];
                int k_f = k_ig + dir_z[m];
                //printf("Immersed.c/ m:%d (i,j,k)=(%d,%d,%d)\n", m, i_f, j_f, k_f);

                Fluid_index.x_index = i_f;
                Fluid_index.y_index = j_f;
                Fluid_index.z_index = k_f;

                Fluid_point.x = xq[i_f];
                Fluid_point.y = yq[j_f];
                Fluid_point.z = zq[k_f];

                ib_node->fluid_point[m] = Fluid_point;
                ib_node->fluid_index[m] = Fluid_index;

            } /* for m */
        }  else if (ib_node->n_fluid == 1) {

            Fluid_index.x_index = i_ig;
            Fluid_index.y_index = j_ig;
            Fluid_index.z_index = k_ig;

            Fluid_point.x = xq[i_ig];
            Fluid_point.y = yq[j_ig];
            Fluid_point.z = zq[k_ig];

            ib_node->fluid_point[0] = Fluid_point;
            ib_node->fluid_index[0] = Fluid_index;

        } /* else if */
    } /* if */else {

        printf("Immersed.c/ Warning. Not a single fluid node was found for the current %c-immersed node at (i,j,k)=(%d,%d,%d) \n", which_quantity, ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index);
        Display_immersed_node(ib_node);
        getchar();
    } /* else */
}
/**************************************************************************************************************/

/* This function sets the correct boundary condition for the immersed nodes of each quantity */
/* Default, u,v: Dirichlet. c: Neumann */
void Immersed_set_q_immersed_boundary_condition(MAC_grid *grid, char which_quantity) {

    Immersed *q_immersed=NULL;

    short int boundary_condition=-1;
    switch (which_quantity) {

    case 'u':
        q_immersed = grid->u_immersed;
        boundary_condition = DIRICHLET;
        break;

    case 'v':
        q_immersed = grid->v_immersed;
        boundary_condition = DIRICHLET;
        break;

    case 'w':
        q_immersed = grid->w_immersed;
        boundary_condition = DIRICHLET;
        break;

    case 'c':
        q_immersed = grid->c_immersed;
        boundary_condition = NEUMANN;
        break;

    default:
        PetscPrintf(PCW, "Immersed.c /Error setting the correct boundary condition for the boundary points. Unknown quantity...\n");
    } /* switch */

    int N = q_immersed->N;
    int g;

    for (g=0; g<N; g++) {

        /* Set the boundary condition fot all the immersed nodes */
        /* Note that for now, all the nodes of the same quantitiy should take the same boundary condition */
        /* One needs to use a more sophisticated approach to impose a combination of those boundary condition */

        q_immersed->ib_nodes[g].boundary_condition = boundary_condition;
    }/* for i */
}
/**************************************************************************************************************/

/* This function returns the i'th immersed node */
ImmersedNode *Immersed_get_ib_node(Immersed *q_immersed, int i) {

    return &(q_immersed->ib_nodes[i]);
}
/**************************************************************************************************************/

/* This function returns the point on the boundary for the current IB node based on the given neighboring fluid node */
PointType *Immersed_find_control_point(ImmersedNode *ib_node, MAC_grid *grid, char which_quantity) {

    double ***q_sdf=NULL;
    switch (which_quantity) {

    case 'u':
        q_sdf = grid->surf->u_sdf;
        break;
    case 'v':
        q_sdf = grid->surf->v_sdf;
        break;
    case 'w':
        q_sdf = grid->surf->w_sdf;
        break;
    case 'c':
        q_sdf = grid->surf->c_sdf;
        break;

    } /* switch */

    /* coordinates of the immersed node */
    double x_im = ib_node->im_point.x;
    double y_im = ib_node->im_point.y;
    double z_im = ib_node->im_point.z;

    /* (i,j,k) of the current ib node */
    int i_im = ib_node->im_index.x_index;
    int j_im = ib_node->im_index.y_index;
    int k_im = ib_node->im_index.z_index;

    /* initial value for the control point */
    double xcont = x_im;
    double ycont = y_im;
    double zcont = z_im;

    double nx = ib_node->n.vx;
    double ny = ib_node->n.vy;
    double nz = ib_node->n.vz;

    double D = Surface_find_Dn_to_surface(q_sdf, grid, i_im, j_im, k_im, &ib_node->n, which_quantity);

    /* Now using the Distance and the vector, find the location of the control point */
    xcont = x_im - nx*D;
    ycont = y_im - ny*D;
    zcont = z_im - nz*D;

    if (which_quantity == 'c') {

        //printf("Immersed.c/ (i,j,k)=(%d,%d,%d) D:%f (nx,ny,nz):(%f,%f,%f) (xim,yim,zim)=(%f,%f,%f)\n", i_im, j_im, k_im, D, nx, ny, nz, x_im, y_im, z_im);
        //printf("Immersed.c/ (xc,yc,zc)=(%f,%f,%f)\n", xcont, ycont, zcont);
        //getchar();
    }

    PointType *control_point = (PointType *)calloc(1, sizeof(PointType));

    control_point->x = xcont;
    control_point->y = ycont;
    control_point->z = zcont;

    return control_point;
}
/**************************************************************************************************************/

/* This function sets the control point for each quantity */
/* By control point, we mean that the point lies on the boundary. It is used for interpolating the value of the immersed node. */
void Immersed_set_q_immersed_control_points(MAC_grid *grid, char which_quantity) {

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    Immersed *q_immersed=NULL;
    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed  = grid->u_immersed;
        break;

    case 'v':
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed  = grid->v_immersed;
        break;

    case 'w':
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed  = grid->w_immersed;
        break;

    case 'c':
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed  = grid->c_immersed;
        break;

    } /* switch */

    /* Total number of the immersed nodes */
    int N = q_immersed->N;

    /* Go through all the nodes in the physical domain and set the physical indices for each node and also the */
    /* And the global index of each immersed node. Global index is simply the number of the ib node */
    int ib;
    for (ib=0; ib<N; ib++) {

        ImmersedNode *ib_node = Immersed_get_ib_node(q_immersed, ib);

        /* Find the control point associated with the current fluid-immersed node */
        PointType *point_c = Immersed_find_control_point(ib_node, grid, which_quantity);

        /* Pass the coordinates of the control points to the immersed node. */
        ib_node->control_point.x = point_c->x;
        ib_node->control_point.y = point_c->y;
        ib_node->control_point.z = point_c->z;

        free(point_c);

    } /* for ib */
}
/**************************************************************************************************************/

/* This function sets the image point for each quantity */
/* By image point, we mean that the mirrored point of the ghost node which is inside the fluid region */
void Immersed_set_q_immersed_image_points(MAC_grid *grid, char which_quantity) {

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    Immersed *q_immersed=NULL;
    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed  = grid->u_immersed;
        break;

    case 'v':
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed  = grid->v_immersed;
        break;

    case 'w':
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed  = grid->w_immersed;
        break;

    case 'c':
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed  = grid->c_immersed;
        break;

    } /* switch */

    /* Total number of the immersed nodes */
    int N = q_immersed->N;

    /* Go through all the nodes in the physical domain and set the physical indices for each node and also the */
    /* And the global index of each immersed node. Global index is simply the number of the ib node */
    int ib;
    for (ib=0; ib<N; ib++) {

        ImmersedNode *ib_node = Immersed_get_ib_node(q_immersed, ib);

        /* Find the image point associated with the current fluid-immersed node */
        PointType *point_i = Immersed_find_image_point(ib_node, grid, which_quantity);

        /* Pass the coordinates of the control points to the immersed node. */
        ib_node->image_point.x = point_i->x;
        ib_node->image_point.y = point_i->y;
        ib_node->image_point.z = point_i->z;

        free(point_i);

    } /* for ib */
}
/**************************************************************************************************************/

/* This function returns the image point of the current IB node based on normal direction to the surface */
PointType *Immersed_find_image_point(ImmersedNode *ib_node, MAC_grid *grid, char which_quantity) {

    PointType *image_point;
    double ***q_sdf;
    double ximage, yimage, zimage;
    double x_im, y_im, z_im;
    int i_im, j_im, k_im;
    double nx, ny, nz;
    double dx, dy, dz;
    double D, ds;
    double epsilon;
    image_point = (PointType *)calloc(1, sizeof(PointType));

    switch (which_quantity) {

    case 'u':
        q_sdf = grid->surf->u_sdf;
        break;
    case 'v':
        q_sdf = grid->surf->v_sdf;
        break;
    case 'w':
        q_sdf = grid->surf->w_sdf;
        break;
    case 'c':
        q_sdf = grid->surf->c_sdf;
        break;

    } /* switch */

    /* coordinates of the immersed node */
    x_im = ib_node->im_point.x;
    y_im = ib_node->im_point.y;
    z_im = ib_node->im_point.z;

    /* (i,j,k) of the current ib node */
    i_im = ib_node->im_index.x_index;
    j_im = ib_node->im_index.y_index;
    k_im = ib_node->im_index.z_index;

    nx = ib_node->n.vx;
    ny = ib_node->n.vy;
    nz = ib_node->n.vz;

    dx = grid->dx_c[i_im];
    dy = grid->dy_c[j_im];
    dz = grid->dz_c[k_im];

    ds = MyMath_get_point_point_distance(&ib_node->im_point, &ib_node->control_point);

    /* This is to ensure that the image point is far from the current surrounding box and lies within the next box */
    D  = 2.0*ds;

    /* First image point: */
    ximage = ib_node->im_point.x + nx*D;
    yimage = ib_node->im_point.y + ny*D;
    zimage = ib_node->im_point.z + nz*D;

    image_point->x = ximage;
    image_point->y = yimage;
    image_point->z = zimage;


    /* Two image points. They could be used for interpolation. They are not used in the immersed boundary interpolation. */
    epsilon = max(0.001*dx, 1e-12);
    double Dtest= min( fabs(dx/nx), fabs(dy/ny)) + epsilon;
    Dtest = min(Dtest, fabs(dz/nz));
//My debug
    Dtest = max(D, Dtest);
    Dtest = max(dx, dy);
    Dtest = max(dz, Dtest);

    /* Assume that the grid is uniform close to the solid surface, d = sqrt(dx^2+dy^2+dz^2) */
    Dtest = sqrt(3.0) * dy;
    ib_node->image_point1.x = ib_node->control_point.x + nx*Dtest;
    ib_node->image_point1.y = ib_node->control_point.y + ny*Dtest;
    ib_node->image_point1.z = ib_node->control_point.z + nz*Dtest;

    ib_node->image_point2.x = ib_node->control_point.x + nx*2.0*Dtest;
    ib_node->image_point2.y = ib_node->control_point.y + ny*2.0*Dtest;
    ib_node->image_point2.z = ib_node->control_point.z + nz*2.0*Dtest;

    //printf("Immersed.c/ (x_im,y_im)=(%f,%f) (xig,yig)=(%f,%f) (nx,ny)=(%f,%f) D:%f\n", x_im, y_im, ximage, yimage, nx, ny, D);
    //getchar();

    return image_point;
}
/**************************************************************************************************************/

/* This function sets the normal unit vector to the surface. It is approximated by the vector at the current immersed node */
void Immersed_find_surface_normals(MAC_grid *grid, char which_quantity) {

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    Immersed *q_immersed=NULL;
    double ***q_sdf=NULL;
    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed  = grid->u_immersed;
        q_sdf = grid->surf->u_sdf;

        break;

    case 'v':
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed  = grid->v_immersed;
        q_sdf = grid->surf->v_sdf;
        break;


    case 'w':
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed  = grid->w_immersed;
        q_sdf = grid->surf->w_sdf;
        break;

    case 'c':
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed  = grid->c_immersed;
        q_sdf = grid->surf->c_sdf;
        break;
    } /* switch */

    ImmersedNode *ib_node;
    /* Total number of the immersed nodes */
    int N = q_immersed->N;
    int ib;
    /* find the normal unit vector at the current location of the immersed node. */
    for (ib=0; ib<N; ib++) {

        ib_node = Immersed_get_ib_node(q_immersed, ib);

        int i_im = ib_node->im_index.x_index;
        int j_im = ib_node->im_index.y_index;
        int k_im = ib_node->im_index.z_index;

        /* n=(nx,ny,nz) vector normal to the surface*/
        ib_node->n = Surface_compute_normal(q_sdf, grid, i_im, j_im, k_im, which_quantity) ;

    } /* for ib */


}
/**************************************************************************************************************/

/* This function returns the index of the back_south_west corner of the node for the box of 8 fluid nodes surrounding one arbitraty node at the given location p */
Indices  Immersed_find_box(ImmersedNode *ib_node, PointType *p, MAC_grid *grid, char which_quantity) {

    static const short int dir_x[8] = {0, 1, 0, 1, 0, 1, 0, 1};
    static const short int dir_y[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    static const short int dir_z[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    Indices Bottom_Left_Back_index;

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;

    switch (which_quantity) {

    case 'u':
        Grid_get_q_status = &Grid_get_u_status;
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        break;
    case 'v':
        Grid_get_q_status = &Grid_get_v_status;
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        break;
    case 'w':
        Grid_get_q_status = &Grid_get_w_status;
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        break;
    case 'c':
        Grid_get_q_status = &Grid_get_c_status;
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        break;
    } /* switch */

    Bottom_Left_Back_index.x_index = NOT_FOUND;
    Bottom_Left_Back_index.y_index = NOT_FOUND;
    Bottom_Left_Back_index.z_index = NOT_FOUND;

    /* coordinates of the current point */
    double x = p->x;
    double y = p->y;
    double z = p->z;

    /* To avoid the problems with the double-precision if (==), we perturb the current image node in 8 direction and see if it lies within a box with 8 fluid (or immersed) nodes at its corners */
    const double tol_x = max(0.005*grid->dx_c[ib_node->im_index.x_index], 1e-12);
    const double tol_y = max(0.005*grid->dy_c[ib_node->im_index.y_index], 1e-12);
    const double tol_z = max(0.005*grid->dz_c[ib_node->im_index.z_index], 1e-12);
    const double pert_x[8]={+tol_x, -tol_x, +tol_x, -tol_x, +tol_x, -tol_x, +tol_x, -tol_x};
    const double pert_y[8]={+tol_y, +tol_y, -tol_y, -tol_y, +tol_y, +tol_y, -tol_y, -tol_y};
    const double pert_z[8]={+tol_z, +tol_z, +tol_z, +tol_z, -tol_z, -tol_z, -tol_z, -tol_z};
    double x_p, y_p, z_p;

    int g, m = 0;
    short int found = NO;
    int i_ig, j_ig, k_ig;
    int i_, j_, k_;

    while ( (!found) && (m<8) ) {

        x_p = x + pert_x[m];
        y_p = y + pert_y[m];
        z_p = z + pert_z[m];

        i_ig = Grid_get_x_index(x_p, grid, which_quantity);
        j_ig = Grid_get_y_index(y_p, grid, which_quantity);
        k_ig = Grid_get_z_index(z_p, grid, which_quantity);

        short int is_good_box = YES;
        /* Check if the node is inside the box with the south-west-back corner node at (i_ig,j_ig,k_ig) */
        /* This box should have all fluid (or immersed) nodes at each corner */
        /* Immersed node can be the current node by itself */
        for (g=0; g<8; g++) {

            i_ = i_ig + dir_x[g];
            j_ = j_ig + dir_y[g];
            k_ = k_ig + dir_z[g];

            /* If one needs to exclude the immersed node from interplation, this if statement should be re-inserted */
            //if ( ( i_ != ib_node->im_index.x_index) ||
            //     ( j_ != ib_node->im_index.y_index) ||
            //     ( k_ != ib_node->im_index.z_index) ) {

            int status = Grid_get_q_status(grid, i_, j_, k_);
            if ( (status != FLUID) && (status != IMMERSED) ) {

                is_good_box = NO; break;
            } /* if */
            // } /* if */
        } /* for g */

        found = is_good_box;

        if (found) {

            if (
                    (x >= xq[i_ig]) && (x <= xq[i_ig+1]) &&
                    (y >= yq[j_ig]) && (y <= yq[j_ig+1]) &&
                    (z >= zq[k_ig]) && (z <= zq[k_ig+1])
               ) {

            } else {
                found = NO;
            }
        }

        m++;
    } /* while */

    if (found) {

        ib_node->n_fluid = 8;
        Bottom_Left_Back_index.x_index = i_ig;
        Bottom_Left_Back_index.y_index = j_ig;
        Bottom_Left_Back_index.z_index = k_ig;
    } else {

        printf("Immersed.c/ Warning!!! A box is not found to surround the immersed node at ig(i,j,k)=(%d,%d,%d) \n", ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index);
        printf("Immersed.c/ Now, assigning only the closest fluid node for the interpolation\n");

        PointType p_;
        double D_ig_p_;
        double D_min = 1e12;
        int g_found = NOT_FOUND;

        i_ig = Grid_get_x_index(x, grid, which_quantity);
        j_ig = Grid_get_y_index(y, grid, which_quantity);
        k_ig = Grid_get_z_index(z, grid, which_quantity);

        i_ig = max(i_ig, 0);
        j_ig = max(j_ig, 0);
        k_ig = max(k_ig, 0);
        i_ig = min(i_ig, grid->NI-1);
        j_ig = min(j_ig, grid->NJ-1);
        k_ig = min(k_ig, grid->NK-1);

        if ( (i_ig == NOT_FOUND) || (j_ig == NOT_FOUND) || (k_ig == NOT_FOUND) ||
             (i_ig == GVG_INF) || (j_ig == GVG_INF) || (k_ig == GVG_INF) ) {

            x = ib_node->im_point.x;
            y = ib_node->im_point.y;
            z = ib_node->im_point.z;

            i_ig = ib_node->im_index.x_index;
            j_ig = ib_node->im_index.y_index;
            k_ig = ib_node->im_index.z_index;
        } /* if */

        /* Find the closest fluid node to the image node and assign its value to the image node. This is not the most accurate thing to do, but could happen in very rare situations such as nodes close to the solid walls */
        for (g=0; g<8; g++) {

            i_ = i_ig + dir_x[g];
            j_ = j_ig + dir_y[g];
            k_ = k_ig + dir_z[g];
            //printf("im_index(i,j,k=(%d,%d,%d) neigh(i,j,k)=(%d,%d,%d) status:%d\n", ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index, i_, j_, k_, Grid_get_q_status(grid, i_, j_, k_));

            if (Grid_get_q_status(grid, i_, j_, k_) == FLUID) {

                p_.x = xq[i_];
                p_.y = yq[j_];
                p_.z = zq[k_];
                D_ig_p_ = MyMath_get_point_point_distance(p, &p_);

                if (D_ig_p_ < D_min) {

                    D_min = D_ig_p_;
                    g_found = g;
                }
                found = YES;
            }
        } /* for g */
        if (found) {

            /* reduce the number of the neighboring fluid nodes found to only "one" for this current immersed node */
            /* This is to avoid the code to crash */
            /* To prevent this to happen, try to avoid very sharp solid surfaces, and also watch for the box boundaries */
            ib_node->n_fluid = 1;
            Bottom_Left_Back_index.x_index = i_ig + dir_x[g_found];
            Bottom_Left_Back_index.y_index = j_ig + dir_y[g_found];
            Bottom_Left_Back_index.z_index = k_ig + dir_z[g_found];
        } /* if */ else {

            printf("Immersed.c/ Warning!!! Could not find even one fluid node neighboring current image node. Possible reason could be due to the mistakes in the definition of the solid surface. \n");
            ib_node->n_fluid = 0;
        }
    }
    return Bottom_Left_Back_index;
}
/**************************************************************************************************************/

void Immersed_is_valid(ImmersedNode *ib_node) {

    double s = 0.0;
    if (ib_node->boundary_condition == NEUMANN) {
        s = -1.0;
    } else if (ib_node->boundary_condition == DIRICHLET) {
        s = 1.0;
    }

    short int bad_node = NO;
    int g;
    for (g=0; g< ib_node->n_fluid; g++) {
        if ( (s*ib_node->fluid_coef[g]) < 0.0) {
            bad_node = YES;
            break;
        }
    }
    if (bad_node) {
        printf("Immersed.c/ Warning. IBM coefficient takes the sign that might cause the linear system go unstable.\n");
        Display_immersed_node(ib_node);
        getchar();
    }
}
