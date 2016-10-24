/* Note: Function names are in alphabetical order (Not exactly) */
/* grid.c
  This file defines the functions to create MAC_grid objects
*/
#include "definitions.h"
#include "DataTypes.h"
#include "Grid.h"
#include "Memory.h"
#include "Surface.h"
#include "Display.h"
#include "Debugger.h"
#include "Communication.h"
#include "Immersed.h"
#include "MyMath.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* This function allocates momory for the grid object */
/* It also holds the object holding the information regarding solid interface and solid object */
MAC_grid *Grid_create(Parameters *params) {

    double *xc, *yc, *zc;
    double *xu, *yu, *zu;
    double *xv, *yv, *zv;
    double *xw, *yw, *zw;
    double *dx_c, *dy_c, *dz_c;
    double *metric_xc, *metric_yc, *metric_zc;
    double *metric_xu, *metric_yv, *metric_zw;
    int NX_cell, NY_cell, NZ_cell;
    int NX, NY, NZ;
    double Lx, Ly, Lz;
    double dx, dy, dz;
    double dEta, dXi, dPhi;
    double **interface_position;
    int **interface_y_index;
    MAC_grid *new_grid;
    int Ghost_Nodes;
    int ierr;
    int Is, Js, Ks;
    int nx, ny, nz;
    int Is_g, Js_g, Ks_g;
    int nx_g, ny_g, nz_g;
    int success_flag;

    new_grid = (MAC_grid *)malloc(sizeof(MAC_grid));
    Memory_check_allocation(new_grid);

    /* Domain length */
    new_grid->Lx = params->Lx;
    new_grid->Ly = params->Ly;
    new_grid->Lz = params->Lz;

    /* Real physical number of cell_grids */
    NX_cell = params->NX;
    NY_cell = params->NY;
    NZ_cell = params->NZ;

    /* Global (world) number of grids. A half cell is added in each direction. Same on all processors */
    NX = NX_cell + 1;
    NY = NY_cell + 1;
    NZ = NZ_cell + 1;

    /* To change this, go to file "definition.h" */
    /* Allow PETSc to perform the domain decomposition */
#ifdef AUTOMATIC_DOMAIN_PARTITIONING
    /*******************************************************************/
    /* Automatic (PETSc) domain decompostion */
    /*******************************************************************/
    /* Create DA for the ghost nodes associated with the ENO scheme (convective terms) */
    Ghost_Nodes = params->ENO_scheme_order;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_STAR, NX, NY, NZ, PETSC_DECIDE, PETSC_DECIDE,
                      PETSC_DECIDE, 1, Ghost_Nodes, PETSC_NULL, PETSC_NULL, PETSC_NULL, &new_grid->DA_3D);

    /* Create DA associated with tagging the nodes (flags: solid, fluid, immersed) */
    /* Note that using Immersed Boundary Method, one has to used BOX stencil */
    Ghost_Nodes = 3;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_BOX, NX, NY, NZ, PETSC_DECIDE, PETSC_DECIDE,
                      PETSC_DECIDE, 1, Ghost_Nodes, PETSC_NULL, PETSC_NULL, PETSC_NULL, &new_grid->DA_3D_BOX);
    PETScErrAct(ierr);

    /* Create DA associate with the linear systems. Uses a 7-point finite difference scheme. */
    Ghost_Nodes = 1;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_STAR, NX, NY, NZ, PETSC_DECIDE, PETSC_DECIDE,
                      PETSC_DECIDE, 1, Ghost_Nodes, PETSC_NULL, PETSC_NULL, PETSC_NULL, &new_grid->DA_3D_Lsys);
#else /* MANUAL_DOMAIN_PARTITIOING */

    /* Number of processors in each direction */
    /* Each processors owns the same number of grid nodes of the subdomain */
    int npx = 4;//params->size/15;
    int npy = 3;
    int npz = 3;

	int nnx[]={46,65,30,30};
	int nny[]={48,48,49};
	int nnz[]={48,48,49};

	if (0) {

    Ghost_Nodes = params->ENO_scheme_order;
    int success1 = Communication_create_DA3D(npx, npy, npz, DA_STENCIL_STAR, Ghost_Nodes, params, &new_grid->DA_3D);

    Ghost_Nodes = 3;
    int success2 = Communication_create_DA3D(npx, npy, npz, DA_STENCIL_BOX, Ghost_Nodes, params, &new_grid->DA_3D_BOX);

    Ghost_Nodes = 1;
    int success3 = Communication_create_DA3D(npx, npy, npz, DA_STENCIL_STAR, Ghost_Nodes, params, &new_grid->DA_3D_Lsys);

    if ( (!success1) || (!success2) || (!success3) ) {

        PetscPrintf(PCW, "Grid.c/ Error!!! Could not perform manual domain decomposition using (npx,npy,npz)=(%d,%d,%d) processors.\n Trying assigning the right number of processors in each direction.\n Or you could allow PETSc do that automatically.\n", npx, npy, npz);
    } /* if */
	}

    Ghost_Nodes = params->ENO_scheme_order;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_STAR, NX, NY, NZ, npx, npy, npz, 1, Ghost_Nodes, nnx, nny, nnz, &new_grid->DA_3D); 

    Ghost_Nodes = 3;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_STAR, NX, NY, NZ, npx, npy, npz, 1, Ghost_Nodes, nnx, nny, nnz, &new_grid->DA_3D_BOX); 

    Ghost_Nodes = 1;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_STAR, NX, NY, NZ, npx, npy, npz, 1, Ghost_Nodes, nnx, nny, nnz, &new_grid->DA_3D_Lsys); 

#endif

    /* Now, get the start corner index and number of nodes on the current processor */
    ierr = DAGetCorners(new_grid->DA_3D, &Is, &Js, &Ks, &nx, &ny, &nz); PETScErrAct(ierr);

    DAView(new_grid->DA_3D, 0);
    Communication_analyze_comm_area(new_grid->DA_3D);
    //getchar();

    /* index start corner (excluding the ghost nodes) */
    new_grid->G_Is = Is;
    new_grid->G_Js = Js;
    new_grid->G_Ks = Ks;

    /* index of end corner(+1) on the current processor (excluding ghost nodes) */
    new_grid->G_Ie = Is + nx;
    new_grid->G_Je = Js + ny;
    new_grid->G_Ke = Ks + nz;


    /* Now, get the start corner index and number of nodes on the current processor including ghost nodes	*/
    ierr = DAGetGhostCorners(new_grid->DA_3D, &Is_g, &Js_g, &Ks_g, &nx_g, &ny_g, &nz_g); PETScErrAct(ierr);

    /* index start corner on the current processor including ghost nodes*/
    new_grid->L_Is = Is_g;
    new_grid->L_Js = Js_g;
    new_grid->L_Ks = Ks_g;

    /* index of end corner(+1) on the current processor (including ghost nodes) */
    new_grid->L_Ie = Is_g + nx_g;
    new_grid->L_Je = Js_g + ny_g;
    new_grid->L_Ke = Ks_g + nz_g;

    /* Create is_solid (global and local) arrays based on the given distributed array */
    /* grid center */
    ierr = DACreateGlobalVector(new_grid->DA_3D_BOX, &new_grid->G_c_status); PETScErrAct(ierr);
    ierr = DACreateLocalVector(new_grid->DA_3D_BOX, &new_grid->L_c_status); PETScErrAct(ierr);

    /* u-grid */
    ierr = VecDuplicate(new_grid->G_c_status, &new_grid->G_u_status); PETScErrAct(ierr);
    ierr = VecDuplicate(new_grid->L_c_status, &new_grid->L_u_status); PETScErrAct(ierr);

    /* v-grid */
    ierr = VecDuplicate(new_grid->G_c_status, &new_grid->G_v_status); PETScErrAct(ierr);
    ierr = VecDuplicate(new_grid->L_c_status, &new_grid->L_v_status); PETScErrAct(ierr);

    /* w-grid */
    ierr = VecDuplicate(new_grid->G_c_status, &new_grid->G_w_status); PETScErrAct(ierr);
    ierr = VecDuplicate(new_grid->L_c_status, &new_grid->L_w_status); PETScErrAct(ierr);

    /* pass the 3d double pointer to the local versions of the q_status arrays.  */
    /* This is done to save a huge amount time by avoiding calling DAVecGetArrays() */
    ierr   = DAVecGetArray(new_grid->DA_3D_BOX, new_grid->L_c_status, (void ***)&new_grid->c_status); PETScErrAct(ierr);
    ierr   = DAVecGetArray(new_grid->DA_3D_BOX, new_grid->L_u_status, (void ***)&new_grid->u_status); PETScErrAct(ierr);
    ierr   = DAVecGetArray(new_grid->DA_3D_BOX, new_grid->L_v_status, (void ***)&new_grid->v_status); PETScErrAct(ierr);
    ierr   = DAVecGetArray(new_grid->DA_3D_BOX, new_grid->L_w_status, (void ***)&new_grid->w_status); PETScErrAct(ierr);

    /* Vector used to update the diagonal part of the linear system matrices for u,v,w and c */
    /* To save memory, only create one and pass the pointer of the others to it */
    ierr = DACreateGlobalVector(new_grid->DA_3D, &new_grid->lsys_diagonal); PETScErrAct(ierr);

    /* create dummy vectors. They can be used anywhere in the code. For now, we use them for transposing the velocities */
    ierr = DACreateGlobalVector(new_grid->DA_3D, &new_grid->G_dummy1); PETScErrAct(ierr);
    ierr = DACreateLocalVector(new_grid->DA_3D, &new_grid->L_dummy1); PETScErrAct(ierr);

    ierr = DACreateGlobalVector(new_grid->DA_3D, &new_grid->G_dummy2); PETScErrAct(ierr);
    ierr = DACreateLocalVector(new_grid->DA_3D, &new_grid->L_dummy2); PETScErrAct(ierr);

    /* For now, create the grid properties on each individual processor since the amount of memory is not
 huge for grid information */
    /**********************************************************/
    /* coordinates of the grid-variables */
    /****************************/
    /* First, cell center coordinates. Add one to the right of each direction */
    xc = Memory_allocate_1D_double(NX);
    yc = Memory_allocate_1D_double(NY);
    zc = Memory_allocate_1D_double(NZ);

    /****************************/
    /* u-grid coordinates. Add one to the right of each direction */
    xu = Memory_allocate_1D_double(NX);
    yu = Memory_allocate_1D_double(NY);
    zu = Memory_allocate_1D_double(NZ);

    /****************************/
    /* v-grid coordinates. Add one to the right of each direction */
    xv = Memory_allocate_1D_double(NX);
    yv = Memory_allocate_1D_double(NY);
    zv = Memory_allocate_1D_double(NZ);

    /****************************/
    /* w-grid coordinates. Add one to the right of each direction */
    xw = Memory_allocate_1D_double(NX);
    yw = Memory_allocate_1D_double(NY);
    zw = Memory_allocate_1D_double(NZ);


    /**********************************************************/
    /* create arrays to hold the values of each cell in each direction. This is done to save some
 computational time */
    dx_c = Memory_allocate_1D_double(NX);
    dy_c = Memory_allocate_1D_double(NY);
    dz_c = Memory_allocate_1D_double(NZ);


    /**********************************************************/
    /* metric coefficients */

    /* (dx/deta)^-1 at cell center (i,j,k) */
    metric_xc = Memory_allocate_1D_double(NX);

    /* (dy/dxi)^-1 at cell center (i,j,k) */
    metric_yc = Memory_allocate_1D_double(NY);

    /* (dz/dphi)^-1 at cell center (i,j,k) */
    metric_zc = Memory_allocate_1D_double(NZ);

    /* (dx/deta)^-1 at cell u_grid (i+1/2,j,k) */
    metric_xu = Memory_allocate_1D_double(NX);

    /* (dy/dxi)^-1 at cell v_grid (i,j+1/2,k) */
    metric_yv = Memory_allocate_1D_double(NY);

    /* (dz/dphi)^-1 at cell w_grid (i,j+1/2,k) */
    metric_zw = Memory_allocate_1D_double(NZ);


    /***********************************/
    /* The physical position of the bottom boundary */
    interface_position = Memory_allocate_2D_double(NX, NZ);

    /* y-index(j) of the (first fluid node) on the bottom geometry */
    interface_y_index  = Memory_allocate_2D_int(NX, NZ);

    Lx = params->Lx;
    Ly = params->Ly;
    Lz = params->Lz;

    /*define step size in xyz and Eta-Xi-Phi coordinate system */
    dEta  = 1.0/(double)NX_cell;
    dXi   = 1.0/(double)NY_cell;
    dPhi  = 1.0/(double)NZ_cell;

    /* Real physical grid dimensions */
    dx = Lx * dEta;
    dy = Ly * dXi;
    dz = Lz * dPhi;

    /* Cell center coordinates */
    new_grid->xc = xc;
    new_grid->yc = yc;
    new_grid->zc = zc;

    /* u-grid coordinates */
    new_grid->xu = xu;
    new_grid->yu = yu;
    new_grid->zu = zu;

    /* v-grid coordinates */
    new_grid->xv = xv;
    new_grid->yv = yv;
    new_grid->zv = zv;

    /* w-grid coordinates */
    new_grid->xw = xw;
    new_grid->yw = yw;
    new_grid->zw = zw;

    /* Grid lengths */
    new_grid->dx_c = dx_c;
    new_grid->dy_c = dy_c;
    new_grid->dz_c = dz_c;


    new_grid->metric_xc = metric_xc;
    new_grid->metric_yc = metric_yc;
    new_grid->metric_zc = metric_zc;
    new_grid->metric_xu = metric_xu;
    new_grid->metric_yv = metric_yv;
    new_grid->metric_zw = metric_zw;

    new_grid->NX    = NX;
    new_grid->NY    = NY;
    new_grid->NZ    = NZ;
    new_grid->NT    = NX*NY*NZ;

    /* Number of grids excluding the last half cell */
    new_grid->NI    = NX-1;
    new_grid->NJ    = NY-1;
    new_grid->NK    = NZ-1;

    /* Length of grid (assuming uniform grid) */
    new_grid->dx   = dx;
    new_grid->dy   = dy;
    new_grid->dz   = dz;

    new_grid->dEta = dEta;
    new_grid->dXi  = dXi;
    new_grid->dPhi = dPhi;

    new_grid->interface_position = interface_position;
    new_grid->interface_y_index  = interface_y_index;

    /* Now, set up the staggered grid coordinates and metric coeficinets */
    /* User can read the grid coordinates from the file: grid.inp. Look at the read_grid_from_file for
 moew info (below) */
    if (params->ImportGridFromFile == YES ) {

        success_flag = Grid_import_grid_from_file(new_grid, params);
        /* Could not find the input files, i.e.
    Grid_x.inp
    Grid_y.inp
    Grid_z.inp

   */
        /* Switch to uniform grid formaulation */
        if (!success_flag) {

            Grid_generate_mesh_coordinates(new_grid, params, 'u');
            Grid_generate_mesh_coordinates(new_grid, params, 'v');
            Grid_generate_mesh_coordinates(new_grid, params, 'w');
            Grid_generate_mesh_coordinates(new_grid, params, 'c');

            params->UniformGridX = YES;
            params->UniformGridY = YES;
            params->UniformGridZ = YES;
        } else {

            PetscPrintf(PCW, "Grid has been imported from file successfully. \n");
        }

    } else {

        Grid_generate_mesh_coordinates(new_grid, params, 'u');
        Grid_generate_mesh_coordinates(new_grid, params, 'v');
        Grid_generate_mesh_coordinates(new_grid, params, 'w');
        Grid_generate_mesh_coordinates(new_grid, params, 'c');
    }

    /* Now, based on the grid, generate the metric coefficients */
    Grid_compute_metric_coefficients(new_grid, 'u');
    Grid_compute_metric_coefficients(new_grid, 'v');
    Grid_compute_metric_coefficients(new_grid, 'w');
    Grid_compute_metric_coefficients(new_grid, 'c');

    //Display_grid(new_grid->metric_xc, new_grid->NY, "metric-xc");
    //Display_grid(new_grid->metric_xu, new_grid->NY, "metric-yu");
    //getchar();
    //Display_grid(new_grid->metric_yc, new_grid->NY, "metric-yc");
    //Display_grid(new_grid->metric_yv, new_grid->NY, "metric-yv");
    //getchar();
    //Display_grid(new_grid->metric_zc, new_grid->NY, "metric-zc");
    //Display_grid(new_grid->metric_zw, new_grid->NY, "metric-zw");
    //getchar();

    /* Create the surface object which represents the location of the surface implicitly using a levelset approach */
    new_grid->surf = Surface_create(new_grid);

    return new_grid;
}
/**********************************************************************************************************/

/* This function releases allocated memory for the variable grid-type */
void Grid_destroy(MAC_grid *grid) {

    int NZ, NY;
    int ierr;

    NZ = grid->NZ;
    NY = grid->NY;

    free(grid->xc);
    free(grid->yc);
    free(grid->zc);

    free(grid->xu);
    free(grid->yu);
    free(grid->zu);

    free(grid->xv);
    free(grid->yv);
    free(grid->zv);

    free(grid->xw);
    free(grid->yw);
    free(grid->zw);

    free(grid->metric_xc);
    free(grid->metric_yc);
    free(grid->metric_zc);

    free(grid->metric_xu);
    free(grid->metric_yv);
    free(grid->metric_zw);

    free(grid->dx_c);
    free(grid->dy_c);
    free(grid->dz_c);

    Memory_free_2D_array(NZ, (void **)grid->interface_position);
    Memory_free_2D_array(NZ, (void **)grid->interface_y_index);

    /* Release the local pointers to the vectors */
    ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->L_c_status, (void ***)&grid->c_status); PETScErrAct(ierr);
    ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->L_u_status, (void ***)&grid->u_status); PETScErrAct(ierr);
    ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->L_v_status, (void ***)&grid->v_status); PETScErrAct(ierr);
    ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->L_w_status, (void ***)&grid->w_status); PETScErrAct(ierr);

    VecDestroy(grid->G_c_status);
    VecDestroy(grid->L_c_status);

    VecDestroy(grid->G_u_status);
    VecDestroy(grid->L_u_status);

    VecDestroy(grid->G_v_status);
    VecDestroy(grid->L_v_status);

    VecDestroy(grid->G_w_status);
    VecDestroy(grid->L_w_status);

    VecDestroy(grid->lsys_diagonal);

    VecDestroy(grid->G_dummy1);
    VecDestroy(grid->G_dummy2);
    VecDestroy(grid->L_dummy1);
    VecDestroy(grid->L_dummy2);

    if ( (grid->surf) != NULL)   {

        Surface_destroy(grid->surf, grid);
    }

    Immersed_destroy(grid->u_immersed, grid);
    Immersed_destroy(grid->v_immersed, grid);
    Immersed_destroy(grid->w_immersed, grid);
    Immersed_destroy(grid->c_immersed, grid);

    DADestroy(grid->DA_3D);
    DADestroy(grid->DA_3D_BOX);
    DADestroy(grid->DA_3D_Lsys);

    free(grid);
}
/**********************************************************************************************************/

/* This function creates the grid coordinates for Uniform mesh!! */
void Grid_generate_mesh_coordinates(MAC_grid *grid, Parameters *params, char which_quantity) {

    int i, j, k;
    double xu_grid, yu_grid, zu_grid;
    double xv_grid, yv_grid, zv_grid;
    double xw_grid, yw_grid, zw_grid;
    double xc_grid, yc_grid, zc_grid;

    /* Half cell is added to the far most place */
    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;

    double dx = grid->dx;
    double dy = grid->dy;
    double dz = grid->dz;

    double xmin = params->xmin;
    double ymin = params->ymin;
    double zmin = params->zmin;

    switch (which_quantity) {

    case 'u':

        xu_grid = xmin;
        for (i=0; i<NX; i++) {

            grid->xu[i]  = xu_grid;
            xu_grid     += dx;
        }

        yu_grid = ymin + 0.5*dy;
        for (j=0; j<NY; j++) {

            grid->yu[j]  = yu_grid;
            yu_grid     += dy;
        }

        zu_grid = zmin + 0.5*dz;
        for (k=0; k<NZ; k++) {

            grid->zu[k]  = zu_grid;
            zu_grid     += dz;
        }

        break;

    case 'v':
        xv_grid = xmin + 0.5*dx;
        for (i=0; i<NX; i++) {

            grid->xv[i]  = xv_grid;
            xv_grid     += dx;
        }

        yv_grid = ymin;
        for (j=0; j<NY; j++) {

            grid->yv[j]  = yv_grid;
            yv_grid     += dy;
        }

        zv_grid = zmin + 0.5*dz;
        for (k=0; k<NZ; k++) {

            grid->zv[k]  = zv_grid;
            zv_grid     += dz;
        }
        break;

    case 'w':
        xw_grid = xmin + 0.5*dx;
        for (i=0; i<NX; i++) {

            grid->xw[i]  = xw_grid;
            xw_grid     += dx;
        }

        yw_grid = ymin + 0.5*dy;
        for (j=0; j<NY; j++) {

            grid->yw[j]  = yw_grid;
            yw_grid     += dy;
        }

        zw_grid = zmin;
        for (k=0; k<NZ; k++) {

            grid->zw[k]  = zw_grid;
            zw_grid     += dz;
        }
        break;

    case 'c':
        xc_grid = xmin + 0.5*dx;
        for (i=0; i<NX; i++) {

            grid->xc[i]  = xc_grid;
            xc_grid     += dx;
        }

        yc_grid = ymin + 0.5*dy;
        for (j=0; j<NY; j++) {

            grid->yc[j]  = yc_grid;
            yc_grid     += dy;
        }

        zc_grid = zmin + 0.5*dz;
        for (k=0; k<NZ; k++) {

            grid->zc[k]  = zc_grid;
            zc_grid     += dz;
        }
        break;

    default:

        PetscPrintf(PCW, "Grid.c/ Error generating Mesh.\nNo '%c' variable is defined\n", which_quantity);
    } /* switch */
}
/**********************************************************************************************************/

/* This function calculates the metric coefficients base on the current coordinates */
void Grid_compute_metric_coefficients(MAC_grid *grid, char which_quantity) {

    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;

    double dEta = grid->dEta;
    double dXi  = grid->dXi;
    double dPhi = grid->dPhi;

    double *xu = grid->xu;
    double *yv = grid->yv;
    double *zw = grid->zw;

    double *xc = grid->xc;
    double *yc = grid->yc;
    double *zc = grid->zc;

    double *dx_c = grid->dx_c;
    double *dy_c = grid->dy_c;
    double *dz_c = grid->dz_c;

    double *metric_xu = grid->metric_xu;
    double *metric_yv = grid->metric_yv;
    double *metric_zw = grid->metric_zw;

    double *metric_xc = grid->metric_xc;
    double *metric_yc = grid->metric_yc;
    double *metric_zc = grid->metric_zc;

    switch (which_quantity) {

    case 'u': {

        int i;
        for (i=1; i<NX-1; i++) {

            metric_xu[i] = dEta / (xc[i] - xc[i-1]);
            dx_c[i]      = xu[i+1] - xu[i];
        }
        /* First and last nodes on the boundaries */
        dx_c[0]         = xu[1] - xu[0];
        dx_c[NX-1]      = xu[NX-1] - xu[NX-2];

        metric_xu[0]    = 0.5*dEta / (xc[0] - xu[0]);
        metric_xu[NX-1] = 0.5*dEta / (xu[NX-1] - xc[NX-2]);
    }
        break;

    case 'v': {

        int j;
        for (j=1; j<NY-1; j++) {

            metric_yv[j] = dXi / (yc[j] - yc[j-1]);
            dy_c[j]      = yv[j+1] - yv[j];
        }
        /* First and last nodes on the boundaries */
        dy_c[0]         = yv[1] - yv[0];
        dy_c[NY-1]      = yv[NY-1] - yv[NY-2];

        metric_yv[0]    = 0.5*dXi / (yc[0] - yv[0]);
        metric_yv[NY-1] = 0.5*dXi / (yv[NY-1] - yc[NY-2]);
    }
        break;

    case 'w': {
        int k;
        for (k=1; k<NZ-1; k++) {

            metric_zw[k] = dPhi / (zc[k] - zc[k-1]);
            dz_c[k]      = zw[k+1] - zw[k];
        }
        /* First and last nodes on the boundaries */
        dz_c[0]         = zw[1] - zw[0];
        dz_c[NZ-1]      = zw[NZ-1] - zw[NZ-2];

        metric_zw[0]    = 0.5*dPhi / (zc[0] - zw[0]);
        metric_zw[NZ-1] = 0.5*dPhi / (zw[NZ-1] - zc[NZ-2]);
    }
        break;


    case 'c': {

        int i, j, k;
        for (i=0; i<NX-1; i++) {

            metric_xc[i] = dEta / (xu[i+1] - xu[i]);
        }
        /* First and last nodes on the boundaries */
        metric_xc[NX-1] = dEta / (xc[NX-1] - xc[NX-2]);

        for (j=0; j<NY-1; j++) {

            metric_yc[j] = dXi / (yv[j+1] - yv[j]);
        }
        /* First and last nodes on the boundaries */
        metric_yc[NY-1] = dXi / (yc[NY-1] - yc[NY-2]);

        for (k=0; k<NZ-1; k++) {

            metric_zc[k] = dPhi / (zc[k+1] - zc[k]);
        }
        /* First and last nodes on the boundaries */
        metric_zc[NZ-1] = dPhi / (zc[NZ-1] - zc[NZ-2]);
    }
        break;

    default:
        PetscPrintf(PCW, "Grid.c/ Error generating Mesh.\nNo '%c' variable is not defined\n", which_quantity);

    } /* switch */
}
/**********************************************************************************************************/

/* This function describes the interface */
void Grid_describe_interface(MAC_grid *grid, Parameters *params) {

    double **interface_position;
    double *xc;
    double *yc;
    double *zc;
    int NX, NZ;
    int i, k;
    double Lx=params->Lx;
    double Ly=params->Ly;
    double Lz=params->Lz;

    int Bottom_Geometry_Shape_Type;
    /*
  0:   Flat surface on the bottom
  1:   Inclinded surface   y=f(x)
  2:   One Errorfucntion ramp at the beginning y=f(x)
  3:   1D Gaussian bump (in y=f(x)
  4:   2D Gaussian bump     y=f(x,z)
  5:   Errorfunction ramp + 1d bump : y=f_1(x) + f_2(x)
  6:   Errorfunction ramp + 2d bump : y =f_1(x) + f_2(x,z)
*/
    double m_ramp;
    double ramp_start_x;
    double ramp_end_x;
    double ramp_start_y;
    double ramp_end_y;

    double max_height;
    double ramp_center;

    double bump_center_x;
    double bump_height;
    double bump_width;

    double bump_center_z;
    double ramp, bump;

    double channel_end_x   = 1.5/2.5*Lx;
    double channel_width_z = 0.1/0.5*Lz;

    /* case 7 */
    double obstacle_center_x;
    double obstacle_center_z;
    double obstacle_height;
    double obstacle_radius_top;
    double obstacle_radius_bottom;


    /* case 8 */
    double *xr, *zr;
    double d, d2, d2min;
    int c;

    int success = NO;
    char filename[100];

    double xstep, step_height;

    /* If the flag is on, then try to import the bottom interface from "Bottom_Topography.inp" file */
    if (params->ImportBottomInterface) {

        sprintf(filename, "Bottom_Topography.inp");
        success = Grid_import_bottom_interface(grid, params, filename) ;
    }
    /* If it was successful, then return. Otherwise, generate the interface based on the default value */
    if (success) {

        PetscPrintf(PCW, "Bottom interface topography has been imported from file successfully.\n");
        return;
    }

    Lx   = params->Lx;
    Ly   = params->Ly;
    Lz   = params->Lz;

    interface_position = grid->interface_position;

    xc = grid->xc;
    yc = grid->yc;
    zc = grid->zc;

    NX = grid->NX;
    NZ = grid->NZ;

    Bottom_Geometry_Shape_Type = 4;
    /* Other options are */
    /*
  0:   Flat surface on the bottom
  1:   Inclinded surface   y=f(x)
  2:   One Errorfucntion ramp at the beginning y=f(x)
  3:   1D Gaussian bump (in y=f(x)
  4:   2D Gaussian bump     y=f(x,z)
  5:   Errorfunction ramp + 1d bump : y=f_1(x) + f_2(x)
  6:   Errorfunction ramp + 2d bump : y =f_1(x) + f_2(x,z)
  7: Al Jaidi's experiment
  11: half cylinder

*/
    switch (Bottom_Geometry_Shape_Type) {

    case 0: /* Flat surface on the bottom */
        break;

    case 1: /* Inclinded surface */

        ramp_start_x = 0.0;
        ramp_end_x   = Lx;
        ramp_start_y = 0.0*Ly;
        ramp_end_y   = 0.0*Ly;
        m_ramp       = (ramp_start_y - ramp_end_y) / (ramp_start_x - ramp_end_x);

        break;

    case 2: /* One Errorfucntion ramp at the beginning */

        max_height  = 0.5*Ly;
        ramp_center = 0.25*Lx;

        break;

    case 3: /* 1D Gaussian bump */

        bump_center_x = 4.0;//0.5*Lx;
        bump_height   = 0.5*Ly;
        bump_width    = 0.5*bump_height; /* not exactly, standard deviation */

        break;

    case 4: /* 2D Gaussian bump */

        //        bump_center_x = 5.50;
        bump_center_x = 1.0;//0.5*Lx;
        bump_center_z = 0.50*Lz;
        bump_height   = 0.5;
        bump_width    = 0.25;

        break;

    case 5: /* Errorfunction ramp + 1d bump */

        max_height  = 0.5*Ly;
        ramp_center = 0.25*Lx;


        bump_center_x = 0.6*Lx;
        bump_height   = 0.5*Ly;
        bump_width    = 0.5*bump_height; /* not exactly, standard deviation */

        break;

    case 6: /* Errorfunction ramp + 2d bump */

        max_height  = 0.5*Ly;
        ramp_center = 0.25*Lx;

        bump_center_x = 0.6*Lx;
        bump_center_z = 0.5*Lz;
        bump_height   = 0.5*Ly;
        bump_width    = 0.41*bump_height; /* not exaclty, standard deviation */

        break;

    case 7: /* Al'Jaidi 's experiment */

        channel_end_x   = 1.5/2.5*Lx;
        channel_width_z = 0.1/0.5*Lz;

        obstacle_center_x      = channel_end_x + 5.0*0.3;
        obstacle_center_z      = 0.0;
        obstacle_height        = 0.048*5.0;
        obstacle_radius_top    = 0.09*5.0;
        obstacle_radius_bottom = 0.15*5.0;

        break;
    case 8: /* Flow in the river */

        /* (xr,zr): coordinates of the centerline of the river in th xz plane */
        zr = Memory_allocate_1D_double(NX);
        xr = grid->xc;

        for (i=0; i<NX; i++) {

            zr[i] = Lz/2.0 + 1.5*erf( (xc[i]-Lx/2.0)/2.0 );
        } /* for */

        break;
    case 9: /* step in x-direction */

        xstep = 0.333*Lx;
        step_height = 0.5*Ly;
        break;

    case 10: /* Run C5 case Kubo's paper */

        break;
    case 11: /* half cylinder */
        break;
    } /* switch */

//My debug
    double offset = grid->dy_c[0]*2.001;
    double h_s;

    for (k=0; k<NZ; k++) {
        for (i=0; i<NX; i++) {

            switch (Bottom_Geometry_Shape_Type) {

            case 0: /* Flat surface on the bottom */

                h_s = 0.0;
                break;

            case 1: /* Inclinded surface */

                if ( (xc[i] > ramp_start_x) && (xc[i] <= ramp_end_x) ){

                    h_s = ramp_start_y + m_ramp*(xc[i] - ramp_start_x);
                }
                else {

                    h_s = 0.0;
                }
                break;

            case 2: /* One Errorfucntion ramp at the beginning */

                h_s = max_height*0.5*(1.0 + 1.0*erf( (ramp_center-xc[i]) / (0.5*max_height)));

                break;

            case 3: /* 1D Gaussian bump */

                h_s = bump_height * exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width));

                break;

            case 4: /* 2D Gaussian bump */

                h_s = bump_height *
                        exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width)) *
                        exp(-(zc[k]-bump_center_z)*(zc[k]-bump_center_z)/(2.0*bump_width*bump_width));
                break;

            case 5: /* Errorfunction ramp + 1d bump */

                ramp = max_height*0.5*(1.0 + 1.0*erf( (ramp_center-xc[i]) / (0.5*max_height)));
                bump = bump_height * exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width));
                h_s = bump + ramp;
                break;

            case 6: /* Errorfunction ramp + 2d bump */

                ramp = max_height*0.5*(1.0 + 1.0*erf( (ramp_center-xc[i]) / (0.5*max_height)));
                bump = bump_height *
                        exp(-(xc[i]-bump_center_x)*(xc[i]-bump_center_x)/(2.0*bump_width*bump_width)) *
                        exp(-(zc[k]-bump_center_z)*(zc[k]-bump_center_z)/(2.0*bump_width*bump_width));

                h_s = bump + ramp;
                break;

            case 7: /* Al'Jaidi 's experiment */

                if (xc[i] <= channel_end_x) {

                    if ( zc[k] <= channel_width_z ) { /* inside channel */

                        h_s = 0.0;
                    } else { /* solid */

                        h_s = Ly;
                    } /* if-else */
                } else {

                    double r = sqrt( (xc[i] - obstacle_center_x)*(xc[i] - obstacle_center_x)
                                     + (zc[k] - obstacle_center_z)*(zc[k] - obstacle_center_z) ) ;

                    if (r <= obstacle_radius_bottom) { //&& (r >= obstacle_radius_top) ) {

                        if (r >= obstacle_radius_top) {

                            h_s = obstacle_height - (r - obstacle_radius_top) / (obstacle_radius_bottom - obstacle_radius_top ) * obstacle_height;
                        } else {

                            h_s = obstacle_height;
                        }
                    } else {

                        h_s = 0.0;
                    } /* else */
                } /* else */
                break;

            case 8: /* Flow in the river */

                d2min = 1.0e10;
                for (c=0; c<NX; c++) {

                    /* Find the distance of the current node from each node on the center line of river */
                    d2 = (xc[i] - xr[c])*(xc[i] - xr[c]) + (zc[k] - zr[c])*(zc[k] - zr[c]);
                    if (d2 < d2min) {

                        d2min = d2;
                    } /* if */
                } /* for */
                d = sqrt(d2min);

                double S = 2.0*GVG_PI;

                h_s = 0.5*Ly*(1.0 - cos(2.0*GVG_PI*d/S)*exp(-0.5*d*d/(0.25*0.25*S*S) ) );
                break;

            case 9: /* step in x-direction */

                if (xc[i] < xstep) {

                    h_s = step_height;
                } else {

                    h_s = 0.0;
                }
                break;

            case 10:
            {
                /* reference length */
                double L_r    = 0.1;
                /* length of initial straigh gate */
                double gate_x = 0.5/L_r;
                /* height of initial gate */
                double height = 0.1/L_r;
                /* end location of ramp */
                double ramp_length = 1.0/L_r;
                /* hump length */
                double h_l = 1.0/L_r;
                double hump_height = 0.036/L_r;

                double x_e1 = gate_x;
                double x_e2 = x_e1 + ramp_length;
                double x_e3 = x_e2 + 0.5*h_l;
                double x_e4 = x_e3 + 0.5*h_l;
                double x_e5 = x_e4 + 0.5*h_l;
                double x_e6 = x_e5 + 0.5*h_l;
                double x_e7 = x_e6 + 0.5*h_l;
                double x_e8 = x_e7 + 0.5*h_l;
                double m_hump = hump_height/(0.5*h_l);

                m_ramp = -height / ramp_length;
                if (xc[i] <= x_e1) {

                    ramp = height;

                } else if (xc[i] <= x_e2) {

                    ramp = height + m_ramp*(xc[i] - x_e1);

                } else if (xc[i] <= x_e3) {

                    ramp = m_hump*(xc[i] - x_e2);

                } else if (xc[i] <= x_e4) {

                    ramp = hump_height - m_hump*(xc[i] - x_e3);

                } else if (xc[i] <= x_e5) {

                    ramp = m_hump*(xc[i] - x_e4);

                } else if (xc[i] <= x_e6) {

                    ramp = hump_height - m_hump*(xc[i] - x_e5);

                } else if (xc[i] <= x_e7) {

                    ramp = m_hump*(xc[i] - x_e6);

                } else if (xc[i] <= x_e8) {

                    ramp = hump_height - m_hump*(xc[i] - x_e7);
                } else {

                    ramp = 0.0;
                }

                h_s = ramp;
            }
                break;

            case 11: {

                double _x_center = params->Lx/3.0;
                double _Radius = 0.5;
                if (fabs(xc[i] - _x_center) <= _Radius) {
                    h_s = sqrt(_Radius*_Radius - (xc[i] - _x_center)*(xc[i] - _x_center));
                } else {
                    h_s = 0.0;
                }

                break;
            }
            default:

                h_s = 0.0;
            } /* switch */

            interface_position[k][i] = h_s + offset;
        } /* for i*/
    }/* for k*/
}
/**********************************************************************************************************/

/* This function returns the state of the grid center node */
/* FLUID, SOLID, IMMERSED, BOUNDARY, OUTSIDE */
int Grid_get_c_status(MAC_grid *grid, int x_index, int y_index, int z_index) {

    static int NI, NJ, NK;

    /* Global number of grids excluding the half cell */
    NI   = grid->NI;
    NJ   = grid->NJ;
    NK   = grid->NK;

    /* Out of regular domain */
    if ( (x_index < 0) || (x_index > NI) ||
         (y_index < 0) || (y_index > NJ) ||
         (z_index < 0) || (z_index > NK) ) {

        return OUTSIDE;
    }
    else {

        return ( (int) grid->c_status[z_index][y_index][x_index] );
    } /*else */
}
/***************************************************************************************************/

/* This function returns the state of the u-grid node */
/* FLUID, SOLID, IMMERSED, BOUNDARY, OUTSIDE */
int Grid_get_u_status(MAC_grid *grid, int x_index, int y_index, int z_index) {

    static int NI, NJ, NK;

    /* Global number of grids excluding the half cell */
    NI   = grid->NI;
    NJ   = grid->NJ;
    NK   = grid->NK;

    /* Out of regular domain */
    if ( (x_index < 0) || (x_index > NI) ||
         (y_index < 0) || (y_index > NJ) ||
         (z_index < 0) || (z_index > NK) ) {

        return OUTSIDE;
    }
    else {

        return ( (int) grid->u_status[z_index][y_index][x_index] );
    } /*else */

}
/***************************************************************************************************/

/* This function returns the state of the v-grid node */
/* FLUID, SOLID, IMMERSED, BOUNDARY, OUTSIDE */
int Grid_get_v_status(MAC_grid *grid, int x_index, int y_index, int z_index) {

    static int NI, NJ, NK;

    /* Global number of grids excluding the half cell */
    NI   = grid->NI;
    NJ   = grid->NJ;
    NK   = grid->NK;

    /* Out of regular domain */
    if ( (x_index < 0) || (x_index > NI) ||
         (y_index < 0) || (y_index > NJ) ||
         (z_index < 0) || (z_index > NK) ) {

        return OUTSIDE;
    }
    else {

        return ( (int) grid->v_status[z_index][y_index][x_index] );
    } /*else */

}
/***************************************************************************************************/

/* This function returns the state of the w-grid node */
/* FLUID, SOLID, IMMERSED, BOUNDARY, OUTSIDE */
int Grid_get_w_status(MAC_grid *grid, int x_index, int y_index, int z_index) {

    static int NI, NJ, NK;

    /* Global number of grids excluding the half cell */
    NI   = grid->NI;
    NJ   = grid->NJ;
    NK   = grid->NK;

    /* Out of regular domain */
    if ( (x_index < 0) || (x_index > NI) ||
         (y_index < 0) || (y_index > NJ) ||
         (z_index < 0) || (z_index > NK) ) {

        return OUTSIDE;
    }
    else {

        return ( (int) grid->w_status[z_index][y_index][x_index] );
    } /*else */

}
/***************************************************************************************************/


/* This function reads the coordinates of the grid in each direction from file. 
Note that you need three files:

1- "Grid_x.inp"
2- "Grid_y.inp"
3- "Grid_z.inp"

In each file (txt files), the data is in the form:
First N nodes are the cell center coordinates, e.g. 
if "NX" is the "number of cells", then the first NX coordinates are corresponding to the xc. 
Then, next NX+1 are the xu coordinates. 
The same rule is applied to the other 2 coordinates. 
If it can not find any of those files, then the grid is assumed to be uniform again */
int Grid_import_grid_from_file(MAC_grid *grid, Parameters *params) {

    int success = YES;
    int i, j, k;
    char filename[100];
    int NX, NY, NZ;
    double xc_, xu_;
    double yc_, yv_;
    double zc_, zw_;
    FILE *file_in;

    NX = params->NX;
    NY = params->NY;
    NZ = params->NZ;

    /* First, read the grid in x-direction */
    sprintf(filename, "Grid_x.inp");
    file_in  = fopen(filename, "r");
    if (file_in == NULL) {

        PetscPrintf(PCW, "Grid.c/ Warning!\nCould not open the grid file \"%s\". Using uniform grid formulation\n", filename);
        success = NO;
    } else {

        for (i=0; i<NX; i++) {

            int n_read = fscanf(file_in, "%lf\n", &xc_);
            if (n_read != 1) {
                PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
            }
            grid->xc[i] = xc_;
        }
        /* Last outside of domain node */
        grid->xc[NX] = xc_ + (grid->xc[NX-1] - grid->xc[NX-2]);

        for (i=NX; i<(2*NX+1); i++) {

            int n_read = fscanf(file_in, "%lf\n", &xu_);
            if (n_read != 1) {
                PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
            }

            grid->xu[i-NX] = xu_;
        }
        fclose(file_in);
    }
    /******************************************/

    /* Read the grid in y-direction */
    sprintf(filename, "Grid_y.inp");
    file_in  = fopen(filename, "r");
    if (file_in == NULL) {

        PetscPrintf(PCW, "Grid.c/ Warning!\nCould not open the grid file \"%s\". Using uniform grid formulation\n", filename);
        success = NO;
    } else {

        for (j=0; j<NY; j++) {

            int n_read = fscanf(file_in, "%lf\n", &yc_);
            if (n_read != 1) {
                PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
            }
            grid->yc[j] = yc_;
        }
        /* Last outside of domain node */
        grid->yc[NY] = yc_ + (grid->yc[NY-1] - grid->yc[NY-2]);

        for (j=NY; j<(2*NY+1); j++) {

            int n_read = fscanf(file_in, "%lf\n", &yv_);
            if (n_read != 1) {
                PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
            }
            grid->yv[j-NY] = yv_;
        }
    }
    fclose(file_in);
    /******************************************/

    /* Read the grid in y-direction */
    sprintf(filename, "Grid_z.inp");
    file_in  = fopen(filename, "r");
    if (file_in == NULL) {

        PetscPrintf(PCW, "Grid.c/ Warning!\nCould not open the grid file \"%s\". Using uniform grid formulation\n", filename);
        success = NO;
    } else {

        for (k=0; k<NZ; k++) {

            int n_read = fscanf(file_in, "%lf\n", &zc_);
            if (n_read != 1) {
                PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
            }
            grid->zc[k] = zc_;
        }
        /* Last outside of domain node */
        grid->zc[NZ] = zc_ + (grid->zc[NZ-1] - grid->zc[NZ-2]);

        for (k=NZ; k<(2*NZ+1); k++) {

            int n_read = fscanf(file_in, "%lf\n", &zw_);
            if (n_read != 1) {
                PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
            }
            grid->zw[k-NZ] = zw_;
        }
    }

    /* Now copy the coordinates for the other grid quantities */
    if (success) {

        for (i=0; i<=NX; i++) {

            grid->xv[i] = grid->xc[i];
            grid->xw[i] = grid->xc[i];
        } /* for i*/

        for (j=0; j<=NY; j++) {

            grid->yu[j] = grid->yc[j];
            grid->yw[j] = grid->yc[j];
        } /* for j*/

        for (k=0; k<=NZ; k++) {

            grid->zu[k] = grid->zc[k];
            grid->zv[k] = grid->zc[k];
        } /* for k*/

    }
    return success;
}
/***************************************************************************************************/

/* This function imports the cell-center bottom interface exact location. This is for the cases where there
is not exact analytical function for describing the bottom interface.
Format would be:  

 x_c z_c y_c
*/
int Grid_import_bottom_interface(MAC_grid *grid, Parameters *params, char *filename) {


    int success = YES;
    FILE *file_in;
    int i, k;
    int NX, NZ;
    double xc_, zc_, y_interface;

    NX = params->NX;
    NZ = params->NZ;

    file_in = fopen(filename, "r");
    if (file_in == NULL) {

        PetscPrintf(PCW, "Warning!\nCould not open the bottom interface grid file \"%s\". Using the default interface\n", filename);
        success = NO;
    } else {

        for (i=0; i<NX; i++) {

            for (k=0; k<NZ; k++) {

                int n_read = fscanf(file_in, "%lf %lf %lf\n", &xc_, &zc_, &y_interface);
                if (n_read != 3) {
                    PetscPrintf(PCW, "Grid.c/ Error reading %s. Probably the file-format is not correct.\n", filename);
                }
                grid->interface_position[k][i] = 30.0*y_interface;
            } /* for k*/
        } /* for i */

        fclose(file_in);
    } /* else */

    return success;
}
/***************************************************************************************************/

/* This function tags the nodes using the given signed-distance-function for each quantity */
/* sdf >  0.0: FLUID */
/* sdf <= 0.0: IMMERSED or SOLID or BUFFER_SOLID */
void Grid_tag_q_nodes_using_surface(MAC_grid *grid, SurfaceType *surf, char which_quantity) {

    /* Number of grid points excluding the half cell added */
    int NI = grid->NI;
    int NJ = grid->NJ;
    int NK = grid->NK;

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;


    int i_start=-1, i_end=-1;
    int j_start=-1, j_end=-1;
    int k_start=-1, k_end=-1;
    void (*Grid_set_q_status)(MAC_grid *, int, int, int, int)=NULL;
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    double ***q_sdf=NULL;

    switch (which_quantity) {

    case 'u':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = max(1, Is);        /* i = 0, BOUNDARY node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = Js;                /* j=0, could be FLUID, IMMERSED or SOLID node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = Ks;      	         /* k=0,  could be FLUID, IMMERSED or SOLID node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* signed distance function: distance from the solid interface */
        q_sdf = surf->u_sdf;

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_u_status;
        Grid_get_q_status = &Grid_get_u_status;
        break;

    case 'v':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = Is;                /* i = 0, BOUNDARY node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = max(1, Js);        /* j=0, BOUNDARY node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = Ks;      	         /* k=0,  could be FLUID, IMMERSED or SOLID node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* signed distance function: distance from the solid interface */
        q_sdf = surf->v_sdf;

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_v_status;
        Grid_get_q_status = &Grid_get_v_status;

        break;

    case 'w':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = Is;        /* i = 0, BOUNDARY node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = Js;                /* j=0, could be FLUID, IMMERSED or SOLID node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = max(1, Ks);        /* k=0,  BOUNDARY node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* signed distance function: distance from the solid interface */
        q_sdf = surf->w_sdf;

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_w_status;
        Grid_get_q_status = &Grid_get_w_status;

        break;

    case 'c':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = Is;                /* i = 0, could be FLUID, IMMERSED or SOLID node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = Js;                /* j=0, could be FLUID, IMMERSED or SOLID node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = Ks;      	         /* k=0,  could be FLUID, IMMERSED or SOLID node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* signed distance function: distance from the solid interface */
        q_sdf = surf->c_sdf;

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_c_status;
        Grid_get_q_status = &Grid_get_c_status;
        break;

    default:
        PetscPrintf(PCW, "Grid.c/ Could not start tagging the nodes. Unknown quantity\n");
    } /* switch */

    double s_p, s_e, s_w, s_n, s_s, s_b, s_f;
    int i, j, k;
    int index ;

    for (k=k_start; k<k_end; k++) {
        for (j=j_start; j<j_end; j++) {
            for (i=i_start; i<i_end; i++) {

                /* current node */
                s_p = q_sdf[k][j][i];

                /* If the node is very close to the boundary, assume it is a fluid node */
                if (s_p >= 0.0) {

                    Grid_set_q_status(grid, i, j, k, FLUID);


                } else { /* Inside the solid region */


                    /* current node */
                    s_p = q_sdf[k][j][i];

                    /* Now, check if at least one of the neighbors are inside fluid. That way, the current node is an IMMERSED node */
                    /* east node */
                    index = min(i+1, NI-1);
                    s_e = q_sdf[k][j][index];

                    /* west node */
                    index = max(i-1, 0);
                    s_w = q_sdf[k][j][index];

                    /* north node */
                    index = min(j+1, NJ-1);
                    s_n = q_sdf[k][index][i];

                    /* south node */
                    index = max(j-1, 0);
                    s_s = q_sdf[k][index][i];

                    /* front node */
                    index = min(k+1, NK-1);
                    s_f = q_sdf[index][j][i];

                    /* back node */
                    index = max(k-1, 0);
                    s_b = sign(q_sdf[index][j][i]);

                    if ( ( (s_p*s_e) < 0.0) ||
                         ( (s_p*s_w) < 0.0) ||
                         ( (s_p*s_n) < 0.0) ||
                         ( (s_p*s_s) < 0.0) ||
                         ( (s_p*s_b) < 0.0) ||
                         ( (s_p*s_f) < 0.0) ) {

                        Grid_set_q_status(grid, i, j, k, IMMERSED);

                    } /* if */ else {

                        Grid_set_q_status(grid, i, j, k, SOLID);
                    } /* if */
                } /* else */

                if (which_quantity == 'c') {

                    //Debugger_check_q_status(i, j, k, grid, params, which_quantity);
                }

            } /* for i */
        } /* for j */
    } /* for k */
}
/***************************************************************************************************/

/* This functions sets the BOUNDARY flag for the box boundaries. The last nodes at (x=Lx, y=Ly, z=Lz) (This is due to the fact that we have addded a half cell to the end of all the grids so they have equal numbers in all directions 
Also, depedning on the grid (u, v or w), it would do the same for x=0, y=0, z=0 */
void Grid_tag_q_box_boundary_nodes(MAC_grid *grid, char which_quantity) {

    /* Total number of grid points in the domain, including the half cell added to the very end */
    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;

    void (*Grid_set_q_status)(MAC_grid *, int, int, int, int)=NULL;
    switch (which_quantity) {
    /* pointer to functions */
    case 'u':
        Grid_set_q_status = &Grid_set_u_status;
        break;

    case 'v':
        Grid_set_q_status = &Grid_set_v_status;
        break;

    case 'w':
        Grid_set_q_status = &Grid_set_w_status;
        break;

    case 'c':
        Grid_set_q_status = &Grid_set_c_status;
        break;

    default:
        PetscPrintf(PCW, "Grid.c/ Could not tag the box boundary nodes. Unknown quantity\n");
    } /* switch */

    int i, j, k;
    /* yz plane */
    for (k=0; k<NZ; k++) {
        for (j=0; j<NY; j++) {

            Grid_set_q_status(grid, NX-1, j, k, BOUNDARY);
            if (which_quantity == 'u') {

                Grid_set_q_status(grid, 0, j, k, BOUNDARY);
            } /* if */
        } /* for j*/
    }/* for k*/

    /* xz plane */
    for (k=0; k<NZ; k++) {
        for (i=0; i<NX; i++) {

            Grid_set_q_status(grid, i, NY-1, k, BOUNDARY);
            if (which_quantity == 'v') {

                Grid_set_q_status(grid, i, 0, k, BOUNDARY);
            } /* if */
        } /* for i */
    } /* for k */

    /* xy plane */
    for (j=0; j<NY; j++) {
        for (i=0; i<NX; i++) {

            Grid_set_q_status(grid, i, j, NZ-1, BOUNDARY);
            if (which_quantity == 'w') {

                Grid_set_q_status(grid, i, j, 0, BOUNDARY);
            } /* if */
        } /* for i */
    } /* for j*/


}
/*****************************************************************************************/

/* This function sets velocity status: FLUID, SOLID, IMMERSED, BOUNDARY */
void Grid_set_u_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status) {

    double ***status;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* if the node is on entire domain (plus the half cell) set the status, else, do nothing */
    if ( (x_index >= Is) && (x_index < Ie) &&
         (y_index >= Js) && (y_index < Je) &&
         (z_index >= Ks) && (z_index < Ke) ) {

        ierr   = DAVecGetArray(grid->DA_3D_BOX, grid->G_u_status, (void ***)&status);
        PETScErrAct(ierr);

        status[z_index][y_index][x_index] = (double) which_status;
        ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->G_u_status, (void ***)&status);
        PETScErrAct(ierr);

    } /* end if */
}
/***************************************************************************************************/

/* This function sets velocity status: FLUID, SOLID, IMMERSED, BOUNDARY */
void Grid_set_v_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status) {

    double ***status;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* if the node is on entire domain (plus the half cell) set the status, else, do nothing */
    if ( (x_index >= Is) && (x_index < Ie) &&
         (y_index >= Js) && (y_index < Je) &&
         (z_index >= Ks) && (z_index < Ke) ) {

        ierr   = DAVecGetArray(grid->DA_3D_BOX, grid->G_v_status, (void ***)&status);
        PETScErrAct(ierr);

        status[z_index][y_index][x_index] = (double) which_status;
        ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->G_v_status, (void ***)&status);
        PETScErrAct(ierr);

    } /* end if */
}
/***************************************************************************************************/

/* This function sets velocity status: FLUID, SOLID, IMMERSED, BOUNDARY */
void Grid_set_w_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status) {

    double ***status;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* if the node is on entire domain (plus the half cell) set the status, else, do nothing */
    if ( (x_index >= Is) && (x_index < Ie) &&
         (y_index >= Js) && (y_index < Je) &&
         (z_index >= Ks) && (z_index < Ke) ) {

        ierr   = DAVecGetArray(grid->DA_3D_BOX, grid->G_w_status, (void ***)&status);
        PETScErrAct(ierr);

        status[z_index][y_index][x_index] = (double) which_status;
        ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->G_w_status, (void ***)&status);
        PETScErrAct(ierr);

    } /* end if */
}
/***************************************************************************************************/

/* This function sets velocity status: FLUID, SOLID, IMMERSED, BOUNDARY */
void Grid_set_c_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status) {

    double ***status;
    int ierr;
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    /* if the node is on entire domain (plus the half cell) set the status, else, do nothing */
    if ( (x_index >= Is) && (x_index < Ie) &&
         (y_index >= Js) && (y_index < Je) &&
         (z_index >= Ks) && (z_index < Ke) ) {

        ierr   = DAVecGetArray(grid->DA_3D_BOX, grid->G_c_status, (void ***)&status);
        PETScErrAct(ierr);

        status[z_index][y_index][x_index] = (double) which_status;
        ierr   = DAVecRestoreArray(grid->DA_3D_BOX, grid->G_c_status, (void ***)&status);
        PETScErrAct(ierr);

    } /* end if */
}
/***************************************************************************************************/

/* This function tags all the nodes in the domain based on the bottom surface */
void Grid_identify_geometry(MAC_grid *grid) {

    /* First, tag the box boundary nodes */
    Grid_tag_q_box_boundary_nodes(grid, 'u') ;
    Grid_tag_q_box_boundary_nodes(grid, 'v') ;
    Grid_tag_q_box_boundary_nodes(grid, 'w') ;
    Grid_tag_q_box_boundary_nodes(grid, 'c') ;

    /* Tag all the nodes based on the location of the inteface for u, v, w and c grid nodes */
    Grid_tag_q_nodes_using_surface(grid, grid->surf, 'u');
    Grid_tag_q_nodes_using_surface(grid, grid->surf, 'v');
    Grid_tag_q_nodes_using_surface(grid, grid->surf, 'w');
    Grid_tag_q_nodes_using_surface(grid, grid->surf, 'c');

    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_u_status, &grid->L_u_status, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_v_status, &grid->L_v_status, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_w_status, &grid->L_w_status, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_c_status, &grid->L_c_status, 'I');

    /* Call this part if you want to solve a BUFFER_SOLID layer inside the solid region
    Note  that the solution inside the solid region do not influence the solution inside the fluid region,
    however, there exist some artificial oscillations in the normal velocity and shear stress which can be significanly reduced
    if we set at least one or two layers of BUFFER_SOLID nodes. */
#ifdef SOLVE_INSIDE_SOLID

    int width = 2;
    Grid_tag_q_buffer_solid_nodes(grid, width, 'u');
    Grid_tag_q_buffer_solid_nodes(grid, width, 'v');
    Grid_tag_q_buffer_solid_nodes(grid, width, 'w');
    //Grid_tag_q_buffer_solid_nodes(grid, width, 'c');

    /* Now again, re-update the ghost node values */
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_u_status, &grid->L_u_status, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_v_status, &grid->L_v_status, 'I');
    Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_w_status, &grid->L_w_status, 'I');
    //Communication_update_ghost_nodes(&grid->DA_3D_BOX, &grid->G_c_status, &grid->L_c_status, 'I');

#endif

    Grid_set_interface_y_index(grid);

    //Display_DA_3D_L_data(grid->L_c_status, grid, params, "c_status", 'c');
    //Display_DA_3D_L_data(grid->L_u_status, grid, params, "u_status", 'u');
    //Display_DA_3D_L_data(grid->L_v_status, grid, params, "v_status", 'v');
    //Display_DA_3D_L_data(grid->L_w_status, grid, params, "w_status", 'w');
    //getchar();

}
/***************************************************************************************************/

/* This function finds the index of the grid for the given x value.
Note that the returned index (grid) is less than the given input position (one grid to the left in the array) */
int Grid_get_x_index(double x, MAC_grid *grid, char which_quantity) {

    double *xq=NULL;
    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
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

    static int NI;
    NI=grid->NI;
    if (x < xq[0]) return NOT_FOUND; /* The given point is outside of the domain, to the left */
    if (x > xq[NI-1]) return GVG_INF; /* The given point is outside of the domain, to the right */

    /* Do binary search to locate the node */
    return MyMath_binary_search(xq, x, 0, NI-1);
}
/***************************************************************************************************/

/* This function finds the index of the grid for the given y value.
Note that the returned index (grid) is less than the given input position */
int Grid_get_y_index(double y, MAC_grid *grid, char which_quantity) {

    double *yq=NULL;
    switch (which_quantity) {

    case 'u':
        yq = grid->yu;
        break;
    case 'v':
        yq = grid->yv;
        break;
    case 'w':
        yq = grid->yw;
        break;
    case 'c':
        yq = grid->yc;
        break;
    } /* switch */

    static int NJ;
    NJ = grid->NJ;

    if (y < yq[0]) return NOT_FOUND; /* The given point is outside of the domain, to the left */
    if (y > yq[NJ-1]) return GVG_INF; /* The given point is outside of the domain, to the right */

    /* Do binary search to locate the node */
    return MyMath_binary_search(yq, y, 0, NJ-1);
}
/***************************************************************************************************/

/* This function finds the index of the grid for the given y value.
Note that the returned index (grid) is less than the given input position */
int Grid_get_z_index(double z, MAC_grid *grid, char which_quantity) {

    double *zq=NULL;
    switch (which_quantity) {

    case 'u':
        zq = grid->zu;
        break;
    case 'v':
        zq = grid->zv;
        break;
    case 'w':
        zq = grid->zw;
        break;
    case 'c':
        zq = grid->zc;
        break;
    } /* switch */

    static int NK;
    NK = grid->NK;

    if (z < zq[0]) return NOT_FOUND; /* The given point is outside of the domain, to the left */
    if (z >= zq[NK-1]) return GVG_INF; /* The given point is outside of the domain, to the right */

    /* Do binary search to locate the node */
    return MyMath_binary_search(zq, z, 0, NK-1);
}
/***************************************************************************************************/

/* This function assumes that we only have a bottom boundry and sets the first y-index correspoding to the IMMERSED node to it */
/* If there is no surface (first interior node is a FLUID node, then it sets j=0. If surface does not lay on the current processor, it sets it to -1 */
void Grid_set_interface_y_index(MAC_grid *grid) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;


    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    for (k=0; k<grid->NZ; k++) {
        for (i=0; i<grid->NX; i++) {

            grid->interface_y_index[k][i] = NOT_FOUND;
        } /* for */
    } /* for */

    for (k=Ks; k<Ke; k++) {
        for (i=Is; i<Ie; i++) {

            /* set the index of first IMMERSED (or FLUID for no solid boundary on the bottom) node */
            for (j=Js; j<Je; j++) {

                if (Grid_get_c_status(grid, i, j, k) == IMMERSED) {

                    grid->interface_y_index[k][i] = j;
                    break;
                } else {
                    if ( (j == 0) && (Grid_get_c_status(grid, i, j, k) == FLUID) ) {

                        grid->interface_y_index[k][i] = j;
                        break;
                    } /* if */
                } /* if */
            } /* for j */
        } /* for i */
    } /* for k */


}

/***************************************************************************************************/

/* This function sets the BUFFER_SOLID nodes for the given width around each IMMERSED node, i.e. 2*width */
/* Set width=1 at least. Do not use a very large value for width, a good number  could be less than 3. */
void Grid_tag_q_buffer_solid_nodes(MAC_grid *grid, int width, char which_quantity) {

    /* Number of grid points excluding the half cell added */
    int NI = grid->NI;
    int NJ = grid->NJ;
    int NK = grid->NK;

    /* Start index of bottom-left-back corner on current processor */


    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;


    int i_start=-1, i_end=-1;
    int j_start=-1, j_end=-1;
    int k_start=-1, k_end=-1;
    void (*Grid_set_q_status)(MAC_grid *, int, int, int, int)=NULL;
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;

    switch (which_quantity) {

    case 'u':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = max(1, Is);        /* i = 0, BOUNDARY node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = Js;                /* j=0, could be FLUID, IMMERSED or SOLID node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = Ks;      	         /* k=0,  could be FLUID, IMMERSED or SOLID node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_u_status;
        Grid_get_q_status = &Grid_get_u_status;
        break;

    case 'v':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = Is;                /* i = 0, BOUNDARY node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = max(1, Js);        /* j=0, BOUNDARY node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = Ks;      	         /* k=0,  could be FLUID, IMMERSED or SOLID node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_v_status;
        Grid_get_q_status = &Grid_get_v_status;

        break;

    case 'w':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = Is;        /* i = 0, BOUNDARY node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = Js;                /* j=0, could be FLUID, IMMERSED or SOLID node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = max(1, Ks);        /* k=0,  BOUNDARY node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_w_status;
        Grid_get_q_status = &Grid_get_w_status;

        break;

    case 'c':

        /* start and end of the for-loops */
        /* The other nodes are the box boudarieds */
        i_start = Is;                /* i = 0, could be FLUID, IMMERSED or SOLID node */
        i_end   = min(Ie, NI);     /* i = NI, BOUDARY node */

        j_start = Js;                /* j=0, could be FLUID, IMMERSED or SOLID node */
        j_end   = min(Je, NJ);     /* j=NJ, BOUNDARY node */

        k_start = Ks;      	         /* k=0,  could be FLUID, IMMERSED or SOLID node */
        k_end   = min(Ke, NK);     /* k=NK, BOUNDARY node */

        /* Pointer to the functions */
        Grid_set_q_status = &Grid_set_c_status;
        Grid_get_q_status = &Grid_get_c_status;
        break;

    default:
        PetscPrintf(PCW, "Grid.c/ Could not start tagging the nodes. Unknown quantity '%c'\n", which_quantity);
    } /* switch */

    int i, j, k;
    int kx, ky, kz;

    /* Go through the domain, and for each IMMERSED node, tag the neighboring SOLID nodes as BUFFER_SOLID node */
    /* The number of neihgboring BUFFER_SOLID is (2*width) */
    for (k=k_start; k<k_end; k++) {
        for (j=j_start; j<j_end; j++) {
            for (i=i_start; i<i_end; i++) {

                if (Grid_get_q_status(grid, i, j, k) == SOLID) {

                    short int done = NO;
                    for (kz=-width; kz<=width && !done; kz++) {
                        for (ky=-width; ky<=width && !done; ky++) {
                            for (kx=-width; kx<=width && !done; kx++) {

                                int i_n = i+kx;
                                int j_n = j+ky;
                                int k_n = k+kz;

                                if (Grid_get_q_status(grid, i_n, j_n, k_n) == IMMERSED) {

                                    Grid_set_q_status(grid, i, j, k, BUFFER_SOLID);
                                    done = YES;
                                } /* if */
                            } /* for kx */
                        } /* for ky */
                    } /* for kz */
                } /* if */
            } /* for i */
        } /* for j */
    } /* for k */
}

/* This function checks wether a point lies within the range of current processor incliding the ghost nodes */
short int Grid_is_local(int i, int j, int k, MAC_grid *grid) {

    /* Check if the current node (with global index i,j,k) lies within the current processor including the ghost nodes */
    if ( (i < grid->L_Is) || (j < grid->L_Js) || (k < grid->L_Ks) ||
            (i >= grid->L_Ie) || (j >= grid->L_Je) || (k >= grid->L_Ke) ) {
        return NO;
    } /*if */
    return YES;
}

