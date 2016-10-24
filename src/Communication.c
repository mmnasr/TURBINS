#include "definitions.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Communication.h"
#include <stdlib.h>
#include <stdio.h>


/* This function reduces (sums) all the global arrays to the W_array. */
/* set status: 
1- status: REDUCE_TO_MASTER : result will be stored only on processor zero 
2- status: REDUCE_TO_ALL    : result will be stores on all processors 
*/
void Communication_reduce_2D_arrays(double **G_send_array, double **W_recv_array, Indices *G_s, Indices *G_e, Indices *W_e, int status, Parameters *params) {

    int i, j;
    int index;
    int W_nt;
    double *recv_buffer=NULL;

    /* Start index [y_s,x_s] */
    /* End  index  [y_e-1, x_e-1] */

    int x_s = G_s->x_index;
    int y_s = G_s->y_index;

    int x_e = G_e->x_index;
    int y_e = G_e->y_index;

    /* Last index (+1) of the 2D world array */
    int X_e = W_e->x_index;
    int Y_e = W_e->y_index;

    /* Total number of grid points in the World array */
    W_nt = X_e * Y_e;


    double *send_buffer = Memory_allocate_1D_double(W_nt);

    /* First pack the data into a 1D array */
    for (j=y_s; j<y_e; j++) {
        for (i=x_s; i<x_e; i++) {

            /* World index of the point in the array */
            index = j*X_e + i;
            send_buffer[index] = G_send_array[j][i];
        } /* for i*/
    } /* for j*/

    if (status == REDUCE_TO_ALL) {

        double *recv_buffer = Memory_allocate_1D_double(W_nt);

        /* This will send to all processors after computing deposit_height */
        (void)MPI_Allreduce( (void *)send_buffer, (void *)recv_buffer, W_nt, MPI_DOUBLE, MPI_SUM, PCW);

    } else if (status == REDUCE_TO_MASTER) {

        if (params->rank == MASTER) {

            recv_buffer = Memory_allocate_1D_double(W_nt);
        } else {

            recv_buffer = NULL;
        } /* else */

        /* Store the final result on processor zero */
        MPI_Reduce( (void *)send_buffer, (void *)recv_buffer, W_nt, MPI_DOUBLE, MPI_SUM, MASTER, PCW);
    } /* if */

    /* Paste back the reduced data to the W_array in the 2D array format */
    if (recv_buffer != NULL ) {

        index = 0;
        /* Now, unpack the 1D array back into the 2D array */
        for (j=0; j<Y_e; j++) {
            for (i=0; i<X_e; i++) {

                W_recv_array[j][i] = recv_buffer[index];
                index++;
            } /* for i*/
        } /* for j*/
    } /* if recv_buffer */
    free(send_buffer);

    if (recv_buffer != NULL) {

        free(recv_buffer);
    }

}
/***************************************************************************************************/

/* This function, communicate between processors and update the ghost nodes on the local vector on current
processor */
void Communication_update_ghost_nodes(DA *DA_array, Vec *global, Vec *local, char status) {

    int ierr;

    /* Replace old values: INSERTVALUES */
    if ( (status == 'I') || (status == 'i') ) {

        ierr = DAGlobalToLocalBegin(*DA_array, *global, INSERT_VALUES, *local);PETScErrAct(ierr);
        ierr = DAGlobalToLocalEnd(*DA_array, *global, INSERT_VALUES, *local);PETScErrAct(ierr);

    } else { /* ADD_VALUES */

        ierr = DAGlobalToLocalBegin(*DA_array, *global, ADD_VALUES, *local);PETScErrAct(ierr);
        ierr = DAGlobalToLocalEnd(*DA_array, *global, ADD_VALUES, *local);PETScErrAct(ierr);
    }
} 
/***************************************************************************************************/

/* Do domain decompostion manually */
/* The user is responsible for the right number of grids on each local processor */
/*******************************************************************/
/* npx: number of processor in x-direction  */
/* npy: number of processor in y-direction  */
/* npz: number of processor in z-direction  */
/* Make sure that: "npx*npy*npz = np" with "np" being the total number of processors */
int Communication_create_DA3D(int npx, int npy, int npz, DAStencilType stencil_type, int ghost_node, Parameters *params, DA *DA_new) {

    /* Array containing the number of nodes on each processor (in x,y, and z-directions) */
    int *nx = Memory_allocate_1D_int(npx);
    int *ny = Memory_allocate_1D_int(npy);;
    int *nz = Memory_allocate_1D_int(npz);;

    if ( ( npx == 0) || (npy == 0) || (npz == 0) ) {

        PetscPrintf(PCW, "Communication.c/ Error performing domain decomposition manually. Cannot set number of processors equal to zero in any direction. \n");
        PetscPrintf(PCW, "Communication.c/ npx:%d npy:%d npz:%d np:%d\n", npx, npy, npz, params->size);

        return NO;
    } /* if */

    /* Total number of processors should be equal to size=npx*npy*npz */
    if (npx*npy*npz != params->size) {

        PetscPrintf(PCW, "Communication.c/ Error performing domain decomposition manually. Total number of processors is not equal to the assigned values in each direction, i.e. size=npx*npy*npz.\n");
        PetscPrintf(PCW, "Communication.c/ npx:%d npy:%d npz:%d np:%d\n", npx, npy, npz, params->size);

        return NO;
    } /* if */

    /* Global number of nodes in x-direction (including the half-cell) */
    int NX = params->NX+1;
    int rem = NX % npx;
    int m;
    for (m=0; m<npx; m++) {

        int temp = NX/npx;
        /* Add one point to the m^th processor if NX/npx is not divisible by zero */
        if (m < rem) {
            temp++;
        }
        nx[m] = temp;
    } /* for m */

    /* Global number of nodes in y-direction (including the half-cell) */
    int NY = params->NY+1;
    rem = NY % npy;
    for (m=0; m<npy; m++) {

        int temp = NY/npy;
        /* Add one point to the m^th processor if NY/npy is not divisible by zero */
        if (m < rem) {
            temp++;
        }
        ny[m] = temp;
    }

    /* Global number of nodes in z-direction (including the half-cell) */
    int NZ = params->NZ+1;
    rem = NZ % npz;
    for (m=0; m<npz; m++) {

        int temp = NZ/npz;
        /* Add one point to the m^th processor if NZ/npz is not divisible by zero */
        if (m < rem) {
            temp++;
        }
        nz[m] = temp;
    }

    int ierr = DACreate3d(PCW, DA_NONPERIODIC, stencil_type, NX, NY, NZ, npx, npy, npz, 1, ghost_node, nx, ny, nz, DA_new);
    PETScErrAct(ierr);

    free(nx);
    free(ny);
    free(nz);

    return YES;
}
/***************************************************************************************************/

/* This function computes the total communication area among processors */
void Communication_analyze_comm_area(DA DA_3D) {

    int NX, NY, NZ;
    int npx, npy, npz;
    int ierr = DAGetInfo(DA_3D, PETSC_NULL, &NX, &NY, &NZ, &npx, &npy, &npz, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

    const int *nx;
    const int *ny;
    const int *nz;

    /* Get the information of the parallel layout and nodes owned on each processor */
    ierr = DAGetOwnershipRanges(DA_3D, &nx, &ny, &nz);

    int area_comm_tot = 0;
    int i, j, k;
    /* Go through all processor and compute the surface area between each processors */
    /* Note that for end walls, we do not add the surface area. This has to be fixed if periodic boundary condition is uded */
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            for (k=0; k<npz; k++) {

                int a_p = 2*(nx[i]*ny[j] + nx[i]*nz[k] + ny[j]*nz[k]);
                if (i == 0) {
                    a_p -= ny[j]*nz[k];
                }
                if (i == npx-1) {
                    a_p -= ny[j]*nz[k];
                }
                if (j == 0) {
                    a_p -= nx[i]*nz[k];
                }
                if (j == npy-1) {
                    a_p -= nx[i]*nz[k];
                }
                if (k == 0) {
                    a_p -= nx[i]*ny[j];
                }
                if (k == npz-1) {
                    a_p -= nx[i]*ny[j];
                }

                area_comm_tot += a_p;
            } /* for k*/
        } /* for j*/
    } /* for k*/

    /* percentage of the total surface area to the total volume */
    double s = (double) area_comm_tot / (double) (NZ*NY*NX) * 100.0;
    PetscPrintf(PCW, "Communication.c/ Total commnucation area (excluding ghost nodes): %d S[comm_area/total_volume](\%%):%f", area_comm_tot, s);
}
