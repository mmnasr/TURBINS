#include "definitions.h"
#include "DataTypes.h"
#include "Communication.h"
#include "Grid.h"
#include "Conc.h"
#include "Velocity.h"
#include "Memory.h"
#include "Surface.h"
#include "Immersed.h"
#include "Writer.h"
#include "Extract.h"
#include "Output.h"
#include "Pressure.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

Output *Output_create(MAC_grid *grid, Parameters *params) {

    Output *new_output;
    int ierr;
    int NX, NY, NZ;
    int Ghost_Nodes;
    int Is, Js, Ks;
    int Is_g, Js_g, Ks_g;
    int nx, ny, nz;
    int nx_g, ny_g, nz_g;
    int nt;

    new_output = (Output *)malloc(sizeof(Output));
    Memory_check_allocation(new_output);

    /* SINGLE_FILE_OUTPUT or PARTITIONED_FILES_OUTPUT */
    new_output->type = params->output_type;

    NX = grid->NX ;
    NY = grid->NY ;
    NZ = grid->NZ ;

    /* Only one level of ghost nodes for output purposes is enough */
    Ghost_Nodes = 1;
    new_output->Ghost_Nodes = Ghost_Nodes;
#ifdef AUTOMATIC_DOMAIN_PARTITIONING /* To change this, goto 'definitions.h' */
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_BOX, NX, NY, NZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, Ghost_Nodes, PETSC_NULL, PETSC_NULL, PETSC_NULL, &new_output->DA_3D);
    PETScErrAct(ierr);
#else /* MANUAL_DOMAIN_PARTITIONING */

   /* Number of processors in each direction */
    /* Each processors owns the same number of grid nodes of the subdomain */
    int npx = 4;//params->size/15;
    int npy = 3;
    int npz = 3;

        int nnx[]={46,65,30,30};
        int nny[]={48,48,49};
        int nnz[]={48,48,49};

    Ghost_Nodes = 1;
    ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_BOX, NX, NY, NZ, npx, npy, npz, 1, Ghost_Nodes, nnx, nny, nnz, &new_output->DA_3D);

	if (0) {

	    int success = Communication_create_DA3D(npx, npy, npz, DA_STENCIL_BOX, Ghost_Nodes, params, &new_output->DA_3D);

    if (!success)  {

        PetscPrintf(PCW, "Output.c/ Error!!! Could not perform manual domain decomposition using (npx,npy,npz)=(%d,%d,%d) processors.\n Trying assigning the right number of processors in each direction.\n Or you could allow PETSc do that automatically.\n", npx, npy, npz);
    } /* if */

	}

#endif

    /* generate global and local 3D data arrays */
    /* Local arrays include the ghost nodes from neighboring processors */
    ierr = DACreateGlobalVector(new_output->DA_3D, &new_output->G_data); PETScErrAct(ierr);
    ierr = DACreateLocalVector(new_output->DA_3D, &new_output->L_data); PETScErrAct(ierr);

    /* get corners of the local part of the array on current processor */
    ierr = DAGetCorners(new_output->DA_3D, &Is, &Js, &Ks, &nx, &ny, &nz); PETScErrAct(ierr);
    /* get corners of the local part of the array including ghost nodes on current processor */
    ierr = DAGetGhostCorners(new_output->DA_3D, &Is_g, &Js_g, &Ks_g, &nx_g, &ny_g, &nz_g);

    /* start and end-1 corners of the global data on current processor */
    new_output->G_Is = Is;
    new_output->G_Js = Js;
    new_output->G_Ks = Ks;

    new_output->G_Ie = Is + nx;;
    new_output->G_Je = Js + ny;
    new_output->G_Ke = Ks + nz;

    /* start and end-1 corners of the global data including ghost nodes on current processor */
    new_output->L_Is = Is_g;
    new_output->L_Js = Js_g;
    new_output->L_Ks = Ks_g;

    new_output->L_Ie = Is_g + nx_g;
    new_output->L_Je = Js_g + ny_g;
    new_output->L_Ke = Ks_g + nz_g;

    /* total number of nodes including ghost nodes on current processor */
    nt = nx_g * ny_g * nz_g;

    /* Set the output precision, GVG_FLOAT or GVG_DOUBLE */
    Output_set_precision(new_output, GVG_FLOAT);

    /* allocate 1d float or double array */
    if (new_output->precision == GVG_FLOAT) {

        new_output->L_data_1D_float  = Memory_allocate_1D_float(nt);
        new_output->L_data_1D = (void *)new_output->L_data_1D_float;

        new_output->xu_float = Memory_allocate_1D_float(nx_g);
        new_output->yu_float = Memory_allocate_1D_float(ny_g);
        new_output->zu_float = Memory_allocate_1D_float(nz_g);

        new_output->xv_float = Memory_allocate_1D_float(nx_g);
        new_output->yv_float = Memory_allocate_1D_float(ny_g);
        new_output->zv_float = Memory_allocate_1D_float(nz_g);

        new_output->xw_float = Memory_allocate_1D_float(nx_g);
        new_output->yw_float = Memory_allocate_1D_float(ny_g);
        new_output->zw_float = Memory_allocate_1D_float(nz_g);

        new_output->xc_float = Memory_allocate_1D_float(nx_g);
        new_output->yc_float = Memory_allocate_1D_float(ny_g);
        new_output->zc_float = Memory_allocate_1D_float(nz_g);

        new_output->xu = (void *)new_output->xu_float;
        new_output->yu = (void *)new_output->yu_float;
        new_output->zu = (void *)new_output->zu_float;

        new_output->xv = (void *)new_output->xv_float;
        new_output->yv = (void *)new_output->yv_float;
        new_output->zv = (void *)new_output->zv_float;

        new_output->xw = (void *)new_output->xw_float;
        new_output->yw = (void *)new_output->yw_float;
        new_output->zw = (void *)new_output->zw_float;

        new_output->xc = (void *)new_output->xc_float;
        new_output->yc = (void *)new_output->yc_float;
        new_output->zc = (void *)new_output->zc_float;

    } else {

        new_output->L_data_1D_float = NULL;
    } /* else float */
    if (new_output->precision == GVG_DOUBLE) {

        new_output->L_data_1D_double  = Memory_allocate_1D_double(nt);
        new_output->L_data_1D = (void *)new_output->L_data_1D_double;

        new_output->xu_double = Memory_allocate_1D_double(nx_g);
        new_output->yu_double = Memory_allocate_1D_double(ny_g);
        new_output->zu_double = Memory_allocate_1D_double(nz_g);

        new_output->xv_double = Memory_allocate_1D_double(nx_g);
        new_output->yv_double = Memory_allocate_1D_double(ny_g);
        new_output->zv_double = Memory_allocate_1D_double(nz_g);

        new_output->xw_double = Memory_allocate_1D_double(nx_g);
        new_output->yw_double = Memory_allocate_1D_double(ny_g);
        new_output->zw_double = Memory_allocate_1D_double(nz_g);

        new_output->xc_double = Memory_allocate_1D_double(nx_g);
        new_output->yc_double = Memory_allocate_1D_double(ny_g);
        new_output->zc_double = Memory_allocate_1D_double(nz_g);

        new_output->xu = (void *)new_output->xu_double;
        new_output->yu = (void *)new_output->yu_double;
        new_output->zu = (void *)new_output->zu_double;

        new_output->xv = (void *)new_output->xv_double;
        new_output->yv = (void *)new_output->yv_double;
        new_output->zv = (void *)new_output->zv_double;

        new_output->xw = (void *)new_output->xw_double;
        new_output->yw = (void *)new_output->yw_double;
        new_output->zw = (void *)new_output->zw_double;

        new_output->xc = (void *)new_output->xc_double;
        new_output->yc = (void *)new_output->yc_double;
        new_output->zc = (void *)new_output->zc_double;

    } else {

        new_output->L_data_1D_double = NULL;
    } /* else double */

    /* Copies the local grid coordinates and applies the conversion (float or double) */
    Output_copy_local_grid_data(new_output, grid, 'u');
    Output_copy_local_grid_data(new_output, grid, 'v');
    Output_copy_local_grid_data(new_output, grid, 'w');
    Output_copy_local_grid_data(new_output, grid, 'c');

    /* create writer for current data layout */
    new_output->writer = Writer_create(nx_g, ny_g, nz_g, new_output->precision);

    /* set local grid coordinates pointer */
    Output_set_grid_data(new_output, 'c') ;

    /* Sets the base and data output directory */
    sprintf(new_output->base_dir, "./output/");
    sprintf(new_output->data_dir, "%sdata/", new_output->base_dir);
    sprintf(new_output->bin_dir, "%sbinary/", new_output->base_dir);
    sprintf(new_output->hdf5_dir, "%shdf5/", new_output->base_dir);
    sprintf(new_output->immersed_dir, "%simmersed/", new_output->base_dir);

    /* create the base_directory folder to store the output files */
    if (params->rank == MASTER) {

        if ( mkdir(new_output->base_dir, 0777) != 0) {
            if (errno == EEXIST) {

                PetscPrintf(PCW, "Output.c/ Warning!!! Output directory %s exists. Data maybe overwritten.\n", new_output->base_dir);
            } else {
                PetscPrintf(PCW, "Output.c/ Error creating the folder: %s", new_output->base_dir);
            }
        } /* if */
        /* Used for output data with .dat extensions */
        if ( mkdir(new_output->data_dir, 0777) != 0) {

            if (errno == EEXIST) {

                PetscPrintf(PCW, "Output.c/ Warning!!! Output directory %s exists. Data maybe overwritten.\n", new_output->data_dir);
            } else {
                PetscPrintf(PCW, "Output.c/ Error creating the folder: %s", new_output->data_dir);
            }
        } /* if */
        if ( mkdir(new_output->bin_dir, 0777) != 0) {

            if (errno == EEXIST) {

                PetscPrintf(PCW, "Output.c/ Warning!!! Output directory %s exists. Data maybe overwritten.\n", new_output->bin_dir);
            } else {
                PetscPrintf(PCW, "Output.c/ Error creating the folder: %s", new_output->bin_dir);
            }
        } /* if */

        if ( mkdir(new_output->hdf5_dir, 0777) != 0) {

            if (errno == EEXIST) {

                PetscPrintf(PCW, "Output.c/ Warning!!! Output directory %s exists. Data maybe overwritten.\n", new_output->hdf5_dir);
            } else {
                PetscPrintf(PCW, "Output.c/ Error creating the folder: %s", new_output->hdf5_dir);
            }
        } /* if */

        if ( mkdir(new_output->immersed_dir, 0777) != 0) {

            if (errno == EEXIST) {

                PetscPrintf(PCW, "Output.c/ Warning!!! Output directory %s exists. Data maybe overwritten.\n", new_output->immersed_dir);
            } else {
                PetscPrintf(PCW, "Output.c/ Error creating the folder: %s", new_output->immersed_dir);
            }
        } /* if */


    } /* if rank */

    /* Since we are performing a hardware operation wait for the operation to be finished. I had to add this, otherwise it would crash */
    MPI_Barrier(PCW);

    new_output->update_last_plane = YES;

    /* set the output type : BINARY, HDF5, or ASCII */
    new_output->type = params->output_type;

    switch (new_output->type) {

    case BINARY:
        sprintf(new_output->data_3d_file_ext, "bin");
        sprintf(new_output->data_3d_dir, "%s", new_output->bin_dir);
        break;

    case ASCII:
        sprintf(new_output->data_3d_file_ext, "dat");
        sprintf(new_output->data_3d_dir, "%s", new_output->data_dir);
        break;

    case HDF5:
        sprintf(new_output->data_3d_file_ext, "h5");
        sprintf(new_output->data_3d_dir, "%s", new_output->hdf5_dir);
        break;

    default:
        PetscPrintf(PCW, "Output.c/ Warning. Unknown output type %d. It should be 1 (BINARY), 2 (ASCII), or 3 (HDF5). Setting it to Binary(1).\n", new_output->type);
        new_output->type = BINARY;
        sprintf(new_output->data_3d_dir, "%s", new_output->bin_dir);
    } /* switch */

    if (params->rank == MASTER) {

        /* Write the grid information to files in binary format */
        Output_grid_binaries(new_output, grid, 'u');
        Output_grid_binaries(new_output, grid, 'v');
        Output_grid_binaries(new_output, grid, 'w');
        Output_grid_binaries(new_output, grid, 'c');
        Output_grid_binaries(new_output, grid, 'p');

    } /* grid */

    Output_grid_q_hdf5(new_output, grid, 'u') ;
    Output_grid_q_hdf5(new_output, grid, 'v') ;
    Output_grid_q_hdf5(new_output, grid, 'w') ;
    Output_grid_q_hdf5(new_output, grid, 'c') ;
    Output_grid_q_hdf5(new_output, grid, 'p') ;


    MPI_Barrier(PCW);
    return (new_output);
}
/*************************************************************/

void Output_destroy(Output *output) {

    int ierr;

    Writer_destroy(output->writer);
    if (output->precision == GVG_FLOAT) {

        free(output->xu_float);
        free(output->yu_float);
        free(output->zu_float);

        free(output->xv_float);
        free(output->yv_float);
        free(output->zv_float);

        free(output->xw_float);
        free(output->yw_float);
        free(output->zw_float);

        free(output->xc_float);
        free(output->yc_float);
        free(output->zc_float);

        free(output->L_data_1D_float);

    } else if (output->precision == GVG_DOUBLE) {

        free(output->xu_double);
        free(output->yu_double);
        free(output->zu_double);

        free(output->xv_double);
        free(output->yv_double);
        free(output->zv_double);

        free(output->xw_double);
        free(output->yw_double);
        free(output->zw_double);

        free(output->xc_double);
        free(output->yc_double);
        free(output->zc_double);

        free(output->L_data_1D_double);
    } /* else */
    ierr = VecDestroy(output->G_data); PETScErrAct(ierr);
    ierr = VecDestroy(output->L_data); PETScErrAct(ierr);

    ierr = DADestroy(output->DA_3D); PETScErrAct(ierr);

    free(output);
}
/*************************************************************/


/* This function converts a 3D double array to a 1D float or double array */
void Output_convert_3D_to_1D(Output *output, double ***data_3D, void *data_1D) {

    /* Start index of bottom-left-back corner on the current processor including ghost nodes*/
    int Is_g = output->L_Is;
    int Js_g = output->L_Js;
    int Ks_g = output->L_Ks;

    /* End index of top-right-front corner on the current processor including ghost nodes */
    int Ie_g = output->L_Ie;
    int Je_g = output->L_Je;
    int Ke_g = output->L_Ke;

    int index = 0;
    int i, j, k;

    /* Depending on the precision, cast the data_1D array */
    if (output->precision == GVG_FLOAT) {

/* cast the void* output 1D array to float* */
        float *data = (float *)data_1D;

        for (k=Ks_g; k<Ke_g; k++) {
            for (j=Js_g; j<Je_g; j++) {
                for (i=Is_g; i<Ie_g; i++) {

                    /* Copy data from the 3D array into 1D-float array */
                    data[index] = (float)data_3D[k][j][i];
                    index++;
                } /* for i*/
            } /* for j*/
        } /* for k*/

    } else if (output->precision == GVG_DOUBLE){

        double *data = (double *)data_1D;
        for (k=Ks_g; k<Ke_g; k++) {
            for (j=Js_g; j<Je_g; j++) {
                for (i=Is_g; i<Ie_g; i++) {

                    /* Copy data from the 3D array into 1D-double array */
                    data[index] = data_3D[k][j][i];
                    index++;
                } /* for i*/
            } /* for j*/
        } /* for k*/
    } /* else */
}
/*************************************************************/

/* This function copies the local grid data for any given quantity (float or double conversion) */
void Output_copy_local_grid_data(Output *output, MAC_grid *grid, char which_quantity) {

    double *xq = NULL;
    double *yq = NULL;
    double *zq = NULL;

    void *x_out = NULL;
    void *y_out = NULL;
    void *z_out = NULL;

    switch (which_quantity) {
    case 'u':

        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        x_out = output->xu;
        y_out = output->yu;
        z_out = output->zu;
        break;
    case 'v':

        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        x_out = output->xv;
        y_out = output->yv;
        z_out = output->zv;
        break;
    case 'w':

        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        x_out = output->xw;
        y_out = output->yw;
        z_out = output->zw;
        break;
    case 'c':

        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        x_out = output->xc;
        y_out = output->yc;
        z_out = output->zc;
        break;
    } /* switch */

    /* Start index of bottom-left-back corner on the current processor including ghost nodes*/
    int Is_g = output->L_Is;
    int Js_g = output->L_Js;
    int Ks_g = output->L_Ks;

    /* End index of top-right-front corner on the current processor including ghost nodes */
    int Ie_g = output->L_Ie;
    int Je_g = output->L_Je;
    int Ke_g = output->L_Ke;

    int i, j, k;
    int index = 0;
    double temp=0.0;
    if (output->precision == GVG_FLOAT) {

        float *x_grid = (float *)x_out;
        index = 0;
        for (i=Is_g; i<Ie_g; i++) {

            /* Copy data from the 1D array into 1D-float array */
            temp = xq[i];

            x_grid[index] = (float)temp;
            index++;
        }

        float *y_grid = (float *)y_out;
        index = 0;
        for (j=Js_g; j<Je_g; j++) {

            temp = yq[j];
            /* Copy data from the 1D array into 1D-float array */
            y_grid[index] = (float)temp;
            index++;
        }

        float *z_grid = (float *)z_out;
        index = 0;
        for (k=Ks_g; k<Ke_g; k++) {

            temp = zq[k];
            /* Copy data from the 1D array into 1D-float array */
            z_grid[index] = (float)temp;
            index++;
        }

    } else if (output->precision == GVG_DOUBLE) {

        double *x_grid = (double *)x_out;
        index = 0;
        for (i=Is_g; i<Ie_g; i++) {

            /* Copy data from the 1D array into 1D-float array */
            temp = xq[i];
            x_grid[index] = temp;
            index++;
        }

        double *y_grid = (double *)y_out;
        index = 0;
        for (j=Js_g; j<Je_g; j++) {

            temp = yq[j];
            y_grid[index] = temp;
            index++;
        }

        double *z_grid = (double *)z_out;
        index = 0;
        for (k=Ks_g; k<Ke_g; k++) {

            temp = zq[k];
            z_grid[index] = temp;
            index++;
        }
    } /* else */
}
/*************************************************************/

/* This function sets the pointers to the given (local) grid coordinates */
void Output_set_grid_data(Output *output, char which_quantity) {

    switch (which_quantity) {
    case 'u':

        output->writer->X1_grid = output->xu;
        output->writer->X2_grid = output->yu;
        output->writer->X3_grid = output->zu;
        break;
    case 'v':

        output->writer->X1_grid = output->xv;
        output->writer->X2_grid = output->yv;
        output->writer->X3_grid = output->zv;
        break;
    case 'w':

        output->writer->X1_grid = output->xw;
        output->writer->X2_grid = output->yw;
        output->writer->X3_grid = output->zw;
        break;
    case 'c':

        output->writer->X1_grid = output->xc;
        output->writer->X2_grid = output->yc;
        output->writer->X3_grid = output->zc;
        break;

    } /* switch */

}
/*************************************************************/

void Output_set_filename(Output *output, char *filename) {

    Writer_set_filename(output->writer, filename);
}
/*************************************************************/

void Output_set_output_data(Output *output, void *data) {

    Writer_set_output_data(output->writer, data);
}
/*************************************************************/

void Output_reset_write_grid_flag(Output *output) {

    Writer_reset_write_grid_flag(output->writer) ;
}
/*************************************************************/

void Output_set_write_grid_flag(Output *output) {

    Writer_set_write_grid_flag(output->writer) ;
}
/*************************************************************/

void Output_set_append_data_flag(Output *output) {

    Writer_set_append_data_flag(output->writer) ;
}
/*************************************************************/

void Output_reset_append_data_flag(Output *output) {

    Writer_reset_append_data_flag(output->writer) ;
}
/*************************************************************/

void Output_write_data(Output *output) {

    /* First, open the file */
    Writer_open_file(output->writer);

    /* write grid information */
    Writer_write_grid(output->writer);

    /* write data to file */
    Writer_write_data(output->writer);

    /* close file */
    Writer_close_file(output->writer);
}
/*************************************************************/

/* This function checks if the file has been opened correctly */
void Output_check_file_open(FILE *file, char *filename) {

    if (file == NULL) {

        PetscPrintf(PETSC_COMM_SELF, "Output.c/ Error opening the file %s\n", filename);
    } /* if */
}
/*************************************************************/

/* This function checks if the data has been written correctly to the binary file */
void Output_check_file_write(int n_items_desired, int n_items_actual) {

    if (n_items_desired != n_items_actual) {

        printf("Output.c/ ******** Warning ********...\n There was an error while writing the data to file\n");
    }
}
/*************************************************************/

/* This function writes the grid info on each processor to file. This is done for post processoring purposes */
void Output_write_grid_partition_info(Output *output, MAC_grid *grid, Parameters *params) {

    FILE *fileout;
    int numprocs;
    int nx_g, ny_g, nz_g;
    int Is_g, Js_g, Ks_g;
    int Ie_g, Je_g, Ke_g;
    int rank;
    char filename[MAX_PATH_LEN];

    numprocs = params->size;
    for(rank=0; rank < numprocs; rank++) {

        if (rank == params->rank) {

            /* append data to the already existing file (Generated by processor zero) */
            sprintf(filename, "%sGridPartition.txt", output->base_dir);
            if (rank > 0) {

                fileout = fopen(filename, "a");

            } else { /* processor zero. Generate the new file */

                fileout = fopen(filename, "w");
                fprintf(fileout, "Number of Processors:\n%d\n", numprocs);
                /* Half cell added added to the end of each direction */
                fprintf(fileout, "Global Number of Grids (x,y,z):\n%d %d %d\n", grid->NX, grid->NY, grid->NZ);
                fprintf(fileout, "Output ghost node levels:\n%d\n", output->Ghost_Nodes);

            }
            Output_check_file_open(fileout, filename);
            /* start index of data on current processor */

            /* Start index of bottom-left-back corner on the current processor including ghost nodes*/
            Is_g = output->L_Is;
            Js_g = output->L_Js;
            Ks_g = output->L_Ks;

            /* End index of top-right-front corner on the current processor including ghost nodes */
            Ie_g = output->L_Ie;
            Je_g = output->L_Je;
            Ke_g = output->L_Ke;

            /* number of grids on current processor */
            nx_g = Ie_g - Is_g;
            ny_g = Je_g - Js_g;
            nz_g = Ke_g - Ks_g;

            fprintf(fileout, "**************************************************\n");
            fprintf(fileout, "Processor ID:\n%d\n(x,y,z)_start index:\n%d %d %d\n", rank, Is_g, Js_g,Ks_g);
            //fprintf(fileout, "Number of Grids (x,y,z):\n%d %d %d\n", output->G_Ie-output->G_Is, output->G_Je-output->G_Js, output->G_Ke-output->G_Ks);
            fprintf(fileout, "Local Number of Grids (x,y,z):\n%d %d %d\n", nx_g, ny_g, nz_g);

            fclose(fileout);
        } /* if */

        /* to make the code run in serial for this part */
        MPI_Barrier(PCW);
    } /* for rank */
}
/*************************************************************/

/* This function writes suspended (total interior concentration mass) to file */
/* It only writes on processor zero */
void Output_write_ascii_suspended_mass(Output *output, Concentration **c, Parameters *params, double time) {

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        int iconc, NConc;
        char filename[MAX_PATH_LEN];
        NConc = params->NConc;

        for (iconc=0; iconc<NConc; iconc++) {


            sprintf(filename, "%sS_susp_mass_%d.dat", output->data_dir, iconc);
            /* open the file and append to the end of it */
            FILE *fileout = fopen(filename, "a");
            fprintf(fileout, "%f %f\n", time, c[iconc]->W_total_susp_mass);
            fclose(fileout);

        } /* for */
    }/* if */
}
/*************************************************************/

/* This function writes average height of the gravity front to file at the given time. 
Note that only processor zero writes to file */
void Output_write_ascii_average_height(Output *output, Concentration **c, MAC_grid *grid, Parameters *params, int timestep) {

    FILE *fileout;
    int i, NX;
    int iconc, NConc;
    int filenumber;
    char filename[50];

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        filenumber = 10000 + timestep;
        NConc      = params->NConc;
        NX         = grid->NX;
        for (iconc=0; iconc<NConc; iconc++) {

            sprintf(filename, "%sS_ave_height_c%d_%d.dat", output->data_dir, iconc, filenumber);
            /* open the file and append to the end of it */
            fileout = fopen(filename, "w");
            /* write the whole 1D array to file. h=h_ave(x) */
            for (i=0; i<NX-1; i++) {

                fprintf(fileout, "%f %f\n", grid->xc[i], c[iconc]->W_ave_height_x[i]);
            }
            fclose(fileout);

        } /* for */
    }/* if */
}
/*************************************************************/

/* This function writes the front location of the current to the file */
void Output_write_ascii_front_location(Output *output, double front_location, Parameters *params, double time) {

    FILE *fileout;
    char filename[MAX_PATH_LEN];

    /* Only write on processor zero */
    if (params->rank == MASTER) {


        sprintf(filename, "%sS_front_location.dat", output->data_dir);
        /* open the file and append to the end of it */
        fileout = fopen(filename, "a");
        /* write the front location at the given time */
        fprintf(fileout, "%f %2.14f\n", time, front_location);

        fclose(fileout);
    }/* if */
}
/*************************************************************/

/* This function writes 2D deposit height to ASCII file for each concentratio field */
void Output_write_ascii_deposit_height(Output *output, Concentration **c, MAC_grid *grid, Parameters *params, int timestep) {

    FILE *fileout;
    int i, k, NI, NK;
    int iconc, NConc;
    int filenumber;
    char filename[MAX_PATH_LEN];

    /* Only processor zero performs writing */
    if (params->rank == MASTER) {

        NConc      = params->NConc;
        filenumber = 10000 + timestep;

        NK = c[0]->NK;
        NI = c[0]->NI;

        for (iconc=0; iconc<NConc; iconc++) {

            /* Only write to file if we have particle, i.e. u_s > 0.0 field */
            if (c[iconc]->Type == PARTICLE) {

                sprintf(filename, "%sS_deposit_height_c%d_%d.dat", output->data_dir, iconc, filenumber);
                /* open the new file */
                fileout = fopen(filename, "w");

                for (k=0; k<NK; k++) {
                    for (i=0; i<NI; i++) {

                        /* write the deposit height at each (z,x) node */
                        fprintf(fileout, "%f %f %2.12f\n", grid->zc[k], grid->xc[i], c[iconc]->W_deposit_height[k][i]);
                    } /* for i */
                } /* for k*/
                fclose(fileout);
            } /* if */

        } /* for */
    }/* if */
}
/*************************************************************/

/* This function writes 
1- Total kinetic energy, 
2- Viscous dissipation 
3- Active, Active_grad and Passive potential energies (of the particles)
to the file. 
*/
void Output_write_ascii_energies(Output *output, Concentration **c, Velocity *u, Parameters *params, double time) {

    FILE *fileout;
    char filename[MAX_PATH_LEN];

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        sprintf(filename, "%sS_kinetic_energy.dat", output->data_dir);
        /* Open the file for kinetic energy */
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */

        /* Kinetic energy of the fluid */
        fprintf(fileout, "%f %f\n", time, u->W_kinetic_energy);
        fclose(fileout);
        /********************************************************************************************/

        sprintf(filename, "%sS_viscous_dissipation_rate.dat", output->data_dir);
        /* Open the file for viscous dissipation */
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */

        /* Kinetic energy of the fluid */
        fprintf(fileout, "%f %2.14f\n", time, u->W_dissipation_rate);
        fclose(fileout);

        /********************************************************************************************/
        /********************************************************************************************/
        /* Open the file: potential energy */
        sprintf(filename, "%sS_potential_energy.dat", output->data_dir);
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */
        Output_check_file_open(fileout, filename);
        fprintf(fileout, "%f ", time);
        int iconc, NConc=params->NConc;
        for (iconc=0; iconc<NConc; iconc++) {

            /* Potential energy associated with the current concentration field */
            fprintf(fileout, "%2.14f ", c[iconc]->W_Ep);
        } /* for iconc */
        fprintf(fileout, "\n");
        fclose(fileout);
        /********************************************************************************************/
        /********************************************************************************************/

        /* Open the file for passive potential energy */
        sprintf(filename, "%sS_passive_potential_energy.dat", output->data_dir);
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */
        Output_check_file_open(fileout, filename);
        fprintf(fileout, "%f ", time);
        for (iconc=0; iconc<NConc; iconc++) {
            fprintf(fileout, "%2.14f ", c[iconc]->W_Ep_passive);
        } /* for iconc */
        fprintf(fileout, "\n");
        fclose(fileout);
        /********************************************************************************************/
        /********************************************************************************************/

        /* Open the file potential energy rate */
        sprintf(filename, "%sS_potential_energy_rate.dat", output->data_dir);
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */
        Output_check_file_open(fileout, filename);

        fprintf(fileout, "%f ", time);
        for (iconc=0; iconc<NConc; iconc++) {
            fprintf(fileout, "%2.14f ", c[iconc]->W_Ep_rate);
        } /* for iconc */
        fprintf(fileout, "\n");
        fclose(fileout);
        /********************************************************************************************/
        /********************************************************************************************/

        /* Open the file: Stokes dissipation rate */
        sprintf(filename, "%sS_Stokes_diss_rate.dat", output->data_dir);
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */
        Output_check_file_open(fileout, filename);

        fprintf(fileout, "%f ", time);
        for (iconc=0; iconc<NConc; iconc++) {

            fprintf(fileout, "%2.14f ", c[iconc]->W_Stokes_dissipation_rate);
        } /* for iconc */
        fprintf(fileout, "\n");
        fclose(fileout);
        /********************************************************************************************/
        /********************************************************************************************/

        /* Open the file: Stokes dissipation rate */
        sprintf(filename, "%sS_kinetic_potential_conversion_rate.dat", output->data_dir);
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */
        Output_check_file_open(fileout, filename);

        fprintf(fileout, "%f ", time);
        for (iconc=0; iconc<NConc; iconc++) {

            fprintf(fileout, "%2.14f ", c[iconc]->W_kinetic_potential_conversion_rate);
        } /* for iconc */
        fprintf(fileout, "\n");
        fclose(fileout);
        /********************************************************************************************/
        /********************************************************************************************/

        /* Open the file: Ep passive rate */
        sprintf(filename, "%sS_passive_potential_energy_rate.dat", output->data_dir);
        if (time > 1e-11) {

            fileout = fopen(filename, "a");
        } else { /* overwrite pre-existing file not from this simulation */

            fileout = fopen(filename, "w");
        } /* else */
        Output_check_file_open(fileout, filename);

        fprintf(fileout, "%f ", time);
        for (iconc=0; iconc<NConc; iconc++) {

            fprintf(fileout, "%2.14f ", c[iconc]->W_Ep_passive_rate);
        } /* for iconc */
        fprintf(fileout, "\n");
        fclose(fileout);
        /********************************************************************************************/
        /********************************************************************************************/

    } /* if MASTER */
}
/*************************************************************/

/* This function writes 
Shear stress on the bottom to file 
*/
void Output_write_ascii_shear_stress_bottom(Output *output, Velocity *u, MAC_grid *grid, Parameters *params, int timestep) {

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        int filenumber = 10000+timestep;
        char filename[MAX_PATH_LEN];
        sprintf(filename, "%sS_shear_stress_%d.dat", output->data_dir, filenumber);

        FILE *fileout;
        fileout = fopen(filename, "w");
        if (fileout == NULL) {

            PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filename);
        }

        double *xc = grid->xc;
        double *zc = grid->zc;
        int NX = grid->NX;
        int NZ = grid->NZ;
        int i, k;
        for (k=0; k<NZ; k++) {
            for(i=0; i<NX; i++) {

                fprintf(fileout, "%f %f %2.14lf\n", (float) zc[k], (float) xc[i], u->W_shear_stress_bottom[k][i]);
            }
        }
        fclose(fileout);

    } /* if MASTER */
}
/*************************************************************/

/* This function writes
Shear stress on the bottom to file
*/
void Output_write_ascii_shear_stress_and_tang_velocity_bottom(Output *output, Velocity *u, MAC_grid *grid, Parameters *params, int timestep) {

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        int filenumber = 10000+timestep;
        char filename[MAX_PATH_LEN];
        sprintf(filename, "%sS_shear_stress_and_tang_velocity_%d.dat", output->data_dir, filenumber);

        FILE *fileout;
        fileout = fopen(filename, "w");
        if (fileout == NULL) {

            PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filename);
        }

        double *xc = grid->xc;
        double *zc = grid->zc;
        int NX = grid->NX;
        int NZ = grid->NZ;
        int i, k;
        for (k=0; k<NZ; k++) {
            for(i=0; i<NX; i++) {

                fprintf(fileout, "%f %f %2.14lf %2.14lf %2.14lf %2.14lf\n", (float) zc[k], (float) xc[i], u->W_shear_stress_bottom[k][i], u->G_u_shear[k][i], u->G_v_shear[k][i], u->G_w_shear[k][i]);
            }
        }
        fclose(fileout);

    } /* if MASTER */
}
/*************************************************************/



/* This function writes the cell centered signed distance function to file */
void Output_write_surface_to_file(Output *output, SurfaceType *surf, MAC_grid *grid, Parameters *params, char which_quantity) {

    FILE *file_c;
    FILE *file_2;
    int i, j, k;
    char filename[MAX_PATH_LEN];
    int Is, Js, Ks;
    int Ie, Je, Ke;

    /* generate the file name */
    sprintf(filename, "%ssurf_%c_P%d.dat", output->data_dir, which_quantity, params->rank);

    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    double ***data=NULL;

    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        data = surf->u_sdf;
        break;
    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        data = surf->v_sdf;
        break;
    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        data = surf->w_sdf;
        break;
    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        data = surf->c_sdf;
        break;
    } /* switch */

    /* open the file */
    file_c = fopen(filename, "w");
    if (file_c == NULL) {

        PetscPrintf(PCW, "Output.c/ Error opening %s file ...\n", filename);
    } /* if */

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    sprintf(filename, "%ssurf_%c_2d_P%d.dat", output->data_dir, which_quantity, params->rank);
    file_2 = fopen(filename, "w");
    if (file_2 == NULL) {

        PetscPrintf(PCW, "Output.c/ Error opening %s file ...\n", filename);
    } /* if */

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                fprintf(file_c, "%f %f %f %f\n", xq[i], yq[j], zq[k], data[k][j][i]);
                if (k == Ke/2) {

                    fprintf(file_2, "%f %f %f\n", xq[i], yq[j], data[Ke/2][j][i]);
                }
            } /* for i */
        } /* for j */
    } /* for k */

    fclose(file_c);
    fclose(file_2);
}
/*************************************************************/

/* This function writes the information for the immersed nodes to the files */
void Output_immersed_info(Output *output, MAC_grid *grid, Parameters *params, char which_quantity) {

    FILE *file;
    int rank;
    char filename[MAX_PATH_LEN];
    ImmersedNode *ib_node=NULL;
    Immersed *q_immersed=NULL;

    /* rank of the current processor */
    rank = params->rank;

    /* filename */
    sprintf(filename, "%c_immersed_P%d.dat", which_quantity, rank);
    file = fopen(filename, "w");
    if (file == NULL) {

        PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filename);
    } /* if */

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
    } /* switch */

    int g;
    for (g=0; g<q_immersed->N; g++) {

        ib_node = Immersed_get_ib_node(q_immersed, g);
        fprintf(file, "*********************************************************************\n");
        fprintf(file, "im(%f,%f,%f) at (%d,%d,%d) ct=(%f,%f,%f) image:(%f,%f,%f)\n", ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z, ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index, ib_node->control_point.x, ib_node->control_point.y, ib_node->control_point.z, ib_node->image_point.x, ib_node->image_point.y, ib_node->image_point.z);
        fprintf(file, "Nfluid:%d norm(%f,%f,%f) rhs_coef:%f\n", ib_node->n_fluid, ib_node->n.vx, ib_node->n.vy, ib_node->n.vz, ib_node->boundary_coef);
        //printf("*********************************************************************\n");
        //printf("im(%f,%f,%f) at (%d,%d,%d) ct=(%f,%f,%f) image:(%f,%f,%f)\n", ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z, ib_node->im_index.x_index, ib_node->im_index.y_index, ib_node->im_index.z_index, ib_node->control_point.x, ib_node->control_point.y, ib_node->control_point.z, ib_node->image_point.x, ib_node->image_point.y, ib_node->image_point.z);
        //printf("Nfluid:%d norm(%f,%f,%f) rhs_coef:%f\n", ib_node->n_fluid, ib_node->n.vx, ib_node->n.vy, ib_node->n.vz, ib_node->boundary_coef);
        //getchar();
        int gg;
        for (gg=0; gg<ib_node->n_fluid; gg++) {

            fprintf(file, "Fluid%d(%f,%f,%f) index(%d,%d,%d) coef:%f\n", gg, ib_node->fluid_point[gg].x, ib_node->fluid_point[gg].y, ib_node->fluid_point[gg].z, ib_node->fluid_index[gg].x_index, ib_node->fluid_index[gg].y_index, ib_node->fluid_index[gg].z_index, ib_node->fluid_coef[gg]);
        }
        fprintf(file, "\n*********************************************************************\n");
    } /* for g */

    fclose(file);
}
/*************************************************************/

/* This function writes the information for the immersed nodes to the files */
void Output_immersed_control(Output *output, MAC_grid *grid, Parameters *params, char which_quantity) {

    FILE *file1, *file2;
    int g;
    int rank;
    char filename[MAX_PATH_LEN];
    ImmersedNode *ib_node;
    Immersed *q_immersed=NULL;

    /* rank of the current processor */
    rank = params->rank;

    /* filename */
    sprintf(filename, "%s%c_immersed_control_P%d.dat", output->immersed_dir, which_quantity, rank);
    file1 = fopen(filename, "w");
    if (file1 == NULL) {

        PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filename);
    } /* if */
    sprintf(filename, "%s%c_immersed_im_P%d.dat", output->immersed_dir, which_quantity, rank);
    file2 = fopen(filename, "w");
    if (file2 == NULL) {

        PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filename);
    } /* if */

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
    } /* switch */

    for (g=0; g<q_immersed->N; g++) {

        ib_node = Immersed_get_ib_node(q_immersed, g);
        fprintf(file1, "%f %f %f %f %f %f\n", ib_node->control_point.x, ib_node->control_point.y, ib_node->control_point.z, ib_node->n.vx, ib_node->n.vy, ib_node->n.vz);
        fprintf(file2, "%f %f %f\n", ib_node->im_point.x, ib_node->im_point.y, ib_node->im_point.z);
    } /* for g */

    fclose(file1);
    fclose(file2);
}
/*************************************************************/


/* This function copies the values of the nodes right before the last plane into the half cell added. This is done just for better visualization 
 * quality since the values of those nodes are always automatically solved as zero */
void Output_update_last_planes(Output *output, Vec G_data, Parameters *params) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    double ***data;
    int ierr;
    int NX, NY, NZ;

    /* total number of cells including the half cell added to the end */
    NX = params->NX+1;
    NY = params->NY+1;
    NZ = params->NZ+1;

    /* get 3D array data's pointer*/
    ierr = DAVecGetArray(output->DA_3D, G_data, (void ***)&data);PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    Is = output->G_Is;
    Js = output->G_Js;
    Ks = output->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = output->G_Ie;
    Je = output->G_Je;
    Ke = output->G_Ke;

    /* if the processor owns the Plane at,  x=Lx */
    if (Ie == NX) {

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                data[k][j][Ie-1] = data[k][j][Ie-2];
            } /* for j */
        } /* for k */

    } /* if */

    /* if the processor owns the Plane at,  y=Ly */
    if (Je == NY) {

        for (k=Ks; k<Ke; k++) {
            for (i=Is; i<Ie; i++) {

                data[k][Je-1][i] = data[k][Je-2][i];
            } /* for j */
        } /* for k */

    } /* if */

    /* if the processor owns the Plane at,  z=Lz */
    if (Ke == NZ) {

        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                data[Ke-1][j][i] = data[Ke-2][j][i];
            } /* for j */
        } /* for k */

    } /* if */

    /* restore the array */
    ierr = DAVecRestoreArray(output->DA_3D, G_data, (void ***)&data);PETScErrAct(ierr);
}
/*************************************************************/

/* This function writes 2D dumped deposit height to ASCII file for each concentration field */
void Output_write_ascii_deposit_height_dumped(Output *output, Concentration **c, MAC_grid *grid, Parameters *params, int timestep) {

    FILE *fileout;
    int i, k, NI, NK;
    int iconc, NConc;
    int filenumber;
    char filename[MAX_PATH_LEN];

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        NConc      = params->NConc;
        filenumber = 10000 + timestep;

        NK = c[0]->NK;
        NI = c[0]->NI;

        for (iconc=0; iconc<NConc; iconc++) {

            /* Only write to file if we have particle, i.e. u_s > 0.0 field */
            if (c[iconc]->Type == PARTICLE) {

                sprintf(filename, "%sS_deposit_height_dumped_c%d_%d.dat", output->data_dir, iconc, filenumber);
                /* open the new file */
                fileout = fopen(filename, "w");

                /* Number of grid points in z and x directions */
                for (k=0; k<NK; k++) {
                    for (i=0; i<NI; i++) {

                        /* write the deposited height */
                        fprintf(fileout, "%f %f %2.12f\n", grid->zc[k], grid->xc[i], c[iconc]->W_deposit_height_dumped[k][i]);
                    } /* for i */
                } /* for k*/
                fclose(fileout);
            } /* if */

        } /* for */
    }/* if */
}
/*************************************************************/

/* This function saved the cell-centered velocity into a raw binary file */
void Output_binary_3D_vorticity_vector(Output *output, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, int timestep) {

    int ierr;
    void *data_1D = output->L_data_1D;

    /* Generate the name of the file based on the current time step */
    int filenumber = 10000 + timestep;

    char bin_filename[MAX_PATH_LEN];
    sprintf(bin_filename, "%sVorticity_P%d_%d.bin", output->bin_dir, params->rank, filenumber);


    /* To avoid allocating extra memory, we compute each component of vorticity vector and store the result in output->G_data */
    /* It has the same parallel layout as the original G_data vector used int grid and velocity vector */
    /* First, x-component of vorticity: w_x */
    /*************************************/
    /*************************************/
    Velocity_compute_vorticity(u, v, w, grid, 'x', &output->G_data);
    /* Copy the last nodes into the half cell added (to avoid zero values in the visualization) */
    Output_update_last_planes(output, output->G_data, params);

    /* Now, update ghost nodes */
    Communication_update_ghost_nodes(&output->DA_3D, &output->G_data, &output->L_data, 'I');

    /* get 3D array data's pointer*/
    double ***data=NULL;
    ierr = DAVecGetArray(output->DA_3D, output->L_data, (void ***)&data);PETScErrAct(ierr);

    /* First convert 3d data including ghost nodes to 1d_float array */
    Output_convert_3D_to_1D(output, data, data_1D);

    /* set file name and appropriate flags */
    Output_set_filename(output, bin_filename);
    /* write grid flag: ON*/
    Output_set_write_grid_flag(output);
    /* Do not append data. Generate new file. */
    Output_reset_append_data_flag(output);
    /* set the data pointer */
    Output_set_output_data(output, data_1D);

    /* write grid and data to file */
    Output_write_data(output);

    /* restore array */
    ierr = DAVecRestoreArray(output->DA_3D, output->L_data, (void ***)&data);PETScErrAct(ierr);

    /*************************************/
    /*************************************/
    /* y-component of vorticity: w_y */
    /*************************************/
    /*************************************/
    /* To avoid extra-storage, we compute the result in the output->G_data vector */
    /* It has the same layout as the original G_data vector */
    Velocity_compute_vorticity(u, v, w, grid, 'y', &output->G_data);

    /* Copy the last nodes into the half cell added (to avoid zero values in the visualization) */
    Output_update_last_planes(output, output->G_data, params);

    /* update the ghost nodes */
    Communication_update_ghost_nodes(&output->DA_3D, &output->G_data, &output->L_data, 'I');

    /* get 3D array data's pointer*/
    ierr = DAVecGetArray(output->DA_3D, output->L_data, (void ***)&data);PETScErrAct(ierr);

    /* First convert 3d data including ghost nodes to 1d_float array */
    Output_convert_3D_to_1D(output, data, data_1D);

    /* write grid flag: OFF*/
    Output_reset_write_grid_flag(output);
    /* Append data. Add to the existing file */
    Output_set_append_data_flag(output);
    /* set the data pointer */
    Output_set_output_data(output, data_1D);

    /* write grid and data to file */
    Output_write_data(output);

    /* restore array */
    ierr = DAVecRestoreArray(output->DA_3D, output->L_data, (void ***)&data);PETScErrAct(ierr);

    /*************************************/
    /*************************************/
    /*************************************/
    /* z-component of vorticity: w_z */
    Velocity_compute_vorticity(u, v, w, grid, 'z', &output->G_data);
    /* Copy the last nodes into the half cell added (to avoid zero values in the visualization) */
    Output_update_last_planes(output, output->G_data, params);

    /* update the ghost nodes */
    Communication_update_ghost_nodes(&output->DA_3D, &output->G_data, &output->L_data, 'I');

    /* get 3D array data's pointer*/
    ierr = DAVecGetArray(output->DA_3D, output->L_data, (void ***)&data);PETScErrAct(ierr);

    /* First convert 3d data including ghost nodes to 1d_float array */
    Output_convert_3D_to_1D(output, data, data_1D);

    /* write grid flag: OFF*/
    Output_reset_write_grid_flag(output);
    /* Append data. Add to the existing file */
    Output_set_append_data_flag(output);
    /* set the data pointer */
    Output_set_output_data(output, data_1D);

    /* write grid and data to file */
    Output_write_data(output);

    /* restore array */
    ierr = DAVecRestoreArray(output->DA_3D, output->L_data, (void ***)&data);PETScErrAct(ierr);

}
/*************************************************************/

/* This function sets the output precision */
void Output_set_precision(Output *output, int precision_type) {

    if ( (precision_type == GVG_FLOAT) || (precision_type == GVG_DOUBLE) ) {

        output->precision = precision_type;
    } else {

        PetscPrintf(PCW, "Output.c/ Warning! Only 'float' or 'double' can be used for binary output. Setting it to 'float'\n");
        output->precision = GVG_FLOAT;
    }
}
/*************************************************************/

/* This function stores the grid information into the corresponding files for all flow quantities in binary format.  */
/* This will be used for data conversion */
/* Note that the order is x, y and z */
void Output_grid_binaries(Output *output, MAC_grid *grid, char which_quantity)  {

    double *xq = NULL;
    double *yq = NULL;
    double *zq = NULL;
    char filename[400];

    switch (which_quantity) {
    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        sprintf(filename, "%su_grid.bin", output->bin_dir);
        break;

    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        sprintf(filename, "%sv_grid.bin", output->bin_dir);
        break;

    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        sprintf(filename, "%sw_grid.bin", output->bin_dir);
        break;

    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        sprintf(filename, "%sc_grid.bin", output->bin_dir);
        break;

    case 'p':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        sprintf(filename, "%sp_grid.bin", output->bin_dir);
        break;

    } /* switch */

    FILE *file = fopen(filename, "wb");
    if (file == NULL) {

        printf("Grid.c/ Warning!! Could not open %s file...\n", filename);
    } /* if */

    fwrite(xq, sizeof(double), grid->NX, file);
    fwrite(yq, sizeof(double), grid->NY, file);
    fwrite(zq, sizeof(double), grid->NZ, file);

    fclose(file);
}
/*************************************************************/

/* This function writes the entire parallel data_vec into one file  (binary format) */
/* This is very useful if the user wants to reload the data for post processing */
void OLDOutput_petsc_vec(Vec *data_vec, char *filename) {

    PetscViewer writer;
    int ierr;

    ierr = PetscViewerCreate(PCW, &writer); PETScErrAct(ierr);
    ierr = PetscViewerSetType(writer, PETSC_VIEWER_BINARY); PETScErrAct(ierr);
    ierr = PetscViewerFileSetMode(writer, FILE_MODE_WRITE); PETScErrAct(ierr);
    ierr = PetscViewerBinarySkipInfo(writer); PETScErrAct(ierr); PETScErrAct(ierr);
    ierr = PetscViewerFileSetName(writer, filename); PETScErrAct(ierr);

    //ierr = PetscViewerBinaryOpen(PCW, filename, FILE_MODE_WRITE, &writer); PETScErrAct(ierr);
    ierr = VecView(*data_vec, writer); PETScErrAct(ierr);
    ierr = PetscViewerDestroy(writer); PETScErrAct(ierr);
}
/*************************************************************/

/* This function writes a PETSc matrix to file in BINARY or ASCII format */
void Output_petsc_mat(Mat *mat, char *filename, short int output_type) {

    PetscViewer writer;
    int ierr;

    if (output_type == BINARY) {

        ierr = PetscViewerBinaryOpen(PCW, filename, FILE_MODE_WRITE, &writer); PETScErrAct(ierr);

    } else if (output_type == ASCII){

        ierr = PetscViewerASCIIOpen(PCW, filename, &writer); PETScErrAct(ierr);
    }

    ierr = MatView(*mat, writer); PETScErrAct(ierr);
    ierr = PetscViewerDestroy(writer); PETScErrAct(ierr);
}

/* Write the sedimentation rate for each concentration field */
void Output_write_ascii_sed_rate(Output *output, Concentration **c, MAC_grid *grid, Parameters *params, double time, int timestep) {

    /* Only write on processor zero */
    if (params->rank == MASTER) {

        FILE *fileout;
        char filename[MAX_PATH_LEN];
        int iconc;
        for (iconc=0; iconc<params->NConc; iconc++) {

            /* Open the file */
            sprintf(filename, "%sS_sedim_rate_c%d.dat", output->data_dir, iconc);
            if (time > 1e-11) {

                fileout = fopen(filename, "a");
            } else { /* overwrite pre-existing file not from this simulation */

                fileout = fopen(filename, "w");
            } /* else */

            /* Kinetic energy of the fluid */
            fprintf(fileout, "%f %f\n", time, c[iconc]->W_sed_rate);
            fclose(fileout);
            /********************************************************************************************/
        } /* for iconc */

        for (iconc=0; iconc<params->NConc; iconc++) {

            if (c[iconc]->Type  == PARTICLE) {
                /* Open the file */
                int filenumber = 10000 + timestep;
                sprintf(filename, "%sS_sedim_rate_bottom_c%d_%d.dat", output->data_dir, iconc, filenumber);

                fileout = fopen(filename, "w");
                if (fileout == NULL) {

                    PetscPrintf(PCW, "Output.c/ Could not open %s file\n", filename);
                }
                double *xc = grid->xc;
                double *zc = grid->zc;
                int NX = grid->NX;
                int NZ = grid->NZ;
                int i, k;
                for (k=0; k<NZ; k++) {
                    for(i=0; i<NX; i++) {

                        fprintf(fileout, "%f %f %2.14lf\n", (float) zc[k], (float) xc[i], c[iconc]->W_sed_rate_bottom[k][i]);
                    } /* for i */
                } /* for k */

                fclose(fileout);
            } /* if PARTICLE */
            /********************************************************************************************/
        }
    } /* if MASTER */
}
/*************************************************************/

/* This function writes mapping between the immersed points and their control points imm[k][i] into control[k][i] */
void Output_immersed_control_new(Output *output, MAC_grid *grid, Parameters *params, char which_quantity) {

    char filename[MAX_PATH_LEN];
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
    } /* switch */

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    double **G_point_x = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **G_point_y = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **G_point_z = Memory_allocate_2D_double(grid->NX, grid->NZ);

    double **W_point_x = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **W_point_y = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **W_point_z = Memory_allocate_2D_double(grid->NX, grid->NZ);

    int i, k;
    for (k=Ks; k<Ke; k++) {
        for (i=Is; i<Ie; i++) {
            int j_s = grid->interface_y_index[k][i];
            if ( (j_s >= Js) && (j_s < Je) && (j_s != -1)) {

                /* Now, go through the fluid nodes, find the global index in the matrix, copy the value and use MatSetValues() instead of MatSetValuesStencil. This is done b/c some of the fluid node might be on the neighboring processor and out of the range of the ghost layers. */
                int ib_index = Immersed_get_ib_global_index(q_immersed, i, j_s, k);
                ImmersedNode *ib_node  = Immersed_get_ib_node(q_immersed, ib_index);

                G_point_x[k][i] = ib_node->control_point.x;
                G_point_y[k][i] = ib_node->control_point.y;
                G_point_z[k][i] = ib_node->control_point.z;
            } /* if */
        } /* for */
    } /* for */

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

    Communication_reduce_2D_arrays(G_point_x, W_point_x, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(G_point_y, W_point_y, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(G_point_z, W_point_z, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);

    if (params->rank == MASTER) {
        /* filename */
        sprintf(filename, "%s%c_immersed_control_point_mapping.dat", output->immersed_dir, which_quantity);
        FILE *file1 = fopen(filename, "w");
        if (file1 == NULL) {

            PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filename);
        } /* if */

        for (k=0; k<grid->NK; k++) {
            for (i=0; i<grid->NI; i++) {
                fprintf(file1, "%lf %lf %lf %lf\n", grid->zc[k], grid->xc[i], W_point_z[k][i], W_point_x[k][i]);
            } /* for i */
        } /* for k */
        fclose(file1);
    } /* if MASTER */

    Memory_free_2D_array(grid->NZ, (void **)G_point_x);
    Memory_free_2D_array(grid->NZ, (void **)G_point_y);
    Memory_free_2D_array(grid->NZ, (void **)G_point_z);
    Memory_free_2D_array(grid->NZ, (void **)W_point_x);
    Memory_free_2D_array(grid->NZ, (void **)W_point_y);
    Memory_free_2D_array(grid->NZ, (void **)W_point_z);
}
/*************************************************************/

/* This function writes
Shear stress on the bottom to file
*/
void Output_write_ascii_shear_stress_bottom_test(Output *output, Velocity *u, MAC_grid *grid, Parameters *params, int timestep) {

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    double **G_point_x = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **G_point_y = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **G_point_z = Memory_allocate_2D_double(grid->NX, grid->NZ);

    double **W_point_x = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **W_point_y = Memory_allocate_2D_double(grid->NX, grid->NZ);
    double **W_point_z = Memory_allocate_2D_double(grid->NX, grid->NZ);

    int i, k;
    for (k=Ks; k<Ke; k++) {
        for (i=Is; i<Ie; i++) {
            int j_s = grid->interface_y_index[k][i];
            if ( (j_s >= Js) && (j_s < Je) && (j_s != NOT_FOUND)) {

                /* Now, go through the fluid nodes, find the global index in the matrix, copy the value and use MatSetValues() instead of MatSetValuesStencil. This is done b/c some of the fluid node might be on the neighboring processor and out of the range of the ghost layers. */
                int ib_index = Immersed_get_ib_global_index(grid->c_immersed, i, j_s, k);
                if (ib_index != NOT_FOUND) {
                    ImmersedNode *ib_node  = Immersed_get_ib_node(grid->c_immersed, ib_index);

                    G_point_x[k][i] = ib_node->control_point.x;
                    G_point_y[k][i] = ib_node->control_point.y;
                    G_point_z[k][i] = ib_node->control_point.z;
                }
            } /* if */
        } /* for */
    } /* for */

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
    Communication_reduce_2D_arrays(G_point_x, W_point_x, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(G_point_y, W_point_y, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);
    Communication_reduce_2D_arrays(G_point_z, W_point_z, &G_s, &G_e, &W_e, REDUCE_TO_MASTER, params);


    /* Only write on processor zero */
    if (params->rank == MASTER) {

        int filenumber = 10000+timestep;
        char filenamen[MAX_PATH_LEN];
        sprintf(filenamen, "%sS_shear_stress_new_%d.dat", output->data_dir, filenumber);

        FILE *fileoutn = fopen(filenamen, "w");
        if (fileoutn == NULL) {

            PetscPrintf(PCW, "Output.c/ Error opening %s file\n", filenamen);
        }

        int NX = grid->NX;
        int NZ = grid->NZ;
        int i, k;
        for (k=0; k<NZ-1; k++) {
            for(i=0; i<NX-1; i++) {

                fprintf(fileoutn, "%lf %lf %2.14lf\n", W_point_z[k][i], W_point_x[k][i], u->W_shear_stress_bottom[k][i]);
            }
        }
        fclose(fileoutn);

    } /* if MASTER */

    Memory_free_2D_array(grid->NZ, (void **)G_point_x);
    Memory_free_2D_array(grid->NZ, (void **)G_point_y);
    Memory_free_2D_array(grid->NZ, (void **)G_point_z);
    Memory_free_2D_array(grid->NZ, (void **)W_point_x);
    Memory_free_2D_array(grid->NZ, (void **)W_point_y);
    Memory_free_2D_array(grid->NZ, (void **)W_point_z);
}
/*************************************************************/

void Output_3D_data(Vec v, MAC_grid *grid, char *filename, char which_quantity) {

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

    double ***data=NULL;
    /* Get regular data array for v vector */
    int ierr = DAVecGetArray(grid->DA_3D, v, (void ***)&data); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    FILE *f=fopen(filename, "w");
    int i, j, k;
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                int status = Grid_get_q_status(grid, i, j, k);
                fprintf(f, "(i,j,k)=(%d,%d,%d) s:%d val:%2.14f\n", i, j, k, status, data[k][j][i]);
            }
        }
    }

    /* Get regular data array for v vector */
    ierr = DAVecRestoreArray(grid->DA_3D, v, (void ***)&data); PETScErrAct(ierr);

    fclose(f);
}
/*************************************************************/

/* This function writes the grid coordinates u or v or w or c or p to file in HDF5 format */
void Output_grid_q_hdf5(Output *output, MAC_grid *grid, char which_quantity) {

    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    switch (which_quantity) {
    case 'u':
        xq = grid->xu; yq = grid->yu; zq = grid->zu; break;
    case 'v':
        xq = grid->xv; yq = grid->yv; zq = grid->zv; break;
    case 'w':
        xq = grid->xw; yq = grid->yw; zq = grid->zw; break;
    case 'c':
        xq = grid->xc; yq = grid->yc; zq = grid->zc; break;
    case 'p':
        xq = grid->xc; yq = grid->yc; zq = grid->zc; break;
    } /* switch */

    char filename[FILENAME_MAX];
    sprintf(filename, "%s/%c_grid.h5", output->hdf5_dir, which_quantity);

    PetscViewer    H5viewer;
    // Create the HDF5 viewer
    //int ierr = PetscViewerHDF5Open(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &H5viewer); PETScErrAct(ierr);

    char g[3];
    sprintf(g, "x");
    Output_grid_hdf5(xq, grid->NX, g, H5viewer) ;

    sprintf(g, "y");
    Output_grid_hdf5(yq, grid->NY, g, H5viewer) ;

    sprintf(g, "z");
    Output_grid_hdf5(zq, grid->NZ, g, H5viewer) ;

    int ierr = PetscViewerDestroy(H5viewer); PETScErrAct(ierr);

}
/*************************************************************/

/* This function writes the grid information to file in HDF5 format. Only processor zero does the entire writing */
void Output_grid_hdf5(double *coord, int N, char *gridname, PetscViewer H5viewer) {

    Vec x;

    int ierr = VecCreate(PETSC_COMM_SELF, &x); PETScErrAct(ierr);
    ierr = VecSetSizes(x, N, N); PETScErrAct(ierr);
    ierr = VecSetFromOptions(x); PETScErrAct(ierr);

    double *vcoord=NULL;
    ierr = VecGetArray(x, &vcoord); PETScErrAct(ierr);

    int i;
    for (i=0; i<N; i++) {
        vcoord[i] = coord[i];
    }

    /* set the name. Important when accessing the data in visualization software */
    ierr = VecRestoreArray(x, (void *)&vcoord); PETScErrAct(ierr);
    ierr = PetscObjectSetName((PetscObject) x, gridname); PETScErrAct(ierr);

    int rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); PETScErrAct(ierr);

    // Write the H5 file only on processor zero
    if (rank == MASTER) {

        ierr = VecView(x, H5viewer); PETScErrAct(ierr);
    } /* if */
    ierr = VecDestroy(x); PETScErrAct(ierr);
}
/*************************************************************/

/* This function writes any parallel vector data to the file. */
/* Given formats are ASCII, BINARY, HDF5 */
void Output_petsc_vec(Output *output, Vec q, char *dataname, char *filename) {

    int ierr;
    PetscViewer writer;

    switch (output->type) {

    case HDF5:
        //ierr = PetscViewerHDF5Open(PCW, filename, FILE_MODE_WRITE, &writer); PETScErrAct(ierr); PETScErrAct(ierr);
        break;

    case BINARY:

        ierr = PetscViewerCreate(PCW, &writer); PETScErrAct(ierr);
        ierr = PetscViewerSetType(writer, PETSC_VIEWER_BINARY); PETScErrAct(ierr);
        ierr = PetscViewerFileSetMode(writer, FILE_MODE_WRITE); PETScErrAct(ierr);
        ierr = PetscViewerBinarySkipInfo(writer); PETScErrAct(ierr); PETScErrAct(ierr);
        ierr = PetscViewerFileSetName(writer, filename); PETScErrAct(ierr);
        break;

    case ASCII:

        ierr = PetscViewerASCIIOpen(PCW, filename, &writer); PETScErrAct(ierr);
        break;

    default:
        PetscPrintf(PCW, "Output.c/ Error! Could not create PetscViewer based on the given output type [%d] (be only BINARY, HDF5 or ASCII (not recommended))\n", output->type);
    }/* switch */

    /* Set the dataname for the current vector. */
    /* It is very IMPORTANT when the output format is HDF5. This name is essentially the SAME name used when accessing the data in any visualization software *.h5:dataname */
    /* Also, when loading the data using and/or VecLoadIntoVector() and HDF5 format, user should set the name the same before calling this command. */
    ierr = PetscObjectSetName((PetscObject) q, dataname); PETScErrAct(ierr);

    /* Write the data to file in the given format */
    ierr = VecView(q, writer); PETScErrAct(ierr);

    /* Close the file and destroy the writer */
    ierr = PetscViewerDestroy(writer); PETScErrAct(ierr);

}

/* This function saves the desired flow properties based on the set flag */
void Output_flow_properties(Output *output, GVG_bag *data_bag) {

    Velocity *u = data_bag->u;
    Velocity *v = data_bag->v;
    Velocity *w = data_bag->w;

    Concentration **c = data_bag->c;
    Pressure *p = data_bag->p;
    MAC_grid *grid = data_bag->grid;
    Parameters *params = data_bag->params;

    int timestep = data_bag->output_iter;
    double time = data_bag->time;

    int NConc = params->NConc;
    int iconc;

    /***********************************************/
    /* First, field data (3D) */
    /***********************************************/
    if (params->u_output) {

        Output_velocity(output, u, timestep, params);
    }
    if (params->v_output) {

        Output_velocity(output, v, timestep, params);
    }
    if (params->w_output) {

        Output_velocity(output, w, timestep, params);
    }

    if (params->conc_output) {

        for (iconc=0; iconc<NConc; iconc++) {

            Output_concentration(output, c[iconc], timestep, params);
        } /* for */
    }

    if (params->div_output) {

        Output_divergence(output, p, timestep, params);
    }

    /********************************************************/
    /* Now, 2-D or 1-D data */
    /********************************************************/
    if (params->susp_mass_output) {

        for (iconc=0; iconc<NConc; iconc++) {

            Conc_compute_total_suspended_mass(c[iconc], grid);
        }
        Output_write_ascii_suspended_mass(output, c, params, time);
        PetscPrintf(PCW, "Output.c/ Suspended mass has been written to file successfully.\n");
    }

    if (params->ave_height_output) {

        /* First compute the average height for each concentration field */
        Conc_compute_ave_height_x(c, grid, params);

        /* Now write them to file for each conc. field vs x*/
        Output_write_ascii_average_height(output, c, grid, params, timestep) ;
        PetscPrintf(PCW, "Output.c/ Average height has been written to file successfully.\n");

    }

    if (params->front_location_output) {

        if (!params->ave_height_output) {

            Conc_compute_ave_height_x(c, grid, params);
        } /* else, already was calculated before, avoid repeated work! */

        /* Limit to capture the front location */
        double front_limit = 0.01;
        double tip = Conc_find_front_location(c, grid, params, front_limit) ;

        Output_write_ascii_front_location(output, tip, params, time) ;
        PetscPrintf(PCW, "Output.c/ Front location has been written to file successfully.\n");

    } /* if */

    if (params->deposit_height_output) {

        Conc_update_world_deposited_height(c, grid, params) ;
        Output_write_ascii_deposit_height(output, c, grid, params, timestep) ;
        PetscPrintf(PCW, "Output.c/ Deposit height has been written to file successfully.\n");

    } /* if */

    if (params->energies_output) {

        for (iconc=0; iconc<NConc; iconc++) {

            Conc_compute_potential_energies(c[iconc], v, grid, params) ;
        } /* for */

        Extract_viscous_dissipation_rate(u, v, w, grid, params);
        Extract_total_kinetic_energy(u, v, w, grid);
        Extract_kinetic_potential_conversion_rate(v, c, grid, params);

        u->W_kinetic_potential_conversion_rate = v->W_kinetic_potential_conversion_rate;
        Output_write_ascii_energies(output, c, u, params, time);

        PetscPrintf(PCW, "Output.c/ Energ(ies) has been written to file successfully.\n");

    }

    if (params->shear_stress_output) {

        Extract_wall_shear(u, v, w, grid, params);
        //Extract_wall_shear_quadratic(u, v, w, grid, params);
        Output_write_ascii_shear_stress_bottom(output, u, grid, params, timestep);
        Output_write_ascii_shear_stress_bottom_test(output, u, grid, params, timestep);

        PetscPrintf(PCW, "Output.c/ Bottom shear stress has been written to file successfully.\n");

    }

    if (params->sedim_rate_output) {

        for (iconc=0; iconc<NConc; iconc++) {

            Conc_compute_sedimentation_rate(c[iconc], grid, params);
            Conc_bottom_sed_rate(c[iconc], grid, params) ;
        }
        Output_write_ascii_sed_rate(output, c, grid, params, time, timestep) ;
        PetscPrintf(PCW, "Output.c/ Sedimentation rate has been written to file successfully.\n");

    }


    if (params->p_output) {

        Pressure_compute_real_pressure(p, grid);
        Output_pressure(output, p, timestep, params);
        /* DO NOT COMMENT THIS LINE OUT! You might get some weired behavior in pressure field after outputing stage */
        //VecSet(p->G_data, 0.0);
        //VecSet(p->G_b, 0.0);
    }

	if (params->flow_over_sphere) {

		Extract_sphere_Drag_coefficient(u, v, w, p, grid, params, time) ;
		VecSet(p->G_data, 0.0); 
        	VecSet(p->G_b, 0.0);
	}

}
/************************************************************************/

/* Write the data for velocity (u or v or w) */
void Output_velocity(Output *output, Velocity *vel, int timestep, Parameters *params) {

    int filenumber = 10000 + timestep;
    char filename[MAX_PATH_LEN];

    /* Generate the filename (including the output path) */
    sprintf(filename, "%sVel_%c_%d.%s", output->data_3d_dir, vel->component, filenumber, output->data_3d_file_ext);

    /* copy data */
    int ierr = VecCopy(vel->G_data, output->G_data); PETScErrAct(ierr);

    /* update the values at the last (extra added point) to avoid bad visualization */
    if (output->update_last_plane) {
        Output_update_last_planes(output, output->G_data, params);
    }

    char dataname[5];
    sprintf(dataname, "%c", vel->component);
    Output_petsc_vec(output, output->G_data, dataname, filename) ;
}
/************************************************************************/

/* Writes the concentration data */
void Output_concentration(Output *output, Concentration *c, int timestep, Parameters *params) {

    int filenumber = 10000 + timestep;
    char filename[MAX_PATH_LEN];

    /* Generate the filename (including the output path) */
    sprintf(filename, "%sConc_%d_%d.%s", output->data_3d_dir, c->conc_index, filenumber, output->data_3d_file_ext);

    /* copy data */
    int ierr = VecCopy(c->G_data, output->G_data); PETScErrAct(ierr);

    /* update the values at the last (extra added point) to avoid bad visualization */
    if (output->update_last_plane) {

        Output_update_last_planes(output, output->G_data, params);
    } /* if */

    char dataname[5];
    sprintf(dataname, "c%d", c->conc_index);
    Output_petsc_vec(output, output->G_data, dataname, filename) ;
}
/************************************************************************/

/* Writes the pressure */
void Output_pressure(Output *output, Pressure *p, int timestep, Parameters *params) {

    int filenumber = 10000 + timestep;
    char filename[MAX_PATH_LEN];

    /* Generate the filename (including the output path) */
    sprintf(filename, "%sPress_%d.%s", output->data_3d_dir, filenumber, output->data_3d_file_ext);

    /* copy data */
    int ierr = VecCopy(p->G_data, output->G_data); PETScErrAct(ierr);

    /* update the values at the last (extra added point) to avoid bad visualization */
    if (output->update_last_plane) {

        Output_update_last_planes(output, output->G_data, params);
    }

    char dataname[5];
    sprintf(dataname, "p");
    Output_petsc_vec(output, output->G_data, dataname, filename) ;

}
/************************************************************************/

/* Writes the velocity divergence */
void Output_divergence(Output *output, Pressure *p, int timestep, Parameters *params) {

    int filenumber = 10000 + timestep;
    char filename[MAX_PATH_LEN];

    /* Generate the filename (including the output path) */
    sprintf(filename, "%sDiv_%d.%s", output->data_3d_dir, filenumber, output->data_3d_file_ext);

    /* copy data */
    int ierr = VecCopy(p->G_divergence, output->G_data); PETScErrAct(ierr);

    /* update the values at the last (extra added point) to avoid bad visualization */
    if (output->update_last_plane) {

        Output_update_last_planes(output, output->G_data, params);
    }

    char dataname[5];
    sprintf(dataname, "div");
    Output_petsc_vec(output, output->G_data, dataname, filename) ;

}
/************************************************************************/

/* This function saves the signed-distance-function to file which represents the interface */
void Output_surface(Output *output, SurfaceType *surf, Parameters *params, int timestep) {

    int filenumber = 10000 + timestep;
    char filename[MAX_PATH_LEN];

    /* Generate the filename (including the output path) */
    sprintf(filename, "%sSurface_%d.%s", output->data_3d_dir, filenumber, output->data_3d_file_ext);

    /* copy data */
    int ierr = VecCopy(surf->G_c_sdf, output->G_data); PETScErrAct(ierr);

    /* update the values at the last (extra added point) to avoid bad visualization */
    if (output->update_last_plane) {

        Output_update_last_planes(output, output->G_data, params);
    }

    char dataname[5];
    sprintf(dataname, "surf");
    Output_petsc_vec(output, output->G_data, dataname, filename) ;
}

