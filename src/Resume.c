#include "definitions.h"
#include "DataTypes.h"
#include "Resume.h"
#include "Memory.h"
#include <stdio.h>
#include <stdlib.h>

/* This function creats a Resumer. This would read previous runs's workspace data and start from that. It also,
saves the workspace data for backup */
Resume *Resume_create(int max_vec) {


    Resume *new_resume;
    new_resume = (Resume *)malloc(sizeof(Resume));
    Memory_check_allocation(new_resume);

    /* Pointer to the parallel workspace vectors */
    new_resume->Vec_data = (Vec **)malloc(max_vec*sizeof(Vec *));
    Memory_check_allocation(new_resume->Vec_data);

    /* Filename for each workspace vector */
    new_resume->IO_filename = (char **)malloc(max_vec*sizeof(char *));
    Memory_check_allocation(new_resume->IO_filename);

    int i;
    for (i=0; i<max_vec; i++) {

        new_resume->IO_filename[i] = (char *)malloc(100*sizeof(char));
        Memory_check_allocation(new_resume->IO_filename[i]);

    } /* for i */

    new_resume->max_vec  = max_vec;

    return (new_resume);
}
/***************************************************************************************************/

void Resume_destroy(Resume *resume) {

    int max_vec;

    max_vec = resume->max_vec;
    Memory_free_2D_array(max_vec, (void **)resume->IO_filename);

    free(resume->Vec_data);
    free(resume);

}
/***************************************************************************************************/

void Resume_set_write_flag(Resume *resume) {

    resume->start_write = YES;
}
/***************************************************************************************************/

void Resume_set_read_flag(Resume *resume) {

    resume->start_read = YES;
}
/***************************************************************************************************/

/* Sets the read/write flag for writing the grid info. Write grid info if the geometry changes in time due to
erosion and deposition */
void Resume_set_IO_grid_flag(Resume *resume) {

    resume->IO_grid_flag = YES;
}
/***************************************************************************************************/

void Resume_reset_write_flag(Resume *resume) {

    resume->start_write = NO;
}
/***************************************************************************************************/

void Resume_reset_read_flag(Resume *resume) {

    resume->start_read = NO;
}
/***************************************************************************************************/

/* Resets the read/write flag for writing the grid info. Write grid info if the geometry changes in time due to
erosion and deposition */
void Resume_reset_IO_grid_flag(Resume *resume) {

    resume->IO_grid_flag = NO;
}
/***************************************************************************************************/

void Resume_write_data(Resume *resume) {

    int ierr;
    int i;
    PetscViewer writer;

    if (resume->IO_type == BINARY) {

        for (i=0; i<resume->max_vec; i++) {

            ierr = PetscViewerBinaryOpen(PCW, resume->IO_filename[i], FILE_MODE_WRITE, &writer);
            ierr = VecView(*(resume->Vec_data[i]), writer);
            ierr = PetscViewerDestroy(writer);

        } /* for i */

    } else if (resume->IO_type == ASCII) {

        for (i=0; i<resume->max_vec; i++) {

            ierr = PetscViewerASCIIOpen(PCW, resume->IO_filename[i], &writer);
            ierr = VecView(*(resume->Vec_data[i]), writer);
            ierr = PetscViewerDestroy(writer);
        } /* for i */
    }
}
/***************************************************************************************************/

void Resume_read_data(Resume *resume) {

    int ierr;
    int i;
    PetscViewer reader;

    if (resume->IO_type == BINARY) {

        for (i=0; i<resume->max_vec; i++) {

            ierr = PetscViewerBinaryOpen(PCW, resume->IO_filename[i], FILE_MODE_READ, &reader);
            /* Use VecLoad only if you want to create the vector */
            //ierr = VecLoad(reader, VECMPI, resume->Vec_data[i]);
            /* load vector into an already existing vector */
            ierr = VecLoadIntoVector(reader, *resume->Vec_data[i]);
            ierr = PetscViewerDestroy(reader);

        } /* for i */

    } else if (resume->IO_type == ASCII) {

        for (i=0; i<resume->max_vec; i++) {

            /* Have to fix this. Since  VecLoad() only reads the binary parallel vector */
            /* Note that generally, ASCII I/O is extremely slow! */

        } /* for i */
    }
}
/***************************************************************************************************/

void Resume_set_IO_filename(Resume *resume, int vec_id, char *filename) {


    strcpy(resume->IO_filename[vec_id], filename);
}
/***************************************************************************************************/

void Resume_set_IO_type(Resume *resume, int IO_type) {

    /* Use BINARY. ASCII reader is not complete yet */
    resume->IO_type = IO_type; /* ASCII or BINARY */
}
/***************************************************************************************************/

/* Set the Vec pointer for the "vec_id" vector */
void Resume_set_vec_data(Resume *resume, int vec_id, Vec *data) {

    resume->Vec_data[vec_id] = data;
}
/***************************************************************************************************/

void Resume_set_data(Resume *resume, GVG_bag *data_bag) {

    char *IO_filename_conc;
    int iconc, NConc;
    int id;

    IO_filename_conc = (char *)malloc(50*sizeof(char));

    NConc = data_bag->params->NConc;
    /*
 Resume_set_max_vectors(resume, int max_vec) {
 */

    id = 0;
    char fname[100];
    /* u,v,w data */
    sprintf(fname, "Resume_u.bin");
    Resume_set_vec_data(resume, id, &(data_bag->u->G_data) );
    Resume_set_IO_filename(resume, id, fname);
    id++;

    sprintf(fname, "Resume_v.bin");
    Resume_set_vec_data(resume, id, &(data_bag->v->G_data) );
    Resume_set_IO_filename(resume, id, fname);
    id++;

    sprintf(fname, "Resume_w.bin");
    Resume_set_vec_data(resume, id, &(data_bag->w->G_data) );
    Resume_set_IO_filename(resume, id, fname);
    id++;

    /* Pressure gradients */
    sprintf(fname, "Resume_dpdx.bin");
    Resume_set_vec_data(resume, id, &(data_bag->p->G_dp_dx) );
    Resume_set_IO_filename(resume, id, fname);
    id++;

    sprintf(fname, "Resume_dpdy.bin");
    Resume_set_vec_data(resume, id, &(data_bag->p->G_dp_dy) );
    Resume_set_IO_filename(resume, id, fname);
    id++;

    sprintf(fname, "Resume_dpdz.bin");
    Resume_set_vec_data(resume, id, &(data_bag->p->G_dp_dz) );
    Resume_set_IO_filename(resume, id, fname);
    id++;

    /* Only read or write grid_c if the flag is on. */
    if (resume->IO_grid_flag) {

        /* Geometry of the cell center */
        sprintf(fname, "Resume_gridc.bin");
        Resume_set_vec_data(resume, id, &(data_bag->grid->G_c_status) );
        Resume_set_IO_filename(resume, id, fname);
        id++;
    }

    /* concentration */
    for (iconc=0; iconc<NConc; iconc++) {

        Resume_set_vec_data(resume, id, &(data_bag->c[iconc]->G_data) );
        //sprintf(IO_filename_conc, "Resume_conc%d.bin", iconc);
        //Resume_set_IO_filename(resume, id, &IO_filename_conc[0]);
        sprintf(resume->IO_filename[id], "Resume_conc%d.bin", iconc);
        id++;
    }

    /* For now, only use binary IO, ASCII case has some problems for parallel vector reading */
    Resume_set_IO_type(resume, BINARY);

    free(IO_filename_conc);
}
/***************************************************************************************************/

/* This function writes the primary essential data to file for future simulations */
void Resume_write_primary_data(Resume *resume, GVG_bag *data_bag) {

    /* Always write in ASCII */
    /* Write only on processor zero */
    /* do nothing */
    if (resume->IO_type) {}

    if (data_bag->params->rank == MASTER) {

        FILE *fileout;
        int iter;
        int output_iter;
        double dt, dt_old;
        double time;

        fileout = fopen("Resume_primary.dat", "w");
        if (fileout != NULL) {

            time        = data_bag->time;
            dt          = data_bag->dt;
            dt_old      = data_bag->dt_old;
            iter        = data_bag->iter;
            output_iter = data_bag->output_iter;

            fprintf(fileout, "time=\n%lf\ndt=\n%lf\ndt_old=\n%lf\niteration=\n%d\noutput_iter=\n%d\n", time, dt, dt_old,iter,output_iter);
            fclose(fileout);
        } else {

            PetscPrintf(PCW, "Resume.c/ Error!!!\nCould not open 'Resume_primary.dat' for writing data...\n");
        }
    }
}
/***************************************************************************************************/

/* This function reads the primary essential data from file */
void Resume_read_primary_data(Resume *resume, GVG_bag *data_bag) {

    /* Always write in ASCII */
    /* Write only on processor zero */

    /* do nothing */
    if (resume->IO_type) {}

    FILE *filein;
    filein = fopen("Resume_primary.dat", "r");
    if (filein != NULL) {

        int iter, output_iter;
        double dt, dt_old;
        double time;
        char string[50];

        int n_read=fscanf(filein, "%s\n%lf\n%s\n%lf\n%s\n%lf\n%s\n%d\n%s\n%d", &string[0], &time, &string[0], &dt, &string[0], &dt_old, &string[0], &iter, &string[0], &output_iter);

        /* do nothing */
        if (n_read) {}

        fclose(filein);

        //printf("time=\n%f\ndt=\n%f\ndt_old=\n%f\niteration=\n%d\noutput_iter=\n%d\n", time, dt, dt_old,iter,output_iter);

        data_bag->time        = time;
        data_bag->dt          = dt;
        data_bag->dt_old      = dt_old;
        data_bag->iter        = iter;
        data_bag->output_iter = output_iter;

    } else {

        PetscPrintf(PCW, "Resume.c/ Error!!!\nCould not open 'Resume_primary.dat' for reading data...\n");
    }
}
/***************************************************************************************************/

/* This function writes any 2D data data[NZ][NX] to file. Note that the array should be all stored on processor zero. */
void Resume_write_2D_data(Resume *resume, double **data, char *filename, int NX, int NZ) {

    int rank;
    int ierr = MPI_Comm_rank(PCW, &rank); PETScErrAct(ierr);

    /* Only processor zero performs writing */
    if (rank == MASTER) {

        int nt = NZ*NX;

        double *data_1D = Memory_allocate_1D_double(nt);
        Memory_convert_2D_double_to_1D(data, data_1D, NX, NZ, GVG_DOUBLE);

        if (resume->IO_type == BINARY) {

            FILE *fout = fopen(filename, "wb");
            if (fout == NULL)  {
                PetscPrintf(PCW, "Resume.c/ Error! Could not open the file:%s\n", filename);
            }
            fwrite(data_1D, sizeof(double), nt, fout);


            fclose(fout);
        } /* if */
        else {

            PetscPrintf(PCW, "Resume.c/ Warning. For now, just use BINARY format to store Resume data.\n");
        } /* else */
        free(data_1D);
    } /* if MASTER */
}

/* This function reads any 2D data data[NZ][NX] to file. Note that the array is read by all processors */
void Resume_read_2D_data(Resume *resume, double **data, char *filename, int NX, int NZ) {

    int nt = NZ*NX;
    double *data_1D = Memory_allocate_1D_double(nt);

    int rank, numprocs;
    int ierr = MPI_Comm_rank(PCW, &rank); PETScErrAct(ierr);
    ierr = MPI_Comm_size(PCW, &numprocs); PETScErrAct(ierr);

    int np;

    /* Read the data serially. This is not the most efficient way, but for now since it is performed only once, it should be ok! */
    for (np=0; np<numprocs; np++) {
        if (np == rank) {

            if (resume->IO_type == BINARY) {

                FILE *fin = fopen(filename, "rb");
                if (fin == NULL)  {
                    PetscPrintf(PCW, "Resume.c/ Error! Could not open the file:%s\n", filename);
                }
                int n_read = fread(data_1D, sizeof(double), nt, fin);
                if (n_read != nt) {

                    printf("Resume.c/ Error reading %s file\n", filename);
                }
                fclose(fin);
            } /* if */
            else {

                PetscPrintf(PCW, "Resume.c/ Warning. For now, just use BINARY format to read Resume data.\n");
            } /* else */

        } /* if np == rank */
        MPI_Barrier(PCW);
    } /* for np */

    /* convert the 1D format to 2D array [NZ][NX] */
    Memory_convert_1D_double_to_2D(data_1D, data, NX, NZ) ;

    free(data_1D);
}
