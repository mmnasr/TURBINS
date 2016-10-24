#include "definitions.h"
#include "DataTypes.h"
#include "Memory.h"
#include <stdio.h>
#include <stdlib.h>


/* This function is to check if the memory allocation was successfull for any arbitrary type */
void Memory_check_allocation(void *var) {

    if (var == NULL) {

        printf("**********************************************************\n");
        printf("!!! Error. Not enough Memory. Terminating the program !!!!\n");
        printf("**********************************************************\n");
        PetscFinalize();
        exit(1);
    }
}
/*********************************************************/

/* This fuction frees the allocate memory for a 3D array */
void Memory_free_3D_array(int NY, int NZ, void ***array) {

    int j, k;

    for (k=0; k<NZ; k++) {
        for (j=0; j<NY; j++) {

            free(array[k][j]);
        }
        free(array[k]);
    }
    free(array);
}
/*********************************************************/

/* This fuction frees the allocate memory for a 2D array */
void Memory_free_2D_array(int NZ, void **array) {

    int k;

    for (k=0; k<NZ; k++) {

        free(array[k]);
    }
    free(array);
}
/*********************************************************/

/* This function converts the 3D double array to 1D float array, i.e. data_3D[NZ][NY][NX] ---> data_1D[NX*NY*NZ]  */
void Memory_convert_3D_double_to_1D_float_array(double ***data_3D, float *data_1D, int NX, int NY, int NZ) {

    int i, j, k;
    int index = 0;

    for (k=0; k<NZ; k++) {
        for (j=0; j<NY; j++) {
            for (i=0; i<NX; i++) {

                /* Copy data from the 3D array into 1D-float array */
                data_1D[index] = (float)data_3D[k][j][i];
                index++;
            } /* for i*/
        } /* for j*/
    } /* for k*/
}
/*********************************************************/

/* This function converts the 2D double array to 1D float array, i.e. data_2D[NY][NX] ---> data_1D[NX*NY] */
void Memory_convert_2D_double_to_1D_float_array(double **data_2D, float *data_1D, int NX, int NY) {

    int i, j;
    int index = 0;

    for (j=0; j<NY; j++) {
        for (i=0; i<NX; i++) {

            /* Copy data from the 3D array into 1D-float array */
            data_1D[index] = (float)data_2D[j][i];
            index++;
        } /* for i*/
    } /* for j*/
}
/*********************************************************/

/* This function converts the 2D double array to 1D (float) or (double) array, i.e. data_2D[NY][NX] ---> data_1D[NX*NY] */
void Memory_convert_2D_double_to_1D(double **data_2D, void *data_1D, int NX, int NY, short int conv_type) {

    int i, j;
    int index = 0;

    if (conv_type == GVG_FLOAT) {

        float *out_1D = (float *)data_1D;
        for (j=0; j<NY; j++) {
            for (i=0; i<NX; i++) {

                /* Copy data from the 3D array into 1D-float array */
                out_1D[index] = (float)data_2D[j][i];
                index++;
            } /* for i*/
        } /* for j*/

    } else if (conv_type == GVG_DOUBLE) {


        double *out_1D = (double *)data_1D;
        for (j=0; j<NY; j++) {
            for (i=0; i<NX; i++) {

                /* Copy data from the 2D array into 1D array */
                out_1D[index] = data_2D[j][i];
                index++;
            } /* for i*/
        } /* for j*/
    } /* else */
}
/*********************************************************/


/* This function converts the 1D double array to 2D (double) array, i.e. data_2D[NY][NX] ---> data_1D[NX*NY] */
void Memory_convert_1D_double_to_2D(double *data_1D, double **data_2D, int NX, int NY) {

    int i, j;
    int index = 0;

    for (j=0; j<NY; j++) {
        for (i=0; i<NX; i++) {

            /* Copy data from the 3D array into 1D-float array */
            data_2D[j][i] = data_1D[index];
            index++;
        } /* for i*/
    } /* for j*/
}
/*********************************************************/


/* This function converts the 1D double array to 1D float array */
void Memory_convert_1D_double_to_1D_float_array(double *data_double, float *data_float, int N) {

    int i;

    for (i=0; i<N; i++) {

        data_float[i] = (float) data_double[i];
    }
}
/*********************************************************/

/* This function allocates a 3D array for double type:: array[NZ][NY][NX] */
double ***Memory_allocate_3D_double(int NX, int NY, int NZ) {

    int j, k;
    double ***array=NULL;

    array = (double ***)calloc(NZ, sizeof(double **));
    //PetscMalloc(NZ*sizeof(double **), &array);
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (double **)calloc(NY, sizeof(double *));
        //PetscMalloc(NY*sizeof(double*), &array[k]);
        Memory_check_allocation(array[k]);

        for (j=0; j<NY; j++) {

            array[k][j] = (double *)calloc(NX, sizeof(double));
            //PetscMalloc(NX*sizeof(double), &array[k][j]);
            Memory_check_allocation(array[k][j]);
        }
    }
    return (array);
}
/*********************************************************/

/* This function allocates a 3D array for float type:: array[NZ][NY][NX] */
float ***Memory_allocate_3D_float(int NX, int NY, int NZ) {

    int j, k;
    float ***array=NULL;

    array = (float ***)calloc(NZ, sizeof(float **));
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (float **)calloc(NY, sizeof(float *));
        Memory_check_allocation(array[k]);

        for (j=0; j<NY; j++) {

            array[k][j] = (float *)calloc(NX, sizeof(float));
            Memory_check_allocation(array[k][j]);
        }
    }
    return (array);
}
/*********************************************************/

/* This function allocates a 3D array for int type:: array[NZ][NY][NX] */
int ***Memory_allocate_3D_int(int NX, int NY, int NZ) {

    int j, k;
    int ***array=NULL;

    array = (int ***)calloc(NZ, sizeof(int **));
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (int **)calloc(NY, sizeof(int *));
        Memory_check_allocation(array[k]);

        for (j=0; j<NY; j++) {

            array[k][j] = (int *)calloc(NX, sizeof(int));
            Memory_check_allocation(array[k][j]);
        }
    }
    return (array);
}
/*********************************************************/

/* This function allocates a 2D array for double:::: array[NZ][NY] */
double **Memory_allocate_2D_double(int NY, int NZ) {

    int k;
    double **array=NULL;

    array = (double **)calloc(NZ, sizeof(double *));
    //PetscMalloc(NZ*sizeof(double *), &array);
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (double *)calloc(NY, sizeof(double));
        //PetscMalloc(NY*sizeof(double), &(array[k]));
        Memory_check_allocation(array[k]);
    }

    return (array);
}
/*********************************************************/

/* This function allocates a 2D array for float:::: array[NZ][NY] */
float **Memory_allocate_2D_float(int NY, int NZ) {

    int k;
    float **array=NULL;

    array = (float **)calloc(NZ, sizeof(float *));
    //PetscMalloc(NZ*sizeof(float *), &array);
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (float *)calloc(NY, sizeof(float));
        //PetscMalloc(NY*sizeof(float), &(array[k]));
        Memory_check_allocation(array[k]);

        /* set the array to zero */
    }
    return (array);
}
/*********************************************************/

/* This function allocates a 2D array for int:::: array[NZ][NY] */
int **Memory_allocate_2D_int(int NY, int NZ) {

    int k;
    int **array=NULL;

    array = (int **)calloc(NZ, sizeof(int *));
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (int *)calloc(NY, sizeof(int));
        Memory_check_allocation(array[k]);
    }

    return (array);
}
/*********************************************************/

/* This function allocates a 2D array for char:::: array[NZ][NY] */
char **Memory_allocate_2D_char(int NY, int NZ) {

    int k;
    char **array=NULL;

    array = (char **)calloc(NZ, sizeof(char *));
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        array[k] = (char *)calloc(NY, sizeof(char));
        Memory_check_allocation(array[k]);
    }

    return (array);
}
/*********************************************************/

/* This function allocates a 1D array for double:::: array[N] */
double *Memory_allocate_1D_double(int N) {

    double *array=NULL;

    array = (double *)calloc(N, sizeof(double));
    Memory_check_allocation(array);
    return (array);
}
/*********************************************************/

/* This function allocates a 1D array for float:::: array[N] */
float *Memory_allocate_1D_float(int N) {

    float *array=NULL;

    array = (float *)calloc(N, sizeof(float));
    Memory_check_allocation(array);
    return (array);
}
/*********************************************************/

/* This function allocates a 1D array for int:::: array[N] */
int *Memory_allocate_1D_int(int N) {

    int *array=NULL;

    array = (int *)calloc(N, sizeof(int));
    Memory_check_allocation(array);
    return (array);
}
/*********************************************************/

/* This function allocates a 1D array for int:::: array[N] */
char *Memory_allocate_1D_char(int N) {

    char *array=NULL;

    array = (char *)calloc(N, sizeof(char));
    Memory_check_allocation(array);
    return (array);
}
/*********************************************************/


/* Using PETSc to allocate 1D array for any type */
void Memory_allocate_1D(int type, int N, void **out) {

    size_t size=0;
    switch (type) {

    case GVG_DOUBLE: size = sizeof(double); break;
    case GVG_FLOAT: size = sizeof(float); break;
    case GVG_INT: size = sizeof(int); break;
    case GVG_CHAR: size = sizeof(char); break;
    }
    PetscMalloc(N*size, out);
}

/* This function allocates a 2D array for int:::: array[NZ][NY] */
int **Memory_allocate_2D_int_new(int NY, int NZ) {

    int k;
    int **array=NULL;

    PetscMalloc(NZ*sizeof(int *), &array);
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        PetscMalloc(NY*sizeof(int), &(array[k]));
        Memory_check_allocation(array[k]);
    }

    return (array);
}
/*********************************************************/

void Memory_allocate_2D(int type, int NY, int NZ, void ***array) {


    size_t size=0;
    switch (type) {

    case GVG_DOUBLE: size = sizeof(double); break;
    case GVG_FLOAT: size = sizeof(float); break;
    case GVG_INT: size = sizeof(int); break;
    case GVG_CHAR: size = sizeof(char); break;
    }

    int k;

    PetscMalloc(NZ*sizeof(void *), array);
    Memory_check_allocation(array);

    for (k=0; k<NZ; k++) {

        PetscMalloc(NY*size, (array[k]));
        //(*array)[k] = malloc(NY*size);
        //Memory_check_allocation(array[k]);
    }

}

