#ifndef MEMORY_H
 #define MEMORY_H

void Memory_check_allocation(void *var);
void Memory_free_3D_array(int NY, int NZ, void ***array) ;
void Memory_free_2D_array(int NZ, void **array) ;
void Memory_convert_3D_double_to_1D_float_array(double ***data_3D, float *data_1D, int NX, int NY, int NZ) ;
void Memory_convert_2D_double_to_1D_float_array(double **data_2D, float *data_1D, int NX, int NY) ;
void Memory_convert_1D_double_to_1D_float_array(double *data_double, float *data_float, int N) ;


double ***Memory_allocate_3D_double(int NX, int NY, int NZ) ;
float ***Memory_allocate_3D_float(int NX, int NY, int NZ) ;
int ***Memory_allocate_3D_int(int NX, int NY, int NZ) ;
double **Memory_allocate_2D_double(int NY, int NZ) ;
float **Memory_allocate_2D_float(int NY, int NZ) ;
int **Memory_allocate_2D_int(int NY, int NZ) ;
char **Memory_allocate_2D_char(int NY, int NZ) ;
double *Memory_allocate_1D_double(int N) ;
float *Memory_allocate_1D_float(int N) ;
int *Memory_allocate_1D_int(int N) ;
char *Memory_allocate_1D_char(int N) ;

void Memory_allocate_2D(int type, int NY, int NZ, void ***array);
void Memory_allocate_1D(int type, int N, void **out) ;

int **Memory_allocate_2D_int_new(int NY, int NZ) ;


void Memory_convert_2D_double_to_1D(double **data_2D, void *data_1D, int NX, int NY, short int conv_type) ;
void Memory_convert_1D_double_to_2D(double *data_1D, double **data_2D, int NX, int NY) ;


#endif
