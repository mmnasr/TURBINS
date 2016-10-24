#ifndef MYMATH_H
 #define MYMATH_H

double MyMath_interpolate_quantity(Scalar *q, Scalar q_int, int N);
double MyMath_random(double max);
Indices MyMath_transform_index(int p_index, int NX, int NY, int NZ);

double MyMath_get_point_point_distance(PointType *p1, PointType *p2) ;
void MyMath_GaussJordan_matrix_inv(MatrixType *mat) ;
void MyMath_matrix_inv(MatrixType *mat) ;
double MyMath_matrix_determ(MatrixType *mat) ;

double MyMath_do_trilinear_inter(double x_int, double y_int, double z_int, double ***data, MAC_grid *grid, char which_quantity) ;
int MyMath_binary_search(double *data, double search_val, int min, int max) ;


#endif
