#include "definitions.h"
#include "DataTypes.h"
#include "Grid.h"
#include "MyMath.h"
#include "gvg.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



/* This function interpolates any quantity based on 2 or 4 neighboring nodes's given values and coordinates */
double MyMath_interpolate_quantity(Scalar *q, Scalar q_int, int N) {

    double q0, q1, q2, q3;
    double a_Dx, a_Dy;
    double Wp0, Wp1;
    double Wp2;
    double m;
    double q_top, q_bottom;
    double tol = 1.0e-8;


    /* NOTE: If the # of points is equal to two, just a linear interpoleation in X1 direction */
    if (N == 2) {

        q0 = q[0].Value;
        q1 = q[1].Value;

        m = (q1 - q0) / (q[1].X1 - q[0].X1);

        /* check if m!= inf */
        if (fabs(q[1].X1 - q[0].X1) < tol ) {

            q_int.Value = q[0].Value; /* constant value everywhere */

        } else {

            /* Linear interpolation: y = m*(x-x1) + y1 */
            q_int.Value = m*(q_int.X1 - q[0].X1) + q0;
        } /* else */

        return (q_int.Value);
    } /* if N =2 */

    /*
  2---------3
         |         |
         |         |
         |         |
  0--------1
     */

    if (N == 4) {

        q0 = q[0].Value;
        q1 = q[1].Value;
        q2 = q[2].Value;
        q3 = q[3].Value;

        a_Dx = 1.0/(q[1].X1 - q[0].X1);
        Wp0  = fabs( (q_int.X1 - q[0].X1) * a_Dx);
        Wp1  = 1.0 - Wp0;

        q_top       = Wp0*q2 + Wp1*q3;
        q_bottom    = Wp0*q0 + Wp1*q1;

        a_Dy = 1.0/(q[2].X2 - q[0].X2);
        Wp0  = fabs( (q_int.X2 - q[0].X2) * a_Dy);
        Wp2  = 1.0 - Wp0;

        q_int.Value = Wp0*q_bottom + Wp2*q_top;
        return (q_int.Value);
    }

    return (-1.0);
}
/**************************************************************************************************************/

/* Retruns a random number between [0, max) */
double MyMath_random(double max) {

    double rnd;

    rnd = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
    return (max*rnd);
}
/**************************************************************************************************************/

/* This function finds (x,y,z) index based on the global index on the node */
Indices MyMath_transform_index(int p_index, int NX, int NY, int NZ) {

    Indices Node_index;

    Node_index.x_index = p_index % NX;
    Node_index.y_index = ( (p_index - Node_index.x_index) % (NZ*NY) ) / NX;
    Node_index.z_index = ( p_index - Node_index.x_index - Node_index.y_index*NX ) / (NX*NY);

    return Node_index;

}
/**************************************************************************************************************/

/* This function returns the distance between the two points */
double MyMath_get_point_point_distance(PointType *p1, PointType *p2) {

    double d2;
    d2 = (p1->x - p2->x)*(p1->x - p2->x) + (p1->y - p2->y)*(p1->y - p2->y) + (p1->z - p2->z)*(p1->z - p2->z);

    return (sqrt(d2));
}
/***************************************************************************************************/

/* This function inverses a matrix using Gauss-Jordan method */
/* Input matrix is in mat->A[][] (size by size) */
/* Inverse matrix will be stored in mat->A_inv[][] (size by size) */
void MyMath_GaussJordan_matrix_inv(MatrixType *mat) {

    int *index_col, *index_row, *index_piv;
    int i, j, k, h, hh;
    double max, temp, a_piv;
    int in_row=-1;
    int in_col=-1;
    double _temp;
    double tol = 1.0e-8;
    int n = mat->size;

    index_col = (int *)calloc(n, sizeof(int));
    index_row = (int *)calloc(n, sizeof(int));
    index_piv = (int *)calloc(n, sizeof(int));

    /* copy the values */
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            mat->A_inv[j][i] = mat->A[j][i];
        } /* for i*/
    } /* for j*/

    /* pivot index */
    for (j=0; j<n; j++) {
        index_piv[j] = -1;
    } /* for j */

    for (i=0; i<n; i++) {
        max=0.0;

        for (j=0; j<n; j++) {
            if (index_piv[j] != 0) {
                for (k=0; k<n; k++) {
                    if (index_piv[k] == -1) {
                        if (fabs(mat->A_inv[j][k]) >= max) {
                            max = fabs(mat->A_inv[j][k]);
                            in_row=j;
                            in_col=k;
                        } /* if */
                    } /* if index_piv */
                } /* for k */
            } /* if index_piv */
        }/* for j */
        ++(index_piv[in_col]);

        if (in_row != in_col) {

            for (h=0; h<n; h++) { /* swap */

                _temp = mat->A_inv[in_row][h];
                mat->A_inv[in_row][h] = mat->A_inv[in_col][h];
                mat->A_inv[in_col][h] = _temp;
            } /* for s */
        } /* if in_row */

        index_row[i] = in_row;
        index_col[i] = in_col;
        if ( fabs(mat->A_inv[in_col][in_col]) <= tol ) {

            printf("MyMath.c/ Error inverting the matrix. Singular matrix...\n");
        } /* if */

        a_piv = 1.0/mat->A_inv[in_col][in_col];
        mat->A_inv[in_col][in_col] = 1.0;

        for (h=0; h<n; h++) mat->A_inv[in_col][h] *= a_piv;

        for (hh=0; hh<n; hh++) {

            if (hh != in_col) {
                temp = mat->A_inv[hh][in_col];
                mat->A_inv[hh][in_col]=0.0;
                for (h=0; h<n; h++) mat->A_inv[hh][h] -= mat->A_inv[in_col][h] * temp;
            }
        } /*for hh */

    } /* for i */

    for (h=n-1; h>=0; h--) {
        if (index_row[h] != index_col[h])
            for (k=0; k<n; k++) { /* swap */

                _temp = mat->A_inv[k][index_row[h]];
                mat->A_inv[k][index_row[h]] = mat->A_inv[k][index_col[h]];
                mat->A_inv[k][index_col[h]] = _temp;

            } /* for k */
    } /* for h */

    free(index_col);
    free(index_row);
    free(index_piv);
}
/***************************************************************************************************/

/* This function does a trilinear interpolation for any given 3D data array and an intepolation point at (x_int,y_int,z_int) */
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
#define __func__ "Mymath_do_trilinear()"
double MyMath_do_trilinear_inter(double x_int, double y_int, double z_int, double ***data, MAC_grid *grid, char which_quantity) {

    double x0, y0, z0;
    double x7, y7, z7;
    int i0, j0, k0;
    int m;
    double alpha, beta, gamma;
    double wc[8];
    static short int dir_x[8] = {0, 1, 0, 1, 0, 1, 0, 1};
    static short int dir_y[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    static short int dir_z[8] = {0, 0, 0, 0, 1, 1, 1, 1};

    /* First the grid coordinates to be used for interpolation */
    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    switch (which_quantity) {

    case 'u':

        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        break;

    case 'v':

        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        break;

    case 'w':

        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        break;

    case 'c':

        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        break;

    case 'p':

        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        break;
    default:

        PetscPrintf(PCW, "MyMath.c/ Error setting coordinates. Unknown quantity %c\n", which_quantity);
    } /* switch */

    /* Now, get the index of the nodes to the left and bottom of the interpolation point */
    i0  = Grid_get_x_index(x_int, grid, which_quantity);
    j0  = Grid_get_y_index(y_int, grid, which_quantity);
    k0  = Grid_get_z_index(z_int, grid, which_quantity);

	
    if ( (i0 == -1) || (j0 == -1) || (k0 == -1) ){

        int rank;
        MPI_Comm_rank(PCW, &rank);
        //GVG_printf("Warning, index -1 for interpolation point (%f,%f,%f). Using extrapolation instead.\n", rank, x_int, y_int, z_int);
        //printf("MyMath.c/ [%d] Warning, index -1 for interpolation point (%f,%f,%f). Using extrapolation instead.\n", rank, x_int, y_int, z_int);
        //PetscPrintf(PCW, "MyMath.c/ Rank:%d Warning, index -1 for interpolation point (%f,%f,%f). Using extrapolation instead.\n", rank, x_int, y_int, z_int);
        i0 = max(i0, 0);
        j0 = max(j0, 0);
        k0 = max(k0, 0);
        //getchar();

    } /* if */
    if ( (i0 == GVG_INF) || (j0 == GVG_INF) || (k0 == GVG_INF) ) {

        int rank;
        MPI_Comm_rank(PCW, &rank);
        //GVG_printf("MyMath.c/ Rank:%d Warning, index INF for interpolation point (%f,%f,%f). Using extrapolation instead.\n", rank, x_int, y_int, z_int);
        //printf("MyMath.c/ [%d] Warning, index INF for interpolation point (%f,%f,%f). Using extrapolation instead.\n", rank, x_int, y_int, z_int);
        //PetscPrintf(PCW, "MyMath.c/ Rank:%d Warning, index INF for interpolation point (%f,%f,%f). Using extrapolation instead.\n", rank, x_int, y_int, z_int);
        i0 = min(i0, grid->NI-2);
        j0 = min(j0, grid->NJ-2);
        k0 = min(k0, grid->NK-2);
    } /* if */

    /* coordinates of the the two corner nodes, i.e. F0 and F7 */
    x0 = xq[i0];
    y0 = yq[j0];
    z0 = zq[k0];

    /* F7 */
    x7 = xq[i0+1];
    y7 = yq[j0+1];
    z7 = zq[k0+1];

    if ( ( !Grid_is_local(i0, j0, k0, grid) ) || (!Grid_is_local(i0+1, j0+1, k0+1, grid)) ){

        //printf(PCW, "MyMath.c/ Warning. The given point (%f,%f,%f) lies outside the range of current processor. Cannot perform interpolation.\n", x_int, y_int, z_int);
        return GVG_INF;
    }
    /* Trilinear interpolation coefficients */
    alpha = (x7 - x_int)/(x7 - x0);
    beta  = (y7 - y_int)/(y7 - y0);
    gamma = (z7 - z_int)/(z7 - z0);


    /* Get the interpolation coefficients to find the value of the image node based on the neigboring fluid nodes */
    wc[0] = alpha * beta * gamma;
    wc[1] = (1.0-alpha) * beta * gamma;
    wc[2] = alpha * (1.0-beta) * gamma;
    wc[3] = (1.0-alpha) * (1.0-beta) * gamma;
    wc[4] = alpha * beta * (1.0-gamma);
    wc[5] = (1.0-alpha) * beta * (1.0-gamma);
    wc[6] = alpha * (1.0-beta) * (1.0-gamma);
    wc[7] = (1.0-alpha)*(1.0-beta)*(1.0-gamma);

    double val_int = 0.0;
    /* Go through all 8 nodes and add the contribution to the interpolated value */
    for (m=0; m<8; m++) {

        val_int += wc[m]*data[ k0+dir_z[m] ][ j0+dir_y[m] ][ i0+dir_x[m] ];
    } /* for m */

    return (val_int);
} 
/***************************************************************************************************/

/* This function uses binary search algorithm to locate a value in a sorted data array */
/* input:
      data: sorted data (ascending order)
      search_val: value to be searched in the data array
      min: index of the lower bound in the array (use 0 if not sure)
      max: index of the upper bound in the array (use size-1 if not sure)
   output:
      the index of the node closest to the search value (to the LEFT!)
      e.g.
        data[0] = 0.7
        data[1] = 2.5
        data[2] = 3.6

        search_val = 2.7

        MyMath_binarysearch(data, 2.7, 0, 2) will return 1;
      if the value is out of range it will return -1;
Inspired from:
http://orcik.net/programming/linear-and-binary-search/
*/
int MyMath_binary_search(double *data, double search_val, int min, int max) {

    while (min <= max) {

        int mid = (min + max) / 2;
        if (search_val >= data[mid]) {

            if (search_val < data[mid+1]) return mid; /* index found */
            min = mid + 1;  // repeat search in top half.

        } else if (search_val < data[mid])
            max = mid - 1;   // repeat search in bottom half.
    } /* while */
    return -1;    // failed, search_val out of data range
}
/***************************************************************************************************/

double MyMath_smooth_intp(Scalar *in, Scalar out, int N) {

    double h;
    double w_n;
    int i;
    double R_max=0.0, R;
    double q=0.0;
    PointType *p_n, *p_int;

    out.Value = 0.0;

    p_n   = (PointType *)calloc(N, sizeof(PointType));
    p_int = (PointType *)calloc(1, sizeof(PointType));

    /* interpolation point */
    p_int->x = out.X1;
    p_int->y = out.X2;

    /* R: most distance from the interpolation point */
    for (i=0; i<N; i++) {

        p_n->x = in[i].X1;
        p_n->y = in[i].X2;

        R = MyMath_get_point_point_distance(p_int, p_n);

        R_max = max(R, R_max);
    } /* for */

    R = 1.05*R_max;
    /* compute q */
    for (i=0; i<N; i++) {

        p_n->x = in[i].X1;
        p_n->y = in[i].X2;

        h = MyMath_get_point_point_distance(p_int, p_n);
        q += (R - h)*(R - h)/(R*R*h*h);
    } /* for i*/

    /* compute the interpolation value */
    for (i=0; i<N; i++) {

        p_n->x = in[i].X1;
        p_n->y = in[i].X2;

        h = MyMath_get_point_point_distance(p_int, p_n);
        w_n = (R - h)*(R - h)/(R*R*h*h) ;
        out.Value += w_n*in[i].Value/q;
        printf("MyMath.c/ i:%d R:%f h:%f w:%f q:%f\n", i, R, h, w_n, q);
        //getchar();
    } /* for i*/

    return out.Value;
}
/***************************************************************************************************/
