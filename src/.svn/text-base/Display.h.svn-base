#ifndef DISPLAY_H
 #define DISPLAY_H

void Display_3D_variable(double ***var, int NX, int NY, int NZ, const char *name);
void Display_3D_variable_point_based(double ***var, int NX, int NY, int NZ, const char *name);
void Display_DA_3D_info(MAC_grid *grid, Parameters *params) ;

void Display_2D_outflow(double **outflow, MAC_grid *grid, Parameters *params, const char *q_name) ;
void Display_2D_inflow(double **inflow, MAC_grid *grid, Parameters *params, const char *q_name) ;


void Display_parameters(Parameters *params) ;
void Display_flag(int flag, const char *name) ;

void Display_point(PointType *p, const char *name) ;
void Display_immersed_node(ImmersedNode *ib_node) ;

void Display_grid(double *data, int N, char *name) ;

void Display_DA_3D_L_data(Vec q, MAC_grid *grid, Parameters *params, char *q_name, char which_quantity) ;

void Display_matrix_info(Mat A, const char *name) ;

void Display_memory_info(void) ;


#endif 
