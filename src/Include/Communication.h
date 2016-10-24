#ifndef COMMUNICATION_H
 #define COMMUNICATION_H

void Communication_reduce_2D_arrays(double **G_send_array, double **W_recv_array, Indices *G_s, Indices *G_e, Indices *W_e, int status, Parameters *params) ;
void Communication_update_ghost_nodes(DA *DA_array, Vec *global, Vec *local, char status) ;
int Communication_create_DA3D(int npx, int npy, int npz, DAStencilType stencil_type, int ghost_node, Parameters *params, DA *DA_new) ;
void Communication_analyze_comm_area(DA DA_3D);


#endif
