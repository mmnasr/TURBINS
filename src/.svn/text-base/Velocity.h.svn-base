#ifndef VELOCITY_H
 #define VELOCITY_H

Velocity *Velocity_create(MAC_grid *grid, Parameters *params, char which_velocity) ;
void Velocity_destroy(Velocity *vel) ;
void Velocity_store_old_data(Velocity *vel) ;

int Velocity_get_lsys_index (Velocity *vel, int i, int j, int k, char which) ;

void Velocity_u_form_LHS_matrix(Velocity *u, MAC_grid *grid, Parameters *params) ;
void Velocity_v_form_LHS_matrix(Velocity *v, MAC_grid *grid, Parameters *params) ;
void Velocity_w_form_LHS_matrix(Velocity *w, MAC_grid *grid, Parameters *params) ;

void Velocity_cell_center( Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) ;
void Velocity_u_set_RHS( Velocity *u, Pressure *p, MAC_grid *grid, Parameters *params, double dt);
void Velocity_v_set_RHS( Velocity *v, Pressure *p, Concentration **c, MAC_grid *grid, Parameters *params, double dt);
void Velocity_w_set_RHS( Velocity *w, Pressure *p, MAC_grid *grid, Parameters *params, double dt);

void Velocity_rk_average(Velocity *vel, double old_frac, double new_frac) ;

void Velocity_u_ENO_ghost_cells(ENO_Scheme *ENO, Velocity *u, MAC_grid *grid, Parameters *params, Indices
start_cell, Indices end_cell, char which_direction) ;
void Velocity_v_ENO_ghost_cells(ENO_Scheme *ENO, Velocity *v, MAC_grid *grid, Parameters *params, Indices
start_cell, Indices end_cell, char which_direction) ;
void Velocity_w_ENO_ghost_cells(ENO_Scheme *ENO, Velocity *w, MAC_grid *grid, Parameters *params, Indices
start_cell, Indices end_cell, char which_direction) ;

void Velocity_display(Velocity *vel, MAC_grid *grid, Parameters *params) ;

void Velocity_modify_diagonal(Velocity *vel, MAC_grid *grid, Parameters *params, double dt, double dt_old); 

void Velocity_nonzero_initialize(Velocity *vel, MAC_grid *grid, Parameters *params); 
void Velocity_debug_flow(Velocity *vel, MAC_grid *grid, Parameters *params) ;

void Velocity_store_old_outflow(Velocity *vel) ;
int Velocity_solve(Velocity *vel) ;

void Velocity_setup_lsys_accounting_geometry(Velocity *vel, MAC_grid *grid, Parameters *params) ;
void Velocity_rk_average_outflow(Velocity *vel, double old_frac, double new_frac) ;

void Velocity_compute_vorticity(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, char which_component, Vec *G_vort) ;

void Velocity_transpose_velocities(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, char which_quantity) ;


#endif
