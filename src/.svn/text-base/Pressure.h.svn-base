#ifndef PRESSURE_H
 #define PRESSURE_H

Pressure *Pressure_create( MAC_grid *grid, Parameters *params) ;
void Pressure_destroy(Pressure *p) ;
void Pressure_setup_lsys_accounting_geometry(Pressure *p, MAC_grid *grid, Parameters *params) ;
int Pressure_get_lsys_index (Pressure *p, int i, int j, int k, char which) ;
void Pressure_form_LHS_matrix(Pressure *p, MAC_grid *grid) ;
void Pressure_set_RHS(Pressure *p, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, double dt) ;
void Pressure_project_velocity(Pressure *p, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, double dt) ;
int Pressure_solve(Pressure *p) ;
double Pressure_compute_velocity_divergence(Pressure *p, Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) ;
void Pressure_compute_real_pressure(Pressure *p, MAC_grid *grid) ;

#endif
