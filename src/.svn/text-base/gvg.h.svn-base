#ifndef GVG_H
 #define GVG_H

void GVG_identify_geometry(MAC_grid *grid, Parameters *params) ;
void GVG_setup_lsys_accounting_geometry(Velocity *u, Velocity *v, Velocity *w, Pressure *p, Concentration **c, MAC_grid
*grid, Parameters *params) ;
void GVG_initialize_primitive_data(Velocity *u, Velocity *v, Velocity *w, Concentration **c, MAC_grid *grid,
Parameters *params) ;
double GVG_new_time_step(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params, double current_time, double next_output_time, double dt_inp) ;
double GVG_compute_dt_applying_cfl_condition(Velocity *u, Velocity *v, Velocity *w, MAC_grid *grid, Parameters *params) ;
void GVG_store_flow_variables_old_data(Velocity *u, Velocity *v, Velocity *w, Concentration **c, Parameters *params) ;
int GVG_tvd_rk(GVG_bag *data_bag, Resume *resume) ;

void GVG_solve_Poisson_equation(Pressure *p, MAC_grid *grid, Parameters *params) ;
//GVG_bag *GVG_create_bag(Velocity *u, Velocity *v, Velocity *w, Pressure *p, Concentration **c, MAC_grid *grid, Parameters *params);
GVG_bag *GVG_create_bag(Parameters *params);


void GVG_integrate_all_the_equations_in_time(Velocity *u, Velocity *v, Velocity *w, Pressure *p, Concentration **c, MAC_grid
*grid, ENO_Scheme *ENO, Parameters *params, double dt, double dt_old, short int which_stage) ;

void GVG_write_workspace_data(GVG_bag *data_bag, Resume *resume) ;
void GVG_rk_average_flow_variables(Velocity *u, Velocity *v, Velocity *w, Concentration **c, Parameters *params, double old_frac, double new_frac) ;
void GVG_read_previous_simulation_data(Resume *resume, GVG_bag *data_bag) ;

void GVG_compute_convective_terms(ENO_Scheme *ENO, Velocity *u, Velocity *v, Velocity *w, Concentration *c, MAC_grid *grid, Parameters *params, char which_quantity);
void GVG_printf(const char *format, ...) ;

void GVG_set_value(Vec vec_data, MAC_grid *grid) ;
void GVG_announce_times(double time) ;

void GVG_write_PETSC_object_to_file(void *obj, char *filename, short int what) ;


#endif
