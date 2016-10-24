#ifndef CONC_H
 #define CONC_H

short int Conc_is_solid(Concentration *c, int i, int j, int k) ;
void Conc_modify_diagonal(Concentration *c, MAC_grid *grid, double dt, double dt_old) ;
void Conc_set_lsys_index( Concentration *c, int i, int j, int k, int lsys_index ) ;
void Conc_set_solid_flag(Concentration *c, int i, int j, int k) ;
void Conc_setup_lsys_accounting_geometry(Concentration *c, MAC_grid *grid, Parameters *params) ;
void Conc_destroy(Concentration *c, Parameters *params, int conc_index) ;
void Conc_identify_geometry(Concentration *c, MAC_grid *grid) ;
void Conc_update_boundary( Concentration *c, MAC_grid *grid, Parameters *params) ;
void Conc_form_LHS_matrix( Concentration *c, MAC_grid *grid,  Parameters *params) ;
double Conc_get_solid_value(Concentration *c, MAC_grid *grid, Parameters *params, int column, int row, char which_neighbor) ;
void Conc_ENO_ghost_cells(ENO_Scheme *ENO, Concentration *c, MAC_grid *grid, Parameters *params, Indices
start_cell, Indices end_cell, char which_direction) ;
void Conc_check_valid(Concentration *c);
void Conc_update_erosion_flux(Concentration *c, MAC_grid *grid, Parameters *params) ;
void Conc_store_old_outflow(Concentration *c) ;
void Conc_rk_average(Concentration *c, double old_frac, double new_frac) ;
void Conc_rk_average_outflow(Concentration *c, double old_frac, double new_frac) ;
void Conc_initialize(Concentration *c, MAC_grid *grid, Parameters *params) ;



Concentration *Conc_create(MAC_grid *grid, Parameters *params, int conc_index) ;
void Conc_set_particle_v_velocity(Concentration *c, Velocity *v, Parameters *params) ;
void Conc_compute_particle_settling_speed(Concentration *c, MAC_grid *grid) ;
double Conc_settling_speed_function(double conc, double phi_max, double V_s0) ;
int Conc_solve(Concentration *c) ;
void Conc_set_RHS(Concentration *c, MAC_grid *grid, Parameters *params, double dt) ;
void Conc_store_old_data(Concentration *c) ;

void Conc_compute_total_suspended_mass(Concentration *c, MAC_grid *grid) ;
void Conc_compute_ave_height_x(Concentration **c, MAC_grid *grid, Parameters *params) ;
double Conc_find_front_location(Concentration **c, MAC_grid *grid, Parameters *params, double front_limit) ;
void Conc_integrate_deposited_height(Concentration *c, MAC_grid *grid, double weight_factor, double dt) ;
void Conc_shrink_deposited_height(Concentration **c, MAC_grid *grid, Parameters *params) ;

void Conc_compute_stokes_dissipation_rate(Concentration *c, MAC_grid *grid) ;
void Conc_update_world_stokes_dissipation_rate( Concentration **c, MAC_grid *grid, Parameters *params) ;

void Conc_update_world_deposited_height(Concentration **c, MAC_grid *grid, Parameters *params) ;
void Conc_compute_total_concentration(Concentration **c, Parameters *params) ;

void Conc_update_world_deposited_height_dumped(Concentration **c, MAC_grid *grid, Parameters *params) ;

void Conc_dump_particles(Concentration **c, MAC_grid *grid, Parameters *params) ;

void Conc_reset_nonlocal_deposit_height(Concentration *c, MAC_grid *grid) ;
void Conc_compute_potential_energies(Concentration *c, Velocity *v, MAC_grid *grid, Parameters *params) ;
void Conc_compute_sedimentation_rate(Concentration *c, MAC_grid *grid, Parameters *params) ;

void Conc_bottom_sed_rate(Concentration *c, MAC_grid *grid, Parameters *params) ;


#endif

