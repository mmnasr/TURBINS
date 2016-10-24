#ifndef OUTFLOW_H
 #define OUTFLOW_H

void Outflow_impose_convective_boundary( Velocity *u, Velocity *v, Velocity *w, Concentration **c, MAC_grid *grid, Parameters *params, double dt) ;
void Outflow_vel_impose_convective_boundary( Velocity *vel, MAC_grid *grid, Parameters *params, double dt) ;
void Outflow_conc_impose_convective_boundary(Concentration *c, MAC_grid *grid, Parameters *params, double dt) ;
void Outflow_update_u_velocity_to_conserve_mass(Velocity *u, MAC_grid *grid) ;

void Outflow_vel_copy_data_into_last_plane(Velocity *vel, MAC_grid *grid, int should_set) ;
void Outflow_concentration_copy_data_into_last_plane(Concentration *c, MAC_grid *grid, int should_set) ;
void Outflow_vel_copy_data_into_2D_array(Velocity *vel, MAC_grid *grid) ;
void Outflow_concentration_copy_data_into_2D_array(Concentration *c, MAC_grid *grid) ;


#endif
