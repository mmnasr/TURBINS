#ifndef INFLOW_H
 #define INFLOW_H

void Inflow_u_velocity_profile (Velocity *u, MAC_grid *grid, Parameters *params, double time) ;
void Inflow_conc_profile(Concentration *c, MAC_grid *grid, Parameters *params, double time);

#endif
