#ifndef DEBUGGER_H
 #define DEBUGGER_H

void Debugger_set_velocity_data(Velocity *vel, MAC_grid *grid) ;
void Debugger_set_velocity_RHS(Velocity *vel, MAC_grid *grid) ;
void Debugger_set_pressure_RHS(Pressure *p, MAC_grid *grid) ;
void Debugger_solve_p_Poisson_equation(Pressure *p, MAC_grid *grid) ;
void Debugger_solve_vel_Poisson_equation(Velocity *vel, MAC_grid *grid) ;
int Debugger_check_nan(Vec q, MAC_grid *grid, Parameters *params, char *name) ;


void Debugger_set_q_rhs(Velocity *vel, Concentration *c, Pressure *p, MAC_grid *grid, Parameters *params, char which_quantity) ;
double Debugger_get_analytical_rhs(double x, double y, double z, Parameters *params, short int boundary) ;
double Debugger_get_analytical_sol(double x, double y, double z, MAC_grid *grid, Parameters *params, short int boundary) ;
double Debugger_get_analytical_neumann(double x, double y, double z, Parameters *params, short int boundary) ;
void Debugger_validate_q(Velocity *vel, Concentration *c, Pressure *p, MAC_grid *grid, Parameters *params, char which_quantity) ;

short int Debugger_check_q_status(int i, int j, int k, MAC_grid *grid, Parameters *params, char which_quantity) ;

void Debugger_q_get_max(Vec q, MAC_grid *grid, Parameters *params, char *name) ;


#endif 
