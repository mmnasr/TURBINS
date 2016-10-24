#ifndef SOLVER_H
 #define SOLVER_H

Mat Solver_create_LHS_matrix(int size);
Vec Solver_create_vector(int size);
KSP Solver_get_concentration_solver(Mat lhs, Parameters *params);
KSP Solver_get_pressure_solver(Mat lhs, Parameters *params);
KSP get_streamfunction_solver(Mat lhs, Parameters *params);
KSP Solver_get_velocity_solver(Mat lhs, Parameters *params);
short int Solver_is_converged(KSP solver, char *lsys_name) ;

#endif
