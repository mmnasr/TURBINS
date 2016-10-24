#ifndef INPUT_H
 #define INPUT_H

Parameters *Input_set_parameters(void);
void Input_destroy_parameters(Parameters *params);
void Input_read_parameters_from_file(Parameters *params);
void Input_set_internal_parameters(Parameters *params);

void Input_read_sub_section(FILE *file_in, double *read_params, int N) ;
void Input_read_sub_section_flag(FILE *file_in, int *flags, int N) ;

void Input_set_sphere_params(Parameters *params, MAC_grid *grid) ;

#endif
