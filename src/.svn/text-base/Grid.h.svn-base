#ifndef GRID_H
 #define GRID_H

MAC_grid *Grid_create(Parameters *params);
void Grid_destroy(MAC_grid *grid);
void Grid_generate_mesh_coordinates(MAC_grid *grid, Parameters *params, char which_quantity);
void Grid_compute_metric_coefficients(MAC_grid *grid, char which_quantity);
void Grid_describe_interface(MAC_grid *grid, Parameters *params);
int Grid_get_c_status(MAC_grid *grid, int x_index, int y_index, int z_index);
int Grid_get_u_status(MAC_grid *grid, int x_index, int y_index, int z_index);
int Grid_get_v_status(MAC_grid *grid, int x_index, int y_index, int z_index);
int Grid_get_w_status(MAC_grid *grid, int x_index, int y_index, int z_index);
int Grid_import_bottom_interface(MAC_grid *grid, Parameters *params, char *filename);
void Grid_tag_q_nodes_using_surface(MAC_grid *grid, SurfaceType *surf, char which_quantity);
void Grid_tag_q_box_boundary_nodes(MAC_grid *grid, char which_quantity);
void Grid_set_u_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status);
void Grid_set_v_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status);
void Grid_set_w_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status);
void Grid_set_c_status(MAC_grid *grid, int x_index, int y_index, int z_index, int which_status);
void Grid_identify_geometry(MAC_grid *grid);

int Grid_get_x_index(double x, MAC_grid *grid, char which_quantity) ;
int Grid_get_y_index(double y, MAC_grid *grid, char which_quantity) ;
int Grid_get_z_index(double z, MAC_grid *grid, char which_quantity) ;

int Grid_import_grid_from_file(MAC_grid *grid, Parameters *params) ;
void Grid_set_interface_y_index(MAC_grid *grid) ;
void Grid_tag_q_buffer_solid_nodes(MAC_grid *grid, int width, char which_quantity) ;

short int Grid_is_local(int i, int j, int k, MAC_grid *grid) ;


#endif
