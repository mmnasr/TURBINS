#ifndef H_IMMERSED
	#define H_IMMERSED 
	
Immersed *Immersed_create(MAC_grid *grid, char which_quantity);
void Immersed_destroy(Immersed *q_immersed, MAC_grid *grid) ;
int Immersed_count_immersed_nodes(MAC_grid *grid, char which_quantity) ;
void Immersed_compute_interpolation_coef(ImmersedNode *ib_node, char which_quantity);
void Immersed_set_q_immersed_indices(MAC_grid *grid, char which_quantity);
void Immersed_set_ib_global_index(Immersed *q_immersed, int x_index, int y_index, int z_index, int global_index);
int Immersed_get_ib_global_index(Immersed *q_immersed, int x_index, int y_index, int z_index);
void Immersed_set_ib_xyz_index(ImmersedNode *ib_node, int x_index, int y_index, int z_index);
void Immersed_set_q_immersed_coordinates(MAC_grid *grid, char which_quantity);
void Immersed_setup_q_immersed_nodes(MAC_grid *grid, char which_quantity);
void Immersed_set_q_immersed_neighboring_fluid_nodes(MAC_grid *grid, char which_quantity);
void Immersed_set_q_interpolation_coef(MAC_grid *grid, char which_quantity);
void Immersed_assign_q_neighboring_fluid_nodes(ImmersedNode *ib_node, MAC_grid *grid, char which_quantity);
void Immersed_set_q_immersed_boundary_condition(MAC_grid *grid, char which_quantity);
ImmersedNode *Immersed_get_ib_node(Immersed *q_immersed, int i);
void Immersed_set_q_immersed_segments(MAC_grid *grid, Parameters *params, char which_quantity);
void Immersed_testing_the_control_points(MAC_grid *grid, Parameters *params);
PointType *Immersed_choose_best_control_point(ImmersedNode *ib_node, MAC_grid *grid, Parameters *params);
void Immersed_write_q_immersed_info(MAC_grid *grid, Parameters *params, char which_quantity);
void OLDImmersed_write_q_immersed_special_cases(MAC_grid *grid, Parameters *params, char which_quantity);
PointType *Immersed_find_control_point(ImmersedNode *ib_node, MAC_grid *grid, char which_quantity);
void Immersed_set_q_immersed_control_points(MAC_grid *grid, char which_quantity);

PointType *Immersed_find_image_point(ImmersedNode *ib_node, MAC_grid *grid, char which_quantity);
void Immersed_set_q_immersed_image_points(MAC_grid *grid, char which_quantity) ;
void Immersed_find_surface_normals(MAC_grid *grid, char which_quantity) ;

Indices Immersed_find_box(ImmersedNode *ib_node, PointType *p, MAC_grid *grid, char which_quantity) ;
void Immersed_is_valid(ImmersedNode *ib_node) ;

#endif
