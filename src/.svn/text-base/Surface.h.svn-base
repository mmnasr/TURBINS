#ifndef H_SURFACE
	#define H_SURFACE
	
SurfaceType *Surface_create(MAC_grid *grid);
void Surface_destroy(SurfaceType *surf, MAC_grid *grid)  ;
double Surface_find_Dx_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity);
void Surface_initialize(SurfaceType *surf, MAC_grid *grid, Parameters *params);
void Surface_find_exact_sdf(SurfaceType *surf, MAC_grid *grid);
void Surface_find_q_sdf(SurfaceType *surf, MAC_grid *grid);
double Surface_find_Dx_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity);
double Surface_find_Dy_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity);
double Surface_find_Dz_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity);
double Surface_find_Dn_to_surface(double ***sdf, MAC_grid *grid, int i, int j, int k, VectorType *n, char which_quantity);

VectorType Surface_compute_normal(double ***sdf, MAC_grid *grid, int i, int j, int k, char which_quantity) ;


#endif
