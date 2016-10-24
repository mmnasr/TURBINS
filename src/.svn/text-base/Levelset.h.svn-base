#ifndef H_LEVELSET
	#define H_LEVELSET

/* This function allocates memory for the levelset function */
LevelsetType *Levelset_create(MAC_grid *grid);
void Levelset_destroy(LevelsetType *phi);
void Levelset_reinitialize(LevelsetType *phi, MAC_grid *grid, int max_iterations);
void Levelset_compute_derivatives(LevelsetType *phi, MAC_grid *grid);
void Levelset_do_Godunov(LevelsetType *phi, MAC_grid *grid, double dt);
int Levelset_is_converged(LevelsetType *phi, MAC_grid *grid, int which_error);
void Levelset_store_old_data(LevelsetType *phi); 
double Levelset_compute_fictitious_time(MAC_grid *grid) ;

void Levelset_update_last_planes(LevelsetType *phi, MAC_grid *grid) ;

#endif
