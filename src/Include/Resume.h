#ifndef H_RESUME 
	#define H_RESUME

Resume *Resume_create(int vec_num) ;
void Resume_destroy(Resume *resume) ;
void Resume_set_write_flag(Resume *resume) ;
void Resume_set_read_flag(Resume *resume) ;
void Resume_set_IO_grid_flag(Resume *resume) ;
void Resume_reset_write_flag(Resume *resume) ;
void Resume_reset_read_flag(Resume *resume) ;
void Resume_reset_IO_grid_flag(Resume *resume) ;
void Resume_write_data(Resume *resume) ;
void Resume_read_data() ;
void Resume_set_IO_filename(Resume *resume, int vec_id, char *filename) ;
void Resume_set_IO_type(Resume *resume, int IO_type) ;
void Resume_set_vec_max(Resume *resume, int max_vec) ;
void Resume_set_vec_data(Resume *resume, int vec_id, Vec *data) ;
void Resume_set_data(Resume *resume, GVG_bag *data_bag) ;
void Resume_write_primary_data(Resume *resume, GVG_bag *data_bag) ;
void Resume_read_primary_data(Resume *resume, GVG_bag *data_bag) ;
void Resume_write_2D_data(Resume *resume, double **data, char *filename, int NX, int NZ) ;
void Resume_read_2D_data(Resume *resume, double **data, char *filename, int NX, int NZ) ;


#endif
