#ifndef WRITER_H
 #define WRITER_H

Writer *Writer_create(int N1, int N2, int N3, int precision) ;
void Writer_set_dimension(Writer *writer, int dimension) ;
void Writer_set_filename(Writer *writer, char *filename) ;
void Writer_set_output_type(Writer *writer, short int Output_Type) ;
void Writer_set_append_data_flag(Writer *writer) ;
void Writer_reset_append_data_flag(Writer *writer) ;
void Writer_set_write_grid_flag(Writer *writer) ;
void Writer_reset_write_grid_flag(Writer *writer) ;
void Writer_set_output_data(Writer *writer, void *data) ;
void Writer_open_file(Writer *writer) ;
void Writer_close_file(Writer *writer) ;
void Writer_write_grid(Writer *writer) ;
void Writer_write_data(Writer *writer) ;
void Writer_check_file_open(FILE *file) ;
void Writer_set_grid(Writer *writer, void *X1_grid, void *X2_grid, void *X3_grid) ;
void Writer_check_file_write(int n_items_desired, int n_items_actual) ;

void Writer_check_file_open(FILE *file) ;
void Writer_check_file_open(FILE *file) ;

void Writer_destroy(Writer *writer) ;

#endif
