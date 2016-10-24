#include "definitions.h"
#include "DataTypes.h"
#include "Memory.h"
#include "Writer.h"
#include <stdlib.h>
#include <stdio.h>


Writer *Writer_create(int N1, int N2, int N3, int precision) {

    Writer *new_writer;

    new_writer = (Writer *)malloc(sizeof(Writer));
    Memory_check_allocation(new_writer);

    /* Number of nodes in each direction */
    new_writer->N1        = N1;
    new_writer->N2        = N2;
    new_writer->N3        = N3;

    new_writer->NT        = N1 * N2 * N3;
    /* Dimension of the data */
    new_writer->Dimension = 1; /* 1D */

    if (N2 > 1) {

        new_writer->Dimension = 2; /* 2D */
    }
    if (N3 > 1) {

        new_writer->Dimension = 3; /* 3D */
    }

    if (precision == GVG_FLOAT) {

        new_writer->precision_size = sizeof(float);

    } else if (precision == GVG_DOUBLE) {

        new_writer->precision_size = sizeof(double);
    }

    /* default */
    new_writer->Output_Type = BINARY;  /* other option ASCII */
    new_writer->Append_Data = NO;

    return (new_writer);
}
/**********************************************************************************/


void Writer_destroy(Writer *writer) {

    free(writer);
}
/**********************************************************************************/

/* This function sets the dimension for the current writer */
void Writer_set_dimension(Writer *writer, int dimension) {

    writer->Dimension = dimension;
}
/**********************************************************************************/

/* This function sets the writer filename */
void Writer_set_filename(Writer *writer, char *filename) {

    writer->Filename = filename;
}
/**********************************************************************************/

/* This function sets the output type flag: BINARY or ASCII */
void Writer_set_output_type(Writer *writer, short int Output_Type) {

    writer->Output_Type = Output_Type; /* BINARY or ASCII */
}
/**********************************************************************************/

/* This function sets the flag for Append data to current file */
void Writer_set_append_data_flag(Writer *writer) {

    writer->Append_Data = YES;
}
/**********************************************************************************/

/* This function sets the flag for Append data to current file */
void Writer_reset_append_data_flag(Writer *writer) {

    writer->Append_Data = NO;
}
/**********************************************************************************/

/* This function sets the grid-writer flag on. Grid data is added to the beginning of the file */
void Writer_set_write_grid_flag(Writer *writer) {

    writer->Write_Grid = YES;
}
/**********************************************************************************/

/* This function sets the grid-writer flag on. Grid data is added to the beginning of the file */
void Writer_reset_write_grid_flag(Writer *writer) {

    writer->Write_Grid = NO;
}
/**********************************************************************************/

/* This function sets the data pointer of current writer to the given data */
void Writer_set_output_data(Writer *writer, void *data) {

    writer->output_data = data;
}
void Writer_open_file(Writer *writer) {

    char *filename;
    FILE *file=NULL;

    filename = writer->Filename;

    switch (writer->Output_Type) {

    case BINARY:

        /* add to the end of current binary file */
        if (writer->Append_Data) {

            file = fopen(filename, "ab");

            /* create an empty binary file */
        } else {

            file = fopen(filename, "wb");
        }
        Writer_check_file_open(file);
        break;

    case ASCII:

        /* add to the end of current ASCII file */
        if (writer->Append_Data) {

            file = fopen(filename, "a");

            /* create an empty ASCII file */
        } else {

            file = fopen(filename, "w");
        }
        Writer_check_file_open(file);
        break;

    }
    writer->file = file;
}
/**********************************************************************************/

/* This function closes the cuurently opened file */
void Writer_close_file(Writer *writer) {

    fclose(writer->file);
}
/**********************************************************************************/

void Writer_write_grid(Writer *writer) {

    int nitems;

    int N1        = writer->N1;
    int N2        = writer->N2;
    int N3        = writer->N3;
    int dimension = writer->Dimension;
    FILE *file      = writer->file;

    size_t write_size = writer->precision_size;

    if (writer->Write_Grid) {

        switch(writer->Output_Type) {

        case BINARY:

            /* first coordinate: x */
            nitems = fwrite(writer->X1_grid, write_size, N1, file);
            Writer_check_file_write(N1, nitems);

            /* second coordinate: y */
            if (dimension > 1) {

                nitems = fwrite(writer->X2_grid, write_size, N2, file);
                Writer_check_file_write(N2, nitems);
            }

            /* third coordinate: z */
            if (dimension > 2) {

                nitems = fwrite(writer->X3_grid, write_size, N3, file);
                Writer_check_file_write(N3, nitems);
            }
            break;

        case ASCII:
            //My debug
            //To be fixed
            /*
            int i;
            for (i=0; i<N1; i++) {

                fprintf(file, "%f\n", X1_grid[i]);
            }

            if (dimension > 1) {

                for (i=0; i<N2; i++) {

                    fprintf(file, "%f\n", X2_grid[i]);
                }
            }
            if (dimension > 2) {

                for (i=0; i<N3; i++) {

                    fprintf(file, "%f\n", X3_grid[i]);
                }
            }
            */
            break;
        } /* switch */

    } /* if */
}
/**********************************************************************************/

void Writer_write_data(Writer *writer) {

    int i;
    int nitems;
    int NT;
    float *data;
    FILE *file;

    NT   = writer->NT;
    file = writer->file;

    data = writer->output_data;

    switch(writer->Output_Type) {

    case BINARY:

        /* first coordinate: x */
        nitems = fwrite(writer->output_data, writer->precision_size, NT, file);
        Writer_check_file_write(NT, nitems);

        break;

    case ASCII:

        for (i=0; i<NT; i++) {

            fprintf(file, "%f\n", data[i]);
        }
        break;
    }
}
/**********************************************************************************/

void Writer_check_file_open(FILE *file) {

    if (file == NULL) {

        PetscPrintf(PCW, "Writer.c/ Error opening file...\n");
        PetscPrintf(PCW, "Writer.c/ Exiting the program.\n \n");
        PetscFinalize();
        exit(1);

    }
}
/**********************************************************************************/

/* This function passes the pointer to grid data to writer */
void Writer_set_grid(Writer *writer, void *X1_grid, void *X2_grid, void *X3_grid) {

    writer->X1_grid = X1_grid;
    writer->X2_grid = X2_grid;
    writer->X3_grid = X3_grid;
}
/**********************************************************************************/

/* This function checks if the data has been written correctly to the binary file */
void Writer_check_file_write(int n_items_desired, int n_items_actual) {

    if (n_items_desired != n_items_actual) {

        PetscPrintf(PCW, "Writer.c/ ******** Warning ********...\n There was an error while writing the data to file\n");
    }
}
/**********************************************************************************/

