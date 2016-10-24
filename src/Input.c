/* input.c - input functions*/
/* This file is used to define all the input parameters */
/* to run the code gvg.c */
#include "definitions.h"
#include "DataTypes.h"
#include "Input.h"
#include "Memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*************************************************************************/
/* set_params is used to easily set default parameters                   */
Parameters *Input_set_parameters(void) {

    Parameters *my_params;

    my_params = (Parameters *)malloc(sizeof(Parameters));

    /***********************************************************************/
    /*          Mesh and physical properties of the domain                 */
    /***********************************************************************/
    my_params->NX = 200; /* number of grid points in x direction */
    my_params->NY = 50;  /* number of grid points in y direction */
    my_params->NZ = 50;  /* number of grid points in z direction */

    my_params->xmin = 0.0; /* x: left position in non-dimensional length*/
    my_params->ymin = 0.0; /* y: bottom position in non-dimensional length*/
    my_params->zmin = 0.0; /* y: bottom position in non-dimensional length*/

    my_params->xmax = 8.0;/* x: right position in non-dimensional length*/
    my_params->ymax = 1.0; /* y: top position in non-dimensional length*/
    my_params->zmax = 1.0;/* z: front position in non-dimensional length*/

    /***********************************************************************/
    /*          Flow input properties                                      */
    /***********************************************************************/
    my_params->Re = 1000.0; /* Reynolds number defined as ?????? */ //1000
    my_params->NConc = 1;

    //my_params->Pe = 1000.0; /* Peclet number defined as ?????? */  //1000
    //my_params->V_s0 = -0.02; /* Settling velocity of particles ignoring "Hindered Settling" effect */
    //my_params->epsilon = 0.1; /*Initial volume fraction of the particles */
    my_params->hindered_settling = NO; /* *Hidered settling effect: Settling speed as a function of concentration and epsilon */

    /***********************************************************************/
    /*          Linear system of equation-solver properties                */
    /***********************************************************************/
    my_params->resmax = 1e-8;/* convergence error ?????? */
    my_params->maxit  = 300;  /* Maximum number of iterations for the solution of linear system ?????? */

    /***********************************************************************/
    /*          Simulation and Run properties                              */
    /***********************************************************************/
    my_params->time_max         = 10.0; /* Maximum time in non-dimensional units */
    my_params->output_time_step = 0.2;  /* Output data saving time step */
    my_params->cfl_method       = 2;    /* Method used to apply cfl condition method */
    my_params->alpha = 0.95;           /* Should be used in checking CFL condition*/
    my_params->default_dt = 0.001; /* it is used for the starting point of the simulation */
    my_params->U = 1.0;     /* Outflow convective velocity */

    /***********************************************************************/
    /*          Lock (initial separator of the two fluids) properties      */
    /***********************************************************************/
    my_params->x_fr = 1.0; /* x: initial location of separating sheet: lock  */
    my_params->y_fr = 3.0; /* y: initial location of separating sheet: lock  */
    my_params->z_fr = 2.0; /* z: initial location of separating sheet: lock  */
    my_params->lock_smooth = YES; /* initial condition of the lock ?????? */
    my_params->smooth_length_x = 0.02; /* ????? */
    my_params->smooth_length_y = 0.04; /* ????? */
    my_params->smooth_length_z = 0.04; /* ????? */

    /***********************************************************************/
    /*          Flags indicating flow properties to be saved to file       */
    /***********************************************************************/
    my_params->conc_output = YES;      /* concentration flag */
    my_params->stream_output = NO;     /* stream function flag */
    my_params->vor_output = NO;       /* vorticity flag */
    my_params->u_output = YES;         /* x-component of velocity flag */
    my_params->v_output = YES;         /* y-componenet of velocity flag */
    my_params->w_output = YES;         /* z-componenet of velocity flag */
    my_params->p_output = YES;         /* pressure flag */
    my_params->symmetry_output = NO;  /* Symmtery of the properties flag */
    my_params->ave_height_output = NO;/* Average height flag */
    my_params->front_location_output = NO; /* front location flag */
    my_params->front_speed_output = NO; /* front speed flag */
    my_params->susp_mass_output = NO; /* Suspended mass flag */
    my_params->sedim_rate_output = NO;/* Sedimentation Rate flag */
    my_params->output_format = OPENDX;  /* Output data format: OPENDX, TECPLOT, PARAVIEW*/

    /***********************************************************************/
    /*          Simulation type and appropriate settings                   */
    /***********************************************************************/
    my_params->sim_type = LOCK_EXCHANGE; /* Other option: DRIVEN_CAVITY */
    my_params->outflow = NO; /* Inflow flag */
    my_params->outflow_type = CONVECTIVE_CONSTANT_VELOCITY;
    //my_params->outflow_type = CONVECTIVE_MAXIMUM_VELOCITY;
    my_params->inflow = NO;  /* Outflow flag */
    my_params->transient_inflow = NO; // Smooth transient inflow
    my_params->time_c = 1.0; // Characteristic time for exp decay in transient inflow
    //my_params->viscosity_type = VARIABLE_VISCOSITY;
    my_params->viscosity_type = CONSTANT_VISCOSITY;


    return (my_params);
}/* End of set_params() */

/*************************************************************************/
/* This function releases all the allocated memory for the parameters */
void Input_destroy_parameters(Parameters *params) {

    free(params->Pe);
    free(params->Conc_alpha);
    free(params->V_s0);
    free(params->Phi);
    free(params->conc_type);

    free(params);

}/* End of free_params() */
/*************************************************************************/

void Input_read_parameters_from_file(Parameters *params) {

    FILE *file_in;
    char string[200];
    int i;

    int *flow_flag;
    int *grid_flag;
    int *output_flag;

    double *flow_param;
    double *geom_param;
    double *runtime_param;
    double *lock_param;
    double **conc_param;

    char resume_flag[5];
    int N_geom;
    int N_grid;
    int N_flow_1;
    int N_flow_2;
    int N_runtime;
    int N_flow_flags;
    int N_output_flags;
    int N_lock;
    int iconc, NConc;
    char filename[100];
    sprintf(filename, "input.inp");
    file_in = fopen(filename, "r");

    if (file_in == NULL) {

        PetscPrintf(PCW, "Input.c/ Could not open input file '%s' to read the parameters...\n", filename);
        PetscPrintf(PCW, "Input.c/ Exiting the program.\n");
    } /* if */


    /* Geometric parameters */
    /****************************************************/
    N_geom     = 9;
    geom_param = Memory_allocate_1D_double(N_geom);
    Input_read_sub_section(file_in, geom_param, N_geom);

    /* Grid parameres */
    /****************************************************/
    N_grid    = 4;
    grid_flag = Memory_allocate_1D_int(N_grid);
    Input_read_sub_section_flag(file_in, grid_flag, N_grid);

    //printf("inside input %d %d %d %d\n", grid_flag[0], grid_flag[1], grid_flag[2], grid_flag[3]);

    /* Flow field parameters */
    /****************************************************/
    N_flow_1   = 2;
    flow_param = Memory_allocate_1D_double(N_flow_1);
    Input_read_sub_section(file_in, flow_param, N_flow_1);

    /* Number of concentration fields */
    N_flow_2   = 4;
    NConc      = (int)flow_param[1];
    conc_param = Memory_allocate_2D_double(NConc, N_flow_2);

    for (i=0; i<N_flow_2; i++) {

        if (fgets (string , 150 , file_in) == NULL) {

            PetscPrintf(PCW, "Input.c/ Error reading '%s' file. Probably the structure of the file is altered.\n", filename);
        }

        for (iconc=0; iconc<NConc; iconc++) {

            if (iconc<(NConc-1)) {

                int n_read = fscanf(file_in, "%lf\t", &conc_param[i][iconc]);
                if (n_read != 1) {
                    PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file.\n");
                }
            }
            else {

                int n_read = fscanf(file_in, "%lf\n", &conc_param[i][iconc]);
                if (n_read != 1) {
                    PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file.\n");
                }
                //PetscPrintf(PCW, "N :%lf:\n", conc_param[i][iconc]);
            } /* if */
        } /* for */
    } /* for i */


    /* Runtime parameters */
    /****************************************************/
    N_runtime     = 2;
    runtime_param = Memory_allocate_1D_double(N_runtime);
    Input_read_sub_section(file_in, runtime_param, N_runtime);

    /* Lock Parameters */
    /****************************************************/
    N_lock     = 6;
    lock_param = Memory_allocate_1D_double(N_lock);
    Input_read_sub_section(file_in, lock_param, N_lock);

    /* Flow flags */
    /****************************************************/
    N_flow_flags = 6;
    flow_flag    = Memory_allocate_1D_int(N_flow_flags);
    Input_read_sub_section_flag(file_in, flow_flag, N_flow_flags);

    /* Output flags */
    /****************************************************/
    N_output_flags = 13;
    output_flag    = Memory_allocate_1D_int(N_output_flags);
    Input_read_sub_section_flag(file_in, output_flag, N_output_flags);


    /****************************************************/
    /* Simulation parameters */
    /* Resume or start from t=0 */
    /* Resume data saving time step. */

    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading '%s' file. Probably the structure of the file is altered.\n", filename);
    }
    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading '%s' file. Probably the structure of the file is altered.\n", filename);
    }
    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading '%s' file. Probably the structure of the file is altered.\n", filename);
    }
    if (fgets (resume_flag , 5 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading '%s' file. Probably the structure of the file is altered.\n", filename);
    }

    if (strncmp (resume_flag, "YES",3) == 0) {

        params->resume = YES;
    }
    else {

        params->resume = NO;
    }

    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading '%s' file. Probably the structure of the file is altered.\n", filename);
    }
    int n_read=fscanf(file_in, "%lf\n", &params->resume_saving_timestep);
    if (n_read != 1) {
        PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. \n");
    }



    /* Done with reading the input file */
    fclose (file_in);


    /****************************************************/
    /****************************************************/
    /****************************************************/
    /* Now, put all the read data into the parameter */
    /****************************************************/
    /****************************************************/
    /****************************************************/

    params->NX   = (int) geom_param[0];
    params->NY   = (int) geom_param[1];
    params->NZ   = (int) geom_param[2];
    params->xmin = geom_param[3];
    params->xmax = geom_param[4];
    params->ymin = geom_param[5];
    params->ymax = geom_param[6];
    params->zmin = geom_param[7];
    params->zmax = geom_param[8];

    /***********************************************************************/
    /*          Flow input properties                                      */
    /***********************************************************************/
    params->Re    = flow_param[0];
    params->NConc = (int)flow_param[1];

    /* Now allocate memory based on the number of concentration fields */
    params->Pe          = Memory_allocate_1D_double(NConc);
    params->Conc_alpha  = Memory_allocate_1D_double(NConc);
    params->V_s0        = Memory_allocate_1D_double(NConc);
    params->Phi         = Memory_allocate_1D_double(NConc);
    params->conc_type   = Memory_allocate_1D_int(NConc);

    for (iconc=0; iconc<NConc; iconc++) {

        params->Pe[iconc]         = conc_param[0][iconc];
        params->Conc_alpha[iconc] = conc_param[1][iconc];
        params->V_s0[iconc]       = conc_param[2][iconc];
        params->Phi[iconc]        = conc_param[3][iconc];
        if (fabs(params->V_s0[iconc]) < 1.0e-10) {

            params->conc_type[iconc]  = SALINITY; /* Type of the concentration field */
        } else {

            params->conc_type[iconc]  = PARTICLE; /* Type of the concentration field */
        }
    }

    /***********************************************************************/
    /*          Grid properties                                            */
    /***********************************************************************/
    params->UniformGridX        = grid_flag[0];
    params->UniformGridY        = grid_flag[1];
    params->UniformGridZ        = grid_flag[2];
    params->ImportGridFromFile  = grid_flag[3];

    /***********************************************************************/
    /*          Simulation and Run properties                              */
    /***********************************************************************/
    params->time_max         = runtime_param[0];
    params->output_time_step = runtime_param[1];

    /***********************************************************************/
    /*          Lock (initial separator of the two fluids) properties      */
    /***********************************************************************/
    params->x_fr = lock_param[0];
    params->y_fr = lock_param[1];
    params->z_fr = lock_param[2];
    params->smooth_length_x = lock_param[3];
    params->smooth_length_y = lock_param[4];
    params->smooth_length_z = lock_param[5];


    /***********************************************************************/
    /*          Simulation type and appropriate settings                   */
    /***********************************************************************/
    if (flow_flag[0] == YES) {

        params->viscosity_type = VARIABLE_VISCOSITY;
    }
    else {

        params->viscosity_type = CONSTANT_VISCOSITY;
    }

    params->sedimentation     = flow_flag[1];
    params->hindered_settling = flow_flag[2];
    params->resuspension      = flow_flag[3];
    params->inflow            = flow_flag[4];
    params->outflow           = flow_flag[5];

    /***********************************************************************/
    /*          Flags indicating flow properties to be saved to file       */
    /***********************************************************************/
    params->conc_output           = output_flag[0];
    params->u_output              = output_flag[1];
    params->v_output              = output_flag[2];
    params->w_output              = output_flag[3];
    params->p_output              = output_flag[4];
    params->stream_output         = output_flag[5];
    params->vor_output            = output_flag[6];
    params->symmetry_output       = output_flag[7];
    params->ave_height_output     = output_flag[8];
    params->front_location_output = output_flag[9];
    params->susp_mass_output      = output_flag[10];
    params->sedim_rate_output     = output_flag[11];
    params->front_speed_output    = output_flag[12];
    //params->deposit_height_output = output_flag[13];

    Memory_free_2D_array(N_flow_2, (void **)conc_param);
    free(geom_param);
    free(flow_param);
    free(grid_flag);
    free(runtime_param);
    free(flow_flag);
    free(output_flag);
    free(lock_param);

}
/*************************************************************************/

/* This function is used to adjust the other parameters which are not read from file. */
void Input_set_internal_parameters(Parameters *params) {

    short int RECALCULATE_MULTIGRAIN_PROPERTIES = YES;

    params->Lx = params->xmax - params->xmin;
    params->Ly = params->ymax - params->ymin;
    params->Lz = params->zmax - params->zmin;
    /***********************************************************************/
    /*          Linear system of equation-solver properties                */
    /***********************************************************************/
    params->resmax = 1e-12;/* convergence error ?????? */
    params->maxit  = 1000;  /* Maximum number of iterations for the solution of linear system ?????? */

    /***********************************************************************/
    /*          Simulation and Run properties                              */
    /***********************************************************************/
    params->cfl_method = 2;         /* Method used to apply cfl condition method */
    params->alpha      = 0.50;      /* Should be used in checking CFL condition*/
    params->default_dt = 0.01;    /* it is used for the starting of the simulation */
    params->U = 1.0;                /* Convective outflow velocity */

    /***********************************************************************/
    /*          Lock (initial separator of the two fluids) properties      */
    /***********************************************************************/
    params->lock_smooth = YES; /* initial condition of the lock ?????? */

    /***********************************************************************/
    /*          Flags indicating flow properties to be saved to file       */
    /***********************************************************************/
    params->output_format = OPENDX;  /* Output data format: OPENDX, TECHPLOT, PARAVIEW*/
    params->div_output    = NO;

    /***********************************************************************/
    /*          Simulation type and appropriate settings                   */
    /***********************************************************************/
    params->sim_type = LOCK_EXCHANGE;
    //params->sim_type = DRIVEN_CAVITY;

    params->outflow_type = CONVECTIVE_MAXIMUM_VELOCITY;
    //params->outflow_type = CONVECTIVE_CONSTANT_VELOCITY;

    params->transient_inflow = NO;
    params->time_c = 0.035; /* Characteristic time for exp decay in transient inflow */

    params->ENO_scheme_order = 3; /* Max: 3 */


    /* Here, based on the given parameters for multiple concentration fields, we are calculting the
 coefficients based on the first conc field */
    if (RECALCULATE_MULTIGRAIN_PROPERTIES) {

        double Et = 0.0;
        double Temp_Conc_alpha0;
        double Ei;
        int iconc;

        for (iconc=0; iconc<params->NConc; iconc++) {

            Ei  = params->Conc_alpha[iconc] / params->Conc_alpha[0];
            Et += Ei;
        }

        Temp_Conc_alpha0 = params->Conc_alpha[0];

        for (iconc=0; iconc<params->NConc; iconc++) {

            Ei  = params->Conc_alpha[iconc] / Temp_Conc_alpha0;
            /* Rescale the "alpha" coefficients (Used to calculate the total concentration term
   in v-mom equation */
            params->Conc_alpha[iconc] = Ei / Et;

            /* Find the physical volume fraction based on the first one: Used for Hindered settling
   effect and also variable viscosity */
            params->Phi[iconc]        = Ei * params->Phi[0];

            /* Find the settling speed based on the first one */
            params->V_s0[iconc]       *= -1.0;
        }
    }

    //params->resume_saving_timestep = 3.0; /* time step which the whole simulation data is saved to future run*/
    //params->resume = YES;
    //params->deposit_height_output = YES;

    /* No change due to erosion or deposition */
    params->constant_geometry = YES;

    /* energies output flag */
    params->energies_output = YES;

    /* flag to compute shear stress on the bottom boundary */
    params->shear_stress_output = YES;
    params->deposit_height_output = YES;
    params->dump_conc = YES;

    /* Flag indicating if the bottom interface should be imported from file */
    params->ImportBottomInterface = NO;

    /* BINARY, HDF5, ASCII */
    params->output_type = BINARY;

    /* Linear system matrix output flag */
    params->output_lsys_matrix = NO;

}
/*************************************************************************/

/* This function reads a subsection of the input file (including N doubles) */
void Input_read_sub_section(FILE *file_in, double *read_params, int N) {

    int i;
    char string[200];

    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
    }
    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
    }

    for (i=0; i<N; i++) {

        if (fgets (string , 150 , file_in) == NULL) {

            PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
        }
        int n_read = fscanf(file_in, "%lf\n", &read_params[i]);
        if (n_read != 1) {
            PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file.\n");
        }
        //puts(string);
        //PetscPrintf(PCW, "%lf:\n", read_params[i]);
    } /* for */
}
/*************************************************************************/

/* This function reads a subsection of the input file (including N YES or NO flags) */
void Input_read_sub_section_flag(FILE *file_in, int *flags, int N) {

    int i;
    char string[200];
    char **read_params;

    /* 2d arrays of characters */
    read_params = Memory_allocate_2D_char(5, N);

    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
    }
    if (fgets (string , 150 , file_in) == NULL) {

        PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
    }

    for (i=0; i<N; i++) {

        if (fgets (string , 150 , file_in) == NULL) {

            PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
        }
        if (fgets (read_params[i] , 5 , file_in) == NULL) {

            PetscPrintf(PCW, "Input.c/ Error reading 'input.inp' file. Probably the structure of the file is altered.\n");
        }

        //puts(string);

        if (strncmp (read_params[i],"YES",3) == 0) {

            flags[i] = YES;
        }
        else {

            flags[i] = NO;
        }
    }

    Memory_free_2D_array(N, (void **)read_params);
}

void Input_set_sphere_params(Parameters *params, MAC_grid *grid) {

	params->flow_over_sphere =  YES; 
	params->sphere_cent_x = 0.34375*params->Lx; /* Lx = 14 */
	params->sphere_cent_y = 0.5*params->Ly;     /* Ly = 14 */
	params->sphere_cent_z = 0.5*params->Lz;     /* Lz = 14 */

	params->sphere_radius = 0.5; 
	params->current_rank_has_sphere = NO; 

	int Is = grid->G_Is;
	int Js = grid->G_Js;
	int Ks = grid->G_Ks;

	/* End index of top-right-front corner on current processor */
	int Ie = grid->G_Ie;
	int Je = grid->G_Je;
	int Ke = grid->G_Ke;

 	if  ( (Is < grid->NI/2) && (Ie > grid->NI/2) &&
			(Js < grid->NJ/2) && (Je > grid->NJ/2) &&
			(Ks < grid->NK/2) && (Ke > grid->NK/2) ) {

		params->current_rank_has_sphere = YES; 
		printf("Processor %d has the sphere. \n", params->rank); 
	}
}
