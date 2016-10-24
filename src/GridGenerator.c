#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Theser are defined for the Mesh Generation */
#define UNIFORM 0
#define NONUNIFORM_ONE_SIDED 1
#define NONUNIFORM_TWO_SIDED 2

/* Stretching function ID */
#define ONE_SIDED_TANH 1
#define ONE_SIDED_TAN  2
#define ONE_SIDED_POL  3
#define TWO_SIDED_TANH 4
#define TWO_SIDED_TAN  5
#define TWO_SIDED_POL  6

/* Slope of the nonuniform region if there is only one region. Other than that, the slope of the nonuniform region will be computed using the uniform region */
#define END_SLOPE 0.35

struct subregion {
	
	double lmin, lmax; /* min and max of the sub region */
	double N; /* Number of grids in the current domain */
	short int type; /* Type of the region */
	short int StretchingFunctionID; /* Stretching function ID */
	/* Stretching function properties */
	double S0, S1;
	double B, A, Q;
};
typedef struct subregion SubRegion; 

struct gridinfo {
	
	int Nsub_region; /* Number of sub regions */
	double Lmin, Lmax; /* Max and minumum coordinate of the domain */
	double N; /* Total number of grid points */
	SubRegion **sub_region; /* array holding all the sub_region properties */
	char  *filename; 
	double *position_c; /* position of the cell center */
	double *position_s; /* position of the side cell, (one more than the number of grids) */
}; 
typedef struct gridinfo GridInfo; 

GridInfo *GridGen_create(int Nsub_region, int N_total, double Lmin, double Lmax) ;
void GridGen_destroy(GridInfo *grid) ;
void GridGen_set_sub_region_type(GridInfo *grid, short int *region_type) ;
void GridGen_set_sub_region_range(GridInfo *grid, double *region_range) ;
void GridGen_set_sub_region_Ngrid(GridInfo *grid, int *regions_N) ;
void GridGen_compute_sub_region_coordinates(GridInfo *grid, int region_id) ;
void GridGen_set_stretching_function_type(SubRegion *sub_region) ;
double GridGen_stretching_function_value(SubRegion *sub_region, double eta) ;
double GridGen_compute_Q(double slope, short int which_equation) ;
void GridGen_set_filename(GridInfo *grid, char *filename) ;
void GridGen_write_grid(GridInfo *grid) ;
void GridGen_setup_first_sub_region_properties(GridInfo *grid) ;
void GridGen_sub_region_boundary_properties(GridInfo *grid, int region_id) ; 



int main()  {
	
	/* User should set these parameters */

/*********************************************************************/
/*********************************************************************/
/* INPUT THE PARAMETERS HERE */
/*********************************************************************/
/*********************************************************************/
	/* Number of sub regions. At least should be 1 */
        int   N_reg = 3;
	
	/* Coordinates of each sub region. */
        double Region_range[]={0.0, 0.2, 1.8, 2.0};
	
	/* Number of grids within each sub region */
        int   Region_Ngrid[]={5, 25, 5};
	/* 
		1- UNIFORM: uniform grid 
		2- NONUNIFORM_TWO_SIDED: non-uniform with the slopes defined at both ends 
		3- NONUNIFORM_ONE_SIDED: non-uniform with the slope  defined at the beginnning 
	        Try to have uniform grid bounding NONUNIFORM_TWO_SIDED sub-region. 
		*/

        short int Region_type[]={UNIFORM, NONUNIFORM_TWO_SIDED, UNIFORM};

	/* Output file name */
	/* Note that the way the result is written to file is:
		1- First, grid side (e.g. x position of u-grid in a MAC-grid)                 (N_total+1)
		2- Second, grid cell center (e.g. x position of center-grid in a MAC-grid)    (N_total)
	    */
	char  Filename[]="grid.dat";
/*********************************************************************/
/*********************************************************************/
/* End of INPUT section */
/*********************************************************************/
/*********************************************************************/

	int   N_total;
	int i; 
	int Nsub_region; 
	GridInfo *grid; 


	N_total = 0; 
	for (i=0; i<N_reg; i++) {
		
		N_total += Region_Ngrid[i]; 
	}
	grid = GridGen_create(N_reg, N_total, Region_range[0], Region_range[N_reg]); 
	Nsub_region = grid->Nsub_region; 

	
	printf("Grid has been initialized successfully\n"); 
	GridGen_set_filename(grid, Filename); 


	GridGen_set_sub_region_type (grid, Region_type);
	GridGen_set_sub_region_range(grid, Region_range);
	GridGen_set_sub_region_Ngrid(grid, Region_Ngrid);


	for (i=0; i<Nsub_region; i++) {
		
		GridGen_sub_region_boundary_properties(grid, i);
	} /* for i*/
	
	/* Compute the grid coordinates for each direction */
	for (i=0; i<Nsub_region; i++) {
		
		GridGen_compute_sub_region_coordinates(grid, i);
	}
	
	GridGen_write_grid(grid);
	printf("Grid has been written to \"%s\" successfully.\n", grid->filename); 
	
	GridGen_destroy(grid); 
}
/**************************************************************************************/

/* This function creates a grid structure based on the input paramaters */
GridInfo *GridGen_create(int Nsub_region, int N_total, double Lmin, double Lmax) {
	
	int i;
	GridInfo *new_grid; 
	
	new_grid = (GridInfo *)malloc(sizeof(GridInfo)); 

	new_grid->Nsub_region = Nsub_region; 

	/* Now, allocate memory for each sub region grids */
	new_grid->sub_region   = (SubRegion **)malloc(new_grid->Nsub_region*sizeof(SubRegion *));
	for (i=0; i<Nsub_region; i++) {
			
		new_grid->sub_region[i] = (SubRegion *)malloc(sizeof(SubRegion));		
	}
	
	/* Max and min of the region */
	new_grid->Lmin = Lmin; 
	new_grid->Lmax = Lmax; 

	/* Total number of cells */
	new_grid->N = N_total; 

	/* Position of the cell center and side cell */
	new_grid->position_c = (double *)calloc(N_total, sizeof(double)); 
	new_grid->position_s = (double *)calloc(N_total+1, sizeof(double)); 
	
	return new_grid; 
}
/**************************************************************************************/

/* This function releases the allocated memory for the grid */
void GridGen_destroy(GridInfo *grid) {

	int i; 
	
	free(grid->position_c);
	free(grid->position_s);

	for (i=0; i<grid->Nsub_region; i++) {
			
		free(grid->sub_region[i]); 
	}
	free(grid->sub_region); 
	free(grid); 
}
/**************************************************************************************/

/* This function sets the number of sub_regions in each dierction. Note that for uniform grid, we only have one
region. For the non-uniform regions, one can set any arbitrary regions > 1 */
void GridGen_set_sub_region_type(GridInfo *grid, short int *region_type) {
	
	int i; 
	
	for (i=0; i<grid->Nsub_region; i++) {
		
		grid->sub_region[i]->type = region_type[i]; 
	}
}
/**************************************************************************************/

/* This function sets the lmax (lmin) of each sub region */
void GridGen_set_sub_region_range(GridInfo *grid, double *region_range) {

	int Nsub_region;
	int i, id; 
	
	/* Number of subregions */
	Nsub_region = grid->Nsub_region;

	for (i=0; i<Nsub_region; i++) {
		
		grid->sub_region[i]->lmax = region_range[i+1];
	}

	grid->sub_region[0]->lmin = region_range[0]; 
	/* Set the next grid lmin to the lmax of the previous region */
	for (i=1; i<Nsub_region; i++) {
		
		grid->sub_region[i]->lmin = grid->sub_region[i-1]->lmax;
		
	} /* for i */
	
	/* Just to make sure, for a one region system, use only xmin and xmax */ 
	if (Nsub_region == 1) {

		grid->sub_region[0]->lmin = grid->Lmin; 
		grid->sub_region[0]->lmax = grid->Lmax; 
	} 
}
/**************************************************************************************/
	
/* This function sets the number of grid points in each sub region */
void GridGen_set_sub_region_Ngrid(GridInfo *grid, int *regions_N) {

	int Nsub_region;
	int i, id; 
	
	/* Number of subregions */
	Nsub_region = grid->Nsub_region;

	for (i=0; i<Nsub_region; i++) {
		
		grid->sub_region[i]->N = regions_N[i];
	}
	
	/* Just to make sure, for a one region system, use only N */
	if (Nsub_region == 1) {

		grid->sub_region[0]->N = grid->N; 
	} 
}
/**************************************************************************************/

/* This function computes the sub_region coordinates */
void GridGen_compute_sub_region_coordinates(GridInfo *grid, int region_id) {
	
	
	int N;
	int i, index; 
	int Nregion; 
	int global_start_index; 
	double dEta;
	double eta;
	double Length;
	double lmin, lmax; 
	double *position_c, *position_s; 
	
	/* global grid coordinates */
	position_s = grid->position_s; 
	position_c = grid->position_c; 
	
	Nregion = grid->sub_region[region_id]->N; /* Number of grid cells in the current sub region */
	dEta    = 1.0/(double) Nregion; /* uniform delta_e in the transfomred coordinates system */
	
	lmin   = grid->sub_region[region_id]->lmin; 
	lmax   = grid->sub_region[region_id]->lmax; 
	Length = lmax - lmin; 

	global_start_index = 0; 

	for (i=0; i<region_id; i++) {
		
		global_start_index += grid->sub_region[i]->N; 
	}
	/* First, coordinate at the cell face */
	eta  = 0.0;
	for (i=0 ; i<=Nregion ; i++) {
		
		index             = i+ global_start_index; 
		position_s[index] = lmin + Length * GridGen_stretching_function_value(grid->sub_region[region_id], eta);
		eta              += dEta;
		//printf("region:%d eta:%f x_u:%f\n", region_id, eta, position_s[index]);
		//getchar();
		
	}

	/* cell center coordinates */
	eta  = 0.5*dEta;
	for (i=0 ; i<Nregion ; i++) {
		
		index             = i+ global_start_index; 
		position_c[index] = lmin + Length * GridGen_stretching_function_value(grid->sub_region[region_id], eta);
		eta              += dEta;
	}

}
/**************************************************************************************/

/* This function sets the stretching function based on the slope given at the boundary nodes */
void GridGen_set_stretching_function_type(SubRegion *sub_region) {
	
	double tol = 0.01;
	double B, S0; 
	short int S_ID=-1; 
	
	B  = sub_region->B; 
	S0 = sub_region->S0; 

	switch (sub_region->type) {
		
		/* If the sub region is uniform */
		case UNIFORM: 
			
			S_ID = UNIFORM;
			break; 
		
		/* If the nonuniform grid boundary condition is given at one side (probably only the left side) */
		case NONUNIFORM_ONE_SIDED: 
			
			if (S0 > 1.0) {

				S_ID = ONE_SIDED_TANH;
			
			} else if (S0 < 1.0) {

				S_ID = ONE_SIDED_TAN;
			}
			if ( fabs(S0 - 1.0) < tol )  {
			
				S_ID = ONE_SIDED_POL;
			}
			break; 
		
		/* If the nonuniform grid boundary condition is given at both ends */
		case NONUNIFORM_TWO_SIDED: 
			
			if (B > 1.0) {

				S_ID = TWO_SIDED_TANH;
			
			} else if (B < 1.0) {

				S_ID = TWO_SIDED_TAN;
			}
			if ( fabs(B - 1.0) < tol )  {
			
				S_ID = TWO_SIDED_POL;
			}
			break; 
			
	} /* switch */

	//printf("Region Type:%d Found S_ID:%d B:%f S0:%f Q:%f\n", sub_region->type, S_ID, B, S0, sub_region->Q); 
	//getchar(); 
	sub_region->StretchingFunctionID = S_ID; 
}
/**************************************************************************************/

/* This function returns the value of the real coordinate based on the given stretching function for the given
subregion */
double GridGen_stretching_function_value(SubRegion *sub_region, double eta) {
	
	double value;
	double Q, A, B, S0;
	
	Q  = sub_region->Q;
	S0 = sub_region->S0;
	A  = sub_region->A;
	B  = sub_region->B;
	
	switch (sub_region->StretchingFunctionID)  {

		case UNIFORM:
			
			value = eta;		
			break;
			
		case ONE_SIDED_TANH:
			
			value = 1.0 + tanh( Q*(eta - 1.0) ) / tanh(Q) ;		
			break;

		case ONE_SIDED_TAN:
			
			value = 1.0 + tan( Q*(eta - 1.0) ) / tan(Q) ;		
			break;
		
		case ONE_SIDED_POL:
			
			value = eta * ( 1.0 - 0.5*(S0-1) * ( 1.0-eta )*( 2.0-eta ) );
			break;
			
		case TWO_SIDED_TANH:
			
			value = 0.5 + 0.5* tanh( Q*(eta - 0.5) ) / tanh(0.5*Q)  ;
			value = value / (A + (1.0 - A)*value); 
			break;

		case TWO_SIDED_TAN:
			
			value = 0.5 + 0.5* tan( Q*(eta - 0.5) ) / tan(0.5*Q)  ;
			value = value / (A + (1.0 - A)*value); 
			break;
		
		case TWO_SIDED_POL:
			
			value = eta * ( 1.0 + 2.0*(B-1.0) * (eta - 0.5)*( 1.0-eta ) );
			value = value / (A + (1.0 - A)*value); 
			break;
		default:
			printf("Invalid function ID\n");
			break;
	}		
	
	//printf("SF_ID:%d B:%f Q:%f S0:%f value:%f\n", sub_region->StretchingFunctionID, B, Q, S0, value); 
	return(value);
}
/**************************************************************************************/

/* This function sets the next sub_region's boundary properties based on the previous one */
void GridGen_sub_region_boundary_properties(GridInfo *grid, int region_id) {
	
	double Q;
	double S0, S1;
	double A, B;
	double L_p, L_c, L_n;
	int   N_p, N_c, N_n;
	int Nsub_region;
	
	
	/* current region's properties */
	L_c  = grid->sub_region[region_id]->lmax - grid->sub_region[region_id]->lmin;
	N_c  = grid->sub_region[region_id]->N;
	
	Nsub_region = grid->Nsub_region;
	
	if (region_id == 0) {

		GridGen_setup_first_sub_region_properties(grid);
		
		S0 = grid->sub_region[region_id]->S0;
		S1 = grid->sub_region[region_id]->S1;
	} else {
	
		if (grid->sub_region[region_id]->type == UNIFORM) {
		
			S0 = 1.0;
			S1 = 1.0;
		
		} else {
	
			L_p = grid->sub_region[region_id-1]->lmax - grid->sub_region[region_id-1]->lmin;
			N_p = grid->sub_region[region_id-1]->N;
		
			/* This ensures that the slope at first grids of the two successive regions are the same */
			S0 = (L_c/L_p) * (double) N_p/(double)N_c * grid->sub_region[region_id-1]->S1;

			if ( (Nsub_region - region_id) > 1 ) {
		
				if ( ( grid->sub_region[region_id]->type == NONUNIFORM_TWO_SIDED ) && (grid->sub_region[region_id+1]->type == UNIFORM) ) {
			
					L_n = grid->sub_region[region_id+1]->lmax - grid->sub_region[region_id+1]->lmin;
					N_n = grid->sub_region[region_id+1]->N;

					/* Two sided function is NOT anti-symmetric */
					S1 = (L_c/L_n) * (double) N_n/(double)N_c ;
					//printf("ID:%d Done\n", region_id);
					//getchar(); 

				} else {

					/* For now, assume that the two sided nonuniform mesh is anti-symmetric */
					S1 = S0;
				}
			} else {
				/* For now, assume that the two sided nonuniform mesh is anti-symmetric */
				S1 = S0;
			}
		} /* else UNIFROM */

	}

	/* This ensures that the slope and first grids of the two successive regions are the same */
	A = sqrt(S0 / S1);
	B = sqrt(S0 * S1);
	grid->sub_region[region_id]->S0 = S0;
	grid->sub_region[region_id]->S1 = S1;
	grid->sub_region[region_id]->A  = A;
	grid->sub_region[region_id]->B  = B;


	GridGen_set_stretching_function_type(grid->sub_region[region_id]);

	switch (grid->sub_region[region_id]->StretchingFunctionID)  {

		case UNIFORM:
			
			break;
			
		case ONE_SIDED_TANH:
			
			Q = 0.5*GridGen_compute_Q(S0, ONE_SIDED_TANH);
			break;

		case ONE_SIDED_TAN:
			
			Q = 0.5*GridGen_compute_Q(S0, ONE_SIDED_TAN);
			break;
		
		/* Have to fix this */
		case ONE_SIDED_POL:
			
			break;
			
		case TWO_SIDED_TANH:
			
			Q = GridGen_compute_Q(B, TWO_SIDED_TANH);
			break;

		case TWO_SIDED_TAN:
			
			Q = GridGen_compute_Q(B, TWO_SIDED_TAN);
			break;

		/* Have to fix this */
		case TWO_SIDED_POL:
			
			break;
		default:
			printf("Invalid function ID\n");
			break;
	} /* switch */
	
	//printf("Region ID:%d S0:%f S1:%f A:%f B:%f Q:%f\n", region_id, S0, S1, A, B, Q); 
	grid->sub_region[region_id]->Q = Q;	
}
/**************************************************************************************/

/* This functtion solves the two possilble equations to find the Parameter Q. */
/* 1- slope = sinh(x)/x;
/* 2- slope = sin(x)/x;
*/
/* See paper: Marcel Vinokur, 1980, NASA report. 
"One one-dimensional Stretching functions for finite difference calculations" */
double GridGen_compute_Q(double slope, short int which_equation) {

	
	double x;
	double y, y_b;
	double v, w;
	double pi = acos(-1.0);
	
	y = slope; 
	
	
	if ( (which_equation == ONE_SIDED_TANH) || (which_equation == TWO_SIDED_TANH) ) {
		
		if (y < 2.7829681) {
			
			y_b = y - 1.0;
			x   = sqrt(6.0*y_b) * (1.0 - 0.15*y_b + 0.057321429*y_b*y_b - 
			      0.024907295*pow(y_b, 3.0) + 0.0077424461*pow(y_b, 4.0) - 0.0010794123 *
			      pow(y_b, 5.0) );  
		}
		else {
			
			v = log(y);
			w = 1.0/y - 0.028527431;
			
			x = v + (1.0 + 1.0/v) * log(2.0*v) - 0.02041793 + 0.24902722*w + 
			    1.9496443*w*w - 2.6294547*pow(w, 3.0) + 8.56795911*pow(w, 4.0);   
		}
	} 

	if ( (which_equation == ONE_SIDED_TAN) || (which_equation == TWO_SIDED_TAN) ) {
		
		if (y < 0.26938972) {
			
			x = pi* ( 1.0 - y + y*y - (1.0+pi*pi/6.0) * pow(y, 3.0) + 6.794732*pow(y, 4.0) 
			    - 13.205501*pow(y, 5.0) + 11.726095*pow(y, 6.0) );
		}
		else { /* y < 1.0 */
			
			y_b = 1.0 - y;
			x   = sqrt(6.0*y_b) * (1.0 + 0.15*y_b + 0.057321429*y_b*y_b + 
			      0.048774238*pow(y_b, 3.0) - 0.053337753*pow(y_b, 4.0) - 0.075845134 *
			      pow(y_b, 5.0) );  
		}
	} 

	return (x);
}
/**************************************************************************************/

/* This function sets the grid output filename */
void GridGen_set_filename(GridInfo *grid, char *filename) {
	
	grid->filename = filename; 
	
}
/**************************************************************************************/

/* This function writes the grid postion to file. First cell center coordinates and then the side-cell
coordinates. */
void GridGen_write_grid(GridInfo *grid) {
	
	FILE *f_out; 
	int i, N; 
	double eta, dEta; 
	
	f_out = fopen(grid->filename, "w"); 

	N    = grid->N; 
	dEta = 1.0/(double) N; 

	eta  = 0.5*dEta; 
	for (i=0; i<N; i++) {
		
		fprintf(f_out, "%2.14lf\n", grid->position_c[i]); 
		eta += dEta; 
	}

	eta = 0.0; 
	for (i=0; i<=N; i++) {
		
		fprintf(f_out, "%2.14lf\n", grid->position_s[i]); 
		eta += dEta; 
	}
	
	fclose(f_out); 
}
/**************************************************************************************/

/* This function sets the firts region next to the boundaries mesh properties */
void GridGen_setup_first_sub_region_properties(GridInfo *grid) {
	
	if (grid->sub_region[0]->type == NONUNIFORM_ONE_SIDED) {
		
		grid->sub_region[0]->S0 = END_SLOPE;	
		grid->sub_region[0]->S1 = END_SLOPE;	
	}

	if (grid->sub_region[0]->type == UNIFORM) {
		
		grid->sub_region[0]->S0 = 1.0;	
		grid->sub_region[0]->S1 = 1.0;	
	}

	if (grid->sub_region[0]->type == NONUNIFORM_TWO_SIDED) {
		
		grid->sub_region[0]->S0 = END_SLOPE;	
		grid->sub_region[0]->S1 = END_SLOPE;	

	}
}
/**********************************************************************************************************/
