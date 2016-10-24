#ifndef DATATYPES_H
#define DATATYPES_H

#include "petsc.h"
#include "petscksp.h"
#include "petscda.h"

struct concentration {

    int NX, NY, NZ;  /*Global number of cells in x, y and z*/
    int NT; /* Total number of conc nodes */

    /* Number of cells excluding the half cell added to the boundaries */
    int NI, NJ, NK;

    double Pe;/* Corresponding Peclet number */
    double Phi; /* Volume fraction */
    double v_settl0; /* Constant settling speed */
    double W_total_susp_mass; /* Total suspended mass within the fluid region */
    double W_front_location; /* x-location of gravity current front */

    Vec G_v_particle; /* Total Particle velocity in y-direction */
    Vec L_v_particle; /* Total Particle velocity in y-direction */

    /* Local values */
    Vec L_data, L_data_old;

    /* Global values */
    Vec G_data, G_data_old;

    /* Global and local total concentration fields (it is generated only if NConc > 1 and stored in c[0] */
    Vec G_c_total, L_c_total;

    /* convective terms: uddx+vddy+wddz */
    Vec G_conv;


    /* Settling Speed y-direction */
    Vec G_v_settl;
    Vec L_v_settl;

    /* Index of the start and end corner on the current processor */
    int G_Is, G_Js, G_Ks;
    int G_Ie, G_Je, G_Ke;

    /* Index of the start and end corner on the current processor including ghost nodes */
    int L_Is, L_Js, L_Ks;
    int L_Ie, L_Je, L_Ke;

    /* kv_G kinematic viscosity function (( nu = nu_0 * kv_G ))*/
    double ***kv_G; /* the value at the cell center.*/
    double ***kv_G_bl; /* the value at the cell bottom-left. */
    double **G_inflow;
    double **G_outflow, **G_outflow_old;
    double **E_s;

    double *G_ave_height_x; /* (concentration) Averaged flow height. Only a function of x */
    double *W_ave_height_x; /* World ave_height */
    double **G_deposit_height; /* deposit height for each concentration field on current processor */
    double **W_deposit_height; /* World, deposit height for each concentration field */

    double **G_deposit_height_dumped; /* deposit height for each concentration field on current processor (after dumping whatever is left in the fluid region) */
    double **W_deposit_height_dumped; /* World, deposit height for each concentration field (after dumping whatever is left in the fluid region) */

    /* instantaneous sedimentation rate on the bottom boundary */
    double **W_sed_rate_bottom;
    double **G_sed_rate_bottom;

    /* local and world stokes dissipation rate */
    double G_Stokes_dissipation_rate;
    double W_Stokes_dissipation_rate;

    /* kinetic to potential energy conversion rate */
    double G_kinetic_potential_conversion_rate, W_kinetic_potential_conversion_rate;

    /* Active and passive potential energies */
    double G_Ep, G_Ep_passive;
    double W_Ep, W_Ep_passive;

    /* potential energy rate of change */
    double G_Ep_rate, W_Ep_rate;

    /* passive potential energy rate of change */
    double G_Ep_passive_rate, W_Ep_passive_rate;

    /* sedimentation rate at any given time from the bottom geometry */
    double W_sed_rate;

    int Type; /* PARTICLE, SALINITY or TEMPERATURE */

    int conc_index; /* for more than 1 concentration field, conc_index shows the index of the field respectively. */
    /* 1D arrays used for computing Convective terms using ENO scheme */
    /* Linear system matirces: [A]{x} = {b}* and Solver */
    Vec G_b, G_x;
    Vec L_b, L_x;
    Mat A;
    Vec *lsys_diagonal; /* array used to update the values of the diagonal part of the matrix A */
    KSP solver;
    PC pc;
};
typedef struct concentration Concentration;
/******************************************************************************************/

/* This structure defines the levelset function */
struct levelset_type {

    int NX, NY, NZ;
    int NI, NJ, NK;

    /* phi(x,y,z) */
    Vec G_data, L_data;
    /* phi0(x,y,z) */
    Vec G_data_0, L_data_0;
    /* phi(x,y,z) at previous time step */
    Vec G_data_old, L_data_old;

    /* dphi_dx, dphi_dy, dphi_dz */
    Vec G_Dx, G_Dy, G_Dz;
    Vec L_Dx, L_Dy, L_Dz;
}; 
typedef struct levelset_type LevelsetType; 
/******************************************************************************************/

/* surface type holds the signed distance functions at u,v,w and c grid locations */
struct surface_type {

    int NX, NY, NZ;

    /* Global and local versions of the signed distance function stored at u,v,w and c grid */
    Vec G_c_sdf, L_c_sdf;
    Vec G_u_sdf, L_u_sdf;
    Vec G_v_sdf, L_v_sdf;
    Vec G_w_sdf, L_w_sdf;

    /* To save time, have the 3D pointers to the L_q_sdf ready */
    double ***u_sdf, ***v_sdf, ***w_sdf, ***c_sdf;

    LevelsetType *sdf_object;
};
typedef struct surface_type SurfaceType; 
/******************************************************************************************/

struct pressure {

    int NX, NY, NZ, NT;
    int NI, NJ, NK;

    /* Index of the start and end corner on the current processor */
    int G_Is, G_Js, G_Ks;
    int G_Ie, G_Je, G_Ke;

    /* Index of the start and end corner on the current processor including ghost nodes */
    int L_Is, L_Js, L_Ks;
    int L_Ie, L_Je, L_Ke;

    Vec G_data, L_data;
    Vec G_divergence, L_divergence;

    /* Pressure gradients, To be used in momentum equations */
    Vec G_dp_dx, L_dp_dx;
    Vec G_dp_dy, L_dp_dy;
    Vec G_dp_dz, L_dp_dz;

    Vec G_dphi_dx, L_dphi_dx;
    Vec G_dphi_dy, L_dphi_dy;
    Vec G_dphi_dz, L_dphi_dz;

    Vec G_b, L_b;
    Mat A;
    KSP solver;
    PC pc;

};
typedef struct pressure Pressure;
/******************************************************************************************/

struct parameters {

    int NX, NY, NZ;
    double Re;
    double *Pe;
    double *Phi;
    double *Conc_alpha; /* Contribution of each conc field in the buoyancy term in v-mom equation */
    double *V_s0; // Settling velocity of partciles: y-direction
    int *conc_type; /* Type of the concentration field: SALINITY, PARTICLE, TEMPERATURE */
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double Lx, Ly, Lz;
    int NConc; /* Number of concentration fields */

    double resmax;
    double time_max;
    int maxit;

    double output_time_step;

    double x_fr, y_fr, z_fr;
    int lock_smooth;
    double smooth_length_x, smooth_length_y, smooth_length_z;

    double resume_saving_timestep;
    short int output_format; /* Output format indicator */
    short int conc_output, stream_output, vor_output, u_output, v_output, w_output, p_output;
    short int div_output;
    short int ave_height_output, susp_mass_output;
    short int deposit_height_output;
    short int symmetry_output, front_location_output, front_speed_output;
    short int energies_output;
    short int shear_stress_output;
    short int sedim_rate_output;
    short int dump_conc;

    short int cfl_method;
    short int hindered_settling;
    short int resuspension;

    double default_dt;
    double alpha, U;

    short int sim_type, outflow, outflow_type, inflow, transient_inflow, viscosity_type;
    short int sedimentation;

    short int UniformGridX, UniformGridY, UniformGridZ;
    short int ImportGridFromFile;
    short int ImportBottomInterface;
    short int constant_geometry; /* Indicates whether the bottom geometry is constant or changing due to erosion and deposition */
    short int resume; /* Flag which indicates if the simulation starts from t=0 or from a previous run */
    double time_c; // Characteristic time for exp. decay in transient inflow

    int ENO_scheme_order;

    /* rank and size of the processors */
    int rank;
    int size;

    /* BINARY, ASCII, HDF5 (for 3D data) */
    short int output_type;

    /* Flag indicating if the linear-system matrices are written to file when generated at the begining */
    short int output_lsys_matrix;
};
typedef struct parameters Parameters; 
/******************************************************************************************/

struct velocity {

    /* Which component of the velocity */
    /* 'u', 'v' or 'w' */
    char component;

    /* Global Number of velocity nodes in x, y and z direction */
    int NX, NY, NZ;

    /* 0Number of nodes excluding the half cell added to the boundaries */
    int NI, NJ, NK;

    /* Global Total number of points */
    int NT;

    /* Index of the start and end corner on the current processor */
    int G_Is, G_Js, G_Ks;
    int G_Ie, G_Je, G_Ke;

    /* Index of the start and end corner on the current processor including ghost nodes */
    int L_Is, L_Js, L_Ks;
    int L_Ie, L_Je, L_Ke;

    /*
        - data    : value of the velocity
        - data_old: value of the velocity of the old time step
        - data_bc : value of the velocity at the cell_center


*/
    /* Local values */
    Vec L_data, L_data_old, L_data_bc;

    /* Global values */
    Vec G_data, G_data_old, G_data_bc;

    /* Cell centerd velocity, the ghost nodes are based on a BOX stencil. This is going to be used in computing shear stress */
    Vec Gb_data_bc, Lb_data_bc;

    /* Local quantities */
    Vec *L_u_transposed; /*hold u_velocity component in this velocity's coordinates */
    Vec *L_v_transposed; /*hold v_velocity component in this velocity's coordinates */
    Vec *L_w_transposed; /*hold w_velocity component in this velocity's coordinates */

    /* Global quantities */
    Vec *G_u_transposed; /*hold u_velocity component in this velocity's coordinates */
    Vec *G_v_transposed; /*hold v_velocity component in this velocity's coordinates */
    Vec *G_w_transposed; /*hold w_velocity component in this velocity's coordinates */

    /* convective terms: uddx+vddy+wddz */
    Vec G_conv;


    /* dummy Local vector based on a box stencil */
    Vec L_data_box;

    /* 2D arrays which hold outflow data */
    double **G_outflow, **G_outflow_old;

    /* 2D array holds the inflow */
    double **G_inflow;

    /* Total influx to the region */
    double G_influx;

    /* (World) influx to the region (same on all processors) */
    double W_influx;

    /* Global and World kinetic energy of the fluid */
    double G_kinetic_energy, W_kinetic_energy;
    /* Global and World viscous dissipation rate */
    double G_dissipation_rate, W_dissipation_rate;

    /* Kinetic to potential energy conversion rate, -\int v*c dVol */
    double G_kinetic_potential_conversion_rate, W_kinetic_potential_conversion_rate;

    /* Total shear stress on the bottom boundary. Store them only in 'u' velocity */
    double **G_shear_stress_bottom, **W_shear_stress_bottom;

    /* Tangential velocity on the bottom boundary */
    double **G_u_shear, **G_v_shear, **G_w_shear;
    double **W_u_shear, **W_v_shear, **W_w_shear;

    long int lsys_count; /*total number of points in the linear system (number of fluid points)*/

    /* Linear solver matirces: [A]{x} = {b}* and Solver */
    /* 1D Distributed array used for the vectors of the lsys */
    Vec G_b, G_x;
    Vec L_b, L_x;
    Vec *lsys_diagonal; /* array used to update the values of the diagonal part of the matrix A */
    Mat A;
    KSP solver;
    PC pc;

};
typedef struct velocity Velocity;
/******************************************************************************************/

/* Eno scheme variables. To be used for computing convective terms */
struct eno_scheme {

    int Nmax; /* Maximum number of nodes */
    int ENO_order;
    /* Used for the nodes that have not enough fluid nodes for the boudary-ghost nodes */
    int Reduced_ENO_order;
    double *D0, *D1, *D2, *D3;
    double *dQ1, *dQ2, *dQ3;
    double *Position, *Velocity;
    double *VdQdX;
    /* reference pointer address. This is done because of the shifting of the D0, Position and Velocity
 arrays on the boundaries */
    double *D0_base;         /* D0[0] */
    double *Position_base;   /* Position[0] */
    double *Velocity_base;   /* Velocity[0] */
};
typedef struct eno_scheme ENO_Scheme;
/******************************************************************************************/

/* S = S(X1, X2, X3) */
struct scalar {

    double X1, X2, X3;
    double Value;
};
typedef struct scalar Scalar;
/******************************************************************************************/

/* Mesh info */
struct meshinfo {

    short int type;
    char quantity;
    char coordinate;
    int N, NI;
    double lmin, lmax;
    double S0, S1;
    double B, A, Q;
    short int StretchingFunctionID;

};
typedef struct meshinfo MeshInfo;
/******************************************************************************************/

/* To keep the 3 index corresponding to variable (u, v, w, c or p) node */
struct indices {

    int x_index;
    int y_index;
    int z_index;
};
typedef struct indices Indices; 
/******************************************************************************************/

struct datawriter {

    short int Output_Type; /* ASCII or BINARY */
    short int Append_Data;
    /* YES: Append to the end of existing file */
    /* NO : Overwrite */
    short int Write_Grid; /* YES: add the grid coordinates to the beginning of the file */
    int Dimension; /* 1,2 or 3 */

    /* sizeof(float) or sizeof(double) */
    size_t precision_size;
    void *output_data; /* In binary case, we always change the 2-3 D array to a 1D array */

    void *X1_grid; /* direction 1 grid position : x */
    void *X2_grid; /* direction 2 grid position : y */
    void *X3_grid; /* direction 3 grid position : z */

    int N1, N2, N3; /* number of data in N1, N2 and N3 direction: NX, NY, NZ */
    int NT; /* total number of points to be writen to file */
    char *Filename; /* file name */
    FILE *file;    /* pointer to the opened file */
};
typedef struct datawriter Writer; 
/******************************************************************************************/

struct outputwriter {

    DA DA_3D;
    Vec G_data, L_data;

    /* Index of the start and end corner on the current processor */
    int G_Is, G_Js, G_Ks;
    int G_Ie, G_Je, G_Ke;

    /* Index of the start and end corner on the current processor including ghost nodes */
    int L_Is, L_Js, L_Ks;
    int L_Ie, L_Je, L_Ke;

    int Ghost_Nodes; /* Number of output ghost nodes. Warning: Different than ENO_order */
    Writer *writer;

    char base_dir[MAX_PATH_LEN]; /* Name of the base directory to-be-created to store the output files */
    char bin_dir[MAX_PATH_LEN];  /* Name of the binary directory to-be-created to store the output files in binary format */
    char data_dir[MAX_PATH_LEN]; /* Name of the data directory to-be-created to store the output files */
    char hdf5_dir[MAX_PATH_LEN]; /* Name of the data directory to-be-created to store the output files in hdf5 format*/
    char immersed_dir[MAX_PATH_LEN]; /* Name of the data directory to-be-created to store immersed information */

    char data_3d_dir[MAX_PATH_LEN]; /* Name of directory to-be-created to store the output (3d) files in hdf5, binary or ascii format */
    char data_3d_file_ext[4];        /* extension of the output file: "dat" for ASCII, "bin" for BINARY, and "h5" for HDF5 */

    int precision; /* the type of output: GVG_FLOAT or GVG_DOUBLE */

    short int type; /* BINARY, ASCII, or HDF5 */
    short int update_last_plane; /* YES or NO. Copy the last interior values into the half cell added at the end */

    void *L_data_1D; /* pointer to the 1D data (float or double) */
    float *L_data_1D_float;  /* memory allocated for 1D array of data (float precision) */
    double *L_data_1D_double;/* memory allocated for 1D array of data (double precision) */

    /* pointers to the local grid coordinates (to be saved either in float or double precision) */
    void *xu, *yu, *zu;
    void *xv, *yv, *zv;
    void *xw, *yw, *zw;
    void *xc, *yc, *zc;

    float *xu_float, *yu_float, *zu_float;
    float *xv_float, *yv_float, *zv_float;
    float *xw_float, *yw_float, *zw_float;
    float *xc_float, *yc_float, *zc_float;

    double *xu_double, *yu_double, *zu_double;
    double *xv_double, *yv_double, *zv_double;
    double *xw_double, *yw_double, *zw_double;
    double *xc_double, *yc_double, *zc_double;

};
typedef struct outputwriter Output;
/******************************************************************************************/

struct resume_reader_writer {

    char **IO_filename; /* Array of filenames corresponding to parallel vector data */
    int IO_type; /* ASCII or BINARY. Use BINARY for now.  */
    int max_vec; /* maximum number of parallel 3D vectors */
    short int start_read;  /* start reading flag */
    short int start_write; /* start writing flag */
    short int IO_grid_flag; /* flag indicating wheter grid_c is stored/read from file. This is useful for
 the simulation where the geometry could change due to erosion or deposition */
    Vec **Vec_data;  /* pointer to vectors of data which have to be written (or read from file) to file (parallel vectors) */

};
typedef struct resume_reader_writer Resume; 
/******************************************************************************************/

/* vector (vx,vy,vz); */
struct vector_type {

    double vx, vy, vz;

};
typedef struct vector_type VectorType;  
/******************************************************************************************/

/* point (x,y,z); */
struct point_type {

    double x, y, z;

};
typedef struct point_type PointType;  
/******************************************************************************************/

/* This sttucture holds the information used for one immersed node */
struct immersed_node {
    /* IBM_MAX: maximum number of "potential" neighboring nodes for each immersed node */
    /* IBM_MAX: 9 EAST, WEST, SOUTH, NORTH, BACK, FRONT, and OTHER1, OTHER2, OTHER3. */
    PointType fluid_point[IBM_MAX];   /* Fluid points chosen for the interpolation */
    PointType control_point;  /* control point (located on the boundary used to enforce the correct boundary conditions */
    PointType im_point;       /* coordinates of the immersed node itself. For the current version, the first node inside the solid region right next to the solid boundary */

    VectorType n;             /* Unit vector normal to the boundary at control point */
    PointType image_point;    /* coordinates of the  image node (inside the fluid region). The mirrored ghost node relative to the solid boundary */

    /* Coordinates of the  image node. The node inside the fluid region mirrored relative to the interface */
    PointType image_point1;
    PointType image_point2;

    Indices neighbor_index[IBM_MAX];  /* indices of the neighboring nodes */
    PointType neighbor_point[IBM_MAX];/* neighboring points coordinates */
    int n_fluid;                      /* Number of the neighboring fluid nodes */
    double fluid_coef[IBM_MAX];       /* Coefficients used to relate the value of the IB node to the neiboring fluid nodes satisfying the correct boundary conditions */
    Indices fluid_index[IBM_MAX];     /* (i,j,k) indices of the chosen fluid nodes for interpolation */
    Indices im_index;                 /* (i,j,k) index of the immersed node */
    short int boundary_condition; /* Type of the boundary condition to be imposed: NEUMANN, DIRICHLET or ROBBIN */
    double boundary_coef;         /* used to impose any nonzero boundary conditions by taking the nonzero term to the rhs of the lsys*/
}; 
typedef struct immersed_node ImmersedNode; 
/******************************************************************************************/

/* Strucutre holding the information for all the immersed node for each quantity */
struct immersed  {

    char quantity;               /* Which flow quantity 'u', 'v', 'w' or 'c' */
    int N;                       /* Total number of the ib_nodes on the current processor */
    /* Array holding the number of IB node in the 3D physical domain */
    /* Note, for non-IMMERSED nodes, we store '-1' and for Immersed nodes, that goes from 0 to N-1 on each processor (locally) */
    Vec G_global_index, L_global_index;
    double ***global_index;      /* To save time, we pass the pointer the global_index to avoid calling some PETSc routines */
    ImmersedNode *ib_nodes;      /* Array of all the ib_nodes */
};
typedef struct immersed Immersed; 
/******************************************************************************************/

/* Structure holding the grid information */
struct mac_grid {

    /* Domain length */
    double Lx, Ly, Lz;

    double *xc, *yc, *zc;
    double *xu, *yu, *zu;
    double *xv, *yv, *zv;
    double *xw, *yw, *zw;
    double *metric_xc, *metric_yc, *metric_zc;
    double *metric_xu, *metric_yv, *metric_zw;
    double *dx_c, *dy_c, *dz_c; /* dx, dy, dz of each cell */
    double dx, dy, dz;
    double dEta, dXi, dPhi;
    int NX, NY, NZ, NT;
    int NI, NJ, NK;
    double **interface_position;
    int   **interface_y_index;

    /* Index of the start and end corner on the current processor */
    int G_Is, G_Js, G_Ks;
    int G_Ie, G_Je, G_Ke;

    /* Index of the start and end corner on the current processor including ghost nodes */
    int L_Is, L_Js, L_Ks;
    int L_Ie, L_Je, L_Ke;

    DA DA_3D; /* Distribued array (parallel) data layout between processors. This takes care of the communications to update the values of the ghost nodes */
    /* Stencil type: STAR (along the grid ghost nodes) [see PETSc reference manual]*/

    DA DA_3D_BOX; /* Distribued array (parallel) data layout between processors. This takes care of the communications to update the values of the ghost nodes */
    /* Stencil type: BOX (along the grid lines and diagonal ghost nodes) [see PETSc reference manual]*/

    DA DA_3D_Lsys; /* Distribued array (parallel) data layout between processors. This takes care of the communications to update the values of the ghost nodes */
    /* Stencil type: STAR (along the grid lines and diagonal ghost nodes) [see PETSc reference manual]*/
    /* This one is solely used to create the lhs matrix and rhs vectors associated with the linear systems for the quantities */


    /* Arrays to store the condition of a grid */
    /* Center of grid: Concentration and Pressure */
    Vec G_c_status, L_c_status;
    /* u-grid */
    Vec G_u_status, L_u_status;
    /* v-grid */
    Vec G_v_status, L_v_status;
    /* w-grid */
    Vec G_w_status, L_w_status;

    /* pointers to the local Vec's of q_status. That is to save time */
    double ***u_status;
    double ***v_status;
    double ***w_status;
    double ***c_status;

    /* Surface type which holds all the data regarding the location of the solid inteface */
    SurfaceType *surf;

    /* flow field quantities holding the required information for the adaptation of the immerse boundary method */
    Immersed *u_immersed, *v_immersed, *w_immersed, *c_immersed;

    Vec lsys_diagonal; /* array used to update the values of the diagonal part of the matrix A */

/* Dummy vectors. They can be used anywhere in the code. Just pass the right pointer address to these vectors */
    Vec G_dummy1, G_dummy2;
    Vec L_dummy1, L_dummy2;

};
typedef struct mac_grid MAC_grid;
/******************************************************************************************/

/* Now, define a structure which holds the pointer to all other primitive variables */
struct gvg_bag {

    Velocity *u, *v, *w;
    Pressure *p;
    Concentration **c;
    MAC_grid *grid;
    Parameters *params;

    double dt, dt_old;
    double time;
    int iter;
    int output_iter;
};
typedef struct gvg_bag GVG_bag;
/******************************************************************************************/

struct matrix_struct {

    /* Max = 10 */
    int size; /* size < Max */
    double det; /* determinent */
    double A[10][10];  /* matrix values */
    double A_inv[10][10]; /* inverse matrix values */

};
typedef struct matrix_struct MatrixType;
/******************************************************************************************/

#endif
