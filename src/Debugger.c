#include "definitions.h"
#include "DataTypes.h"
#include "Pressure.h"
#include "Velocity.h"
#include "Conc.h"
#include "Grid.h"
#include "Debugger.h"
#include "Immersed.h"
#include <math.h>

struct field {
    double u, v, w;
};
typedef struct field Field;

/* Just a junk (ordered) velocity field */
void Debugger_set_velocity_data(Velocity *vel, MAC_grid *grid) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    double ***data;
    int ierr;

    ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&data); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {


                data[k][j][i] = 1.0+(double)i*j*k;

            } /* for i */
        } /* for j */
    } /* for k*/

    ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&data); PETScErrAct(ierr);
}
/**************************************************************************************************************/

/* Just a junk (ordered) rhs for the velocity field */
void Debugger_set_velocity_RHS(Velocity *vel, MAC_grid *grid) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    double ***rhs_vec;
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    int ierr;
    int index = 0;
    ierr = DAVecGetArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    double *xq = NULL;
    double *yq = NULL;
    switch (vel->component) {

    case 'u': Grid_get_q_status = Grid_get_u_status; xq = grid->xu; yq = grid->yu; break;
    case 'v': Grid_get_q_status = Grid_get_v_status; xq = grid->xv; yq = grid->yv; break;
    case 'w': Grid_get_q_status = Grid_get_w_status; xq = grid->xw; yq = grid->yw; break;
    } /* switch */

    double a_Lx = 1.0/grid->Lx;
    double a_Ly = 1.0/grid->Ly;

    double coef = (2.0*PI/grid->Lx)*(2.0*PI/grid->Lx);

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (Grid_get_q_status(grid, i, j, k) == FLUID) {

                    double tx = 2.0*PI*xq[i]*a_Lx;
                    double ty = 2.0*PI*yq[j]*a_Ly;

                    if (vel->component == 'u') {

                        rhs_vec[k][j][i] = coef*2.0*(sin(tx)*cos(ty));
                        //rhs_vec[k][j][i] = 2.0*coef*(sin(tx)*sin(ty));

                    } else if (vel->component == 'v') {

                        rhs_vec[k][j][i] = -coef*2.0*(cos(tx)*sin(ty));
                        //rhs_vec[k][j][i] = 2.0*coef*(sin(tx)*sin(ty));

                    } else {

                        rhs_vec[k][j][i] = 0.0;
                    }

                    index++;

                } else {

                    rhs_vec[k][j][i] = 0.0;
                }
                //printf("Debugger.c/ status:%d (i,j,k)=(%d,%d,%d) rhs:%f x:%f y:%f\n", Grid_get_q_status(grid, i, j, k), i, j, k, rhs_vec[k][j][i], xq[i], yq[j]);

            } /* for i */
        } /* for j */
    } /* for k*/

    ierr = DAVecRestoreArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
}
/**************************************************************************************************************/

/* Just a junk (ordered) rhs for the pressure */
void Debugger_set_pressure_RHS(Pressure *p, MAC_grid *grid) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    double ***rhs_vec;
    int ierr;
    int index = 0;

    ierr = DAVecGetArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (Grid_get_c_status(grid, i, j, k) == FLUID) {

                    //rhs_vec[k][j][i] = 1.0+(double)i*j*k;
                    rhs_vec[k][j][i] = (double)index;

                    index++;

                } else {

                    rhs_vec[k][j][i] = 0.0;
                }

            } /* for i */
        } /* for j */
    } /* for k*/

    ierr = DAVecRestoreArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
}
/**************************************************************************************************************/

/* for debugging purposes */
void Debugger_solve_p_Poisson_equation(Pressure *p, MAC_grid *grid) {

    int ierr;
    double norm1A, norm2A, normIA;
    double norm1b, norm2b, normIb;
    double norm1x, norm2x, normIx;
    PetscViewer viewer;

    Debugger_set_pressure_RHS(p, grid) ;
    Pressure_solve(p);

    /* Prints the data to file */
    ierr = PetscViewerASCIIOpen(PCW, "pResults.txt", &viewer); PETScErrAct(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "Pressure LHS matrix\n"); PETScErrAct(ierr);
    ierr = MatView(p->A, viewer); PETScErrAct(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "Pressure RHS\n"); PETScErrAct(ierr);
    ierr = VecView(p->G_b, viewer); PETScErrAct(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "Pressure Solution\n"); PETScErrAct(ierr);
    ierr = VecView(p->G_data, viewer); PETScErrAct(ierr);


    /* Compute norms */
    ierr = MatNorm(p->A, NORM_1, &norm1A); PETScErrAct(ierr);
    ierr = MatNorm(p->A, NORM_FROBENIUS, &norm2A); PETScErrAct(ierr);
    ierr = MatNorm(p->A, NORM_INFINITY, &normIA); PETScErrAct(ierr);

    ierr = VecNorm(p->G_b, NORM_1, &norm1b); PETScErrAct(ierr);
    ierr = VecNorm(p->G_b, NORM_2, &norm2b); PETScErrAct(ierr);
    ierr = VecNorm(p->G_b, NORM_INFINITY, &normIb); PETScErrAct(ierr);

    ierr = VecNorm(p->G_data, NORM_1, &norm1x); PETScErrAct(ierr);
    ierr = VecNorm(p->G_data, NORM_2, &norm2x); PETScErrAct(ierr);
    ierr = VecNorm(p->G_data, NORM_INFINITY, &normIx); PETScErrAct(ierr);

    ierr = PetscPrintf(PCW, "Pressure Poisson equation with Neumann B.C. is solved...\n"); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, "A  : norm1:%f norm2:%f normI:%f\n", norm1A, norm2A, normIA); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, "rhs: norm1:%f norm2:%f normI:%f\n", norm1b, norm2b, normIb); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, "sol: norm1:%f norm2:%f normI:%f\n", norm1x, norm2x, normIx); PETScErrAct(ierr);

    //Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');

    //MatView(p->A, 0);
    //VecView(p->G_b, 0);

    //Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');
    //getchar();

    ierr = PetscViewerDestroy(viewer); PETScErrAct(ierr);
}
/**************************************************************************************************************/

/* for debugging purposes */
void Debugger_solve_vel_Poisson_equation(Velocity *vel, MAC_grid *grid) {

    int ierr;
    double norm1A, norm2A, normIA;
    double norm1b, norm2b, normIb;
    double norm1x, norm2x, normIx;
    PetscViewer viewer;

    Debugger_set_velocity_RHS(vel, grid) ;
    Velocity_solve(vel);

    switch (vel->component) {

    case 'u': ierr = PetscViewerASCIIOpen(PCW, "uResults.txt", &viewer); PETScErrAct(ierr); break;
    case 'v': ierr = PetscViewerASCIIOpen(PCW, "vResults.txt", &viewer); PETScErrAct(ierr); break;
    case 'w': ierr = PetscViewerASCIIOpen(PCW, "wResults.txt", &viewer); PETScErrAct(ierr); break;
    } /* switch */

    /* Prints the data to file */
    ierr = PetscViewerASCIIPrintf(viewer, "%c-Velocity LHS matrix\n", vel->component); PETScErrAct(ierr);
    ierr = MatView(vel->A, viewer); PETScErrAct(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "%c-Velocity RHS\n", vel->component); PETScErrAct(ierr);
    ierr = VecView(vel->G_b, viewer); PETScErrAct(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "%c-Velocity Solution\n", vel->component); PETScErrAct(ierr);
    ierr = VecView(vel->G_data, viewer); PETScErrAct(ierr);


    /* Compute norms */
    ierr = MatNorm(vel->A, NORM_1, &norm1A); PETScErrAct(ierr);
    ierr = MatNorm(vel->A, NORM_FROBENIUS, &norm2A); PETScErrAct(ierr);
    ierr = MatNorm(vel->A, NORM_INFINITY, &normIA); PETScErrAct(ierr);

    ierr = VecNorm(vel->G_b, NORM_1, &norm1b); PETScErrAct(ierr);
    ierr = VecNorm(vel->G_b, NORM_2, &norm2b); PETScErrAct(ierr);
    ierr = VecNorm(vel->G_b, NORM_INFINITY, &normIb); PETScErrAct(ierr);

    ierr = VecNorm(vel->G_data, NORM_1, &norm1x); PETScErrAct(ierr);
    ierr = VecNorm(vel->G_data, NORM_2, &norm2x); PETScErrAct(ierr);
    ierr = VecNorm(vel->G_data, NORM_INFINITY, &normIx); PETScErrAct(ierr);

    ierr = PetscPrintf(PCW, "%c-Velocity Poisson equation with Dirichlet B.C. is solved...\n", vel->component); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, "A  : norm1:%f norm2:%f normI:%f\n", norm1A, norm2A, normIA); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, "rhs: norm1:%f norm2:%f normI:%f\n", norm1b, norm2b, normIb); PETScErrAct(ierr);
    ierr = PetscPrintf(PCW, "sol: norm1:%f norm2:%f normI:%f\n", norm1x, norm2x, normIx); PETScErrAct(ierr);

    //Display_DA_3D_data(vel->G_data, grid, params, "velocity", vel->component);

    //MatView(vel->A, 0);
    //VecView(vel->G_b, 0);

    //Display_DA_3D_data(vel->G_data, grid, params, "velocity", vel->component);
    //getchar();

    ierr = PetscViewerDestroy(viewer); PETScErrAct(ierr);

}
/**************************************************************************************************************/

/* This function checks if the current data array vector q has any NAN value. */
int Debugger_check_nan(Vec q, MAC_grid *grid, Parameters *params, char *name) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    int ierr;
    int found = NO;
    double ***data;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    ierr = DAVecGetArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if ( (isnan(data[k][j][i])) || (isinf(data[k][j][i]) ) ) {

                    printf("Rank:%d !!!Warning!!!... NAN or INF at (i,j,k)=(%d,%d,%d) found in data array %s\n", params->rank, i, j, k, name);
                    found = YES;
                }
                //}
            } /* for i*/
        } /* for j*/
    } /* for k*/

    ierr = DAVecRestoreArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

    return found;
}
/******************************************************************************************************************************************/

void Debugger_set_q_rhs(Velocity *vel, Concentration *c, Pressure *p, MAC_grid *grid, Parameters *params, char which_quantity) {

    double x, y, z;
    double ***rhs_vec;
    int i, j, k;
    int ib_index;
    ImmersedNode *ib_node;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int ierr;

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    Immersed *q_immersed=NULL;
    short int boundary=-1;
    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        ierr = DAVecGetArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed = grid->u_immersed;
        boundary = DIRICHLET;
        break;
    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        ierr = DAVecGetArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed = grid->v_immersed;
        boundary = DIRICHLET;
        break;
    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        ierr = DAVecGetArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed = grid->w_immersed;
        boundary = DIRICHLET;
        break;
    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        ierr = DAVecGetArray(grid->DA_3D, c->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed = grid->c_immersed;
        boundary = NEUMANN;
        break;
    case 'p':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        ierr = DAVecGetArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed = grid->c_immersed;
        boundary = NEUMANN;
        break;

    } /* switch */

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                x = xq[i];
                y = yq[j];
                z = zq[k];

                if (Grid_get_q_status(grid, i, j, k) == FLUID) {

                    rhs_vec[k][j][i] = Debugger_get_analytical_rhs(x, y, z, params, boundary);
                } else {

                    if (Grid_get_q_status(grid, i, j, k) == IMMERSED) {

                        ib_index = Immersed_get_ib_global_index(q_immersed, i, j, k);
                        ib_node  = Immersed_get_ib_node(q_immersed, ib_index);

                        rhs_vec[k][j][i] = ib_node->boundary_coef * Debugger_get_analytical_sol(ib_node->control_point.x, ib_node->control_point.y, ib_node->control_point.z, grid, params, boundary);

                        if (which_quantity == 'c') {

                            rhs_vec[k][j][i] = ib_node->boundary_coef * Debugger_get_analytical_neumann(ib_node->control_point.x, ib_node->control_point.y, ib_node->control_point.z, params, boundary);
                        } /* if */
                    } else {

                        rhs_vec[k][j][i] = 0.0;
                    } /* else */
                } /* else */
                if ( (which_quantity == 'c') && (Grid_get_q_status(grid, i, j, k) == IMMERSED) ){

                    printf("Debugger.c/ status:%d (i,j,k)=(%d,%d,%d) (x,y,z)=(%f,%f,%f) rhs:%f\n", Grid_get_q_status(grid, i, j, k), i, j, k, x, y, z, rhs_vec[k][j][i]);
                } /* if */
            } /* for i*/
        } /* for j*/
    } /* for k*/

    switch (which_quantity) {

    case 'u':
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'v':
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'w':
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'c':
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'p':
        ierr = DAVecRestoreArray(grid->DA_3D, p->G_b, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    } /* switch */
}
/**************************************************************************************************************/

/* This function returns the rhs value for the laplacian of any quantity */
double Debugger_get_analytical_rhs(double x, double y, double z, Parameters *params, short int boundary) {

    double Ly = params->Ly;
    double Lx = params->Lx;
    double Lz = params->Lz;
    double sol = 0.0;

    if (boundary == DIRICHLET) {

        sol = sin(2.0*PI*x/Lx)*sin(2.0*PI*y/Ly)*sin(2.0*PI*z/Lz) * ( (2.0*PI/Lx)*(2.0*PI/Lx) + (2.0*PI/Ly)*(2.0*PI/Ly) + (2.0*PI/Lz)*(2.0*PI/Lz));
    } else {

        sol = -cos(2.0*PI*x/Lx)*cos(2.0*PI*y/Ly)*cos(2.0*PI*z/Lz) * ( (2.0*PI/Lx)*(2.0*PI/Lx) + (2.0*PI/Ly)*(2.0*PI/Ly) + (2.0*PI/Lz)*(2.0*PI/Lz) );

        double Px = (2.0*x*x*x - 3.0*x*x);
        double Py = (2.0*y*y*y - 3.0*y*y);
        double Pz = (2.0*z*z*z - 3.0*z*z);
        sol = +(12.0*x - 6.0)*Py*Pz +(12.0*y - 6.0)*Px*Pz +(12.0*z - 6.0)*Px*Py;

        sol *= -1.0;

        Px = cos(2.0*PI*x/Lx);
        Py = cos(2.0*PI*y/Ly);
        Pz = cos(2.0*PI*z/Lz);

        sol = Px*Py*Pz*( (2.0*PI/Lx)*(2.0*PI/Lx) + (2.0*PI/Ly)*(2.0*PI/Ly) + (2.0*PI/Lz)*(2.0*PI/Lz) );

    } /* else */


    return sol;


}
/**************************************************************************************************************/

/* This function returns the rhs value for analytical function */
double Debugger_get_analytical_sol(double x, double y, double z, MAC_grid *grid, Parameters *params, short int boundary) {

    double Lx = params->Lx;
    double Ly = params->Ly;
    double Lz = params->Lz;
    double sol=0.0;

    /* Compatibility condition imposed at the top-right-front corner */
    double xx = Lx + 0.5*grid->dx;
    double yy = Ly - 0.5*grid->dy;
    double zz = Lz - 0.5*grid->dz;

    double Px = cos(2.0*PI*xx/Lx);
    double Py = cos(2.0*PI*yy/Ly);
    double Pz = cos(2.0*PI*zz/Lz);

    double K2 = -Px*Py*Pz;

    xx = grid->xc[grid->NI-1]+grid->dx;
    yy = grid->yc[grid->NJ-1];
    zz = grid->zc[grid->NK-1];

    //Px = (2.0*xx*xx*xx - 3.0*xx*xx);
    //Py = (2.0*yy*yy*yy - 3.0*yy*yy);
    //Pz = (2.0*zz*zz*zz - 3.0*zz*zz);

    if (boundary == DIRICHLET) {

        sol = sin(2.0*PI*x/Lx)*sin(2.0*PI*y/Ly)*sin(2.0*PI*z/Lz);
    } else {


        //sol = cos(2.0*PI*x/Lx)*cos(2.0*PI*y/Ly)*cos(2.0*PI*z/Lz) + K1;

        //Px = (2.0*x*x*x - 3.0*x*x);
        //Py = (2.0*y*y*y - 3.0*y*y);
        //Pz = (2.0*z*z*z - 3.0*z*z);

        Px = cos(2.0*PI*x/Lx);
        Py = cos(2.0*PI*y/Ly);
        Pz = cos(2.0*PI*z/Lz);

        sol = Px*Py*Pz;// + K2;

    }

    return sol;
}
/**************************************************************************************************************/

/* This function returns the rhs value for grad_phi . n */
double Debugger_get_analytical_neumann(double x, double y, double z, Parameters *params, short int boundary) {

    double Lx = params->Lx;
    double Ly = params->Ly;
    double Lz = params->Lz;
    double sol=0.0;
    double nx, ny, nz;
    double x_cent = Lx/2.0;
    double y_cent = Ly/2.0;
    double z_cent = Lz/2.0;
    double R = 0.25;

    if (boundary == DIRICHLET) {

        sol = 0.0;
    } else {


        /* flat plate inclined */
        ny = 2.0/sqrt(5.0);
        nx = 1.0/sqrt(5.0);
        nz = 0.0/sqrt(2.0);

        /* Cylinder */
        ny = ( y - y_cent ) / R;
        nz = ( z - z_cent ) / R;
        nx = 0.0;

        /* flat plate inclined */
        //ny = 1.0/sqrt(1.0+2*0.2*0.2);
        //nx = 0.2/sqrt(1.0+2*0.2*0.2);
        //nz = 0.2/sqrt(1.0+2*0.2*0.2);

        double Px = cos(2.0*PI*x/Lx);
        double Py = cos(2.0*PI*y/Ly);
        double Pz = cos(2.0*PI*z/Lz);

        //Px = (2.0*x*x*x - 3.0*x*x);
        //Py = (2.0*y*y*y - 3.0*y*y);
        //Pz = (2.0*z*z*z - 3.0*z*z);

        /* flat plate */
        ny = 1.0;
        nx = 0.0;
        nz = 0.0;

        /* Sphere */
        nx = ( x - x_cent ) / R;
        ny = ( y - y_cent ) / R;
        nz = ( z - z_cent ) / R;

        sol = nx*(6.0*x*x - 6.0*x)*Py*Pz + ny*(6.0*y*y - 6.0*y)*Px*Pz + nz*(6.0*z*z - 6.0*z)*Px*Py;
        /*
  sol = nx*Px1*Py*Pz +
        ny*Px*Py1*Pz +
        nz*Px*Py*Pz1;
*/
    }



    return sol;
}
/**************************************************************************************************************/


void Debugger_validate_q(Velocity *vel, Concentration *c, Pressure *p, MAC_grid *grid, Parameters *params, char which_quantity) {

    double x, y, z;
    int i, j, k;
    double E1 = 0.0;
    double E2 = 0.0;
    double E_inf = 0.0;
    double W_E1 = 0.0;
    double W_E2 = 0.0;
    double W_Einf = 0.0;
    double error;
    double an_sol, num_sol;
    double ***rhs_vec;
    double dx = grid->dx;
    double dy = grid->dy;
    double dz = grid->dz;
    int i_max=-1;
    int j_max=-1;
    int k_max=-1;
    int ierr;
    Immersed *q_immersed;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int mpi_flag;


    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    short int boundary=-1;
    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed = grid->u_immersed;
        boundary = DIRICHLET;
        break;
    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed = grid->v_immersed;
        boundary = DIRICHLET;
        break;
    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed = grid->w_immersed;
        boundary = DIRICHLET;
        break;
    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed = grid->c_immersed;
        boundary = NEUMANN;
        break;
    case 'p':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        ierr = DAVecGetArray(grid->DA_3D, p->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed = grid->c_immersed;
        boundary = NEUMANN;
        break;

    } /* switch */

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    double shift = rhs_vec[0][0][0] - Debugger_get_analytical_sol(xq[0], yq[0], zq[0], grid, params, boundary);

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                x = xq[i];
                y = yq[j];
                z = zq[k];

                if (Grid_get_q_status(grid, i, j, k) == FLUID) {

                    //an_sol  = Debugger_get_analytical_sol(x, y, z, grid, params, boundary);
                    //num_sol = rhs_vec[k][j][i] - shift;
                    double tx = 2.0*PI*x/params->Lx;
                    double ty = 2.0*PI*y/params->Ly;
                    if (which_quantity == 'u') {

                        an_sol  = sin(tx)*cos(ty);
                    } else if (which_quantity == 'v') {

                        an_sol  = -cos(tx)*sin(ty);
                    }
                    num_sol = rhs_vec[k][j][i];

                    //printf("Debugger.c/ validate (i,j,k)=(%d,%d,%d) num:%f an:%f x:%f y:%f\n", i, j, k, num_sol, an_sol, x, y, z);

                } else {

                    an_sol  = 0.0;
                    num_sol = 0.0;
                } /* else */

                error = fabs(an_sol - num_sol);

                E1 += error;
                E2 += error*error;

                if (error > E_inf) {
                    i_max = i; j_max = j; k_max = k;
                    E_inf = error;
                }
            } /* for i*/
        } /* for j*/
    } /* for k*/

    E1 *= dx*dy*dz;
    E2 *= dx*dy*dz;

    switch (which_quantity) {

    case 'u':
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'v':
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'w':
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'c':
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    case 'p':
        ierr = DAVecRestoreArray(grid->DA_3D, p->G_data, (void ***)&rhs_vec); PETScErrAct(ierr);
        break;
    } /* switch */
    /* NORM1 */
    mpi_flag = MPI_Allreduce( &E1, &W_E1, 1, MPI_DOUBLE, MPI_SUM, PCW);
    mpi_flag = MPI_Allreduce( &E2, &W_E2, 1, MPI_DOUBLE, MPI_SUM, PCW);
    mpi_flag = MPI_Allreduce( &E_inf, &W_Einf, 1, MPI_DOUBLE, MPI_MAX, PCW);
    W_E2 = sqrt(W_E2);

    PetscPrintf(PCW, "Debugger.c/ Analyzing %c quantity. \n E1:%2.14f E2:%2.14f Einf:%2.14f at (i,j,k)=(%d,%d,%d)\n", which_quantity, W_E1, W_E2, W_Einf, i_max, j_max, k_max);

}

short int Debugger_check_q_status(int i, int j, int k, MAC_grid *grid, Parameters *params, char which_quantity) {

    double Lx = params->Lx;
    double Ly = params->Ly;
    double Lz = params->Lz;
    double x, y, z;

    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    double *xq=NULL;
    double *yq=NULL;
    double *zq=NULL;
    Immersed *q_immersed=NULL;
    short int in_solid=-1;
    switch (which_quantity) {

    case 'u':
        xq = grid->xu;
        yq = grid->yu;
        zq = grid->zu;
        Grid_get_q_status = &Grid_get_u_status;
        q_immersed = grid->u_immersed;
        break;
    case 'v':
        xq = grid->xv;
        yq = grid->yv;
        zq = grid->zv;
        Grid_get_q_status = &Grid_get_v_status;
        q_immersed = grid->v_immersed;
        break;
    case 'w':
        xq = grid->xw;
        yq = grid->yw;
        zq = grid->zw;
        Grid_get_q_status = &Grid_get_w_status;
        q_immersed = grid->w_immersed;
        break;
    case 'c':
        xq = grid->xc;
        yq = grid->yc;
        zq = grid->zc;
        Grid_get_q_status = &Grid_get_c_status;
        q_immersed = grid->c_immersed;
        break;
    } /* switch */


    x = xq[i];
    y = yq[j];
    z = zq[k];
    /* For instance, for a sphere located at (Lx/2, Ly/2, Lz/2) and Radius of 0.5, we have */
    double d = sqrt( (x - 0.5*Lx)*(x - 0.5*Lx) + (y - 0.5*Ly)*(y - 0.5*Ly) + (z - 0.5*Lz)*(z - 0.5*Lz) ) - 0.25;

    if (d < 0.0) {

        in_solid = YES;
    } else {
        in_solid = NO;
    }

    int a= Grid_get_q_status(grid, i, j, k);
    if ( (a == IMMERSED) || (a == SOLID) ){
        a = YES;
    } else {
        a = NO;

    }

    if (a - in_solid != 0) {

        printf("Debugger.c/ Kooni, n_:%d a_:%d diff:%d\n", a, in_solid, (a-in_solid));
        getchar();
    }

    return in_solid;
}

/* This function returns the maximum value of a PETSC vector, processor rank owning it and the index (i,j,k) to file */
void Debugger_q_get_max(Vec q, MAC_grid *grid, Parameters *params, char *name) {


    /* Get the diagonal vector array */
    double ***data=NULL;
    int ierr = DAVecGetArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

    /* Adds constant 'mod' to the diagonal part of matrix A*/
    /* Start index of bottom-left-back corner on current processor */
    int Is = grid->G_Is;
    int Js = grid->G_Js;
    int Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie = grid->G_Ie;
    int Je = grid->G_Je;
    int Ke = grid->G_Ke;

    /* Go through the entire domain and only insert the delta_t term for the nodes which are not immersed nodes */
    int i, j, k;
    double max_q = 0.0;
    int i_m = -1;
    int j_m = -1;
    int k_m = -1;

    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                if (fabs(data[k][j][i]) > max_q) {

                    max_q = fabs(data[k][j][i]);
                    i_m = i;
                    j_m = j;
                    k_m = k;
                }
            }
        }
    }

    double W_max_q = 0.0;
    (void)MPI_Allreduce (&max_q, &W_max_q, 1, MPI_DOUBLE, MPI_MAX, PCW);

    /* current processor owns the max value */
    if (fabs(W_max_q - max_q) < 1e-9) {

        FILE *f = fopen(name, "aa");
        fprintf(f, "Maximum value of %s is: %2.16lf\n", name, W_max_q);
        fprintf(f, "Happening on processor %d at node (i,j,k)=(%d,%d,%d)\n", params->rank, i_m, j_m, k_m);
        fclose(f);
    } /* if */

    ierr = DAVecRestoreArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);


}

void Debugger_testi_dof(MAC_grid *grid) {

    int NX = grid->NX;
    int NY = grid->NY;
    int NZ = grid->NZ;
    int Ghost_Nodes = 1;
    int dof = 3;
    DA _DA_3D;
    Vec G_data, G_temp;

    int ierr = DACreate3d(PCW, DA_NONPERIODIC, DA_STENCIL_STAR, NX, NY, NZ, PETSC_DECIDE, PETSC_DECIDE,
                      PETSC_DECIDE, dof, Ghost_Nodes, PETSC_NULL, PETSC_NULL, PETSC_NULL, &_DA_3D); PETScErrAct(ierr);

    ierr = DACreateGlobalVector(_DA_3D, &G_data); PETScErrAct(ierr);
    ierr = DACreateGlobalVector(_DA_3D, &G_temp); PETScErrAct(ierr);

    int Is, Js, Ks;
    int nx, ny, nz;
    ierr = DAGetCorners(_DA_3D, &Is, &Js, &Ks, &nx, &ny, &nz); PETScErrAct(ierr);

    int Ie=Is+nx;
    int Je=Js+ny;
    int Ke=Ks+nz;

    Field ***data;
    ierr = DAVecGetArray(_DA_3D, G_data, &data); PETScErrAct(ierr);

    int i, j, k;
    for (k=Ks; k<Ke; k++) {
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                data[k][j][i].u = sin(2.0*PI*grid->xc[i])* sin(2.0*PI*grid->yc[j]) * sin(2.0*PI*grid->zc[k]);
                data[k][j][i].v = sin(2.0*PI*grid->xc[i])* sin(2.0*PI*grid->yc[j]) * sin(2.0*PI*grid->zc[k]);
                data[k][j][i].w = sin(2.0*PI*grid->xc[i])* sin(2.0*PI*grid->yc[j]) * sin(2.0*PI*grid->zc[k]);
            }
        }
    }

    ierr = DAVecRestoreArray(_DA_3D, G_data, &data);

    PetscViewer writer;
    char filename[]="test_10000.bin";
    ierr = PetscViewerCreate(PCW, &writer); PETScErrAct(ierr);
    ierr = PetscViewerSetType(writer, PETSC_VIEWER_BINARY); PETScErrAct(ierr);
    ierr = PetscViewerFileSetMode(writer, FILE_MODE_WRITE); PETScErrAct(ierr);
    ierr = PetscViewerBinarySkipInfo(writer); PETScErrAct(ierr); PETScErrAct(ierr);
    ierr = PetscViewerFileSetName(writer, filename); PETScErrAct(ierr);

    ierr = VecView(G_data, writer);
    ierr = PetscViewerDestroy(writer); PETScErrAct(ierr);

    //VecSet(G_data, -0.0);

    PetscViewer reader;
    ierr = PetscViewerCreate(PCW, &reader); PETScErrAct(ierr);
    ierr = PetscViewerSetType(reader, PETSC_VIEWER_BINARY); PETScErrAct(ierr);
    ierr = PetscViewerFileSetMode(reader, FILE_MODE_READ); PETScErrAct(ierr);
    ierr = PetscViewerBinarySkipInfo(reader); PETScErrAct(ierr); PETScErrAct(ierr);
    ierr = PetscViewerFileSetName(reader, filename); PETScErrAct(ierr);

    ierr = VecLoadIntoVector(reader, G_temp);
    ierr = PetscViewerDestroy(reader); PETScErrAct(ierr);

    VecAYPX(G_temp,-1.0,G_data);
    double error[3];// = 1000.0;
    int loc[3];//=-1;
    VecMax(G_temp, loc, error);
    PetscPrintf(PCW, "Debegger.c/ (e,l)=(%f,%d) (%f,%d) (%f,%d)\n", error[0], loc[0], error[1], loc[1], error[2], loc[2]);
    VecMin(G_temp, loc, error);
    PetscPrintf(PCW, "Debegger.c/ (e,l)=(%f,%d) (%f,%d) (%f,%d)\n", error[0], loc[0], error[1], loc[1], error[2], loc[2]);


    VecDestroy(G_data);
    VecDestroy(G_temp);
    DADestroy(_DA_3D);
}
