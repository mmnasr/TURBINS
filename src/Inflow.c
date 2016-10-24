#include "definitions.h"
#include "DataTypes.h"
#include "Inflow.h"
#include "Velocity.h"
#include "Conc.h"
#include "Grid.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* This function defines the inflow profile in two cases. Constant or time dependent inflow profile */
void Inflow_u_velocity_profile (Velocity *u, MAC_grid *grid, Parameters *params, double time) {

    int j, k;
    int Js, Ks;
    int Je, Ke;
    int first_i_index;
    short int status;
    double **u_inflow;
    double *yv, *yc, *zw, *zc, *yu, *zu;
    double Ly, Lz;
    double dy, dz;
    double dA;
    double influx, W_influx;
    double u_0, y_0, y_1, A, B, pi, Re;
    double time_factor;
    double S_Re, y_solid, y_top, d;
    double z_max, z_min, y_max;

    /* influx on the current processor */
    influx = 0.0;

    /* Update only on the processors which have the first yz plane */
    if (grid->G_Is == 0) {


        pi = 4.0*atan(1.0);

        Ly = params->Ly;
        Lz = params->Lz;

        yc = grid->yc;
        zc = grid->zc;
        yu = grid->yu;
        zu = grid->zu;
        yv = grid->yv;
        zw = grid->zw;

        Re = params->Re;

        /* Parameters defining the inflow profile. */
        u_0 = 0.625;
        y_0 = 0.7*Ly;
        y_1 = 0.8*(0.5*Ly);
        B   = pi / (2.0*y_0);
        A   = u_0 / sin(B*y_0);


        /* solid start at y=y_solid */
        /* Smooth increase from u=0 to u u_0. And then again from u_0 to zero (on the top wall) */
        u_0     = 1.0;
        y_solid = 0.5*params->Ly;
        y_top   = 0.85*params->Ly;
        d       = params->Ly - y_solid;
        S_Re    = sqrt(4.0*Re);

        u_inflow = u->G_inflow;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        first_i_index = 1;

        /* If we only have on slot of inflow */
        z_min = 0.3*params->Lz;
        z_max = 0.7*params->Lz;
        y_max = 0.5*params->Ly;

        for (k=Ks; k<Ke; k++) {

            for (j=Js; j<Je; j++) {

                status = Grid_get_u_status(grid, first_i_index, j, k);
                if (status == FLUID) {


                    if ( (yu[j] <= y_max) && (zu[k] >= z_min) && (zu[k] < z_max) ) {

                        u_inflow[k][j] = 1.0;
                    } else {

                        u_inflow[k][j] = 0.0;
                    }
                    u_inflow[k][j] = 1.0;//0.25;
                    //u_max = 1.0;
                    //yy = yu[j]/params->Ly;
                    //u_inflow[k][j] = 4.0*u_max*(yy - yy*yy);


                    /* */
                    /*
     if ( yu[j] < y_solid ) {

      u_inflow[k][j] = 0.0;
     }
     else if ( (yu[j] >= y_solid) && (yu[j] < y_top) ){

      u_inflow[k][j] = u_0*( erf( (yc[j+1] - y_solid)*S_Re) );;

     } else if (yu[j] >= y_top) {

      u_inflow[k][j] = u_0*( erf( (d - yc[j+1] + y_solid)*S_Re) );
     }
*/

                    dy  = yv[j+1] - yv[j];
                    dz  = zw[k+1] - zw[k];
                    dA  = dy*dz;


                    if (params->transient_inflow) {

                        /* Temporal factor : u(t) = u0 * (1.0 - exp (-t/t_c) ) */
                        time_factor = (1.0 - exp ( -(time + 0.01) / (params->time_c) ) );
                    } else {

                        time_factor = 1.0;
                    }

                    /* Multuply the profile by a time factor */
                    u_inflow[k][j] *= time_factor;

                    /* Calculate influx, i.e. u)*Area */
                    influx       += u_inflow[k][j]*dA;
                    //printf("inflow.c/ (j,k)=(%d,%d) inflow:%f dA:%f\n", j, k, u_inflow[k][j], dA);
                } /* if */
            } /* for j*/
        } /* for k*/
    }  /* if */

    /* Now, send the local influx part and sum them up to get the total influx: W_influx */
    (void)MPI_Allreduce ( (void *)&influx, (void *)&W_influx, 1, MPI_DOUBLE, MPI_SUM, PCW);

    /* Local part of influx */
    u->G_influx = influx;

    /* Total influx (send to all processors) */
    u->W_influx = W_influx;

    printf("inflow.c/ l_influx:%f w_influx:%f\n", influx, W_influx);
    //getchar();
}
/********************************************************************************/

/* This function defines the inflow profile for the concentration field */
void Inflow_conc_profile(Concentration *c, MAC_grid *grid, Parameters *params, double time) {

    int j, k;
    int Js, Ks;
    int Je, Ke;
    int first_i_index;
    double **c_inflow;
    double *yu, *yv, *yc, *zw, *zc;
    double Re;
    double A, B;
    double u_0;
    double Ly, Lz;
    double y_0, y_1;
    double time_factor;
    double pi;
    double z_max, z_min, y_max;

    /* Update only on the processors which have the yz plane at x=0*/
    if (grid->G_Is == 0) {

        pi = 4.0*atan(1.0);

        Ly = params->Ly;
        Lz = params->Lz;

        yc = grid->yc;
        zc = grid->zc;
        yu = grid->yu;
        yv = grid->yv;
        zw = grid->zw;

        Re  = params->Re;

        /* Parameters defining the inflow profile. */
        u_0 = 0.625;
        y_0 = 0.7*Ly;
        y_1 = 0.8*(0.5*Ly);
        B = pi / (2.0*y_0);
        A = u_0 / sin(B*y_0);

        c_inflow = c->G_inflow;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        first_i_index = 0;

        /* If we only have on slot of inflow */
        z_min = 0.3*params->Lz;
        z_max = 0.7*params->Lz;
        y_max = 0.5*params->Ly;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                if (Grid_get_c_status(grid, first_i_index, j, k) != SOLID) {


                    if (yc[j] < 0.5*Ly) {

                        c_inflow[k][j]      = 1.0;//0.5*(1.0 + erf( (y_1 - yc[j])*sqrt(Re)) );
                    }


                    if ( (yc[j] <= y_max) && (zc[k] >= z_min) && (zc[k] < z_max) ) {

                        c_inflow[k][j] = 1.0;
                    } else {

                        c_inflow[k][j] = 0.0;
                    }


                    c_inflow[k][j] = 1.0;
                    if (params->transient_inflow) {

                        /* Temporal factor : c(t) = c0 * (1.0 - exp (-t/t_c) ) */
                        time_factor = (1.0 - exp ( -(time + 0.01) / (params->time_c) ) );
                    } else {

                        time_factor = 1.0;
                    }



                    /* Multiply the profile by a time factor */
                    c_inflow[k][j] *= time_factor;

                } /* if */
            } /* for j*/
        } /* for k*/
    } /* if */
}

