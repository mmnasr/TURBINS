#include "definitions.h"
#include "DataTypes.h"
#include "Outflow.h"
#include "Velocity.h"
#include "Conc.h"
#include "Grid.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* This function imposes the convective B.C. on all velocitites and concentration fields */
void Outflow_impose_convective_boundary( Velocity *u, Velocity *v, Velocity *w, Concentration **c, MAC_grid *grid, Parameters *params, double dt) {

    int NConc, iconc;

    NConc = params->NConc;

    Outflow_vel_impose_convective_boundary(u, grid, params, dt);
    Outflow_vel_impose_convective_boundary(v, grid, params, dt);
    Outflow_vel_impose_convective_boundary(w, grid, params, dt);

    for (iconc=0; iconc<NConc; iconc++) {

        Outflow_conc_impose_convective_boundary(c[iconc], grid, params, dt);
    } /* for */
}
/********************************************************************************/

/* This function imposes the convective boundary condition at the outflow. It also checks if the total mass is conserved */
void Outflow_vel_impose_convective_boundary( Velocity *vel, MAC_grid *grid, Parameters *params, double dt) {

    int ierr;
    int Js, Ks;
    int Je, Ke;
    int j, k;
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    double ***vel_data;
    double **outflow;
    double U, coef, dx;
    int NI;
    int last_i_index;
    short int Should_Conserve_Mass = YES;


    /* Update only on the processors which have the last yz plane */
    if (grid->G_Ie == grid->NX) {

        ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&vel_data);PETScErrAct(ierr);
        outflow = vel->G_outflow;

        /* First, assign the value for the convective velocity */
        /* Either constant or the max velocity in the domain */
        U  = params->U;
        NI = grid->NI; /* Number of grids excluding the half cell */
        dx = (grid->xu[NI] - grid->xu[NI-1]);

        coef =  dt*U/dx;

        /* i(x)index of the last plane before the outflow region */
        last_i_index = NI-1;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        switch (vel->component) {

        case 'u': Grid_get_q_status = &Grid_get_u_status; break;
        case 'v': Grid_get_q_status = &Grid_get_v_status; break;
        case 'w': Grid_get_q_status = &Grid_get_w_status; break;
        } /* switch */

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                /* Check if the node to the left of the outflow node is fluid */
                int status = Grid_get_q_status(grid, last_i_index, j, k);
                if ( ( status == FLUID) || (status == IMMERSED) ){

                    outflow[k][j] = (1.0 - coef) * outflow[k][j] + coef * vel_data[k][j][last_i_index];
                } /* if fluid */

            } /* for j*/
        } /* for k */

        /* restore the velocity data array */
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&vel_data);PETScErrAct(ierr);

    } /* if */

    /* Now, update u to conserve global mass. Call this on all processors since we need the influx too. */
    if ( (vel->component == 'u') && (Should_Conserve_Mass) ) {

        Outflow_update_u_velocity_to_conserve_mass(vel, grid) ;
    } /* if conserve mass */

}
/********************************************************************************/

/* This function imposes the convective boundary condition at the outflow for 1 concentration field */
void Outflow_conc_impose_convective_boundary(Concentration *c, MAC_grid *grid, Parameters *params, double dt) {

    int ierr;
    int Js, Ks;
    int Je, Ke;
    int j, k;
    double ***conc_data;
    double **outflow;
    double U, dx, coef;
    int NI;
    int last_i_index;


    /* Update only on the processors which have the last yz plane */
    if (grid->G_Ie == grid->NX) {

        ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc_data);PETScErrAct(ierr);
        outflow = c->G_outflow;

        /* First, assign the value for the convective velocity */
        /* Either constant or the max velocity in the domain */
        U  = params->U;
        NI = grid->NI; /* Number of grids excluding the half cell */
        dx = (grid->xu[NI] - grid->xu[NI-1]);

        coef =  dt*U/dx;

        /* i(x)index of the last plane before the outflow region */
        last_i_index = NI-1;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                /* Check if the node to the left of the outflow node is fluid */
                int status = Grid_get_c_status(grid, last_i_index, j, k);
                if ( (status == FLUID) || (status == IMMERSED) ){

                    outflow[k][j] = (1.0 - coef) * outflow[k][j] + coef * conc_data[k][j][last_i_index];
                } /* if !solid */

            } /* for j*/
        } /* for k */

        /* restore the velocity data array */
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc_data);PETScErrAct(ierr);
    } /* if */
}
/********************************************************************************/

/* This function checks if total mass is conserved in the whole domain (u-inflow*Area)in - u_outflow*Area)out = 0.0 */
void Outflow_update_u_velocity_to_conserve_mass(Velocity *u, MAC_grid *grid) {

    int j, k;
    int NI=grid->NI;
    double **outflow=NULL;
    double *yv=NULL;
    double *zw=NULL;
    double total_u_outflux;
    double outflow_area;
    double W_outflow_area;
    double W_total_u_outflux;
    double W_total_u_influx;
    double velocity_correction;
    double dy, dz;
    double dA;
    double send_data[2], W_sum_data[2];
    double temp_flux = 0.0;
    double ***data;
    int ierr;

    ierr = DAVecGetArray(grid->DA_3D, u->G_data, (void ***)&data);PETScErrAct(ierr);

    /* Update only on the processors which have the last yz plane */
    /* First, integrate over the outflow area to find the u_outflux */
    total_u_outflux = 0.0;
    outflow_area    = 0.0;

    if (grid->G_Ie == grid->NX) {

        outflow = u->G_outflow;

        yv = grid->yv;
        zw = grid->zw;

        /* i(x)index of the last plane before the outflow region */
        int last_i_index = NI-1;

        /* start index on current processor */
        int Js = grid->G_Js;
        int Ks = grid->G_Ks;

        /* end index on current processor */
        int Je = grid->G_Je;
        int Ke = grid->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                /* Check if the node to the left of the outflow node is fluid */
                int status = Grid_get_u_status(grid, last_i_index, j, k);
                if ( (status == FLUID) || (status == IMMERSED) ) {

                    dy = yv[j+1] - yv[j];
                    dz = zw[k+1] - zw[k];


                    dA            = dy*dz;
                    /* Total outflow area. This is done if there is a geometry (solid) right before the outflow boundary */
                    outflow_area += dA;

                    /* integral over each mesh: u*(dy*dz) */
                    total_u_outflux += dA*outflow[k][j];

                    /* one row before */
                    temp_flux += dA*data[k][j][last_i_index];
                } /* if */
            } /* for j */
        } /* for k*/

        //printf("Outflow.c/ here total outflux:%f\n", total_u_outflux);
    } /* if */

    /* Send two variables to all processors, sum them up and send them back in W_recv_data array */
    send_data[0] = total_u_outflux;
    send_data[1] = outflow_area;

    /* Now, add all the total_u_outflux and outflow area from all processors and get the sum back on all processors */
    /* Better way is to use only processors which have the last yz plane */
    (void)MPI_Allreduce ( (void *)send_data, (void *)W_sum_data, 2, MPI_DOUBLE, MPI_SUM, PCW);


    //printf("Recv: Total outflux:%f Total area:%f \n", W_sum_data[0], W_sum_data[1]);
    //getchar();

    /* Update only on the processors which have the last yz plane */
    /* Now, add or substract the outflow u-velocity correction */
    if (grid->G_Ie == grid->NX) {

        /* (World) Total influx, outflux and outflow area on all processors */
        W_total_u_influx  = u->W_influx;
        W_total_u_outflux = W_sum_data[0];
        W_outflow_area    = W_sum_data[1];


        //printf("Total outflux:%f Total influx:%f Total area:%f \n", W_total_u_outflux, W_total_u_influx, W_outflow_area);
        //getchar();
        /* Now, calculate the correction amount to be added (substracted) to all the u-nodes at the outflow region */
        velocity_correction = (W_total_u_influx - W_total_u_outflux) / W_outflow_area;
        //printf("Outlfow.c/ last_row_flux:%f outflux:%f area:%f correction:%f\n", temp_flux, W_total_u_outflux, W_outflow_area, velocity_correction);
        //getchar();

        /* i(x)index of the last plane before the outflow region */
        int last_i_index = NI-1;

        /* start index on current processor */
        int Js = grid->G_Js;
        int Ks = grid->G_Ks;

        /* end index on current processor */
        int Je = grid->G_Je;
        int Ke = grid->G_Ke;

        double flux_new = 0.0;
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                /* Check if the node to the left of the outflow node is fluid */
                int status = Grid_get_u_status(grid, last_i_index, j, k);
                if ( (status == FLUID) || (status == IMMERSED) ) {

                    //double temp = outflow[k][j];
                    dy = yv[j+1] - yv[j];
                    dz = zw[k+1] - zw[k];

                    dA            = dy*dz;
                    outflow[k][j]   += velocity_correction;
                    flux_new += outflow[k][j]*dA;
                    //printf("Outflow.c/ (i,j,k)=(%d,%d,%d) u-out before:%f after correction:%f flux_old:%2.12f flux_new:%2.12f\n", last_i_index, j, k, temp, outflow[k][j], W_total_u_outflux, flux_new);
                    //getchar();
                    //after_u_flux      += outflow[k][j]*dA;

                } /* if */
            } /* for j */
        } /* for k*/
        //getchar();
    } /* if */

    ierr = DAVecRestoreArray(grid->DA_3D, u->G_data, (void ***)&data);PETScErrAct(ierr);
    //printf("\n Netflux before: %.4e, Netflux after correction: %.4e\n", (flux_correction*Outflow_Area), (total_u_influx - after_u_flux));
}
/********************************************************************************/

/* This function sets/resets data in the last plane of the 3D data array */
void Outflow_vel_copy_data_into_last_plane(Velocity *vel, MAC_grid *grid, int should_set) {

    int j, k;
    int ierr;
    int Js, Ks;
    int Je, Ke;
    int NI;
    double ***vel_data;
    double **outflow;

    /* Copy only on the processors which have the last yz plane */
    if (grid->G_Ie == grid->NX) {

        ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&vel_data);PETScErrAct(ierr);
        outflow = vel->G_outflow;

        /* i(x)index of the last plane */
        NI = grid->NI;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        if (should_set) {

            for (k=Ks; k<Ke; k++) {
                for (j=Js; j<Je; j++) {

                    vel_data[k][j][NI] = outflow[k][j];

                } /* for j*/
            } /* for k */
        } else {

            for (k=Ks; k<Ke; k++) {
                for (j=Js; j<Je; j++) {

                    vel_data[k][j][NI] = 0.0;

                } /* for j*/
            } /* for k */
        }

        /* restore the velocity data array */
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&vel_data);PETScErrAct(ierr);

    } /* if */

}
/********************************************************************************/

/* This function sets/resets data in the last plane of the 3D data array of concentration */
void Outflow_concentration_copy_data_into_last_plane(Concentration *c, MAC_grid *grid, int should_set) {

    int j, k;
    int ierr;
    int Js, Ks;
    int Je, Ke;
    int NI;
    double ***conc_data;
    double **outflow;

    /* Copy only on the processors which have the last yz plane */
    if (grid->G_Ie == grid->NX) {

        ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc_data);PETScErrAct(ierr);
        outflow = c->G_outflow;

        /* i(x)index of the last plane */
        NI = grid->NI;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        if (should_set) {

            for (k=Ks; k<Ke; k++) {
                for (j=Js; j<Je; j++) {

                    conc_data[k][j][NI] = outflow[k][j];

                } /* for j*/
            } /* for k */
        } else { /* reset, i.e. to zero */

            for (k=Ks; k<Ke; k++) {
                for (j=Js; j<Je; j++) {

                    conc_data[k][j][NI] = 0.0;

                } /* for j*/
            } /* for k */
        }

        /* restore the velocity data array */
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc_data);PETScErrAct(ierr);

    } /* if */

}
/********************************************************************************/

/* This function sets/resets data in the last plane of the 3D data array */
void Outflow_vel_copy_data_into_2D_array(Velocity *vel, MAC_grid *grid) {

    int j, k;
    int ierr;
    int Js, Ks;
    int Je, Ke;
    int NI;
    double ***vel_data;
    double **outflow;

    /* Copy only on the processors which have the last yz plane */
    if (grid->G_Ie == grid->NX) {

        ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&vel_data);PETScErrAct(ierr);
        outflow = vel->G_outflow;

        /* i(x)index of the last plane */
        NI = grid->NI;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                outflow[k][j] = vel_data[k][j][NI];

            } /* for j*/
        } /* for k */

        /* restore the velocity data array */
        ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&vel_data);PETScErrAct(ierr);

    } /* if */

}
/********************************************************************************/

/* This function sets/resets data in the last plane of the 3D data array of concentration */
void Outflow_concentration_copy_data_into_2D_array(Concentration *c, MAC_grid *grid) {

    int j, k;
    int ierr;
    int Js, Ks;
    int Je, Ke;
    int NI;
    double ***conc_data;
    double **outflow;

    /* Copy only on the processors which have the last yz plane */
    if (grid->G_Ie == grid->NX) {

        ierr = DAVecGetArray(grid->DA_3D, c->G_data, (void ***)&conc_data);PETScErrAct(ierr);
        outflow = c->G_outflow;

        /* i(x)index of the last plane */
        NI = grid->NI;

        /* start index on current processor */
        Js = grid->G_Js;
        Ks = grid->G_Ks;

        /* end index on current processor */
        Je = grid->G_Je;
        Ke = grid->G_Ke;

        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                outflow[k][j] = conc_data[k][j][NI];

            } /* for j*/
        } /* for k */

        /* restore the velocity data array */
        ierr = DAVecRestoreArray(grid->DA_3D, c->G_data, (void ***)&conc_data);PETScErrAct(ierr);

    } /* if */

}

