#include "definitions.h"
#include "DataTypes.h"
#include "Grid.h"
#include "Display.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void Display_3D_variable(double ***var, int NX, int NY, int NZ, const char *name) {

    int i, j, k;
    char message[50] = "\n3D data for variable: ";

    strcat(message, name);
    puts(message);

    for (k=0; k<NZ; k++) {
        printf ("z:k=%d (y:j,x:i) \n\n", k);
        for (j=0; j<NY; j++) {
            for (i=0; i<NX; i++) {

                printf("q(%d,%d)=%f  ", j, i, var[k][j][i]);
            }
            printf("\n");
        }
        printf("********************************\n");
    }
    printf("**************************************************************\n");
    printf("**************************************************************\n\n");
}
/***************************************************************************************************/

void Display_3D_variable_point_based(double ***var, int NX, int NY, int NZ, const char *name) {

    int i, j, k;
    char message[50] = "\n3D data for variable: ";

    strcat(message, name);
    puts(message);

    printf("q(z,y,x)\n");
    for (k=0; k<NZ; k++) {
        for (j=0; j<NY; j++) {
            for (i=0; i<NX; i++) {

                printf("q(%d,%d,%d)=%f  ", k, j, i, var[k][j][i]);
            }
            printf("\n");
        }
        printf("********************************\n");
        getchar();
    }
    printf("**************************************************************\n");
    printf("**************************************************************\n\n");
}
/***************************************************************************************************/


/* This function prints a 3D Distributed array vector data between processors */
void Display_DA_3D_data(Vec q, MAC_grid *grid, Parameters *params, char *q_name, char which_quantity) {

    int Is, Js, Ks;
    int Ie, Je, Ke;
    int i, j, k;
    double ***data;
    int ierr;

    /* Pointer to the function */
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    switch (which_quantity) {

    case 'u':

        Grid_get_q_status = Grid_get_u_status;
        break;
    case 'v':

        Grid_get_q_status = Grid_get_v_status;
        break;
    case 'w':

        Grid_get_q_status = Grid_get_w_status;
        break;

    case 'c':

        Grid_get_q_status = Grid_get_c_status;
        break;

    case 'p':

        Grid_get_q_status = Grid_get_c_status;
        break;
    } /* swicth */

    ierr = DAVecGetArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    PetscSynchronizedPrintf(PCW, "***********************************\n");
    PetscSynchronizedPrintf(PCW, "Rank:%d Printing \"%s\"\n", params->rank, q_name);
    PetscSynchronizedPrintf(PCW, "StartIndex(k:%d,j:%d,i:%d) EndIndex(k:%d,j:%d,i:%d)\n", Ks, Js, Is, Ke, Je, Ie);

    for (k=Ks; k<Ke; k++) {

        PetscSynchronizedPrintf (PCW, "z:k=%d (y:j,x:i) \n\n", k);
        for (j=Js; j<Je; j++) {
            for (i=Is; i<Ie; i++) {

                //if (!quantity_is_solid(grid, i, j, k)) { /*only include point if it is fluid*/

                PetscSynchronizedPrintf(PCW, "q(%d,%d)=%1.6f  ", j, i, data[k][j][i]);
                //}
            } /* for i*/
            PetscSynchronizedPrintf(PCW,"\n");
        } /* for j*/
    } /* for k*/
    PetscSynchronizedFlush(PCW);

    ierr = DAVecRestoreArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

}
/***************************************************************************************************/


/* This function prints a 3D Distributed array vector data between processors (including ghost nodes) */
void Display_DA_3D_L_data(Vec q, MAC_grid *grid, Parameters *params, char *q_name, char which_quantity) {

    int i, j, k;
    int ierr;

    /* Pointer to the function */
    int (*Grid_get_q_status)(MAC_grid *, int, int, int)=NULL;
    switch (which_quantity) {

    case 'u':

        Grid_get_q_status = &Grid_get_u_status;
        break;
    case 'v':

        Grid_get_q_status = &Grid_get_v_status;
        break;
    case 'w':

        Grid_get_q_status = &Grid_get_w_status;
        break;

    case 'c':

        Grid_get_q_status = &Grid_get_c_status;
        break;

    case 'p':

        Grid_get_q_status = &Grid_get_c_status;
        break;
    } /* swicth */

    double ***data=NULL;
    ierr = DAVecGetArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

    /* Start index of bottom-left-back corner on current processor */
    int Is_g = grid->L_Is;
    int Js_g = grid->L_Js;
    int Ks_g = grid->L_Ks;

    /* End index of top-right-front corner on current processor */
    int Ie_g = grid->L_Ie;
    int Je_g = grid->L_Je;
    int Ke_g = grid->L_Ke;

    PetscSynchronizedPrintf(PCW, "***********************************\n");
    PetscSynchronizedPrintf(PCW, "Rank:%d Printing \"%s\"\n", params->rank, q_name);
    PetscSynchronizedPrintf(PCW, "StartIndex(k:%d,j:%d,i:%d) EndIndex(k:%d,j:%d,i:%d)\n", Ks_g, Js_g, Is_g, Ke_g, Je_g, Ie_g);

    for (k=Ks_g; k<Ke_g; k++) {

        PetscSynchronizedPrintf (PCW, "z:k=%d (y:j,x:i) \n\n", k);
        for (j=Js_g; j<Je_g; j++) {
            for (i=Is_g; i<Ie_g; i++) {

                //if (!quantity_is_solid(grid, i, j, k)) { /*only include point if it is fluid*/

                PetscSynchronizedPrintf(PCW, "q(%d,%d)=%1.6f  ", j, i, data[k][j][i]);
                //}
            } /* for i*/
            PetscSynchronizedPrintf(PCW,"\n");
        } /* for j*/
    } /* for k*/
    PetscSynchronizedFlush(PCW);

    ierr = DAVecRestoreArray(grid->DA_3D, q, (void ***)&data); PETScErrAct(ierr);

}
/***************************************************************************************************/


void Display_DA_3D_info(MAC_grid *grid, Parameters *params) {


    int ierr;
    int Is_g, Js_g, Ks_g;
    int Ie_g, Je_g, Ke_g;
    int Is, Js, Ks;
    int Ie, Je, Ke;
    int rank;

    /* Start index of bottom-left-back corner on current processor including ghost nodes */
    Is_g = grid->L_Is;
    Js_g = grid->L_Js;
    Ks_g = grid->L_Ks;

    /* End index of top-right-front corner on current processor including ghost nodes */
    Ie_g = grid->L_Ie;
    Je_g = grid->L_Je;
    Ke_g = grid->L_Ke;

    /* Start index of bottom-left-back corner on current processor */
    Is = grid->G_Is;
    Js = grid->G_Js;
    Ks = grid->G_Ks;

    /* End index of top-right-front corner on current processor */
    Ie = grid->G_Ie;
    Je = grid->G_Je;
    Ke = grid->G_Ke;

    rank = params->rank;

    ierr = PetscSynchronizedPrintf(PCW, "\nProcessor Rank:%d Is:%d Js:%d Ks:%d\n", params->rank, Is, Js, Ks);
    PETScErrAct(ierr);
    ierr = PetscSynchronizedPrintf(PCW, "Processor Rank:%d Ie:%d Je:%d Ke:%d\n", params->rank, Ie, Je, Ke);
    PETScErrAct(ierr);
    ierr = PetscSynchronizedPrintf(PCW, "Processor Rank:%d Is_g:%d Js_g:%d Ks_g:%d\n", params->rank, Is_g, Js_g, Ks_g);
    PETScErrAct(ierr);
    ierr = PetscSynchronizedPrintf(PCW, "Processor Rank:%d Ie_g:%d Je_g:%d Ke_g:%d\n\n", params->rank, Ie_g, Je_g, Ke_g);
    PETScErrAct(ierr);

    ierr = PetscSynchronizedFlush(PCW); PETScErrAct(ierr);
}	
/***************************************************************************************************/

void Display_2D_outflow(double **outflow, MAC_grid *grid, Parameters *params, const char *q_name) {


    if (grid->G_Ie == grid->NX) {


        /* start index on current processor */
        int Js = grid->G_Js;
        int Ks = grid->G_Ks;

        /* end index on current processor */
        int Je = grid->G_Je;
        int Ke = grid->G_Ke;


        PetscSynchronizedPrintf(PCW, "***********************************\n");
        PetscSynchronizedPrintf(PCW, "Rank:%d Printing \"%s\"\n", params->rank, q_name);
        PetscSynchronizedPrintf(PCW, "StartIndex(k:%d,j:%d) EndIndex(k:%d,j:%d)\n", Ks, Js, Ke, Je);
        PetscSynchronizedPrintf(PCW, "q(j:y,k:z)\n");

        int j, k;
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                PetscSynchronizedPrintf(PCW, "q(%d,%d)=%f ", j, k, outflow[k][j]);

            } /* for j */
            PetscSynchronizedPrintf(PCW, "\n");
        } /* for k*/
        PetscSynchronizedFlush(PCW);
    }
}
/***************************************************************************************************/

void Display_2D_inflow(double **inflow, MAC_grid *grid, Parameters *params, const char *q_name) {


    if (grid->G_Is == 0) {

        /* start index on current processor */
        int Js = grid->G_Js;
        int Ks = grid->G_Ks;

        /* end index on current processor */
        int Je = grid->G_Je;
        int Ke = grid->G_Ke;

        PetscSynchronizedPrintf(PCW, "***********************************\n");
        PetscSynchronizedPrintf(PCW, "Rank:%d Printing \"%s\"\n", params->rank, q_name);
        PetscSynchronizedPrintf(PCW, "StartIndex(k:%d,j:%d) EndIndex(k:%d,j:%d)\n", Ks, Js, Ke, Je);
        PetscSynchronizedPrintf(PCW, "q(j:y,k:z)\n");

        int j, k;
        for (k=Ks; k<Ke; k++) {
            for (j=Js; j<Je; j++) {

                PetscSynchronizedPrintf(PCW, "q(%d,%d)=%f ", j, k, inflow[k][j]);

            } /* for j */
            PetscSynchronizedPrintf(PCW, "\n");
        } /* for k*/
        PetscSynchronizedFlush(PCW);
    }
}
/***************************************************************************************************/

/* This function displays the parameters */
void Display_parameters(Parameters *params) {


    int iconc, NConc;


    printf("********************************************************\n");
    printf("********************************************************\n");
    printf("Printing parameters on processor %d\n", params->rank);

    printf("NX:%d NY:%d NZ:%d\n", params->NX, params->NY, params->NZ);
    printf("xmin:%f xmax:%f\n", params->xmin, params->xmax);
    printf("ymin:%f ymax:%f\n", params->ymin, params->ymax);
    printf("zmin:%f zmax:%f\n\n", params->zmin, params->zmax);


    Display_flag(params->UniformGridX, "Uniform Grid X");
    Display_flag(params->UniformGridY, "Uniform Grid Y");
    Display_flag(params->UniformGridZ, "Uniform Grid Z");
    Display_flag(params->ImportGridFromFile, "Import Grid From File");

    printf("\nTime max:%f\n", params->time_max);
    printf("Output time step:%f\n", params->output_time_step);
    printf("Resume time step:%f\n", params->resume_saving_timestep);


    printf("\nRe:%f\n\n", params->Re);
    NConc = params->NConc;
    for (iconc=0; iconc<NConc; iconc++) {

        printf("iconc:%d\nPe:%f V_s0:%f Conc_alpha:%f Phi:%f\n", iconc, params->Pe[iconc], params->V_s0[iconc], params->Conc_alpha[iconc], params->Phi[iconc]);
    } /* iconc */


    printf("\nLock front x:%f\n", params->x_fr);
    printf("Lock front y:%f\n", params->y_fr);
    printf("Lock front z:%f\n", params->z_fr);

    if (params->viscosity_type == VARIABLE_VISCOSITY) {

        printf("\nVariable viscosity flow.\n");
    } else {

        printf("\nConstant viscosity flow.\n");
    }

    Display_flag(params->sedimentation, "sedimentation");
    Display_flag(params->hindered_settling, "hindered_settling");
    Display_flag(params->resuspension, "resuspension");
    Display_flag(params->inflow, "inflow");
    Display_flag(params->outflow, "outflow");

    /***********************************************************************/
    /*          Flags indicating flow properties to be saved to file       */
    /***********************************************************************/
    printf("\n");
    Display_flag(params->conc_output, "conc output");
    Display_flag(params->u_output, "u output");
    Display_flag(params->v_output, "v output");
    Display_flag(params->w_output, "w output");
    Display_flag(params->p_output, "p output");
    Display_flag(params->ave_height_output, "conc average height output");
    Display_flag(params->front_location_output, "front location output");
    Display_flag(params->susp_mass_output, "suspended mass output");
    Display_flag(params->deposit_height_output, "deposit height output");

    printf("********************************************************\n");
    printf("********************************************************\n");
    printf("\n");

}
/***************************************************************************************************/

/* This function printd the status of a flag */
void Display_flag(int flag, const char *name) {

    if (flag) {

        printf("%s: YES\n", name);
    } else {

        printf("%s: NO\n", name);
    }
}
/***************************************************************************************************/

void Display_point(PointType *p, const char *name) {

    PetscPrintf(PCW, "Printing point:%s ", name);
    PetscPrintf(PCW, "P(x,y,z)=(%f,%f,%f)\n", p->x, p->y, p->z);
}
/***************************************************************************************************/

/* This function displays the information for the immersed node */
void Display_immersed_node(ImmersedNode *ib_node) {

    char s[100];
    int g;
    PetscPrintf(PCW, "Printing information for the immersed node\n");
    Display_point(&ib_node->im_point, "immersed");
    Display_point(&ib_node->control_point, "control point");
    Display_point(&ib_node->image_point, "image point");
    for (g=0; g<ib_node->n_fluid; g++) {

        sprintf(s, "Fluid_%d", g);
        Display_point(&ib_node->fluid_point[g], s);
        printf("<--coef:%lf:", ib_node->fluid_coef[g]);
    } /* for g*/

    PetscPrintf(PCW, "Number of Fluid: %d\n*********************************\n", ib_node->n_fluid);

}
/***************************************************************************************************/

/* Displays any grid info */
void Display_grid(double *data, int N, char *name) {

    int i;
    printf("Display.c/ *************************************/ \n");
    for (i=0; i<N; i++) {

        printf("Display.c/ gridinfo-%s(%d)=%f\n", name, i, data[i]);
    } /* for i */
}
/***************************************************************************************************/

/* Displays basic information for a matrix (PETSc) */
void Display_matrix_info(Mat A, const char *name) {

    MatInfo info;
    MatGetInfo(A, MAT_GLOBAL_SUM, &info);

    int mal  = (int) info.mallocs;
    int memory = (int) info.memory;
    int nz_allocated = (int) info.nz_allocated;
    int nz_used = (int) info.nz_allocated;
    int nz_unneeded = (int) info.nz_unneeded;

    PetscPrintf(PCW, "Display.c/ printing information for the %s matrix. \n", name);
    PetscPrintf(PCW, "Display.c/ allocated non-zeros:%d  used non-zeros:%d unneeded non-zeros:%d\n", nz_allocated, nz_used, nz_unneeded);
    PetscPrintf(PCW, "Display.c/ memory (bytes):%d  number of extra mallocs:%d \n", memory, mal);

}
/***************************************************************************************************/

#ifdef MEMORY_PROFILING
/*  This function prints allocated memory for the primitive variables. It is not 100% accurate */
void Display_memory_info(void) {

    double tot = u_mem + v_mem + w_mem + c_mem + p_mem;
    PetscPrintf(PCW, "Display.c/ ************************************************************\n");
    PetscPrintf(PCW, "Display.c/ ************************************************************\n");
    PetscPrintf(PCW, "Display.c/ **************Memory allocated******************************\n");
    PetscPrintf(PCW, "Display.c/ ************************************************************\n");
    PetscPrintf(PCW, "Display.c/ grid allocated (bytes): %d\n", grid_mem);
    PetscPrintf(PCW, "Display.c/ u-allocated (bytes): %d (percentage:%f\%)\n", u_mem, u_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/    u-solver (bytes): %d (percentage:%f\%)\n", u_solver_mem, u_solver_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/ v-allocated (bytes): %d (percentage:%f\%)\n", v_mem, v_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/    v-solver (bytes): %d (percentage:%f\%)\n", v_solver_mem, v_solver_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/ w-allocated (bytes): %d (percentage:%f\%)\n", w_mem, w_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/    w-solver (bytes): %d (percentage:%f\%)\n", w_solver_mem, w_solver_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/ p-allocated (bytes): %d (percentage:%f\%)\n", p_mem, p_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/    p-solver (bytes): %d (percentage:%f\%)\n", p_solver_mem, p_solver_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/ c-allocated (bytes): %d (percentage:%f\%)\n", c_mem, c_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/    c-solver (bytes): %d (percentage:%f\%)\n", c_solver_mem, c_solver_mem/tot*100.0);
    PetscPrintf(PCW, "Display.c/ ************************************************************\n");
    PetscPrintf(PCW, "Display.c/ ************************************************************\n");

}
#endif
