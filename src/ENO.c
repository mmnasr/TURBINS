#include "definitions.h"
#include "DataTypes.h"
#include "ENO.h"
#include "Memory.h"
#include <stdio.h>
#include <stdlib.h>


/* This function create ENO variable. This is done to avoid runtime memory allocation. */
ENO_Scheme *ENO_create(Parameters *params) {


    ENO_Scheme *new_ENO;
    double *D0=NULL, *D1=NULL, *D2=NULL, *D3=NULL;
    double *dQ1=NULL, *dQ2=NULL, *dQ3=NULL;
    double *VdQdX=NULL;
    double *Position=NULL, *Velocity=NULL;
    int NX, NY, NZ;
    int ENO_order;
    int Nmax, Nmax1;

    new_ENO = (ENO_Scheme *)malloc(sizeof(ENO_Scheme));
    Memory_check_allocation(new_ENO);

    NX = params->NX+1;
    NY = params->NY+1;
    NZ = params->NZ+1;
    ENO_order = params->ENO_scheme_order;

    Nmax1 = max(NX, NY);
    Nmax  = max(Nmax1, NZ);

    Nmax += 2*ENO_order;

    D0  = Memory_allocate_1D_double(Nmax); /* First order divided difference  */
    D1  = Memory_allocate_1D_double(Nmax); /* First order divided difference  */
    dQ1 = Memory_allocate_1D_double(Nmax); /* Derivative of first order ploynomial  */

    if (ENO_order > 1) {

        D2  = Memory_allocate_1D_double(Nmax); /* Second order divided difference */
        dQ2 = Memory_allocate_1D_double(Nmax); /* Derivative of second order ploynomial */
    }

    if (ENO_order > 2) {

        D3  = Memory_allocate_1D_double(Nmax); /* Thirds order divided difference */
        dQ3 = Memory_allocate_1D_double(Nmax); /* Derivative of third order ploynomial  */
    }

    Position = Memory_allocate_1D_double(Nmax); /* Position of each cell */
    Velocity = Memory_allocate_1D_double(Nmax); /* Velocity of each cell */
    VdQdX    = Memory_allocate_1D_double(Nmax); /* Convective term */

    new_ENO->D0 = D0;
    new_ENO->D1 = D1;
    new_ENO->D2 = D2;
    new_ENO->D3 = D3;

    new_ENO->dQ1 = dQ1;
    new_ENO->dQ2 = dQ2;
    new_ENO->dQ3 = dQ3;

    new_ENO->Position = Position;
    new_ENO->Velocity = Velocity;
    new_ENO->VdQdX    = VdQdX;

    new_ENO->ENO_order = ENO_order;
    new_ENO->Nmax      = Nmax;

    return (new_ENO);
}
/***************************************************************************************************/

/* This function releases the memory allocated for the ENO variables */
void ENO_destroy(ENO_Scheme *ENO) {

    free(ENO->D0);
    free(ENO->D1);
    free(ENO->dQ1);

    if (ENO->ENO_order > 1) {

        free(ENO->D2);
        free(ENO->dQ2);
    }

    if (ENO->ENO_order > 2) {

        free(ENO->D3);
        free(ENO->dQ3);
    }

    free(ENO->Position);
    free(ENO->Velocity);
    free(ENO->VdQdX);

    free(ENO);
}
/***************************************************************************************************/

/* This function uses ENO scheme to find convective derivative in the form of udQdx. */
void ENO_compute_convective_terms(ENO_Scheme *ENO, int Nmax_local) {

    double *D0, *D1, *D2, *D3;
    double *dQ1, *dQ2, *dQ3;
    double *Position, *Velocity, *VdQdX;
    double dQ;
    int i;
    int k, kstar;
    double c, cstar;
    double x_i, x_k, x_k1, x_ks1, x_ks2, x_ks;
    int ENO_order;

    D0 = ENO->D0;
    D1 = ENO->D1;
    D2 = ENO->D2;
    D3 = ENO->D3;

    dQ1 = ENO->dQ1;
    dQ2 = ENO->dQ2;
    dQ3 = ENO->dQ3;

    Position = ENO->Position;
    Velocity = ENO->Velocity;
    VdQdX    = ENO->VdQdX;

    ENO_order = ENO->ENO_order;

    /*********************************************/
    /*********************************************/
    /* Start the ENO scheme now */

    /* To reduce the order of ENO scheme close to the boundaries. */
    /* First node next to the boundary: 1st order ENO */
    /* Second node next to the boundary: 2nd order ENO */
    /*
 for (i=0; i<ENO_order+1; i++) {

  right_index = Nmax_local-1-i;

  D1[i]           = 0.0;
  D1[right_index] = 0.0;

  if (ENO_order > 1) {

   D2[i]           = 0.0;
   D2[right_index] = 0.0;
  }

  if (ENO_order > 2) {

   D3[i]           = 0.0;
   D3[right_index] = 0.0;
  }
 }
*/
    /*first divided difference*/
    for (i = 0; i < Nmax_local-1; i++) {
        //for (i = ENO_order-1; i < Nmax_local-ENO_order; i++) {

        D1[i] = (D0[i+1] - D0[i]) / (Position[i+1] - Position[i]);
    }

    if (ENO_order > 1) {
        /*second divided difference*/
        for (i = 1; i < Nmax_local-1; i++) {
            //for (i = ENO_order; i < Nmax_local-ENO_order; i++) {

            D2[i] = ( D1[i] - D1[i-1]) / (Position[i+1] - Position[i-1]);
        }
    }

    if (ENO_order > 2) {
        /*third divided difference*/
        for (i = 1; i < Nmax_local-2; i++) {
            //for (i = ENO_order; i < Nmax_local-1-ENO_order; i++) {

            D3[i] = ( D2[i+1] - D2[i]) / (Position[i+2] - Position[i-1]);
        }
    }

    /*find convective term VdQdX */
    for (i=ENO_order; i<Nmax_local-ENO_order; i++) {

        if ( Velocity[i] > 0.0 ) {

            k = i-1;
        }
        else {
            k = i;
        }

        dQ1[i] = D1[k];
        dQ     = dQ1[i];

        //printf("i: %d   D0(i): %2.10f  D0(i+1):%2.10f \n", i, D0[i], D0[i+1]);
        //printf("i: %d   dQ1: %2.10f: \n", i, dQ1[i]);
        //getchar();
        if (ENO_order > 1) {

            if ( fabs( D2[k] ) <= fabs( D2[k+1] ) ){

                c = D2[k];
                kstar = k-1;
            } else {
                c = D2[k+1];
                kstar = k;
            }

            x_i    = Position[i  ];
            x_k    = Position[k  ];
            x_k1   = Position[k+1];
            dQ2[i] = c * ( (x_i - x_k) + (x_i - x_k1) );

            dQ    += dQ2[i];
            //printf("i: %d   dQ2: %2.10f: \n", i, dQ2[i]);
            //getchar();
        }

        if (ENO_order > 2) {

            if ( fabs( D3[kstar] ) <= fabs( D3[kstar+1] )) {

                cstar = D3[ kstar ];
            } else {

                cstar = D3[kstar+1];
            }

            x_i    = Position[i      ];
            x_ks   = Position[kstar  ];
            x_ks1  = Position[kstar+1];
            x_ks2  = Position[kstar+2];
            dQ3[i] = cstar * ( (x_i-x_ks) * (x_i - x_ks1)  +
                               (x_i-x_ks1)* (x_i - x_ks2)  +
                               (x_i-x_ks) * (x_i - x_ks2) );
            dQ += dQ3[i];

        }

        VdQdX[i] = Velocity[i] * dQ;
        //printf("Index: %d, k: %d, kstar:%d, dQ1:%2.16f dQ2:%2.16f dQ3:%2.16f Vel:%2.16f VdQdX: %2.16f\n", i, k, kstar, dQ1[i], dQ2[i], dQ3[i], Velocity[i], VdQdX[i]);
    } /* for i*/
}
/***************************************************************************************************/

/* This function uses centrail differnce compute convective derivative in the form of udQdx. */
/* Use it only for u, v, w */
void ENO_compute_convective_terms_central(ENO_Scheme *ENO, int Nmax_local) {

        double *D0 = ENO->D0;
        double *Position = ENO->Position;
        double *Velocity = ENO->Velocity;
        double *VdQdX    = ENO->VdQdX;

        int ENO_order = ENO->ENO_order;
        int i;
        for (i=ENO_order; i<Nmax_local-ENO_order; i++) {

                VdQdX[i] = (Velocity[i+1]*D0[i+1]-Velocity[i-1]*D0[i-1]) / (Position[i+1] - Position[i-1]);
        } /* for i */
}
