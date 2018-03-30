/*
 * File: culvert.h
 * Author : Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for culvert.c
 * including all the variables for culverts
 *
 */
#ifndef _culvert_h
#define _culvert_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _culvertT {
REAL *top, // the rigid top for culvert cells, Culvertheight=10000 for non culvert cells Culvertheight[cells]
     *pressure, // store the exact pressure field, transfer results to phys->h before new time step calculation
     *pressure2, // store the outer iteration value for h
     *pressure3, // store the outer iteration value for h
     *pressure4, // store the minimum residual of h
     *Qcoef,   // the coefficients decided by outer loop value for free surface solver 
     *condition, // 1. first store the last guess for h during iteration, 2. store the residual term for converge check
     *condition2, // outer iteration convergence check
     *toparea,  // the wet area for culvert top, no subgrid it will be 0 (no culvert part) or Ac (culvert part), with subgrid it will be 0 (no culvert part) or Aceff(Culverttop)  
     *source1,  // store the initial htmp with iteration
     *source2,  // store hold
     *source3,  // store htmp2
     *source4,  // store htmp3
     *htmp;     // temporary h result during iteration
REAL CdC,eps,sum;  // culvert drag coefficient
int constant; // whether there is constant culvert everywhere or different culvert height,
}culvertT;

culvertT *culvert;

void FreeCulvert(int myproc);
void CheckCulvertCondition(gridT *grid,physT *phys, propT *prop, int myproc);
void CulvertHCoefficients(REAL *coef, REAL *fcoef,gridT *grid, physT *phys, propT *prop, int myproc);
void StoreCulvertPressure(REAL *h, int Nc, int no, int myproc);
void CulvertIterationSource(gridT *grid, physT *phys,  propT *prop, REAL theta, REAL dt, int myproc);
void CulvertInitIteration(gridT *grid, physT *phys, propT *prop, int no, int myproc);
void UpdateCulvertQcoef(gridT *grid, propT *prop,int no, int myproc);
void SetupCulvertmodel(gridT *grid, physT *phys, propT *prop, int myproc);
void SetCulvertDragCoefficient(gridT *grid, physT *phys, int myproc);

#endif
