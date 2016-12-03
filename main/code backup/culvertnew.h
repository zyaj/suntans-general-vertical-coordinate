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

REAL *Culverttop, // the rigid top for culvert cells, Culvertheight=10000 for non culvert cells Culvertheight[cells]
     *Culvertpressure, // store the exact pressure field, transfer results to phys->h before new time step calculation
     *Culvertcondition, // 1. first store the last guess for h during iteration, 2. store the residual term for converge check
     *Culvertsource1,  //store the initial htmp with iteration
     *Culvertsource2,  // store hold
     *Culvertsource3,  // store htmp2
     *Culvertsource4;  // store htmp3

REAL Cdculvert,Culverteps,Culvertsum;  // culvert drag coefficient

int ConstantCulvert; // whether there is constant culvert everywhere or different culvert height,


void ReadCulvertProperties(int myproc);
void AllocateCulvert(gridT *grid,int myproc);
void InitializeCulvert(gridT *grid, physT *phys,int myproc);
void FreeCulvert(int myproc);
void CheckCulvertCondition(REAL *h, gridT *grid, int myproc);
void CulvertHCoefficients(REAL *coef, REAL *fcoef,gridT *grid, physT *phys, propT *prop, int myproc);
void StoreCulvertPressure(REAL *h, int Nc, int no, int myproc);
void CulvertIterationSource(gridT *grid, physT *phys,  propT *prop, REAL theta, REAL dt, int myproc);
void CulvertInitIteration(gridT *grid, physT *phys, propT *prop, int no, int myproc);
#endif
