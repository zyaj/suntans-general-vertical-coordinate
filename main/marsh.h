/*
 * File: marsh.h
 * Author : Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for marsh.c
 * including all the variables for marsh model
 *
 */
#ifndef _marsh_h
#define _marsh_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _marshT {
REAL cdV, Alphav, Hmarsh, Na, Rv; // marsh basic parameters

int IntCdV, Inthmarsh; // switch for marsh model IntCdV, Inthmarsh: 1 read file interpolate 2 read edgecenter data 3 constant 0 use user define function in initialization.c

REAL *CdV, // drag coefficient for each cell [cell]
     *hmarshleft, // the marsh top cell hmarsh left [cell]
     *hmarsh, // the total height for each cell [cell]
     *hmarshcenter,
     *CdVcenter;

int *marshtop; // the upmost layer for with hmarshleft!=0 [cell] 
} marshT;

marshT *marsh;

void SetMarshTop(gridT *grid, physT *phys, int myproc);
void MarshExplicitTerm(gridT *grid, physT *phys, propT *prop, int edgep, REAL theta, REAL dt, int myproc);
double MarshImplicitTerm(gridT *grid, physT *phys, propT *prop, int edgep, int layer, REAL theta, REAL dt, int myproc);
void SetupMarshmodel(gridT *grid,physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm);
void FreeMarsh(gridT *grid, int myproc);
#endif

