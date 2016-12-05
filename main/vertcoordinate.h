/*
 * File: vertcoordinate.h
 * Author : Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for vertcoordinate.c
 * including all the variables for the general vertical coordinate
 *
 */
#ifndef _vertcoordinate_h
#define __vertcoordinate_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _vertT {
REAL **uf,**wf,**vf; // the u v w at the layer center of each edge [Nj][Nke]
REAL **omegaf; // the vertical contravariant flux at the face [Ne][Nke]
REAL **ul,**vl; // the u v at the top and bottom face of each layer at each edge [Nc][Nk+1]
REAL **omega; // the vertical contravariant flux [Nc][Nk+1]
REAL **omegac; // the cell-centered vertical contravariant flux [Nc][Nk]
REAL **zc,**zcold; // the cell center vertical location in the Cartesian coordinate [Nc][Nk]
REAL **dvdx, **dudy, **dwdx, **dwdy, **dzdx, **dzdy; // the cell-centered averaged gradient of different variables
REAL *dsigma; // the dsigma for the sigma coordinate to define the vertical coordinate density [Nkmax]
int vertcoord; // the switch for different choice of vertical coordinates 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma

} vertT;

vertT *vert;

void AllocateVertCoordinate(gridT *grid, int myproc);
void UpdateLayerThickness(gridT *grid, propT *prop, physT *phys,int myproc);
void ComputeUf(gridT *grid, propT *prop, physT *phys, int myproc);
void LayerAveragedContinuity(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeZc(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeCellAveragedHorizontalGradient(REAL **gradient, int direction, REAL **scalar, gridT *grid, propT *prop, physT *phys, int myproc);
void VariationalVertCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
REAL InterpToLayerTopFace(int i, int k, REAL **phi, gridT *grid);
#endif