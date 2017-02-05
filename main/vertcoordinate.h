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
REAL **zf; // the current vertical location z of each edge layer center [Nj][Nke]
REAL **omegaf; // the vertical contravariant flux at the face [Ne][Nke]
REAL **ul,**vl; // the u v at the top and bottom face of each layer at each edge [Nc][Nk+1]
REAL **omega; // the vertical contravariant flux [Nc][Nk+1]
REAL **M; // monitor function value [Nc][Nk]
REAL **dzztmp; //temporary layer thickness [Nc][Nk]
REAL *Msum; // calculate the sum of M value for each cell column [Nc]
//REAL **boundary_omega; // the boundary value of omega, used for momentum advection in Wpredictor [grid->edgedist[5]-grid->edgedist[2]][Nk+1]
REAL **omega_im; // omega_im=fac1*omega+fac2*omega^n+fac3*omega^n-1 [Nc][Nk+1]
REAL **omega_old; // omega^(n) [Nc][Nk+1]
REAL **omega_old2; // omega^(n-1) [Nc][Nk+1]
REAL **omega_star; // the vertical contravariant flux [Nc][Nk+1] for hydrostatic calculation and prepare for scalar transport 
REAL **omegac; // the cell-centered vertical contravariant flux [Nc][Nk]
REAL **U3,**U3_old,**U3_old2; // w-udzdx-vdzdy [Nc][Nk+1] 
REAL *n1;  // the x component of the outpointing vector from cell center to its face center [Nc*maxfaces*Nc]
REAL *n2; // the y component of the outpointing vector from cell center to its face center [Nc*maxfaces*Nc]
REAL **zc,**zcold; // the cell center vertical location in the Cartesian coordinate [Nc][Nk]
REAL **f_r; // the cell center relative vorticity dvdx-dudy [Nc][Nk]
REAL **dvdx, **dudy, **dvdy, **dudx, **dwdx, **dwdy, **dzdx, **dzdy; // the cell-centered averaged gradient of different variables
REAL *dsigma; // the dsigma for the sigma coordinate to define the vertical coordinate density [Nkmax]
REAL *tmp; //temporary array for output
int *Nkeb; // store the layer index for each edge when ze is higher than buffer height [Ne]
REAL *zfb; // store the sum of dzf from bottom layer to Nkeb [Ne] 
//int vertcoord; // the switch for different choice of vertical coordinates 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma
FILE *zcFID, *dzzFID, *omegaFID;
int modifydzf;
} vertT;

vertT *vert;

void AllocateVertCoordinate(gridT *grid, propT *prop, int myproc);
void UpdateLayerThickness(gridT *grid, propT *prop, physT *phys,int index, int myproc);
void InitializeLayerThickness(gridT *grid, propT *prop, physT *phys,int myproc);
void ComputeUf(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeUl(gridT *grid, propT *prop, physT *phys, int myproc);
void LayerAveragedContinuity(REAL **omega, gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int index, int myproc);
void ComputeZc(gridT *grid, propT *prop, physT *phys, int index, int myproc);
void VertCoordinateHorizontalSource(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void ComputeCellAveragedHorizontalGradient(REAL **gradient, int direction, REAL **scalar, gridT *grid, propT *prop, physT *phys, int myproc);
void VariationalVertCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
REAL InterpToLayerTopFace(int i, int k, REAL **phi, gridT *grid);
void VertCoordinateBasic(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeNormalVector(gridT *grid, physT *phys, int myproc);
void OpenVertCoordinateFiles(gridT *grid,int mergeArrays, int myproc);
void OutputVertCoordinate(gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void StoreVertVariables(gridT *grid, physT *phys);
void FindBottomLayer(gridT *grid, propT *prop, physT *phys, int myproc);
void UpdateCellcenteredFreeSurface(gridT *grid, propT *prop, physT *phys, int myproc);
void VerifyFluxHeight(gridT *grid, propT *prop, physT *phys, int myproc);
#endif