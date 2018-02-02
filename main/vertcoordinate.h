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
#define _vertcoordinate_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _vertT {
REAL **uf,**vf, **qcf; // the u v at the layer center of each edge [Nj][Nke]
REAL **wf; // the w at the layer top of each edge [Nj][Nke], the bottom vertical velocity is assumed as 0.
REAL **zf; // the current vertical location z of each edge layer center [Nj][Nke]
REAL **omegaf; // the vertical contravariant flux at the face [Ne][Nke]
REAL **ul,**vl; // the u v at the top and bottom face of each layer at each edge [Nc][Nk+1]
REAL **omega; // the vertical contravariant flux [Nc][Nk+1]
REAL **Mc; // monitor function value at cell center [Nc][Nk]
REAL **dzztmp; //temporary layer thickness [Nc][Nk]
REAL *Msum; // calculate the sum of M value for each cell column [Nc] at cell center
REAL **Me_l; // monitor function value at each edge top face of each layer [Ne][Nke+1]
//REAL **boundary_omega; // the boundary value of omega, used for momentum advection in Wpredictor [grid->edgedist[5]-grid->edgedist[2]][Nk+1]
REAL **omega_im; // omega_im=fac1*omega+fac2*omega^n+fac3*omega^n-1 [Nc][Nk+1]
REAL **omega_old; // omega^(n) [Nc][Nk+1]
REAL **omega_old2; // omega^(n-1) [Nc][Nk+1]
REAL **omegac; // the cell-centered vertical contravariant flux [Nc][Nk]
REAL **U3,**U3_old,**U3_old2; // w-udzdx-vdzdy [Nc][Nk+1] 
REAL *n1;  // the x component of the outpointing vector from cell center to its face center [Nc*maxfaces*Nc]
REAL *n2; // the y component of the outpointing vector from cell center to its face center [Nc*maxfaces*Nc]
REAL **zc,**zcold; // the cell center vertical location in the Cartesian coordinate [Nc][Nk]
REAL **zl; // the cell center vertical location of cell top face in the cartesian coordinate [Nc][Nk+1]
REAL **f_r; // the cell center relative vorticity dvdx-dudy [Nc][Nk]
REAL **f_re; // the edge center relative vorticity dvdx-dudy [Ne][Nke]
REAL **f_rp; // the nodal relative vorticity dvdx-dudy [Np][Nkp]
REAL *Ap; // the nodal area covered by node [Np]
int *Nkpmin; // the minimum number of layers among the edges that include a specific node [Np]
int *typep; // the type of a node, if all the edges include a node have mark=0 or 5, type=1, else -1
REAL *CCNpe; // the inner product between nodal counter-clockwise diction to the positive velocity of a specific edge [2*Ne]
REAL **dzdx, **dzdy, **dqdx, **dqdy; // the cell-centered averaged gradient of different variables
REAL **dvdx, **dudy, **dvdy, **dudx, **dwdx, **dwdy;
REAL *dsigma; // the dsigma for the sigma coordinate to define the vertical coordinate density [Nkmax]
REAL *tmp; //temporary array for output
REAL *tmp_nc; // temporary array for zl
int *Nkeb; // store the layer index for each edge when ze is higher than buffer height [Ne]
REAL *zfb; // store the sum of dzf from bottom layer to Nkeb [Ne] 
FILE *zcFID, *dzzFID, *omegaFID;
int modifydzf;  // whether recalculate the u^im to check whether the flux height is upwind
int dJdtmeth; // how to treat the u/JdJdt term 0 implicit and 1 explicit
int dzfmeth; // how to calculate flux height 1 upwind, 2 lax wendroff, 3 superbee, 4 vanleer
REAL thetaT,vertdzmin;
} vertT;

vertT *vert;

void AllocateVertCoordinate(gridT *grid, propT *prop, int myproc);
void UpdateLayerThickness(gridT *grid, propT *prop, physT *phys,int index, int myproc,int numprocs, MPI_Comm comm);
void InitializeLayerThickness(gridT *grid, propT *prop, physT *phys,int myproc);
void ComputeUf(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeUl(gridT *grid, propT *prop, physT *phys, int myproc);
void LayerAveragedContinuity(REAL **omega, gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int index, int myproc);
void ComputeZc(gridT *grid, propT *prop, physT *phys, int index, int myproc);
void ComputeDSigma(gridT *grid, physT *phys, int myproc);
void VertCoordinateHorizontalSource(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void ComputeCellAveragedHorizontalGradient(REAL **gradient, int direction, REAL **scalar, gridT *grid, propT *prop, physT *phys, int myproc);
void VariationalVertCoordinate(gridT *grid, propT *prop, physT *phys, int myproc,int numprocs, MPI_Comm comm);
void VariationalVertCoordinateAverageMethod(gridT *grid, propT *prop, physT *phys, int myproc);
REAL InterpToLayerTopFace(int i, int k, REAL **phi, gridT *grid);
void VertCoordinateBasic(gridT *grid, propT *prop, physT *phys, int myproc);
void VertCoordinateBasicRestart(gridT *grid, propT *prop, physT *phys, int myproc);
void ComputeNormalVector(gridT *grid, physT *phys, int myproc);
void OpenVertCoordinateFiles(gridT *grid,int mergeArrays, int myproc);
void OutputVertCoordinate(gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void StoreVertVariables(gridT *grid, physT *phys);
void FindBottomLayer(gridT *grid, propT *prop, physT *phys, int myproc);
void UpdateCellcenteredFreeSurface(gridT *grid, propT *prop, physT *phys, int myproc);
void VerifyFluxHeight(gridT *grid, propT *prop, physT *phys, int myproc);
void TvdFluxHeight(gridT *grid, physT *phys, propT *prop, int TVD, MPI_Comm comm, int myproc);
void ComputeNodalData(gridT *grid,int myproc);
void ComputeRelativeVorticity(gridT *grid, physT *phys, propT *prop, int myproc);
#endif