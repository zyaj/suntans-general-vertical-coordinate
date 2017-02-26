/*
 * File: uservertcoordinate.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file include a function to user defined vertical coordinate
 * 
 */

#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "initialization.h"
#include "boundaries.h"
#include "util.h"
#include "tvd.h"
#include "mympi.h"
#include "scalars.h"
#include "vertcoordinate.h"
#include "physio.h"
#include "subgrid.h"

/*
 * Function: UserDefinedVerticalCoordinate
 * User define vertical coordinate 
 * basically it is a user-defined function to calculate the layer thickness based on 
 * different criterion
 * ----------------------------------------------------
 * the original code has already include 1 z-level, 2 isopycnal, 3 sigma, 4 variational 
 */
void UserDefinedVerticalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  for(i=0;i<grid->Nc;i++)
  { 
    k=grid->ctop[i];
    grid->dzz[i][k]=7.0+phys->h[i];
    for(k=grid->ctop[i]+1;k<10;k++)
      grid->dzz[i][k]=7.0;
    for(k=10;k<50;k++)
      grid->dzz[i][k]=1.5;
    for(k=50;k<grid->Nk[i];k++)
      grid->dzz[i][k]=(grid->dv[i]-130.0)/(grid->Nk[i]-50);
  }
}

/*
 * Function: InitializeVerticalCoordinate
 * to setup the initial condition of dzz for user defined vertical coordinate
 * ----------------------------------------------------
 */
void InitializeVerticalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  for(i=0;i<grid->Nc;i++)
  { 
    k=grid->ctop[i];
    grid->dzz[i][k]=7.0;
    for(k=grid->ctop[i]+1;k<10;k++)
      grid->dzz[i][k]=7.0;
    for(k=10;k<50;k++)
      grid->dzz[i][k]=1.5;
    for(k=50;k<grid->Nk[i];k++)
      grid->dzz[i][k]=(grid->dv[i]-130.0)/(grid->Nk[i]-50);
  }	
}

/*
 * Function: InitializeIsopycnalCoordinate
 * User define isopycnal coordinate 
 * define the initial dzz for each cell under isopycnal coordinate
 * ----------------------------------------------------
 */
void InitializeIsopycnalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  REAL ratio=1.0/grid->Nkmax,a=30;
  REAL L_rho=2000,eta,h1=100;

  for(i=0;i<grid->Nc;i++)
  {
    eta=-a*exp(-(grid->xv[i]-60000)*(grid->xv[i]-60000)/L_rho/L_rho);
    for(k=grid->ctop[i];k<grid->Nk[i]/2;k++)   
      grid->dzz[i][k]=(h1-eta)/grid->Nk[i]*2;
    for(k=grid->Nk[i]/2;k<grid->Nk[i];k++)   
      grid->dzz[i][k]=((grid->dv[i]+phys->h[i])-h1+eta)/grid->Nk[i]*2;
  }
}

/*
 * Function: InitializeVariationalCoordinate
 * Initialize dzz for variational vertical coordinate
 * --------zz--------------------------------------------
 */
void InitializeVariationalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  REAL ratio=1.0/grid->Nkmax,a=30;
  REAL L_rho=1000,eta,h1=100;

  for(i=0;i<grid->Nc;i++)
  {
    eta=-a*exp(-(grid->xv[i]-60000)*(grid->xv[i]-60000)/L_rho/L_rho);
    grid->dzz[i][0]=h1-eta;
    grid->dzz[i][1]=(grid->dv[i]+phys->h[i])-grid->dzz[i][0];
  }
}

/*
 * Function: UserDefinedSigmaCoordinate
 * User define sigma coordinate 
 * basically to define the dsigma for each layer
 * ----------------------------------------------------
 */
void InitializeSigmaCoordinate(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k,sum=0;
  for(k=0;k<10;k++)
    vert->dsigma[k]=8.5/300;
  for(k=10;k<30;k++)
    vert->dsigma[k]=1.5/300;
  for(k=30;k<grid->Nkmax;k++)
    vert->dsigma[k]=9.25/300;

  for(i=0;i<grid->Nc;i++)
  {
  	for(k=grid->ctop[i];k<grid->Nk[i];k++)
  	{
  	  grid->dzz[i][k]=vert->dsigma[k]*(grid->dv[i]+phys->h[i]);
  	  grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }
}

/*
 * Function: MonitorFunctionForVariationalMethod
 * calculate the value of monitor function for the variational approach
 * to update layer thickness when nonlinear==4
 * ----------------------------------------------------
 * Mii=sqrt(1-alphaM*(drhodz)^2)
 */
void MonitorFunctionForAverageMethod(gridT *grid, propT *prop, physT *phys, int myproc)
{
}

/*
 * Function: MonitorFunctionForVariationalMethod
 * calculate the value of monitor function for the variational approach
 * to update layer thickness when nonlinear==4
 * solve the elliptic equation using iteration method
 * ----------------------------------------------------
 * Mii=sqrt(1-alphaM*(drhodz)^2)
 */
void MonitorFunctionForVariationalMethod(gridT *grid, propT *prop, physT *phys, int myproc, int iter, int numprocs, MPI_Comm comm)
{
  int i,k,j,nf,neigh,ne,kk,nc1,nc2;
  REAL alphaH=1, alphaV=160, minM=0.15,max,tmp,max_gradient_v=0,max_gradient_h=0,max_gradient_h_global,H1,H2;
}