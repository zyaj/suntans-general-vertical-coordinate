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
	// one for other update scheme
	
}

/*
 * Function: InitializeVerticalCoordinate
 * to setup the initial condition of dzz for user defined vertical coordinate
 * ----------------------------------------------------
 */
void InitializeVerticalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
	// one for other update scheme
	
}

/*
 * Function: InitializeIsopycnalCoordinate
 * User define isopycnal coordinate 
 * define the initial dzz for each cell under isopycnal coordinate
 * ----------------------------------------------------
 */
void InitializeIsopycnalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k,Nkmax=grid->Nkmax;
  REAL ratio=1.0/Nkmax;
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      grid->dzz[i][k]=ratio*(phys->h[i]+grid->dv[i]);
}

/*
 * Function: InitializeVariationalCoordinate
 * Initialize dzz for variational vertical coordinate
 * --------zz--------------------------------------------
 */
void InitializeVariationalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  REAL ratio=1.0/grid->Nkmax;

  for(i=0;i<grid->Nc;i++)
  {
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
    {
      grid->dzz[i][k]=ratio*(grid->dv[i]+phys->h[i]);
      grid->dzzold[i][k]=grid->dzz[i][k];
    }
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
  int i,k;
  for(k=0;k<grid->Nkmax;k++){

  	vert->dsigma[k]=1.0/grid->Nkmax;
  }

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
void MonitorFunctionForVariationalMethod(gridT *grid, propT *prop, physT *phys, int myproc)
{
   int i,k;
   REAL alphaM=160,minM=0.02,max;

   for(i=0;i<grid->Nc;i++)
   {
     max=0;
     vert->Msum[i]=0;
     for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++){
       vert->M[i][k]=1000*(phys->rho[i][k-1]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k-1]+grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
       if(fabs(vert->M[i][k])>max)
         max=fabs(vert->M[i][k]);
     }
     
     // top boundary
     k=grid->ctop[i];
     vert->M[i][k]=1000*(phys->rho[i][k]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
     if(fabs(vert->M[i][k])>max)
       max=fabs(vert->M[i][k]);   
     // bottom boundary
     k=grid->Nk[i]-1;
     vert->M[i][k]=1000*(phys->rho[i][k-1]-phys->rho[i][k])/(0.5*grid->dzz[i][k-1]+0.5*grid->dzz[i][k]);
     if(fabs(vert->M[i][k])>max)
       max=fabs(vert->M[i][k]);   
     if(max<1)
       max=1;
     
     for(k=grid->ctop[i];k<grid->Nk[i];k++){ 
       vert->M[i][k]=1/sqrt(1+alphaM*vert->M[i][k]/max*vert->M[i][k]/max);
       if(vert->M[i][k]<minM)
         vert->M[i][k]=minM;     
       vert->Msum[i]+=vert->M[i][k];
     }

   }
}