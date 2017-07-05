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
  int i,k,Nk1=5,Nk2=6,Nk3=3, Nk4=grid->Nkmax-Nk1-Nk2-Nk3;
  REAL ratio=1.0/grid->Nkmax,a=15,h3=75.0,delta=60;
  REAL L0=381.3850*3,eta,h1=50;
  REAL r,v_sech;
  for(i=0;i<grid->Nc;i++)
  {
    r=(grid->xv[i])/2/L0;
    //r*=r/4;
    v_sech=2/(exp(r)+exp(-r));
    eta=-2*a*v_sech*v_sech;
    // for two layer system
    //for(k=grid->ctop[i];k<grid->Nk[i]/2;k++)   
    //  grid->dzz[i][k]=(h1-eta)/grid->Nk[i]*2;
    //for(k=grid->Nk[i]/2;k<grid->Nk[i];k++)   
      //grid->dzz[i][k]=((grid->dv[i]+phys->h[i])-h1+eta)/grid->Nk[i]*2;
    for(k=0;k<Nk1;k++)
      grid->dzz[i][k]=(h1-eta-delta/2)/Nk1;
    for(k=Nk1;k<(Nk1+Nk2);k++)
      grid->dzz[i][k]=(delta)/Nk2;
    for(k=(Nk1+Nk2);k<grid->Nk[i];k++)
      grid->dzz[i][k]=(grid->dv[i]+phys->h[i]-h1+eta-delta/2)/(grid->Nk[i]-Nk1-Nk2);
    for(k=0;k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
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
  REAL ratio=1.0/grid->Nkmax,a=120;
  REAL L0=6000/2,eta,h1=460;
  REAL r,v_sech;
  for(i=0;i<grid->Nc;i++)
  {
    r=-grid->xv[i]*grid->xv[i]/4/L0/L0;
    v_sech=2/(exp(r)+exp(-r));
    eta=-a*v_sech*v_sech;
    for(k=grid->ctop[i];k<grid->Nk[i]/2;k++)   
      grid->dzz[i][k]=(h1-eta)/grid->Nk[i]*2;
    for(k=grid->Nk[i]/2;k<grid->Nk[i];k++)   
      grid->dzz[i][k]=((grid->dv[i]+phys->h[i])-h1+eta)/grid->Nk[i]*2;
    for(k=0;k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
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
  //for(k=0;k<10;k++)
    //vert->dsigma[k]=8.5/300;
  //for(k=10;k<30;k++)
    //vert->dsigma[k]=1.5/300;
  //for(k=30;k<grid->Nkmax;k++)
    //vert->dsigma[k]=9.25/300;
  
  for(i=0;i<grid->Nc;i++)
  {
  	for(k=grid->ctop[i];k<grid->Nk[i];k++)
  	{
      vert->dsigma[k]=1.0/grid->Nkmax;
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
 * alpha_H define how much horizontal diffusion 
 * alphaH define how much horizontal density gradient 
 * alphaV define how much vertical density gradient
 */
void MonitorFunctionForVariationalMethod(gridT *grid, propT *prop, physT *phys, int myproc, int numprocs, MPI_Comm comm)
{
}