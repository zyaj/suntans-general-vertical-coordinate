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
  int i,k,Nkmax=grid->Nkmax,Nk_noiso=6;
  REAL alpha_s=0.99,delta=4,H=16,a=0.1,L=10,pi=3.14159265358979323846, rho_diff=0.01, drho,rho;
  REAL zbot,ztop, DH;
  drho=rho_diff/(grid->Nkmax-Nk_noiso),ztop,zbot;

  for(i=0;i<grid->Nc;i++)
  {
    DH=(0-(delta/2-H/2+a*cos(pi/L*grid->xv[i])))/Nk_noiso*2;
    for(k=0;k<Nk_noiso/2;k++)
      grid->dzz[i][k]=DH;
    ztop=(delta/2-H/2+a*cos(pi/L*grid->xv[i]));
    for(k=Nk_noiso/2;k<grid->Nk[i]-Nk_noiso/2-1;k++)
    {
       rho=-rho_diff/2+drho*(k-Nk_noiso/2+1);
       zbot=delta/2/atanh(alpha_s)*atanh(-2/rho_diff*rho)-H/2+a*cos(pi/L*grid->xv[i]);
       grid->dzz[i][k]=ztop-zbot;
       ztop=zbot;
    }
    grid->dzz[i][grid->Nk[i]-Nk_noiso/2-1]=ztop-(-delta/2-H/2+a*cos(pi/L*grid->xv[i]));
    DH=(-delta/2+H/2+a*cos(pi/L*grid->xv[i]))/Nk_noiso*2;
    for(k=grid->Nk[i]-Nk_noiso/2;k<grid->Nk[i];k++)
      grid->dzz[i][k]=DH;
  }
}

/*
 * Function: InitializeVariationalCoordinate
 * Initialize dzz for variational vertical coordinate
 * ----------------------------------------------------
 */
void InitializeVariationalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
	// one for prop->n==1
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
  for(k=0;k<grid->Nkmax-2;k++){

  	vert->dsigma[k]=1.0/grid->Nkmax;
  }
  vert->dsigma[grid->Nkmax-1]=1e-12;
  vert->dsigma[grid->Nkmax-2]=2.0/grid->Nkmax-1e-12;

  for(i=0;i<grid->Nc;i++)
  {
  	for(k=grid->ctop[i];k<grid->Nk[i];k++)
  	{
  	  grid->dzz[i][k]=vert->dsigma[k]*(grid->dv[i]+phys->h[i]);
  	  grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }
}