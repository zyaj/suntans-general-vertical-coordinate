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
  int i,k,k0,ktop,Nkmax=grid->Nkmax,Nk_noiso=2,Nk_mixed,Nk_pycnocline,Nk_lower;  
  REAL a0=15, a;
  REAL L0=381.3850*3,h1=75;
  REAL r,v_sech;
  REAL alpha_s=0.99,delta=75;
  REAL rho_diff=0.01,beta=1e-3,eta,*zptr,z;
  
  // Mixed layer is 75-37.5 = 37.5 m
  // Pycnocline is 75 m
  Nk_mixed = 3;
  Nk_pycnocline = 7;
  Nk_lower = Nkmax-Nk_mixed-Nk_pycnocline;

  a=a0;
  for(i=-1;i<grid->Nc;i++) {
    if(i==-1) {
      zptr=grid->dz;
      eta=0;
      ktop=0;
    } else {
      zptr=grid->dzz[i];
      eta=-2*a*pow(1/cosh(grid->xv[i]/2/L0),2);
      ktop=grid->ctop[i];
    }

    for(k=ktop;k<Nk_mixed;k++) 
      zptr[k]=(h1-eta-delta/2)/Nk_mixed;
    for(k=Nk_mixed;k<Nk_mixed+Nk_pycnocline;k++) 
      zptr[k]=delta/Nk_pycnocline;
    for(k=Nk_mixed+Nk_pycnocline;k<Nkmax;k++) 
      zptr[k]=(0*grid->dv[i]+600-h1+eta-delta/2)/Nk_lower;
  }

  for(i=0;i<grid->Nc;i++) {
    z=0;
    for(k=0;k<Nkmax;k++) {
      z-=grid->dzz[i][k];
      if(z<=-grid->dv[i]) {
	grid->dzz[i][k]=z+grid->dzz[i][k]+grid->dv[i];
	if(grid->dzz[i][k]<vert->vertdzmin)
	  grid->dzz[i][k]=vert->vertdzmin;
	break;
      }
    }
    for(k0=k+1;k0<Nkmax;k0++)
      grid->dzz[i][k0]=vert->vertdzmin;
  }

  for(i=0;i<grid->Nc;i++) {
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
  REAL ratio=1.0/grid->Nkmax,a=250;
  REAL L_rho=15000,eta,h1=250,L1=300000,delta=0.000001,Nk1=2;

  for(i=0;i<grid->Nc;i++)
  {
    eta=a*exp(-(grid->xv[i]-L1)*(grid->xv[i]-L1)/L_rho/L_rho);
    for(k=0;k<Nk1;k++)
      grid->dzz[i][k]=(h1+eta-delta/2)/Nk1;
    for(k=Nk1;k<grid->Nk[i]-Nk1;k++)
      grid->dzz[i][k]=delta/(grid->Nk[i]-2*Nk1);
    for(k=grid->Nk[i]-Nk1;k<grid->Nk[i];k++)
      grid->dzz[i][k]=(grid->dv[i]-h1-eta-delta/2)/Nk1;
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
   int i,k;
   REAL alphaM=0,minM=0.15,max;
   // nonlinear=1 or 5 stable with alpham=320
   // nonlinear=2 stable with alpham=60
   // nonlinear=4 stable with alpham=60
 
   for(i=0;i<grid->Nc;i++)
   {
     max=0;
     vert->Msum[i]=0;
     for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++){
       vert->Mc[i][k]=1000*(phys->rho[i][k-1]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k-1]+grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
       if(fabs(vert->Mc[i][k])>max)
         max=fabs(vert->Mc[i][k]);
     }
     
     // top boundary
     k=grid->ctop[i];
     vert->Mc[i][k]=1000*(phys->rho[i][k]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
     if(fabs(vert->Mc[i][k])>max)
       max=fabs(vert->Mc[i][k]);   
     // bottom boundary
     k=grid->Nk[i]-1;
     vert->Mc[i][k]=1000*(phys->rho[i][k-1]-phys->rho[i][k])/(0.5*grid->dzz[i][k-1]+0.5*grid->dzz[i][k]);
     if(fabs(vert->Mc[i][k])>max)
       max=fabs(vert->Mc[i][k]);   
     if(max<1)
       max=1;
     
     for(k=grid->ctop[i];k<grid->Nk[i];k++){ 
       vert->Mc[i][k]=sqrt(1+alphaM*vert->Mc[i][k]/max*vert->Mc[i][k]/max);
       if(vert->Mc[i][k]<minM)
         vert->Mc[i][k]=minM;     
       vert->Msum[i]+=1/vert->Mc[i][k];
     }
   }
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
  int i,k,j,nf,neigh,ne,kk,nc1,nc2;
  REAL normal,alpha_H=0, alphaH=0, alphaV=10, minM=0.15, maxM=100000,max,tmp;
  REAL max_gradient_v,max_gradient_h=0,max_gradient_h_global,H1,H2,rho1,rho2;
  // initialize everything zero
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      vert->Mc[i][k]=0;
  
  for(j=0;j<grid->Ne;j++)
    for(k=0;k<grid->Nke[j]+1;k++)
      vert->Me_l[j][k]=0;

  // calculate Mc first
  for(i=0;i<grid->Nc;i++)
  {
    // calculate gradient
    max_gradient_v=0;
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++){
      vert->Mc[i][k]=RHO0*(phys->rho[i][k-1]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k-1]+grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
      if(fabs(vert->Mc[i][k])>max_gradient_v)
        max_gradient_v=fabs(vert->Mc[i][k]);
    }

    // top boundary
    k=grid->ctop[i];
    vert->Mc[i][k]=RHO0*(phys->rho[i][k]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
    if(fabs(vert->Mc[i][k])>max_gradient_v)
      max_gradient_v=fabs(vert->Mc[i][k]);   

    // bottom boundary
    k=grid->Nk[i]-1;
    vert->Mc[i][k]=RHO0*(phys->rho[i][k-1]-phys->rho[i][k])/(0.5*grid->dzz[i][k-1]+0.5*grid->dzz[i][k]);
    if(fabs(vert->Mc[i][k])>max_gradient_v)
      max_gradient_v=fabs(vert->Mc[i][k]);   
    if(max_gradient_v<1)
      max_gradient_v=1;  
    // calculate monitor function value
    for(k=grid->ctop[i];k<grid->Nk[i];k++){ 
      if(alphaV!=0){
        if(vert->Mc[i][k]/max_gradient_v>(maxM-1)/sqrt(alphaV))
          vert->Mc[i][k]=maxM; 
        else
          vert->Mc[i][k]=sqrt(1+alphaV*vert->Mc[i][k]/max_gradient_v*vert->Mc[i][k]/max_gradient_v);
      } else 
        vert->Mc[i][k]=1;
    }
  }

  // calculate Me_l 
  // calculate gradient
  for(j=0;j<grid->Ne;j++)
  {
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;   
    
    // interior layer
    for(k=grid->etop[j]+1;k<grid->Nke[j];k++)
    {
      rho1=grid->dzz[nc1][k-1]/(grid->dzz[nc1][k]+grid->dzz[nc1][k-1])*phys->rho[nc1][k]+
        grid->dzz[nc1][k]/(grid->dzz[nc1][k]+grid->dzz[nc1][k-1])*phys->rho[nc1][k-1];
      rho2=grid->dzz[nc2][k-1]/(grid->dzz[nc2][k]+grid->dzz[nc2][k-1])*phys->rho[nc2][k]+
        grid->dzz[nc2][k]/(grid->dzz[nc2][k]+grid->dzz[nc2][k-1])*phys->rho[nc2][k-1];
      vert->Me_l[j][k]=RHO0*(rho1-rho2)/grid->dg[j];
      if(fabs(vert->Me_l[j][k])>max_gradient_h)
        max_gradient_h=fabs(vert->Me_l[j][k]);
    }

    // top and bottom
    k=grid->etop[j];
    vert->Me_l[j][k]=RHO0*(phys->rho[nc1][k]-phys->rho[nc2][k])/grid->dg[j];
    if(fabs(vert->Me_l[j][k])>max_gradient_h)
      max_gradient_h=fabs(vert->Me_l[j][k]); 
    k=grid->Nke[j];
    vert->Me_l[j][k]=RHO0*(phys->rho[nc1][k-1]-phys->rho[nc2][k-1])/grid->dg[j];
    if(fabs(vert->Me_l[j][k])>max_gradient_h)
      max_gradient_h=fabs(vert->Me_l[j][k]);
  }

  // find max_global and normalize
  MPI_Reduce(&max_gradient_h,&max_gradient_h_global,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Bcast(&max_gradient_h_global,1,MPI_DOUBLE,0,comm);
  if(max_gradient_h_global<1){
    max_gradient_h_global=1.0;
  }

  for(j=0;j<grid->Ne;j++)
    for(k=grid->etop[j];k<=grid->Nke[j];k++){
      vert->Me_l[j][k]=alpha_H*sqrt(1+alphaH*vert->Me_l[j][k]/max_gradient_h_global*
        vert->Me_l[j][k]/max_gradient_h_global);
    }
}
