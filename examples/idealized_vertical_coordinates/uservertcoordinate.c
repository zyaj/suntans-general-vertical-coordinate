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
  int i,k,j,nf,ne;
  REAL fac1,fac2,fac3,Ac,sum=0,sum1,flux;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(i=0;i<grid->Nc;i++)
  { 
    sum1=0;
    sum=0;
    Ac=grid->Ac[i];
    for(k=25;k<grid->Nk[i];k++){
      grid->dzz[i][k]=(grid->dv[i]-15)/10;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(k<grid->Nke[ne]){
          sum+=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+
            fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
        }
      }
    } 

    Ac=grid->Ac[i];
    k=24;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      if(k<grid->Nke[ne]){
        grid->dzz[i][k]-=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+
          fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
      }
    }
      grid->dzz[i][k]-=sum;

    for(k=6;k<24;k++){
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(k<grid->Nke[ne]){
          grid->dzz[i][k]-=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+
          fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
        }
      }
    }

    for(k=grid->Nk[i]-1;k>=6;k--)
      sum1+=grid->dzz[i][k];

    grid->dzz[i][5]=grid->dv[i]-5-sum1;

    grid->dzz[i][0]=1+phys->h[i];
    for(k=1;k<5;k++)
      grid->dzz[i][k]=1;

    sum1=0;
    for(k=0;k<grid->Nk[i];k++)
      sum1+=grid->dzz[i][k];
    /*if(i==0)
      for(k=0;k<grid->Nk[i];k++)
        printf("dv %e sum1 %e sum %e n %d k %d x %e dzz %e vert %e\n",grid->dv[i]+phys->h[i],sum1,sum,prop->n,k,grid->xv[i],grid->dzz[i][k],vert->omega[i][k]); 
    */
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
  REAL DH,sum;
  REAL alpha_s=0.99,delta=8,a=0.5,L=50,pi=3.14159265358979323846, rho_diff=0.01, drho,rho,A;

  for(i=0;i<grid->Nc;i++)
  { 
    sum=0;
    // first 10 layers zlevel
    grid->dzz[i][0]=phys->h[i]+1;
    for(k=1;k<5;k++)
      grid->dzz[i][k]=1;
    // after 20 layers iso
    A=a*cos(pi/L*grid->xv[i]);
    grid->dzz[i][5]=5-A-delta/2;
    for(k=6;k<24;k++)
      grid->dzz[i][k]=delta/18;
    grid->dzz[i][24]=5+A-delta/2;
    // last 10 layers sigma   
     for(k=25;k<grid->Nk[i];k++){
       grid->dzz[i][k]=(grid->dv[i]-15)/10;    
     }
    for(k=0;k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
    /*if(i==0)
      for(k=0;k<grid->Nk[i];k++)
      {
        sum+=grid->dzz[i][k];
        printf("k %d dz %e sum %e H %e\n",k,grid->dzz[i][k],sum, phys->h[i]+grid->dv[i]);
      }*/  
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
  int i,k,Nkmax=grid->Nkmax,Nk_noiso=20;
  REAL alpha_s=0.99,delta=1,H=16,a=0.1,L=10,pi=3.14159265358979323846, rho_diff=0.01, drho,rho;
  REAL zbot,ztop, DH;
  drho=rho_diff/(grid->Nkmax-Nk_noiso),ztop,zbot;

  for(i=0;i<grid->Nc;i++)
  {
    DH=(0-(delta-H/2+a*cos(pi/L*grid->xv[i])))/Nk_noiso*2;
    for(k=0;k<Nk_noiso/2;k++)
      grid->dzz[i][k]=DH;
    /*ztop=(delta/2-H/2+a*cos(pi/L*grid->xv[i]));
    for(k=Nk_noiso/2;k<grid->Nk[i]-Nk_noiso/2-1;k++)
    {
       rho=-rho_diff/2+drho*(k-Nk_noiso/2+1);
       zbot=delta/2/atanh(alpha_s)*atanh(-2/rho_diff*rho)-H/2+a*cos(pi/L*grid->xv[i]);
       grid->dzz[i][k]=ztop-zbot;
       ztop=zbot;
    }
    grid->dzz[i][grid->Nk[i]-Nk_noiso/2-1]=ztop-(-delta/2-H/2+a*cos(pi/L*grid->xv[i]));
    */
    for(k=Nk_noiso/2;k<grid->Nk[i]-Nk_noiso/2;k++)
    {
      grid->dzz[i][k]=2*delta/(grid->Nk[i]-Nk_noiso);
    }

    DH=(-delta+H/2+a*cos(pi/L*grid->xv[i]))/Nk_noiso*2;
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