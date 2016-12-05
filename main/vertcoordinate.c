/*
 * File: vertcoordinate.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file include all the functions for the new general 
 * vertical coordinate
 * 
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
#include "uservertcoordinate.h"
#include "physio.h"
#include "subgrid.h"

/*
 * Function: AllocateVertCoordinate
 * Allocate space for the general vertical coordinate
 * ----------------------------------------------------
 *
 */
void AllocateandInitializeVertCoordinate(gridT *grid, int myproc)
{ 
  int i,j,k;

  // velocity field at layer center for each edge
  vert->uf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->wf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omegaf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  for(j=0;j<grid->Ne;j++)
  {
    vert->uf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->vf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->wf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->omegaf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    for(k=0;k<grid->Nkc[j];k++)
    {
      vert->uf[j][k]=0;
      vert->vf[j][k]=0;
      vert->wf[j][k]=0;
      vert->omegaf[j][k]=0;
    }
  }

  // velocity field at the top and bottom face of each layer
  vert->ul=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vl=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zc=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zcold=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dvdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dudy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");

  for(i=0;i<grid->Nc;i++)
  {
    vert->ul[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->vl[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omegac[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zc[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zcold[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dvdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dudy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");

    for(k=0;k<grid->Nk[i]+1;k++)
    {
      vert->ul[i][k]=0;
      vert->vl[i][k]=0;
      vert->omega[i][k]=0;
    }
    for(k=0;k<grid->Nk[i];k++)
    {
      vert->omegac[i][k]=0;
      vert->zc[i][k]=0;
      vert->zcold[i][k]=0;
      vert->dvdx[i][k]=0;
      vert->dudy[i][k]=0;
      vert->dwdx[i][k]=0;
      vert->dwdy[i][k]=0;
      vert->dzdx[i][k]=0;      
      vert->dzdy[i][k]=0;
    }
  }

  // read vertical coordinate switch
  vert->vertcoord = MPI_GetValue(DATAFILE,"vertcoord","AllocateVertCoordinate",myproc);
  if(vert->vertcoord==3){
    vert->dsigma=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"AllocateVertCoordinate");
    // initialize dsigma
    UserDefinedSigmaCoordinate(grid, myproc);
  }
}

/*
 * Function: UpdateLayerThickness
 * Calculate dz based on the vertical coordinate chosen
 * ----------------------------------------------------
 * The switch is vertcoord 
 * 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma
 * need to modify grid.c to make sure ctop=0 and Nk=Nkmax for all cells 
 */
void UpdateLayerThickness(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k,j,nf,ne;
  REAL fac1,fac2,fac3,Ac;
  if(prop->n==1 || prop->im==0 || prop->wetdry)
  {
    fac1=prop->theta;
    fac2=1-prop->theta;
    fac3=0;
  } else if(prop->im==1){
    fac1=3.0/4.0;
    fac2=0;
    fac3=1.0/4.0;
  } else {
    fac1=5.0/4.0;
    fac2=-1;
    fac3=3.0/4.0;
  }

  if(vert->vertcoord!=1) {
    for(j=0;j<grid->Ne;j++)
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<grid->Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }

  switch(vert->vertcoord)
  {
    case 0:
      UserDefinedVerticalCoordinate(grid,prop,phys,myproc);
      break;
    case 2:
      for(i=0;i<grid->Nc;i++)
      {
        Ac=grid->Ac[i];
        for(k=grid->ctop[i];k<grid->Nk[i];k++){
          for(nf=0;nf<grid->nfaces[i];nf++) {
            ne = grid->face[i*grid->maxfaces+nf];
            if(k<grid->Nke[ne])
              grid->dzz[i][k]-=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
          }
        }
      }
      break;
    case 3:
      for(i=0;i<grid->Nc;i++)
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          grid->dzz[i][k]=vert->dsigma[k]*(phys->h[i]+grid->dv[i]);
      break;
    case 4:
      VariationalVertCoordinate(grid,prop,phys,myproc);
      break;
  }
}


/*
 * Function: VariationalVertCoordinate
 * use the variational moving mesh approach (Tang & Tang, 2003) and (Koltakov & Fringer, 2013)
 * ----------------------------------------------------
 *
 */
void VariationalVertCoordinate(gridT *grid, propT *prop, physT *phys, int myproc)
{

}


/*
 * Function: ComputeUf
 * Calculate the u v w omegaf at the edge center of each layer and 
 * u v at the top and bottom face of each layer
 * ----------------------------------------------------
 * compute uf vf wf and ul vl
 */
void ComputeUf(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,j,k;
  // uf vf and wf
  for(j=0;j<grid->Ne;j++)
    for(k=grid->etop[j];k<grid->Nke[j];k++)
    {
      vert->uf[j][k]=InterpToFace(j,k,phys->uc,phys->u,grid);
      vert->vf[j][k]=InterpToFace(j,k,phys->vc,phys->u,grid);
      vert->wf[j][k]=InterpToFace(j,k,phys->wc,phys->u,grid);
      vert->omegaf[j][k]=InterpToFace(j,k,vert->omegac,phys->u,grid);
    }

  // ul and vl
  for(i=0;i<grid->Nc;i++)
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
       vert->ul[i][k]=InterpToLayerTopFace(i, k, phys->uc, grid);
       vert->vl[i][k]=InterpToLayerTopFace(i, k, phys->vc, grid);
    }
}

/*
 * Function: InterpToLayerTopBotFace
 * Usage: uface = InterpToLayerTopBotFace(j,k,phys->uc,u,grid);
 * -------------------------------------------------
 * Linear interpolation of a Voronoi-centered value to the top and bottom face of each layer
 * 
 * phi_top = (dz1*u2 + dz2*u1)/(dz1+dz2);
 *
 */
REAL InterpToLayerTopFace(int i, int k, REAL **phi, gridT *grid) {
  int nc1, nc2;
  REAL def1, def2, Dj;
  if(k>grid->ctop[i])
  {
    Dj = grid->dzz[i][k]+grid->dzz[i][k-1];
    def1=grid->dzz[i][k-1];
    def2=grid->dzz[i][k];
  }else{
    Dj = grid->dzz[i][k];
    def1=0;
    def2=grid->dzz[i][k];    
  }

  if(def2==0 && k<grid->Nk[i]-1) {
    // means the interior layer is dry
    return phi[i][k+1];
  } else if(def2==0 && k==grid->Nk[i]-1){
    // means the bottom layer is dry
    return 0.0;
  }else {
    return (phi[i][k-1]*def2+phi[i][k]*def1)/Dj;
  }
}

/*
 * Function: LayerAveragedContinuity
 * Compute omega for for the top and bottom face of each layer
 * from the layer-averaged continuity equation
 * ----------------------------------------------------
 * 
 */
void LayerAveragedContinuity(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k,nf,ne;
  REAL fac1,fac2,fac3,flux,Ac;
  if(prop->n==1 || prop->im==0 || prop->wetdry)
  {
    fac1=prop->theta;
    fac2=1-prop->theta;
    fac3=0;
  } else if(prop->im==1){
    fac1=3.0/4.0;
    fac2=0;
    fac3=1.0/4.0;
  } else {
    fac1=5.0/4.0;
    fac2=-1;
    fac3=3.0/4.0;
  }

  for(i=0;i<grid->Nc;i++)
  {
    Ac=grid->Ac[i];
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--){
      flux=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(k<grid->Nke[ne])
          flux+=(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
      }
      vert->omega[i][k]=vert->omega[i][k+1]-flux-(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;
    }
  }
}

/*
 * Function: ComputeZc
 * Compute the vertical location of each cell center
 * ----------------------------------------------------
 */
void ComputeZc(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k;
  REAL z;

  for(i=0;i<grid->Nc;i++)
  {
    for(k=0;k<grid->Nk[i];k++)
      vert->zcold[i][k]=vert->zc[i][k];
    z=-grid->dv[i];
    vert->zc[i][grid->Nk[i]-1]=z+grid->dzz[i][grid->Nk[i]-1]/2;
    for(k=grid->Nk[i]-2;k>=grid->ctop[i];k--)
      vert->zc[i][k]=vert->zc[i][k+1]+grid->dzz[i][k+1]/2+grid->dzz[i][k]/2;
    if(prop->n==1)
      for(k=0;k<grid->Nk[i];k++)
        vert->zcold[i][k]=vert->zc[i][k];
  }  
}

/*
 * Function: ComputeOmega
 * Compute omega from the definition omega=w-udzdx-vdzdy-wg
 * ----------------------------------------------------
 */
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k;
  for(i=0;i<grid->Nc;i++)
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
      vert->omega[i][k]=phys->w[i][k]-vert->ul[i][k]*InterpToLayerTopFace(i,k,vert->dzdx,grid)-\
      vert->vl[i][k]*InterpToLayerTopFace(i,k,vert->dzdy,grid)-(vert->zc[i][k]+grid->dzz[i][k]/2-\
        vert->zcold[i][k]-grid->dzzold[i][k]/2)/prop->dt;
    }
}

/*
 * Function: ComputeCellAveragedGradient
 * Compute the cell averaged gradient for scalars
 * ----------------------------------------------------
 * direction=0 x gradient, direction=1 y gradient
 */
void ComputeCellAveragedHorizontalGradient(REAL **gradient, int direction, REAL **scalar, gridT *grid, propT *prop, 
	physT *phys, int myproc)
{
  int i,k,nf,ne;
  REAL vec;
  for(i=0;i<grid->Nc;i++)
  {
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
    {
      gradient[i][k]=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        if(direction==0)
          vec=grid->n1[ne];
        else
          vec=grid->n2[ne];
        ne = grid->face[i*grid->maxfaces+nf];
        gradient[i][k]+=scalar[ne][k]*grid->df[ne]*vec;
      }
      gradient[i][k]/=grid->Ac[i];
    }
  }
}
