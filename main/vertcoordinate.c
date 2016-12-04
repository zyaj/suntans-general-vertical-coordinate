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
void AllocateVertCoordinate(gridT *grid, int myproc)
{ 
  int i,j;

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
  }

  // velocity field at the top and bottom face of each layer
  vert->ul=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vl=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zc=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dvdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dudy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");

  for(i=0;i<grid->Nc;i++)
  {
    vert->ul[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->vl[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->zc[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dvdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dudy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
  }

  // read vertical coordinate switch
  vert->vertcoord = MPI_GetValue(DATAFILE,"vertcoord","AllocateVertCoordinate",myproc);
  if(vert->vertcoord==3)
    vert->dsigma=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"AllocateVertCoordinate");
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

}

/*
 * Function: ComputeOmega
 * Compute omega from the definition omega=w-udzdx-vdzdy-wg
 * ----------------------------------------------------
 */
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int myproc)
{

}

/*
 * Function: ComputeZc
 * Compute the vertical location of each cell center
 * ----------------------------------------------------
 */
void ComputeZc(gridT *grid, propT *prop, physT *phys, int myproc)
{

}

/*
 * Function: ComputeCellAveragedGradient
 * Compute the cell averaged gradient for scalars
 * ----------------------------------------------------
 */
void ComputeCellAveragedGradient(REAL **gradient, REAL **scalar, gridT *grid, propT *prop, 
	physT *phys, int myproc)
{

}

