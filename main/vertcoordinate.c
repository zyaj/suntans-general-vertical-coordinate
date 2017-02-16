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
void AllocateandInitializeVertCoordinate(gridT *grid, propT *prop, int myproc)
{ 
  int i,j,k;

  // velocity field at layer center for each edge
  vert->uf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->wf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omegaf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->Nkeb=(int *)SunMalloc(grid->Ne*sizeof(int),"AllocateVertCoordinate");
  vert->zfb=(REAL *)SunMalloc(grid->Ne*sizeof(REAL),"AllocateVertCoordinate");
  vert->modifydzf = MPI_GetValue(DATAFILE,"modifydzf","AllocateVertCoordinate",myproc);
  vert->dJdtmeth = MPI_GetValue(DATAFILE,"dJdtmeth","AllocateVertCoordinate",myproc);

  for(j=0;j<grid->Ne;j++)
  {
    vert->Nkeb[j]=0;
    vert->zfb[j]=0.0;
    vert->uf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->vf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->wf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->omegaf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    for(k=0;k<grid->Nkc[j];k++)
    {
      vert->uf[j][k]=0;
      vert->zf[j][k]=0;
      vert->vf[j][k]=0;
      vert->wf[j][k]=0;
      vert->omegaf[j][k]=0;
    }
  }

  // normal vector
  vert->n1=(REAL *)SunMalloc(grid->Nc*grid->maxfaces*sizeof(REAL),"AllocateVertCoordinate");
  vert->n2=(REAL *)SunMalloc(grid->Nc*grid->maxfaces*sizeof(REAL),"AllocateVertCoordinate");
  // velocity field at the top and bottom face of each layer
  vert->ul=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vl=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega_im=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega_old=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega_old2=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  // U3=w-udzdx-udzdy
  vert->U3=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->U3_old=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->U3_old2=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omegac=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zc=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zcold=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->f_r=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dvdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dudy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dudx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dvdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->M=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzztmp=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->Msum=(REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateVertCoordinate");

  for(i=0;i<grid->Nc;i++)
  {
    vert->ul[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->vl[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega_im[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega_old[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega_old2[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->U3[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->U3_old[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->U3_old2[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omegac[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zc[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zcold[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dvdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dvdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->f_r[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dudy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dudx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->M[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzztmp[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    
    // initialize normal vector
    for(k=0;k<grid->maxfaces;k++)
    {
      vert->n1[i*grid->maxfaces+k]=0;
      vert->n2[i*grid->maxfaces+k]=0;
    }

    for(k=0;k<grid->Nk[i]+1;k++)
    {
      vert->ul[i][k]=0;
      vert->vl[i][k]=0;
      vert->omega[i][k]=0;
      vert->omega_im[i][k]=0;
      vert->omega_old[i][k]=0;
      vert->omega_old2[i][k]=0;
      vert->U3[i][k]=0;
      vert->U3_old[i][k]=0;
      vert->U3_old2[i][k]=0;
    }
    for(k=0;k<grid->Nk[i];k++)
    {
      vert->omegac[i][k]=0;
      vert->zc[i][k]=0;
      vert->zcold[i][k]=0;
      vert->f_r[i][k]=0;
      vert->dvdx[i][k]=0;
      vert->dudy[i][k]=0;
      vert->dudx[i][k]=0;
      vert->dvdy[i][k]=0;
      vert->dwdx[i][k]=0;
      vert->dwdy[i][k]=0;
      vert->dzdx[i][k]=0;      
      vert->dzdy[i][k]=0;
      vert->M[i][k]=0;
      vert->dzztmp[i][k]=0;
    }
    vert->Msum[i]=0;
  }

  // temporary array for output 
  vert->tmp = (REAL *)SunMalloc(10*grid->Nc*sizeof(REAL),"AllocatePhysicalVariables");

  // read vertical coordinate switch
  if(prop->vertcoord==3)
    vert->dsigma=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"AllocateVertCoordinate");
}

/*
 * Function: InitializeLayerThickness
 * Calculate dz for the first time step based on the vertical coordinate chosen
 * ----------------------------------------------------
 * The switch is vertcoord 
 * 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma
 * need to modify grid.c to make sure ctop=0 and Nk=Nkmax for all cells 
 */
void InitializeLayerThickness(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,j,k;
  if(prop->vertcoord!=5){
    for(i=0;i<grid->Nc;i++)
    {
      grid->ctopold[i]=grid->ctop[i]=0;
      if(grid->Nk[i]!=grid->Nkmax){
        printf("There is something wrong in setting the Nkmax for each cell\n");
        exit(1);
      }
    }

    for(j=0;j<grid->Ne;j++)
      grid->etopold[j]=grid->etop[j]=0;
  }

  switch(prop->vertcoord)
  {
    case 0: // user defined
      InitializeVerticalCoordinate(grid,prop,phys,myproc);  
      break;    
    case 2:
      InitializeIsopycnalCoordinate(grid,prop,phys,myproc);
      break;
    case 3:
      InitializeSigmaCoordinate(grid, prop, phys, myproc);
      break;
    case 4:
      InitializeVariationalCoordinate(grid, prop,phys, myproc);
      break;
  } 

  // set the old value as the initial condition
  for(i=0;i<grid->Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
}


/*
 * Function: UpdateLayerThickness
 * Calculate dz based on the vertical coordinate chosen
 * ----------------------------------------------------
 * The switch is vertcoord 
 * 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma
 * need to modify grid.c to make sure ctop=0 and Nk=Nkmax for all cells 
 */
void UpdateLayerThickness(gridT *grid, propT *prop, physT *phys, int index, int myproc)
{
  int i,k,j,nf,ne;
  REAL fac1,fac2,fac3,Ac,sum,sum_old;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  // ctop and etop is always 0 and Nke and Nk is always Nkmax
  if(!index)
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++){
        grid->dzzold[i][k]=grid->dzz[i][k];
        vert->zcold[i][k]=vert->zc[i][k];
      }
    }

  switch(prop->vertcoord)
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
            if(k<grid->Nke[ne]){
              grid->dzz[i][k]-=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+
                fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
            }
              // the function is not right for subgrid module
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
    case 5:
      UpdateDZ(grid,phys,prop, 0); 
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
   int i,k,nf,Neigh;
   REAL thetaT=0.15,thetaS=0.25,sum;
   MonitorFunctionForVariationalMethod(grid, prop, phys, myproc);
   // use the monitor function first
   for(i=0;i<grid->Nc;i++)
   {
     sum=0;
     Neigh=0;
     for(nf=0;nf<grid->nfaces[i];nf++)
       if(grid->neigh[i*grid->maxfaces+nf]!=-1)
        Neigh++;
     for(k=0;k<grid->Nk[i];k++){
        vert->dzztmp[i][k]=0;
        // use the monitor function
        grid->dzz[i][k]=(phys->h[i]+grid->dv[i])/vert->Msum[i]*vert->M[i][k];
        // time averaging with the old step to slow down the time varying
        grid->dzz[i][k]=thetaT*grid->dzz[i][k]+(1-thetaT)*grid->dzzold[i][k];
        grid->dzz[i][k]=(phys->h[i]+grid->dv[i])/((1-thetaT)*phys->h_old[i]+thetaT*phys->h[i]+grid->dv[i])*grid->dzz[i][k];
        // horizontal averaging
        for(nf=0;nf<grid->nfaces[i];nf++){
          if(grid->neigh[i*grid->maxfaces+nf]!=-1)
            vert->dzztmp[i][k]+=grid->dzz[grid->neigh[i*grid->maxfaces+nf]][k]/Neigh;
        }
        grid->dzz[i][k]=grid->dzz[i][k]*(1-thetaS)+vert->dzztmp[i][k]*thetaS;
        sum+=grid->dzz[i][k];
     }
     for(k=0;k<grid->Nk[i];k++){
       grid->dzz[i][k]=(phys->h[i]+grid->dv[i])/sum*grid->dzz[i][k];

     }
   

   }

}

/*
 * Function: FindBottomLayer
 * Calculate Nkeb and zfb for each edge based on dzf
 * ----------------------------------------------------
 * these two values will be further used in SetDragCoefficients and calculate drag terms
 */
 void FindBottomLayer(gridT *grid, propT *prop, physT *phys, int myproc)
 {
    int j,k,Nkeb,nc1,nc2;
    REAL zfb;
    // find 
    for(j=0;j<grid->Ne;j++)
    {
      nc1=grid->grad[2*j];
      nc2=grid->grad[2*j+1];
      if(nc1==-1)
        nc1=nc2;
      if(nc2==-1)
        nc2=nc1;
      // if constant drag coefficient 
      // use the bottom layer as the drag layer
      //if(prop->z0B!=0)
      //{

        zfb=0;
        for(k=grid->Nke[j]-1;k>=grid->etop[j];k--)
        {
           zfb+=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
           // find the layer when zfb>buffer height
           if(zfb>BUFFERHEIGHT){
             vert->Nkeb[j]=k;
             vert->zfb[j]=zfb;
             break;
           }
           // if not find at the top layer, it means the total dzf is smaller than bufferheight=dry
           // assume the drag layer at top layer
           if(k==0)
           {
             vert->Nkeb[j]=k;
             vert->zfb[j]=zfb;
             break;           
           } 
           zfb+=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
        }
      //} else {
        //vert->Nkeb[j]=grid->Nke[j]-1;
        //vert->zfb[j]=grid->dzf[j][grid->Nke[j]-1]/2;
      //}
        if(prop->z0B==0){
          vert->Nkeb[j]=grid->Nke[j]-1;
          vert->zfb[j]=0.5*(grid->dzz[nc1][grid->Nke[j]-1]+grid->dzz[nc2][grid->Nke[j]-1]);
        }
    }
 }

/*
 * Function: ComputeUl
 * Calculate u v at the top and bottom face of each layer
 * ----------------------------------------------------
 * compute ul vl
 */
void ComputeUl(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,j,k;
  for(i=0;i<grid->Nc;i++)
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
       vert->ul[i][k]=InterpToLayerTopFace(i, k, phys->uc, grid);
       vert->vl[i][k]=InterpToLayerTopFace(i, k, phys->vc, grid);
    }
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
  // compute omegac 
  for(i=0;i<grid->Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      phys->wc[i][k]=(phys->w[i][k]+phys->w[i][k+1])/2;
      vert->omegac[i][k]=(vert->omega[i][k]+vert->omega[i][k+1])/2;
    }

  // uf vf and wf
  for(j=0;j<grid->Ne;j++)
    for(k=grid->etop[j];k<grid->Nke[j];k++)
    {
      vert->uf[j][k]=InterpToFace(j,k,phys->uc,phys->u,grid);
      vert->vf[j][k]=InterpToFace(j,k,phys->vc,phys->u,grid);
      vert->wf[j][k]=InterpToFace(j,k,phys->wc,phys->u,grid);
      vert->omegaf[j][k]=InterpToFace(j,k,vert->omegac,phys->u,grid);
    }
}

/*
 * Function: InterpToLayerTopFace
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
    Dj = grid->dzz[i][k]/2+grid->dzz[i][k+1]/2;
    def1=0;
    def2=grid->dzz[i][k]/2;    
  }

  if(def2==0 && k<grid->Nk[i]-1 && k!=grid->ctop[i]) {
    // means the interior layer is dry
    return phi[i][k+1];
  } else if(def2==0 && k==grid->Nk[i]-1){
    // means the bottom layer is dry
    return phi[i][k];
  } else if(k==grid->ctop[i]) {
    return phi[i][k]+def2*(phi[i][k]-phi[i][k+1])/Dj;
  } else {
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
void LayerAveragedContinuity(REAL **omega, gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k,nf,ne;
  REAL fac1,fac2,fac3,flux,Ac;
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(i=0;i<grid->Nc;i++)
  {
    // compute omega from bottom layer using the omega_bot=0
    // worry about the numerical accuracy to ensure omega_top==0
    Ac=grid->Ac[i];
    omega[i][grid->Nk[i]]=0;
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
      flux=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(k<grid->Nke[ne])
          flux+=(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
      }
      if(!prop->subgrid)
        omega[i][k]=(fac1*omega[i][k+1]+fac2*vert->omega_old[i][k+1]+fac3*vert->omega_old2[i][k+1])-
        (fac2*vert->omega_old[i][k]+fac3*vert->omega_old2[i][k])-flux-
        (grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;
      else{
        omega[i][k]=(fac1*omega[i][k+1]*subgrid->Acveff[i][k+1]+
          fac2*vert->omega_old[i][k+1]*subgrid->Acveffold[i][k+1]+
          fac3*vert->omega_old2[i][k+1]*subgrid->Acveffold2[i][k+1])-
        (fac2*vert->omega_old[i][k]*subgrid->Acveffold[i][k]+
          fac3*vert->omega_old2[i][k]*subgrid->Acveffold2[i][k])-flux*Ac-
        (subgrid->Acceff[i][k]*grid->dzz[i][k]-subgrid->Acceffold[i][k]*grid->dzzold[i][k])/prop->dt;        
        omega[i][k]/=subgrid->Acveff[i][k];
      }
      omega[i][k]/=fac1;
    }
    //compute omega_im for updatescalar function
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      if(!prop->subgrid)
        vert->omega_im[i][k]=fac1*omega[i][k]+fac2*vert->omega_old[i][k]+fac3*vert->omega_old2[i][k];
      else
        vert->omega_im[i][k]=(fac2*vert->omega_old[i][k]*subgrid->Acveffold[i][k]+
                fac3*vert->omega_old2[i][k]*subgrid->Acveffold2[i][k]+
                fac1*omega[i][k]*subgrid->Acveff[i][k])/subgrid->Acveff[i][k];
    }
    vert->omega_im[i][grid->Nk[i]]=0;
  }
}

/*
 * Function: ComputeVerticalVelocity
 * Compute w for for the top and bottom face of each layer
 * from the layer-averaged continuity equation
 * ----------------------------------------------------
 * 
 */
void ComputeVerticalVelocity(REAL **omega, REAL **zc, gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k,nf,ne;
  REAL fac1,fac2,fac3,flux,Ac;
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(i=0;i<grid->Nc;i++)
  {
    Ac=grid->Ac[i];
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--){
      flux=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(k<grid->Nke[ne])
          flux+=(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
      }
      omega[i][k]=vert->omega[i][k+1]-flux-(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;
    }
  }
}

/*
 * Function: ComputeZc
 * Compute the vertical location of each cell center
 * also compute the vertical location of each layer center at edge face
 * ----------------------------------------------------
 */
void ComputeZc(gridT *grid, propT *prop, physT *phys, int index, int myproc)
{
  int i,j,k,nc1,nc2;
  REAL z,def1,def2,Dj;
  for(i=0;i<grid->Nc;i++)
  {
    if(!index)
      for(k=0;k<grid->Nk[i];k++)
        vert->zcold[i][k]=vert->zc[i][k];
    z=-grid->dv[i];
    vert->zc[i][grid->Nk[i]-1]=z+grid->dzz[i][grid->Nk[i]-1]/2;
    for(k=grid->Nk[i]-2;k>=grid->ctop[i];k--){
      vert->zc[i][k]=vert->zc[i][k+1]+grid->dzz[i][k+1]/2+grid->dzz[i][k]/2;
    }
    if(prop->n==0)
      for(k=0;k<grid->Nk[i];k++)
        vert->zcold[i][k]=vert->zc[i][k];
  } 
  for(j=0;j<grid->Ne;j++)
  {
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;
    Dj = grid->dg[j];
    def1 = sqrt(pow((grid->xv[nc1]-grid->xe[j]),2)+pow((grid->yv[nc1]-grid->ye[j]),2));
    def2 = Dj-def1;
    for(k=0;k<grid->Nke[j];k++){
      vert->zf[j][k]=vert->zc[nc1][k]*def2/Dj+vert->zc[nc2][k]*def1/Dj;
    }

  }
}

/*
 * Function: ComputeOmega
 * Compute omega from the definition omega=w-udzdx-vdzdy-wg and U3=w-udzdx-vdzdy
 * ----------------------------------------------------
 * index=1, w->omega index=0 omega->w index=-1 w->U3
 * index=1 is never used in phys.c
 * omega is always calculated from continuity
 */
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int index, int myproc)
{
  int i,k;
  REAL alpha=1;
  if(prop->vertcoord==5)
    alpha=0;
  if(index==1)
    for(i=0;i<grid->Nc;i++)
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
      {
        vert->omega[i][k]=phys->w[i][k]-alpha*vert->ul[i][k]*InterpToLayerTopFace(i,k,vert->dzdx,grid)-\
        alpha*vert->vl[i][k]*InterpToLayerTopFace(i,k,vert->dzdy,grid)-((vert->zc[i][k]-vert->zcold[i][k])+(grid->dzz[i][k]-grid->dzzold[i][k])/2)/prop->dt;
      }
  else if(index==0)
    for(i=0;i<grid->Nc;i++)
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
      {
        phys->w[i][k]=vert->omega[i][k]+alpha*vert->ul[i][k]*InterpToLayerTopFace(i,k,vert->dzdx,grid)+\
        alpha*vert->vl[i][k]*InterpToLayerTopFace(i,k,vert->dzdy,grid)+((vert->zc[i][k]-vert->zcold[i][k])+(grid->dzz[i][k]-grid->dzzold[i][k])/2)/prop->dt;
      }
  else
    for(i=0;i<grid->Nc;i++)
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
      {
        vert->U3[i][k]=phys->w[i][k]-alpha*vert->ul[i][k]*InterpToLayerTopFace(i,k,vert->dzdx,grid)-\
        alpha*vert->vl[i][k]*InterpToLayerTopFace(i,k,vert->dzdy,grid);
      }
}

/*
 * Function: VerticalCoordinateHorizontalSource
 * Compute all the preparation due to the new general vertical
 * coordinate for the horizontal source of momentum conservation
 * ----------------------------------------------------
 * 
 */
void VertCoordinateHorizontalSource(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, MPI_Comm comm)
{
   int i,k;
   // compute the relative vorticity
   // 1. compute v and u at each edge center
   ComputeUf(grid, prop, phys,myproc);

   // comment out due to use conservative method
   // 2. compute dvdx and dudy
   //ComputeCellAveragedHorizontalGradient(vert->dvdx, 0, vert->vf, grid, prop, phys, myproc);
   //ComputeCellAveragedHorizontalGradient(vert->dudy, 1, vert->uf, grid, prop, phys, myproc);

   // 3. compute f_r
   //for(i=0;i<grid->Nc;i++)
    //for(k=0;k<grid->Nk[i];k++)
      //vert->f_r[i][k]=vert->dvdx[i][k]-vert->dudy[i][k];
  
   /*if(prop->nonhydrostatic){
     ComputeCellAveragedHorizontalGradient(vert->dwdy, 1, vert->wf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dwdx, 0, vert->wf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dvdy, 1, vert->vf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dudx, 0, vert->uf, grid, prop, phys, myproc);
   }*/
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
        ne = grid->face[i*grid->maxfaces+nf];
        if(direction==0)
          vec=vert->n1[i*grid->maxfaces+nf];
        else
          vec=vert->n2[i*grid->maxfaces+nf];
        gradient[i][k]+=scalar[ne][k]*grid->df[ne]*vec;
      }
      gradient[i][k]/=grid->Ac[i];
    }
  }
}
/*
 * Function: ComputeNormalVector
 * Compute the outpointing normalized vector from cell center to each face center
 * ----------------------------------------------------
 * n1 is the x-component and n2 is y-component
 */
void ComputeNormalVector(gridT *grid, physT *phys, int myproc)
{ 
  int i,nf,ne;
  REAL D,Dx,Dy;
  for(i=0;i<grid->Nc;i++)
    for(nf=0;nf<grid->nfaces[i];nf++)
    {
      ne=grid->face[i*grid->maxfaces+nf];
      D=(grid->xv[i]-grid->xe[ne])*(grid->xv[i]-grid->xe[ne])+
      (grid->yv[i]-grid->ye[ne])*(grid->yv[i]-grid->ye[ne]);
      D=sqrt(D);
      Dx=grid->xe[ne]-grid->xv[i];
      Dy=grid->ye[ne]-grid->yv[i];
      vert->n1[i*grid->maxfaces+nf]=Dx/D;
      vert->n2[i*grid->maxfaces+nf]=Dy/D;
    }
}

/* 
 * Function: OpenVertCoordinateFiles
 * Usage: OpenFiles(myproc);
 * ------------------------------
 * Open all of the files used for i/o to store the file pointers.
 * include dzz, zc, omega for each cell at different time steps
 * 
 */
void OpenVertCoordinateFiles(gridT *grid,int mergeArrays, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  MPI_GetFile(filename,DATAFILE,"zcFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  vert->zcFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"dzzFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  vert->dzzFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"omegaFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  vert->omegaFID = MPI_FOpen(str,"w","OpenFiles",myproc);
}

/*
 * Function: OutputVertcoordinate
 * Usage: OutputVertcoordinate(grid,phys,prop,myproc,numprocs,comm);
 * ---------------------------------------------------------------------------
 * Output the data every ntout steps as specified in suntans.dat
 * now include zc,dzz,omega
 *
 */
void OutputVertCoordinate(gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, j, jptr, k, nwritten, arraySize, writeProc,nc1,nc2;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart) { 
    Write3DData(vert->zc,vert->tmp,prop->mergeArrays,vert->zcFID,
    "Error outputting uc-data!\n",grid,numprocs,myproc,comm);
    Write3DData(grid->dzz,vert->tmp,prop->mergeArrays,vert->dzzFID,
    "Error outputting uc-data!\n",grid,numprocs,myproc,comm);
    Write3DData(vert->omega,vert->tmp,prop->mergeArrays,vert->omegaFID,
    "Error outputting uc-data!\n",grid,numprocs,myproc,comm);
  }

  if(prop->n==prop->nsteps+prop->nstart && myproc==0) {
    fclose(vert->zcFID);
    fclose(vert->dzzFID);
    fclose(vert->omegaFID);
  }
}

/*
 * Function: VertCoordinateBasic
 * Setup the basics and initial condition for the vertical coordinate
 * ----------------------------------------------------
 * 
 */
void VertCoordinateBasic(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k;
  // allocate subgrid struture first
  vert=(vertT *)SunMalloc(sizeof(vertT),"VertCoordinateBasic");

  // allocate necessary variable
  AllocateandInitializeVertCoordinate(grid,prop, myproc);

  // output 
  OpenVertCoordinateFiles(grid,prop->mergeArrays,myproc);

  // if not z-level setup dzz
  InitializeLayerThickness(grid, prop, phys,myproc);

  // compute zc
  ComputeZc(grid,prop,phys,0,myproc);

  // compute normal vector
  ComputeNormalVector(grid,phys,myproc);
  
  // compute cell centered dzdx and dzdy
  ComputeCellAveragedHorizontalGradient(vert->dzdx, 0, vert->zf, grid, prop, phys, myproc);

  ComputeCellAveragedHorizontalGradient(vert->dzdy, 1, vert->zf, grid, prop, phys, myproc);
}

/*
 * Function: StoreVariables
 * Usage: StoreVariables(grid,phys);
 * ---------------------------------
 * Store the old values of s, u, and w into stmp3, u_old, and wtmp2,
 * respectively.
 *
 */
void StoreVertVariables(gridT *grid, physT *phys) {
  int i, j, k, iptr, jptr;

  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) {
      // store omega^n-1 and omega^n
      vert->omega_old2[i][k]=vert->omega_old[i][k];
      vert->U3_old2[i][k]=vert->U3_old[i][k];
      vert->omega_old[i][k]=vert->omega[i][k];
      vert->U3_old[i][k]=vert->U3[i][k];
    }
}

/*
 * Function: UpdateCellCenteredFreeSurface
 * Usage: UpdateCellCenteredFreeSurface(grid,prop,phys,myproc);
 * ---------------------------------
 * update the free surface height after the nonhydrostatic pressure is solved
 *
 */
void UpdateCellcenteredFreeSurface(gridT *grid, propT *prop, physT *phys, int myproc)
{
   int i,iptr,k,nf,ne;
   REAL sum,normal,tmp;
   for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    sum = 0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
      {  
        sum+=(prop->imfac1*phys->u[ne][k]+prop->imfac2*phys->u_old[ne][k]+prop->imfac3*phys->u_old2[ne][k])*
          grid->df[ne]*normal*grid->dzf[ne][k];
      }
    }
    tmp=phys->h_old[i]-sum*prop->dt/grid->Ac[i];
    phys->h[i]=tmp;
  }
}

/*
 * Function: UpdateFluxHeight
 * Usage: UpdateFluxHeight(grid,prop,phys,myproc);
 * ---------------------------------
 * Update the flux height use the upwind the scheme with u*, u^n-1, u^n-2
 * The reason is to ensure a positive layer height for the isopycnal coordinate
 * The flux height set by the SetFluxHeight is calculated by u^n-1 which need to be verified by
 * u^im which is used to calculate the flux into/out of a grid cell.
 * The new flux height will be used to further modified free surface height/ layer thickness
 * only work when vertcoord==2
 */
void VerifyFluxHeight(gridT *grid, propT *prop, physT *phys, int myproc)
{  
  int i, j, k, nc1, nc2;
  REAL dz_bottom, dzsmall=grid->dzsmall,h_uw,utmp,tmp;

  for(j=0;j<grid->Ne;j++) 
  {
    grid->hf[j]=0;
    //for(k=0; k<grid->Nkc[j];k++)
      //grid->dzf[j][k]=0;
  }

  for(j=0;j<grid->Ne;j++) {
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    for(k=0;k<grid->etop[j];k++)
      grid->dzf[j][k]=0;
 
    for(k=grid->etop[j];k<grid->Nke[j];k++){
      utmp=prop->imfac1*phys->u[j][k]+prop->imfac2*phys->u_old[j][k]+prop->imfac3*phys->u_old2[j][k];
      tmp=UpWind(utmp,grid->dzz[nc1][k],grid->dzz[nc2][k]);
      grid->dzf[j][k]=tmp;
    }

    /* This works with Wet/dry but not with cylinder case...*/
    if(grid->etop[j]==grid->Nke[j]-1) {
      // added part
     if(grid->mark[j]==2 && grid->dzf[j][k]<=0.01)
        grid->dzf[j][k]=0.01;
    }

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(grid->dzf[j][k]<=DRYCELLHEIGHT)
        grid->dzf[j][k]=0;
      
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      grid->hf[j]+=grid->dzf[j][k];
  }
}

