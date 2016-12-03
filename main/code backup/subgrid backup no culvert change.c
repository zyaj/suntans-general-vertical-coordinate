/*
 * File: subgrid.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file contains the functions to use subgrid scale model method
 * based on Casulli paper (2009)
 * Now it can only work with quad and triangle
 * Now it can only work when Nk=1
 * 08/01/2014
 */
#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "initialization.h"
#include "util.h"
#include "culvert.h"
#include "mympi.h"
#include "subgrid.h"

void ReadSubgridProperties(int myproc);
void AllocateandInitializeSubgrid(gridT *grid, int myproc);
void CalculateSubgridXY(gridT *grid, int myproc);
void CalculateSubgridCellp(gridT *grid, int myproc);
void InterpolateSubgridDepth(gridT *grid,physT *phys, int myproc);
void CalculateVolumeProfile(gridT *grid, int myproc);
void CalculateHProfile(gridT *grid, int myproc);
void CalculateAcProfile(gridT *grid, int myproc);
void CalculateVolumeProfile(gridT *grid, int myproc);
void CalculateFluxHeightProfile(gridT *grid, int myproc);
void CalculateWetperimeterProfile(gridT *grid, int myproc);
void SubgridCheck(gridT *grid, int myproc);
REAL CalculateWetArea(int nc, int nfaces, int maxfaces, REAL h);
REAL CalculateWetperimeter(int ne, REAL h);
REAL CalculateFluxHeight(int ne, REAL h);
REAL UpdateFluxHeight(int ne, REAL h);
REAL UpdateWetperi(int ne, REAL h);
REAL UpdateAceff(int nc, REAL h);
REAL UpdateVeff(int nc, REAL h);

void OutputSubgrid(gridT *grid, physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm);

/*
 * Function: ReadSubgridProperties
 * Usage: Read Subgrid Properties
 * ----------------------------------------------------
 * Based on suntans.dat, load in the important parameters for
 * subgrid model
 */
void ReadSubgridProperties(int myproc) {
  subgrid->segN = MPI_GetValue(DATAFILE,"segN","ReadSubgridProperties",myproc); 
  subgrid->disN = MPI_GetValue(DATAFILE,"disN","ReadSubgridProperties",myproc); 
  subgrid->meth = MPI_GetValue(DATAFILE,"subgridmeth","ReadSubgridProperties",myproc); 
  subgrid->dpint =MPI_GetValue(DATAFILE,"subgriddpint","ReadSubgridProperties",myproc); 
}

/*
 * Function: AllocateandInitializeSubgrid
 * Usage: Allocate Subgrid variable
 * ----------------------------------------------------
 * allocate space and initialize(set all zero) for all subgrid variables
 *
 */
void AllocateandInitializeSubgrid(gridT *grid, int myproc) {
  int i,N,Ntotal,k;

  N=grid->Nc*(grid->maxfaces-2);
  // allocate xp yp dp for subgrid scal
  Ntotal=N*(subgrid->segN+1)*(subgrid->segN+2)/2;
  subgrid->xp = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->yp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->dp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  for(i=0;i<Ntotal;i++)
  {
    subgrid->xp[i]=0.0;
    subgrid->yp[i]=0.0;
    subgrid->dp[i]=0.0;
  } 
 
  // allocate xpe ype dpe for each edge
  Ntotal=grid->Ne*(subgrid->segN+1);
  subgrid->xpe = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->ype = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->dpe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  for(i=0;i<Ntotal;i++)
  {
    subgrid->xpe[i]=0.0;
    subgrid->ype[i]=0.0;
    subgrid->dpe[i]=0.0;
  }  

  // allocate subcellp
  Ntotal=N*(pow(subgrid->segN,2))*3;
  subgrid->cellp = (int *)SunMalloc(Ntotal*sizeof(int), "AllocateSubgrid");
  for(i=0;i<Ntotal;i++)
    subgrid->cellp[i]=0;

  // allocate profile for effective area and volume for each cell
  Ntotal=grid->Nc*(subgrid->disN+1);
  subgrid->hprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->Vprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->Acprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->hmin = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSubgrid");

  for(i=0;i<Ntotal;i++)
  {
    subgrid->hprof[i]=0.0;
    subgrid->Vprof[i]=0.0;
    subgrid->Acprof[i]=0.0;
  }
 
  for(i=0;i<grid->Nc;i++)
    subgrid->hmin[i]=0.0;


  // allocate profile for flux height profile
  Ntotal=grid->Ne*(subgrid->disN+1);
  subgrid->fluxhprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->hprofe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->Wetperiprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->dzboteff=(REAL *)SunMalloc(grid->Ne*sizeof(REAL),"AllocateSubgrid");
  
  for(i=0;i<Ntotal;i++){
    subgrid->fluxhprof[i]=0.0;
    subgrid->hprofe[i]=0.0;
    subgrid->Wetperiprof[i]=0.0;
  }
  for(i=0;i<grid->Ne;i++)
    subgrid->dzboteff[i]=0.0;  

  // allocate the backup of original Ac and effective wet Ac
  subgrid->Acbackup = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Acwet = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Aceff = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Aceffold = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Veff = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Acceff = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acceffold = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acveffold = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acveff = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");

  for(i=0;i<grid->Nc;i++)
  {
    subgrid->Acbackup[i]=grid->Ac[i];
    subgrid->Acwet[i]=grid->Ac[i];
    subgrid->Aceff[i]=grid->Ac[i];
    subgrid->Aceffold[i]=grid->Ac[i];
    subgrid->Veff[i]=0;
    subgrid->Acceff[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acceffold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acveffold[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acveff[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateSubgrid");
    for(k=0;k<grid->Nk[i];k++)
    {
      subgrid->Acceff[i][k]=0;
      subgrid->Acceffold[i][k]=0;
      subgrid->Acveff[i][k]=0;
      subgrid->Acveffold[i][k]=0;   
    }
    subgrid->Acveff[i][grid->Nk[i]]=0;
    subgrid->Acveffold[i][grid->Nk[i]]=0;
  }
}

/*
 * Function:  CalculateSubgridXY
 * Usage: calculate all the points XY for each subcell
 * ----------------------------------------------------
 * the index is from the bottom to top from left to right 
 *
 */
void CalculateSubgridXY(gridT *grid, int myproc)
{
  int i,j,k,l,m,base;
  REAL dx1,dx2,dy1,dy2;

  // calculate x y for sub cell
  for(i=0;i<grid->Nc;i++)
  {
    m=0;
    for(j=0;j<(grid->nfaces[i]-2);j++)
    {
      if(j==0)
      {
        base=grid->cells[i*grid->maxfaces];
        dx1=(grid->xp[grid->cells[i*grid->maxfaces+1]]-grid->xp[base])/subgrid->segN;
        dx2=(grid->xp[grid->cells[i*grid->maxfaces+2]]-grid->xp[base])/subgrid->segN;
        dy1=(grid->yp[grid->cells[i*grid->maxfaces+1]]-grid->yp[base])/subgrid->segN;
        dy2=(grid->yp[grid->cells[i*grid->maxfaces+2]]-grid->yp[base])/subgrid->segN;        
      } else {
        base=grid->cells[i*grid->maxfaces+2];
        dx1=(grid->xp[grid->cells[i*grid->maxfaces+3]]-grid->xp[base])/subgrid->segN;
        dx2=(grid->xp[grid->cells[i*grid->maxfaces]]-grid->xp[base])/subgrid->segN;
        dy1=(grid->yp[grid->cells[i*grid->maxfaces+3]]-grid->yp[base])/subgrid->segN;
        dy2=(grid->yp[grid->cells[i*grid->maxfaces]]-grid->yp[base])/subgrid->segN; 
      }
      for(k=0;k<=subgrid->segN;k++)
      {
        for(l=k;l<=subgrid->segN;l++)
        {
          subgrid->xp[(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2*i+m]=grid->xp[base]+dx1*(l-k)+dx2*(k);
          subgrid->yp[(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2*i+m]=grid->yp[base]+dy1*(l-k)+dy2*(k);
          m++;
        }  
      }    
    }
    // check m correct
    if(m!=((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2)){
      printf("Error computing subgrid->xp/yp in function CalculateSubgridXY\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }  
    /*if(i==0){
      printf("%f %f %f %f \n",subgrid->xp[0],subgrid->yp[0],subgrid->xp[1],subgrid->yp[1]); 
      printf("%f %f %f %f \n",subgrid->xp[2],subgrid->yp[2],subgrid->xp[3],subgrid->yp[3]); 
      printf("%f %f %f %f \n",subgrid->xp[4],subgrid->yp[4],subgrid->xp[5],subgrid->yp[5]); 
      printf("%d %d %d %d \n",grid->cells[0],grid->cells[1],grid->cells[2],grid->cells[3]);
    }*/
  }

  
  // calculate X Y for subedge
  for(i=0;i<grid->Ne;i++)
  {
    base=grid->edges[NUMEDGECOLUMNS*i];
    dx1=(grid->xp[grid->edges[NUMEDGECOLUMNS*i+1]]-grid->xp[base])/subgrid->segN;
    dy1=(grid->yp[grid->edges[NUMEDGECOLUMNS*i+1]]-grid->yp[base])/subgrid->segN;
    for(j=0;j<=subgrid->segN;j++)
    {
      subgrid->xpe[i*(subgrid->segN+1)+j]=grid->xp[base]+dx1*j;
      subgrid->ype[i*(subgrid->segN+1)+j]=grid->yp[base]+dy1*j;      
    }
  }
}

/*
 * Function: CalculateSubgridCellp
 * Usage: calculate all the point index group for each subcell
 * ----------------------------------------------------
 * from bottom to top
 *
 */
void CalculateSubgridCellp(gridT *grid, int myproc)
{
  int nc,nf,i,j,kk,N,base;
  N=subgrid->segN;
  for(nc=0;nc<grid->Nc;nc++)
  {
    kk=0;
    for(nf=0;nf<(grid->nfaces[nc]-2);nf++)
    {
      for(i=1;i<=N;i++)
      {
        // calculate the index for the starting point for each layer 
        if(i==1)
          base=(nc*(grid->maxfaces-2)+nf)*(N+2)*(N+1)/2;
        else
          base+=N+3-i; 
 
        for(j=1;j<=(2*N-2*i+1);j++)
        {
          if(j<=(N+1-i))
          {
            subgrid->cellp[nc*(grid->maxfaces-2)*N*N*3+kk*3]=base+j-1;
            subgrid->cellp[nc*(grid->maxfaces-2)*N*N*3+kk*3+1]=base+j;
            subgrid->cellp[nc*(grid->maxfaces-2)*N*N*3+kk*3+2]=base+N+1-i+j;
            kk++;
          } else {
            subgrid->cellp[nc*(grid->maxfaces-2)*N*N*3+kk*3]=base+j-(N+1-i);
            subgrid->cellp[nc*(grid->maxfaces-2)*N*N*3+kk*3+1]=base+j;
            subgrid->cellp[nc*(grid->maxfaces-2)*N*N*3+kk*3+2]=base+j+1;
            kk++;
          }
        }
      }
    }
    if(kk!=(grid->nfaces[nc]-2)*N*N)
    {
      printf("Error computing subgrid->cellp in function CalculateSubgridCellp\n");
      MPI_Finalize();
      exit(EXIT_FAILURE); 
    }
  }
}


/*
 * Function: InterpolateSubgridDepth
 * Usage: calculate depth for all subpoints
 * ----------------------------------------------------
 * input the same file as in grid.c to interpolate 
 * depth for all sub points
 *
 */
void InterpolateSubgridDepth(gridT *grid, physT *phys, int myproc)
{
  int i, j, n, nc,nc1,nc2,ne, Nd, scaledepth,ncount,base,k,intdd;
  REAL *xd, *yd, *d, scaledepthfactor,depthelev,intd,min;
  char str[BUFFERLENGTH];
  FILE *ifile;

  if(subgrid->dpint==1)
  {
    scaledepth=(int)MPI_GetValue(DATAFILE,"scaledepth","InterpolateSubgridDepth",myproc);
    scaledepthfactor=MPI_GetValue(DATAFILE,"scaledepthfactor","InterpolateSubgridDepth",myproc);
    depthelev=MPI_GetValue(DATAFILE,"depthelev","InterpolateSubgridDepth",myproc);

    Nd = MPI_GetSize(INPUTDEPTHFILE,"InterpolateSubgridDepth",myproc);
    xd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridDepth");
    yd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridDepth");
    d = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridDepth");

    ifile = MPI_FOpen(INPUTDEPTHFILE,"r","InterpolateSubgridDepth",myproc);
    for(n=0;n<Nd;n++) {
      xd[n]=getfield(ifile,str);
      yd[n]=getfield(ifile,str);
      d[n]=getfield(ifile,str);//fabs(getfield(ifile,str));
      if(scaledepth)
        d[n]*=scaledepthfactor;
      d[n]+=depthelev;    
    }

    // interpolate subcell depth
    ncount=grid->Nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
    Interp(xd,yd,d,Nd,&(subgrid->xp[0]), &(subgrid->yp[0]),&(subgrid->dp[0]),ncount,grid->maxfaces);

    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        intd=subgrid->dp[base+n];
        k=0;

        // adjust accuracy as the output grid.c
        if(intd>=1)
        {
          while(floor(intd)!=0)
          {
            intd=intd/10.0;
            k++;
          }
          intdd=floor(subgrid->dp[base+n]*1e6/pow(10,k-1));
          if((subgrid->dp[base+n]*1e6/pow(10,k-1)-intdd)>=0.5)
            intdd+=1;

          subgrid->dp[base+n]=intdd/1e6*pow(10,k-1);
        } else if(intd<1 && intd!=0) {

          while(floor(intd)==0)
          {
            intd=intd*10.0;
            k++;
          }
          intdd=floor(subgrid->dp[base+n]*1e6*pow(10,k));
          if((subgrid->dp[base+n]*1e6*pow(10,k)-intdd)>=0.5)
            intdd+=1;

          subgrid->dp[base+n]=intdd/1e6/pow(10,k);         
                    
        }
      }
    }


    // get subgrid->hmin
    for(i=0;i<grid->Nc;i++)
    {
      for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
      {
        if(subgrid->hmin[i]>-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmin[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
      }
    }

    // interpolate subedge depth
    ncount=grid->Ne*(subgrid->segN+1);
    Interp(xd,yd,d,Nd,&(subgrid->xpe[0]), &(subgrid->ype[0]),&(subgrid->dpe[0]),ncount,grid->maxfaces);

    for(n=0;n<ncount;n++)
    {
      intd=subgrid->dpe[n];
      k=0;
      if(intd>=1)
      {
        while(floor(intd)!=0)
        {
          intd=intd/10.0;
          k++;
        }
        intdd=floor(subgrid->dpe[n]*1e6/pow(10,k-1));
        if((subgrid->dpe[n]*1e6/pow(10,k-1)-intdd)>=0.5)
          intdd+=1;

        subgrid->dpe[n]=intdd/1e6*pow(10,k-1);

      } else if(intd<1 && intd!=0){
        while(floor(intd)==0)
        {
          intd=intd*10.0;
          k++;
        }
        intdd=floor(subgrid->dpe[n]*1e6*pow(10,k));
        if((subgrid->dpe[n]*1e6*pow(10,k)-intdd)>=0.5)
          intdd+=1;
 
        subgrid->dpe[n]=intdd/1e6/pow(10,k);  
      }      
    }

  } else {

    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        subgrid->dp[base+n]=ReturnSubgridPointDepth(subgrid->xp[base+n],subgrid->yp[base+n],grid->xv[nc],grid->yv[nc]);
        // added
        //if(subgrid->dp[base+n]>grid->dv[nc])
          //subgrid->dp[base+n]=grid->dv[nc];
        intd=subgrid->dp[base+n];
        k=0;

        // adjust accuracy as the output grid.c
        if(intd>=1)
        {
          while(floor(intd)!=0)
          {
            intd=intd/10.0;
            k++;
          }
          intdd=floor(subgrid->dp[base+n]*1e6/pow(10,k-1));
          if((subgrid->dp[base+n]*1e6/pow(10,k-1)-intdd)>=0.5)
            intdd+=1;

          subgrid->dp[base+n]=intdd/1e6*pow(10,k-1);
        } else if(intd<1 && intd!=0) {

          while(floor(intd)==0)
          {
            intd=intd*10.0;
            k++;
          }
          intdd=floor(subgrid->dp[base+n]*1e6*pow(10,k));
          if((subgrid->dp[base+n]*1e6*pow(10,k)-intdd)>=0.5)
            intdd+=1;

          subgrid->dp[base+n]=intdd/1e6/pow(10,k);         
                    
        }

      }
    }

    // get subgrid->hmin
    for(i=0;i<grid->Nc;i++)
    {
      for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
      {
        if(subgrid->hmin[i]>-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmin[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
      }
    }


    ncount=grid->Ne*(subgrid->segN+1);
    for(n=0;n<ncount;n++)
    {
      /*ne=floor(n/(subgrid->segN+1));
      nc1=grid->grad[ne*2];
      nc2=grid->grad[ne*2+1];
      if(nc1==-1)
        nc1=nc2;
      else if(nc2==-1)
        nc2=nc1;

      
      min=grid->dv[nc1];
      if(grid->dv[nc2]<min)
        min=grid->dv[nc2];
       */
      subgrid->dpe[n]=ReturnSubgridPointeDepth(subgrid->xpe[n],subgrid->ype[n]);
   
      //if(subgrid->dpe[n]>min)
       // subgrid->dpe[n]=min;

      intd=subgrid->dpe[n];
      k=0;
      if(intd>=1)
      {
        while(floor(intd)!=0)
        {
          intd=intd/10.0;
          k++;
        }
        intdd=floor(subgrid->dpe[n]*1e6/pow(10,k-1));
        if((subgrid->dpe[n]*1e6/pow(10,k-1)-intdd)>=0.5)
          intdd+=1;

        subgrid->dpe[n]=intdd/1e6*pow(10,k-1);

      } else if(intd<1 && intd!=0){
        while(floor(intd)==0)
        {
          intd=intd*10.0;
          k++;
        }
        intdd=floor(subgrid->dpe[n]*1e6*pow(10,k));
        if((subgrid->dpe[n]*1e6*pow(10,k)-intdd)>=0.5)
          intdd+=1;
 
        subgrid->dpe[n]=intdd/1e6/pow(10,k);  
      }      

    }

  }
}

/*
 * Function: CalculateWetArea
 * Usage: calculate total wet area for different h
 * ----------------------------------------------------
 * meth=1 use regular subgrid method (linear interpolation)
 * meth=2 use cell center depth (average for all three points) to 
 * meth=3 use user defined function 
 *
 */
REAL CalculateWetArea(int nc, int nfaces, int maxfaces, REAL h)
{
  REAL ac1,Ac,ac,x1,y1,x2,y2,x3,y3,d1,d2,d3,x,y,d,r;
  int i,ntotal,N; 
  N=subgrid->segN;
  ntotal=(nfaces-2)*(pow(N,2));
  ac1=subgrid->Acbackup[nc]/ntotal;
  Ac=0;
  for(i=0;i<ntotal;i++)
  { 
    if(subgrid->meth==1){
      // get x y d for all three points of each subcell
      d1=subgrid->dp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3]];
      d2=subgrid->dp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+1]];
      d3=subgrid->dp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+2]];

      // make them in order by d d1<=d2<=d3
      if(d1>d3)
      {
        d=d1;
        d1=d3;
        d3=d;
      } 
      if(d2>d3)
      {
        d=d2;
        d2=d3;
        d3=d;
      } 
      if(d2<d1)
      {
        d=d2;
        d2=d1;
        d1=d;
      } 
      // for different condition use different method to calculate ac
      if(d1==d2 && d2==d3 )//|| fabs(d3-d1)<=0.01)
      {
        if(h<-d3)
          ac=0;
        else
          ac=ac1;
      
      } else if(d1==d2){
        if(h>=-d1)
          ac=ac1;
        else if(h<-d1 && h>=-d3){
          r=pow(((h+d3)/(d3-d1)),2);
          ac=ac1*r;      
        } else {
          ac=0;
        }
      } else if(d2==d3){
        if(h>=-d1)
          ac=ac1;
        else if(h<-d1 && h>=-d2){
          r=pow((-(h+d1)/(d3-d1)),2);
          ac=ac1*(1-r);      
        } else {
          ac=0;
        }
      } else {
        if(h>=-d1)
          ac=ac1;
        else if(h<-d1 && h>=-d2)
        {
          r=pow((-(h+d1)/(d2-d1)),2);
          ac=ac1*(d3-d2)/(d3-d1)+ac1*(d2-d1)/(d3-d1)*(1-r);
        } else if(h<-d2 && h>=-d3) {
          r=pow(((h+d3)/(d3-d2)),2);
          ac=ac1*(d3-d2)/(d3-d1)*r;
        } else { 
          ac=0;
        }
      }
    } else if(subgrid->meth==2){
      for(i=0;i<ntotal;i++)
      {
        // get x y d for all three points of each subcell
        d1=subgrid->dp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3]];
        d2=subgrid->dp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+1]];
        d3=subgrid->dp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+2]];
        if(h>-(d1+d2+d3)/3)
          ac=ac1;
        else
          ac=0;
      }
    } else if(subgrid->meth==3){
      x1=subgrid->xp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3]];
      x2=subgrid->xp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+1]];
      x3=subgrid->xp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+2]];
      y1=subgrid->yp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3]];
      y2=subgrid->yp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+1]];
      y3=subgrid->yp[subgrid->cellp[nc*(maxfaces-2)*N*N*3+i*3+2]];
      ac=ReturnSubCellArea(x1,y1,x2,y2,x3,y3,h);   
    }
    Ac+=ac;
  }
  return Ac;
}



/*
 * Function: CalculateHProfile
 * Usage: the distribution of h for Ac prof
 * ----------------------------------------------------
 * calculate the dmax and dmin to distribute 
 *
 */
void CalculateHProfile(gridT *grid, int myproc)
{
   int i,nc,ntotal,ne,base;
   REAL min,max,dmax,dmin,dh;   

   // get h profile for each cell
   for(nc=0;nc<grid->Nc;nc++)
   {
     ntotal=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
     base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
     max=subgrid->dp[base];
     min=subgrid->dp[base];

     for(i=1;i<ntotal;i++)
     {
       if(subgrid->dp[base+i]>max)
       {
         max=subgrid->dp[base+i];
       }
       if(subgrid->dp[base+i]<min)
       {
         min=subgrid->dp[base+i];
       }
     }
     if(max!=min)
     {
       dh=(max-min)/subgrid->disN;
     } else {
       dh=0;
     }
     for(i=0;i<=(subgrid->disN);i++)
     {
       subgrid->hprof[nc*(subgrid->disN+1)+i]=-max+dh*i;
     }
   } 

   // get h profile for each edge
   for(ne=0;ne<grid->Ne;ne++)
   {
     base=ne*(subgrid->segN+1);
     max=subgrid->dpe[base];
     min=subgrid->dpe[base];
     for(i=1;i<(subgrid->segN+1);i++)
     {
       if(subgrid->dpe[base+i]>max)
       {
         max=subgrid->dpe[base+i];
       }
       if(subgrid->dpe[base+i]<min)
       {
         min=subgrid->dpe[base+i];
       } 
     }

     if(max!=min)
     {
       dh=(max-min)/subgrid->disN;
     } else {
       dh=0;
     }

     for(i=0;i<=(subgrid->disN);i++)
     {
       subgrid->hprofe[ne*(subgrid->disN+1)+i]=-max+dh*i;
     }
   }  
}


/*
 * Function: CalculateAcProfile
 * Usage: calculate Ac and h profile for each cell 
 * ----------------------------------------------------
 * use CalculateWetArea to get the total wet area for each h 
 *
 */
void CalculateAcProfile(gridT *grid, int myproc)
{
  int nc,i,base;
  
  // get the profile of Ac for each cell   
  for(nc=0;nc<grid->Nc;nc++)
  {
    base=nc*(subgrid->disN+1);
    subgrid->Acprof[base+subgrid->disN]=subgrid->Acbackup[nc];
    
    for(i=0;i<subgrid->disN;i++)
    {
      subgrid->Acprof[base+i]=CalculateWetArea(nc, grid->nfaces[nc], grid->maxfaces, subgrid->hprof[base+i]);
    }
  } 
}

/*
 * Function: CalculateFluxHeight
 * Usage: calculate flux height profile for each edge 
 * ----------------------------------------------------
 * calculate crossectional wet area first and then divide by grid->df[ne] 
 *
 */
REAL CalculateFluxHeight(int ne, REAL h)
{
  int i,base;
  REAL r,r1,R,d1,d2,d,x1,x2,y1,y2;
  R=0;
  if(subgrid->meth==1){
    for(i=0;i<subgrid->segN;i++)
    {
      base=ne*(subgrid->segN+1);
      d1=subgrid->dpe[base+i];
      d2=subgrid->dpe[base+i+1];
      if(d1>d2)
      {
        d=d1;
        d1=d2;
        d2=d;
      }
      if(d1==d2)
      { 
        if(h<=-d1)
          r=0;
        else
          r=(h+d1)/subgrid->segN;
      } else {
        if(h<=-d2)
          r=0;
        else if(h>-d2 && h<=-d1)
        { 
          r1=(h+d2)/(d2-d1);
          r=r1*(h+d2)/2/subgrid->segN;
        } else {
          r=((h+d1)+(d2-d1)/2)/subgrid->segN;
        }
      }
      R+=r;
    }
  } else if(subgrid->meth==2) {
    for(i=0;i<subgrid->segN;i++)
    {
      base=ne*(subgrid->segN+1);
      d1=subgrid->dpe[base+i];
      d2=subgrid->dpe[base+i+1];
      if(h>-(d1+d2)/2)
        r=(h+(d1+d2)/2)/subgrid->segN;
      else
        r=0;
      R+=r;
    }
  } else if(subgrid->meth==3){
    for(i=0;i<subgrid->segN;i++)
    {
      base=ne*(subgrid->segN+1);
      x1=subgrid->xpe[base+i];
      x2=subgrid->xpe[base+i+1];
      y1=subgrid->ype[base+i];
      y2=subgrid->ype[base+i+1];       
      r=ReturnFluxHeight(x1,y1,x2,y2,h)/subgrid->segN;
    }
    R+=r;
  }

  return R;
}

/*
 * Function: CalculateFluxHeightProfile
 * Usage: calculate Ac and h profile for each cell 
 * ----------------------------------------------------
 * use CalculateWetArea to get the total wet area for each h 
 *
 */
void CalculateFluxHeightProfile(gridT *grid, int myproc)
{
  int i,j,ne,base;
  REAL k,d1,d2,d;
  for(ne=0;ne<grid->Ne;ne++)
  {
    base=ne*(subgrid->disN+1);
    subgrid->fluxhprof[base]=0;
    for(i=1;i<=subgrid->disN;i++)
    {
      subgrid->fluxhprof[base+i]=CalculateFluxHeight(ne,subgrid->hprofe[base+i]);
    } 
  }
}

/*
 * Function: CalculateWetperimeter 
 * Usage: calculate the wet perimeter for the bottom layer 
 * ----------------------------------------------------
 * calculate the wet perimeter based on depth for each subgrid point 
 *
 */
REAL CalculateWetperimeter(int ne, REAL h)
{
  int i,base;
  REAL r,r1,dis,R,d1,d2,d,x1,x2,y1,y2;
  R=0;
  if(subgrid->meth==1){
    for(i=0;i<subgrid->segN;i++)
    {
      base=ne*(subgrid->segN+1);
      d1=subgrid->dpe[base+i];
      d2=subgrid->dpe[base+i+1];
      dis=pow((subgrid->xpe[base+i]-subgrid->xpe[base+i+1]),2)+pow((subgrid->ype[base+i]-subgrid->ype[base+i+1]),2);

      if(d1>d2)
      {
        d=d1;
        d1=d2;
        d2=d;
      }
      if(d1==d2)
      { 
        if(h<=-d1)
          r=0;
        else
          r=pow(dis,0.5);
      } else {
        if(h<=-d2)
          r=0;
        else if(h>-d2 && h<=-d1)
        { 
          r1=(h+d2)/(d2-d1);
          r=r1*pow((pow((d1-d2),2)+dis),0.5);
        } else {
          r=pow((pow((d1-d2),2)+dis),0.5);
        }
      }
      R+=r;
    }
  } else if(subgrid->meth==2) {
    for(i=0;i<subgrid->segN;i++)
    {
      base=ne*(subgrid->segN+1);
      dis=pow((subgrid->xpe[base+i]-subgrid->xpe[base+i+1]),2)+pow((subgrid->ype[base+i]-subgrid->ype[base+i+1]),2);
      d1=subgrid->dpe[base+i];
      d2=subgrid->dpe[base+i+1];
      if(h>-(d1+d2)/2)
        r=pow((pow((d2-d1),2)+dis),0.5);
      else
        r=0;
      R+=r;
    }
  }

  return R;
}

/*
 * Function: CalculateWetperimeterProfile
 * Usage: calculate Ac and h profile for each cell 
 * ----------------------------------------------------
 * use CalculateWetArea to get the total wet area for each h 
 *
 */
void CalculateWetperimeterProfile(gridT *grid, int myproc)
{
  int i,j,ne,base;
  REAL k,d1,d2,d;
  for(ne=0;ne<grid->Ne;ne++)
  {
    base=ne*(subgrid->disN+1);
    subgrid->Wetperiprof[base]=0;
    for(i=1;i<=subgrid->disN;i++)
    {
      subgrid->Wetperiprof[base+i]=CalculateWetperimeter(ne,subgrid->hprofe[base+i]);
    } 
  }
}

/*
 * Function: CalculateVolumeProfile
 * Usage: calculate Volume profile for each cell 
 * ----------------------------------------------------
 * use Acprofile to get the total volume of water for each h 
 * will be used in UpdateAc to change grid->Ac
 *
 */
void CalculateVolumeProfile(gridT *grid, int myproc)
{
  int i,nc,base;
  REAL dh;
  for(nc=0;nc<grid->Nc;nc++)
  {
    base=nc*(subgrid->disN+1);
    subgrid->Vprof[base]=0;
    dh=subgrid->hprof[base+1]-subgrid->hprof[base];
    for(i=1;i<=subgrid->disN;i++)
    {
      subgrid->Vprof[base+i]=subgrid->Vprof[base+i-1]+dh*(subgrid->Acprof[base+i]+subgrid->Acprof[base+i-1])/2;
    }
  }
}

/*
 * Function: SubgridCheck
 * Usage: check whether nan happens for subgrid method
 * ----------------------------------------------------
 * for Acprof Vprof hprof hprofe fluxhprof
 * 
 */
void SubgridCheck(gridT *grid, int myproc)
{
  int i, k, j, base;
  // for hprof, acprof, vprof
  for(i=0;i<grid->Nc;i++)
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
    {
      if(subgrid->hprof[base+j]!=subgrid->hprof[base+j])
        printf("i %d cell j %dth hprof is %f\n",i,j,subgrid->hprof[base+j]);
      if(subgrid->Acprof[base+j]!=subgrid->Acprof[base+j])
        printf("i %d cell j %dth Acprof is %f\n",i,j,subgrid->Acprof[base+j]);
      if(subgrid->Vprof[base+j]!=subgrid->Vprof[base+j])
        printf("i %d cell j %dth Vprof is %f\n",i,j,subgrid->Vprof[base+j]);
    }
  }
  // for hprofe, fluxhprof
  for(i=0;i<grid->Ne;i++)
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
    {
      if(subgrid->hprofe[base+j]!=subgrid->hprofe[base+j])
        printf("i %d cell j %dth hprofe is %f\n",i,j,subgrid->hprof[base+j]);
      if(subgrid->fluxhprof[base+j]!=subgrid->fluxhprof[base+j])
        printf("i %d cell j %dth fluxhprofe is %f\n",i,j,subgrid->fluxhprof[base+j]);
      if(subgrid->Wetperiprof[base+j]!=subgrid->Wetperiprof[base+j])
        printf("i %d cell j %dth wetperiprofe is %f\n",i,j,subgrid->Wetperiprof[base+j]);
    }
  }
 
  // for Acceff Acveff
  for(i=0;i<grid->Nc;i++)
  {
    for(k=0;k<grid->Nk[i];k++)
    { 
      if(subgrid->Acceff[i][k]!=subgrid->Acceff[i][k])
        printf("i %d cell k %dth layer Acceff is %f\n",i,k,subgrid->Acceff[i][k]);
      if(subgrid->Acveff[i][k]!=subgrid->Acveff[i][k])
        printf("i %d cell k %dth layer top Acveff is %f\n",i,k,subgrid->Acceff[i][k]);      
    }
 
    if(subgrid->Acveff[i][grid->Nk[i]]!=subgrid->Acveff[i][grid->Nk[i]])
        printf("i %d cell bottom layer bottom Acveff is %f\n",i,subgrid->Acveff[i][grid->Nk[i]]);      
  }   
}

/*
 * Function: SubgridBasic
 * Usage: setup everything for subgrid method
 * ----------------------------------------------------
 * used in Phys.c at first time step to setup subgrid method
 *
 */
void SubgridBasic(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  // allocate subgrid struture first
  subgrid=(subgridT *)SunMalloc(sizeof(subgridT),"SubgridBasic");

  // read subgrid properties
  ReadSubgridProperties(myproc);

  // allocate space for subgrid method
  AllocateandInitializeSubgrid(grid, myproc);

  // calculate all the x y for subcell and subedege
  CalculateSubgridXY(grid, myproc);

  // calculate the index for each subcell
  CalculateSubgridCellp(grid, myproc);

  // calculate all the for all the points in subgrid
  InterpolateSubgridDepth(grid, phys,myproc);

  // calculate the max and min depth to get a reasonable range for h
  CalculateHProfile(grid, myproc);

  // calculate the wet area profile for future use
  CalculateAcProfile(grid, myproc);

  // calculate the cell volume profile for future use
  CalculateVolumeProfile(grid, myproc);

  // calculate the edge flux height profile for future use
  CalculateFluxHeightProfile(grid, myproc);

  // calculate wet perimeter profile for the last layer
  CalculateWetperimeterProfile(grid,myproc);

  // check subgrid
  SubgridCheck(grid,myproc);
  
  // output subgrid for test
  OutputSubgrid(grid,phys,prop,myproc,numprocs,comm);

  printf("finish subgrid basic module \n");
}

/*
 * Function: UpdateSubgridVeff
 * Usage: update subgrid->Veff to represent the change of water volume
 * ----------------------------------------------------
 * use Vprof to get water volume;
 *
 */
void UpdateSubgridVeff(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int nc, i,base;
  REAL dh,h,V,dV;
  for(nc=0;nc<grid->Nc;nc++)
  {
    h=phys->h[nc];
    if(prop->culvertmodel && h>Culverttop[nc])
      h=Culverttop[nc];
    subgrid->Veff[nc]=UpdateVeff(nc,h);
  }
}

/*
 * Function: UpdateVeff
 * Usage: basic function to calculate Veff to represent the change of water volume
 * ----------------------------------------------------
 * use Vprof to get water volume;
 *
 */
REAL UpdateVeff(int nc, REAL h)
{
  int  i,base;
  REAL dh,V,dV;

  base=nc*(subgrid->disN+1); 
  if(h<subgrid->hprof[base])
    V=0;
  else if(h>=subgrid->hprof[base+subgrid->disN])
  {
    V=subgrid->Vprof[base+subgrid->disN]+(h-subgrid->hprof[base+subgrid->disN])*subgrid->Acbackup[nc];
  }
  else {
    i=1;
    while(i<=subgrid->disN)
    {
      if(h<subgrid->hprof[base+i])
        break;
      i++;
    }
    dh=subgrid->hprof[base+i]-subgrid->hprof[base+i-1];  
    dV=subgrid->Vprof[base+i]-subgrid->Vprof[base+i-1]; 
    V=subgrid->Vprof[base+i-1]+(h-subgrid->hprof[base+i-1])/dh*dV;
  }
  return V;
}


/*
 * Function: UpdateSubgridAceff
 * Usage: update subgrid->Aceff at free surface to represent the change of wet area
 * ----------------------------------------------------
 * use Acprof to calculate effective wet area for the iteration of h
 * to solve mass conservation problem for free surface
 *
 */
void UpdateSubgridAceff(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int nc, i,base;
  REAL dh,h,Ac,dac;
  for(nc=0;nc<grid->Nc;nc++)
  {
    h=phys->h[nc];
    subgrid->Aceffold[nc]=subgrid->Aceff[nc];
    subgrid->Aceff[nc]=UpdateAceff(nc,h);
    if(prop->culvertmodel && h>=Culverttop[nc])
      subgrid->Aceff[nc]=0;
    // may be not correct
    //if(h<subgrid->hprof[nc*(subgrid->disN+1)+1])
      //if(subgrid->Acprof[nc*(subgrid->disN+1)]!=0)
        //subgrid->Aceff[nc]=subgrid->Acprof[nc*(subgrid->disN+1)];
      //else
        //subgrid->Aceff[nc]=subgrid->Acprof[nc*(subgrid->disN+1)+1];
  }
}


/*
 * Function: UpdateAceff
 * Usage:  basic function to calculate wet area based on h
 * ----------------------------------------------------
 * use Acprof to calculate effective wet area for the iteration of h
 * to solve mass conservation problem
 *
 */
REAL UpdateAceff(int nc, REAL h)
{
  int i, base;
  REAL dh,Ac,dac;
   
  base=nc*(subgrid->disN+1);
    
  if(h<subgrid->hprof[base])
    Ac=0;
  else if(h>=subgrid->hprof[base+subgrid->disN])
    Ac=subgrid->Acbackup[nc];
  else {
    i=1;
    while(i<=subgrid->disN)
    {
      if(h<subgrid->hprof[base+i])
        break;
      i++;
    }
    dh=subgrid->hprof[base+i]-subgrid->hprof[base+i-1];  
    dac=subgrid->Acprof[base+i]-subgrid->Acprof[base+i-1]; 
    Ac=subgrid->Acprof[base+i-1]+(h-subgrid->hprof[base+i-1])/dh*dac;
  } 
  
  return Ac;
}

/*
 * Function: UpdateSubgridVerticalAceff
 * Usage: update vertical Aceff for each cell and layer to represent the change of wet area
 * ----------------------------------------------------
 * use to solve scalar transport problem
 *
 */
void UpdateSubgridVerticalAceff(gridT *grid, physT *phys, propT *prop, int option, int myproc)
{
  // if option == 0 means this function is for scalar preparation, the only difference is when ctop[i]<=ctopold[i], we only calculate Acceff and Acveff at ctopold[i]
  // if option == 1 mean this function is used after updatescalars, then to calculate the exact Acceff and Acveff for all layer from ctop[i] to grid->Nk[i]-1
  REAL v,ac,ztop,zbot,dznew;
  int i,k,ktop;
   
  // prepare for Updatescalars
  if(option==0){
    for(i=0;i<grid->Nc;i++)
    {
      // store the old values
      for(k=0;k<grid->Nk[i];k++) 
        subgrid->Acceffold[i][k]=subgrid->Acceff[i][k];
 
      for(k=0;k<=grid->Nk[i];k++)
        subgrid->Acveffold[i][k]=subgrid->Acveff[i][k];


      if(grid->ctop[i]>=grid->ctopold[i]) {
        ktop=grid->ctop[i];
        dznew=grid->dzz[i][ktop];
      } else {
        ktop=grid->ctopold[i];
        dznew=0;
        for(k=grid->ctop[i];k<=grid->ctopold[i];k++) 
          dznew+=grid->dzz[i][k];      
      }

      // beyond ctop the ac is zero
      for(k=0;k<ktop;k++)
      {
        subgrid->Acceff[i][k]=0;
        subgrid->Acveff[i][k]=0;
      }
    
      // from bottom to top
      zbot=-grid->dv[i];
      ztop=zbot+grid->dzz[i][grid->Nk[i]-1];

      // update bottom wet area for the bottom vertical layer 
      subgrid->Acveff[i][grid->Nk[i]]=UpdateAceff(i,zbot);

      for(k=grid->Nk[i]-1;k>ktop;k--)
      {  
        // update mean cell center wet area based on volume
        if(grid->dzz[i][k]>DRYCELLHEIGHT)
          subgrid->Acceff[i][k]=(UpdateVeff(i,ztop)-UpdateVeff(i,zbot))/grid->dzz[i][k];
        else
          subgrid->Acceff[i][k]=0;

        subgrid->Acveff[i][k]=UpdateAceff(i,ztop);

        if(k!=ktop+1){
          zbot=ztop;
          ztop=zbot+grid->dzz[i][k-1];
        }
      } 
      if(ktop!=(grid->Nk[i]-1))
      {
        zbot=ztop;
        ztop=zbot+dznew;
      } else {
        zbot=-grid->dv[i];
        ztop=zbot+dznew;
      }

      if(dznew>DRYCELLHEIGHT)
        subgrid->Acceff[i][ktop]=(UpdateVeff(i,ztop)-UpdateVeff(i,zbot))/dznew;
      else
        subgrid->Acceff[i][ktop]=0;

      subgrid->Acveff[i][ktop]=UpdateAceff(i,ztop);
    }
  } else {
    for(i=0;i<grid->Nc;i++)
    {
      if(grid->ctop[i]<grid->ctopold[i])
      {   
        // from bottom to top
        zbot=-grid->dv[i];
        ztop=zbot+grid->dzz[i][grid->Nk[i]-1];

        for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
        {  
          // update mean cell center wet area based on volume
          if(grid->dzz[i][k]>DRYCELLHEIGHT)
            subgrid->Acceff[i][k]=(UpdateVeff(i,ztop)-UpdateVeff(i,zbot))/grid->dzz[i][k];
          else
            subgrid->Acceff[i][k]=0;

          subgrid->Acveff[i][k]=UpdateAceff(i,ztop);

          if(k!=grid->ctop[i])
          {
            zbot=ztop;
            ztop=zbot+grid->dzz[i][k-1];
          }
        }
      }
    } 
  }
}


/*
 * Function: UpdateSubgridFluxHeight
 * Usage: update grid->dzf to get new flux height 
 * ----------------------------------------------------
 * use fluxhprof to calculate new effective grid->dzf
 *
 */
void UpdateSubgridFluxHeight(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int ne, nc,nc1,nc2, i,base,k,j;
  REAL dh,h,fh,fhr,dfhr,dfh,culverttop,u,hbot,hc,wetperi;
  for(ne=0;ne<grid->Ne;ne++)
  {
    nc1 = grid->grad[2*ne];
    nc2 = grid->grad[2*ne+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;
    
    u=phys->u[ne][grid->etop[ne]];

    //hc is for upwind free surface height for edge ne
    hc=UpWind(u,phys->h[nc1],phys->h[nc2]);

    if(hc==phys->h[nc1])
      nc=nc1;
    else
      nc=nc2;

    // h is the free surface height used to calculate flux height can be upwind or central differencing
    h=UpWind(u,phys->h[nc1],phys->h[nc2]);

    if(prop->culvertmodel)
    {
      culverttop=UpWind(phys->u[ne][0],Culverttop[nc1],Culverttop[nc2]);
      if(h>culverttop)
        h=culverttop;
    }

    for(k=0;k<grid->etop[ne];k++)
      grid->dzf[ne][k]=0;
    
    hbot=hc;    
    for(k=grid->ctop[nc];k<grid->etop[ne];k++)
      hbot-=grid->dzz[nc][k];

    for(k=grid->etop[ne];k<grid->Nke[ne]-1;k++)
    { 
      hbot-=grid->dzf[ne][k]; // still old dzf from the original dzf
      grid->dzf[ne][k]=UpdateFluxHeight(ne,h)-UpdateFluxHeight(ne,hbot);
      h=hbot;
    }
    grid->dzf[ne][grid->Nke[ne]-1]=UpdateFluxHeight(ne,h);

    if(((grid->dv[nc]+phys->h[nc])<=1.001*DRYCELLHEIGHT))
    {
      grid->dzf[ne][grid->Nke[ne]-1]=0;
    }

    wetperi=UpdateWetperi(ne,h);
    if(grid->dzf[ne][k]==0)
      wetperi=grid->df[ne];

    // added part
    k=grid->Nke[ne]-1;
    
    if(grid->etop[ne]==k && grid->mark[ne]==2 && grid->dzf[ne][k]<=0.01)
    {  
       grid->dzf[ne][k]=0.01;
       wetperi=grid->df[ne];
    }
    
    for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
      if(grid->dzf[ne][k]<=DRYCELLHEIGHT)
        grid->dzf[ne][k]=0;

    subgrid->dzboteff[ne]=grid->dzf[ne][grid->Nke[ne]-1]*grid->df[ne]/wetperi;
  }   
}

/*
 * Function: UpdateFluxHeight
 * Usage: calculate average flux height based on fluxhprof
 * ----------------------------------------------------
 * used in UpdateSubgridFluxheight
 *
 */
REAL UpdateFluxHeight(int ne, REAL h)
{
  int base,i;
  REAL fh,dh,dfh;
  base=ne*(subgrid->disN+1); 
   
  if(h<subgrid->hprofe[base])
  {
    fh=0;
  }
  else if(h>=subgrid->hprofe[base+subgrid->disN])
  {
    dh=h-subgrid->hprofe[base+subgrid->disN];
    fh=subgrid->fluxhprof[base+subgrid->disN]+dh;
  } else {
    i=1;
    while(i<=subgrid->disN)
    {
      if(h<subgrid->hprofe[base+i])
        break;
      i++;
    }
    dh=subgrid->hprofe[base+i]-subgrid->hprofe[base+i-1];  
    dfh=subgrid->fluxhprof[base+i]-subgrid->fluxhprof[base+i-1];
    fh=subgrid->fluxhprof[base+i-1]+(h-subgrid->hprofe[base+i-1])/dh*dfh;
  }
  return fh;
}

/*
 * Function: UpdateWetperi
 * Usage: calculate wet perimeter based on wetperiprof
 * ----------------------------------------------------
 * used in UpdateSubgridFluxheight
 *
 */
REAL UpdateWetperi(int ne, REAL h)
{
  int base,i;
  REAL fh,dh,dfh;
  base=ne*(subgrid->disN+1); 
   
  if(h<subgrid->hprofe[base])
  {
    fh=0;
  }
  else if(h>=subgrid->hprofe[base+subgrid->disN])
  {
    fh=subgrid->Wetperiprof[base+subgrid->disN];
  } else {
    i=1;
    while(i<=subgrid->disN)
    {
      if(h<subgrid->hprofe[base+i])
        break;
      i++;
    }
    dh=subgrid->hprofe[base+i]-subgrid->hprofe[base+i-1];  
    dfh=subgrid->Wetperiprof[base+i]-subgrid->Wetperiprof[base+i-1];
    fh=subgrid->Wetperiprof[base+i-1]+(h-subgrid->hprofe[base+i-1])/dh*dfh;
  }
  return fh;
}

/*
 * Function: OutputSubgrid
 * Usage: output hprof acprof Vprof hprofe fluxhprof 
 * ----------------------------------------------------
 * use fluxhprof to calculate new effective grid->dzf
 *
 */
void OutputSubgrid(gridT *grid, physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm)
{
  int i,j,base,N;
  char str1[BUFFERLENGTH];
  FILE *ofile; 

  // output hprof for each cell
  MPI_GetFile(str1,DATAFILE,"HprofFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
      fprintf(ofile,"%e ",subgrid->hprof[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output hprofe for each edge
  MPI_GetFile(str1,DATAFILE,"HprofeFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Ne;i++) 
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
      fprintf(ofile,"%e ",subgrid->hprofe[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output Acprof for each cell
  MPI_GetFile(str1,DATAFILE,"AcprofFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
      fprintf(ofile,"%e ",subgrid->Acprof[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);
  
  // output Vprof for each cell
  MPI_GetFile(str1,DATAFILE,"VprofFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
      fprintf(ofile,"%e ",subgrid->Vprof[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output fluxhprof for each edge
  MPI_GetFile(str1,DATAFILE,"FluxhprofFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Ne;i++) 
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
      fprintf(ofile,"%e ",subgrid->fluxhprof[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output Wetperiprof for each edge
  MPI_GetFile(str1,DATAFILE,"WetperiprofFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Ne;i++) 
  {
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++)
      fprintf(ofile,"%e ",subgrid->Wetperiprof[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);


  // output cellp for each cell
  MPI_GetFile(str1,DATAFILE,"cellpFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*(pow(subgrid->segN,2)*3)*(grid->maxfaces-2);
    for(j=0;j<((grid->nfaces[i]-2)*(pow(subgrid->segN,2)*3));j++)
      fprintf(ofile,"%d ",subgrid->cellp[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output subgrid->hmin for each cell
  MPI_GetFile(str1,DATAFILE,"subgridDFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  
  for(i=0;i<grid->Nc;i++) 
    fprintf(ofile,"%e %e %e\n",grid->xv[i],grid->yv[i],(-subgrid->hmin[i]));
  fclose(ofile);


  // output xp yp dp for each cell
  MPI_GetFile(str1,DATAFILE,"subpointFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->maxfaces-2);
    for(j=0;j<(((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->nfaces[i]-2));j++){
      fprintf(ofile,"%d %e %e %e\n",i,subgrid->xp[base+j],subgrid->yp[base+j],subgrid->dp[base+j]);
    }
    //fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output xpe ype dpe for each edge
  MPI_GetFile(str1,DATAFILE,"subpointeFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Ne;i++)
  { 
    base=i*(subgrid->segN+1);
    for(j=0;j<=subgrid->segN;j++)
      fprintf(ofile,"%d %e %e %e\n",i,subgrid->xpe[base+j],subgrid->ype[base+j],subgrid->dpe[base+j]);    
    //fprintf(ofile,"\n");  
  }
  fclose(ofile);
 
  // output subgrid->hmin
  MPI_GetFile(str1,DATAFILE,"subdmaxFile","OutputSubgrid",myproc);
  ofile = MPI_FOpen(str1,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Nc;i++)
  { 
      fprintf(ofile,"%e %e %e\n",grid->xv[i],grid->yv[i],(-subgrid->hmin[i]));    
    //fprintf(ofile,"\n");  
  }
  fclose(ofile);  
}

/*
 * Function: UpdateSubgridDz
 * Usage: new UpdateDz for subgrid in the consideration when
 * subgrid->hmin<-grid->dv
 * ----------------------------------------------------
 * use to arrange grid->dzz
 * cannot change grid->Nk due to phys variable based on grid Nk
 */

void UpdateSubgridDZ(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int i, j, k, ne1, ne2, Nc=grid->Nc, Ne=grid->Ne, flag, nc1, nc2;
  REAL z, dzz1, dzz2;

  // don't need to recompute for linearized FS
  if(prop->linearFS) {
    return;
  }

  // If this is not an initial call then set dzzold to store the old value of dzz
  // and also set the etopold and ctopold pointers to store the top indices of
  // the grid.
  for(j=0;j<Ne;j++)
    grid->etopold[j]=grid->etop[j];
  for(i=0;i<Nc;i++) 
  {
    grid->ctopold[i]=grid->ctop[i];
    for(k=0;k<grid->ctop[i];k++)
      grid->dzzold[i][k]=0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
  }

  // First set the thickness of the bottom grid layer.  If this is a partial-step
  // grid then the dzz will vary over the horizontal at the bottom layer.  Otherwise,
  // the dzz at the bottom will be equal to dz at the bottom.
  for(i=0;i<Nc;i++) {
    z = 0;
    for(k=0;k<grid->Nk[i];k++)
      z-=grid->dz[k];
    grid->dzz[i][grid->Nk[i]-1]=grid->dz[grid->Nk[i]-1]-subgrid->hmin[i]+z;
  }

  // Loop through and set the vertical grid thickness when the free surface cuts through 
  // a particular cell.
  if(grid->Nkmax>1) {
    for(i=0;i<Nc;i++) {
      z = 0;
      flag = 0;

      for(k=0;k<grid->Nk[i];k++) {
        z-=grid->dz[k];
        if(phys->h[i]>=z) 
          if(!flag) {
            if(k==grid->Nk[i]-1) {
              grid->dzz[i][k]=phys->h[i]-subgrid->hmin[i];
              grid->ctop[i]=k;
            } else if(phys->h[i]==z) {
              grid->dzz[i][k]=0;
              grid->ctop[i]=k+1;
            } else {
              grid->dzz[i][k]=phys->h[i]-z;
              grid->ctop[i]=k;
            }
            flag=1;
          } else {
            if(k==grid->Nk[i]-1) 
              grid->dzz[i][k]=grid->dz[k]-subgrid->hmin[i]+z;
            else 
              if(z<subgrid->hmin[i])
                grid->dzz[i][k]=0;
              else 
                grid->dzz[i][k]=grid->dz[k];
          } 
        else 
          grid->dzz[i][k]=0;
      }
    }
  } else 
    for(i=0;i<Nc;i++) 
      grid->dzz[i][0]=-subgrid->hmin[i]+phys->h[i];


  // Now set grid->etop and ctop which store the index of the top cell  
  for(j=0;j<grid->Ne;j++) {
    ne1 = grid->grad[2*j];
    ne2 = grid->grad[2*j+1];
    if(ne1 == -1)
      grid->etop[j]=grid->ctop[ne2];
    else if(ne2 == -1)
      grid->etop[j]=grid->ctop[ne1];
    else if(grid->ctop[ne1]<grid->ctop[ne2])
      grid->etop[j]=grid->ctop[ne1];
    else
      grid->etop[j]=grid->ctop[ne2];
  }

}





