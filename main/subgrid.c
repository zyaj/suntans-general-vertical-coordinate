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
#include "boundaries.h"
#include "sediments.h"
#include "subgrid.h"
#include "util.h"
#include "tvd.h"
#include "mympi.h"
#include "subgrid.h"
#include "scalars.h"
#include "physio.h"
#include "wave.h"
#include "marsh.h"
#include "culvert.h"
#include "vertcoordinate.h"

void ReadSubgridProperties(propT *prop, int myproc);
void AllocateandInitializeSubgrid(gridT *grid, propT *prop, int myproc);
void CalculateSubgridXY(gridT *grid, int myproc);
void CalculateSubgridCellp(gridT *grid, int myproc);
void InterpolateSubgridDepth(gridT *grid,physT *phys, propT *prop,int myproc,int numprocs);
void InterpolateSubgridHmarsh(gridT *grid,physT *phys, int myproc, int numprocs);
void InterpolateSubgridCdV(gridT *grid,physT *phys, int myproc, int numprocs);
void CalculateVolumeProfile(gridT *grid, int myproc);
void CalculateHProfile(gridT *grid, int myproc);
void CalculateAcProfile(gridT *grid, int myproc);
void CalculateAreaRatio(gridT *grid, int myproc);
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
REAL UpdateFreeSurface(int nc, REAL V);
void OutputSubgrid(gridT *grid, physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm);

/*
 * Function: ReadSubgridProperties
 * Usage: Read Subgrid Properties
 * ----------------------------------------------------
 * Based on suntans.dat, load in the important parameters for
 * subgrid model
 */
void ReadSubgridProperties(propT *prop, int myproc) {
  subgrid->segN = MPI_GetValue(DATAFILE,"segN","ReadSubgridProperties",myproc); 
  subgrid->disN = MPI_GetValue(DATAFILE,"disN","ReadSubgridProperties",myproc); 
  subgrid->meth = MPI_GetValue(DATAFILE,"subgridmeth","ReadSubgridProperties",myproc); 
  subgrid->dzfmeth = MPI_GetValue(DATAFILE,"subgriddzfmeth","ReadSubgridProperties",myproc); 
  subgrid->dpint =MPI_GetValue(DATAFILE,"subgriddpint","ReadSubgridProperties",myproc); 
  subgrid->dragpara =MPI_GetValue(DATAFILE,"subgriddragpara","ReadSubgridProperties",myproc); 
  subgrid->eps =MPI_GetValue(DATAFILE,"subgrideps","ReadSubgridProperties",myproc); 
  
  if(prop->marshmodel)
  {
    subgrid->hmarshint =MPI_GetValue(DATAFILE,"subgridhmarshint","ReadSubgridProperties",myproc); 
    subgrid->cdvint =MPI_GetValue(DATAFILE,"subgridcdvint","ReadSubgridProperties",myproc); 

  }
  if(prop->computeSediments)
    subgrid->erosionpara= MPI_GetValue(DATAFILE,"subgriderosionpara","ReadSubgridProperties",myproc); 
}

/*
 * Function: AllocateandInitializeSubgrid
 * Usage: Allocate Subgrid variable
 * ----------------------------------------------------
 * allocate space and initialize(set all zero) for all subgrid variables
 *
 */
void AllocateandInitializeSubgrid(gridT *grid, propT *prop, int myproc) {
  int i,N,Ntotal,k;
  REAL tmp;
  N=grid->Nc*(grid->maxfaces-2);
  // allocate xp yp dp hmarshp for subgrid scal
  Ntotal=N*(subgrid->segN+1)*(subgrid->segN+2)/2;
  subgrid->xp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->yp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->dp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  if(prop->marshmodel)
  {
    subgrid->hmarshp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
    subgrid->cdvp = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  }
  for(i=0;i<Ntotal;i++)
  {
    subgrid->xp[i]=0.0;
    subgrid->yp[i]=0.0;
    subgrid->dp[i]=0.0;
    if(prop->marshmodel)
    {
      subgrid->hmarshp[i]=0.0;
      subgrid->cdvp[i]=0.0;

    }
  } 

  // allocate dc hmarshc for subgrid cell
  Ntotal=N*(pow(subgrid->segN,2));
  subgrid->dc = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->d1 = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->d2 = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->d3 = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  // only work with erosionpara
  subgrid->H_sub = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->A_sub = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");

  if(subgrid->erosionpara)
  {
    subgrid->Cdratio = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
    subgrid->taubsub = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  }

  if(prop->marshmodel)
  {
      subgrid->hmarshc = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
      subgrid->cdvc = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  }

  for(i=0;i<Ntotal;i++)
  {
    subgrid->dc[i]=0.0;
    subgrid->d1[i]=0.0;
    subgrid->d2[i]=0.0;
    subgrid->d3[i]=0.0;
    subgrid->A_sub[i]=0.0;
    subgrid->H_sub[i]=0.0;
    if(subgrid->erosionpara)
      subgrid->Cdratio[i]=0.0;
    if(prop->marshmodel)
    {
      subgrid->hmarshc[i]=0.0;
      subgrid->cdvc[i]=0.0;
    }
    if(subgrid->erosionpara)
      subgrid->taubsub[i]=0.0;
  }

  // allocate xpe ype dpe hmarshpe for each edge
  Ntotal=grid->Ne*(subgrid->segN+1);
  subgrid->xpe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->ype = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  subgrid->dpe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  
  if(prop->marshmodel)
  {
    subgrid->hmarshpe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
    subgrid->cdvpe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");

  }
  
  for(i=0;i<Ntotal;i++)
  {
    subgrid->xpe[i]=0.0;
    subgrid->ype[i]=0.0;
    subgrid->dpe[i]=0.0;
    if(prop->marshmodel)
    {
      subgrid->hmarshpe[i]=0.0;
      subgrid->cdvpe[i]=0.0;
    }
  }  

  // depth and marsh height at subedge center 
  Ntotal=grid->Ne*(subgrid->segN);
  subgrid->de = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  if(prop->marshmodel)
  {
    subgrid->hmarshe = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
    subgrid->cdve = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  }
  for(i=0;i<Ntotal;i++)
  {
    subgrid->de[i]=0.0;
    if(prop->marshmodel){
      subgrid->hmarshe[i]=0.0;
      subgrid->cdve[i]=0.0;
    }
  }    

  // allocate subcellp
  Ntotal=N*(pow(subgrid->segN,2))*3;
  subgrid->cellp = (int *)SunMalloc(Ntotal*sizeof(int), "AllocateSubgrid");
  for(i=0;i<Ntotal;i++)
    subgrid->cellp[i]=0;

  // allocate he
  Ntotal=grid->Ne;
  subgrid->he = (REAL *)SunMalloc(Ntotal*sizeof(REAL),"AllocateSubgrid");
  for(i=0;i<Ntotal;i++)
    subgrid->he[i]=0;

  // allocate profile for effective area and volume for each cell
  Ntotal=grid->Nc*(subgrid->disN+1);
  subgrid->hprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->Vprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->Acprof = (REAL *)SunMalloc(Ntotal*sizeof(REAL), "AllocateSubgrid");
  subgrid->hmin = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSubgrid");
  subgrid->hmax = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSubgrid");
  subgrid->hiter = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSubgrid");
  subgrid->hiter_min = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSubgrid");

  for(i=0;i<Ntotal;i++)
  {
    subgrid->hprof[i]=0.0;
    subgrid->Vprof[i]=0.0;
    subgrid->Acprof[i]=0.0;
  }
 
  for(i=0;i<grid->Nc;i++)
  {
    subgrid->hmin[i]=0.0;
    subgrid->hmax[i]=0.0;
    subgrid->hiter[i]=0.0;
    subgrid->hiter_min[i]=0.0;
  }

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
  subgrid->varNc = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->rhs = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->residual = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Verr = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Aceff = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Acratio = (REAL *)SunMalloc((grid->maxfaces-2)*grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Aceffold = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Heff = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Heffold = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Veff = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Veffold = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  subgrid->Acceff = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acceffold = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acveffold = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acveffold2 = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->Acveff = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");

  // sediment transport for 1d-s and 2d-s subgrid model 
  // for deposition
  if(prop->computeSediments){
    subgrid->Asedi = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
    subgrid->Vsedi = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
    tmp=MPI_GetValue(DATAFILE,"subgriddelta","AllocateSubgrid",myproc); 
    subgrid->delta = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
    subgrid->Cdmean = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateSubgrid");
  }


  for(i=0;i<grid->Nc;i++)
  {
    if(prop->computeSediments){
      subgrid->Cdmean[i]=0.0;
      subgrid->Asedi[i]=0.0;
      subgrid->Vsedi[i]=0.0;
      subgrid->delta[i]=tmp;
    }
    subgrid->Acbackup[i]=grid->Ac[i];
    subgrid->Acwet[i]=grid->Ac[i];
    subgrid->Aceff[i]=grid->Ac[i];
    subgrid->Aceffold[i]=grid->Ac[i];
    subgrid->varNc[i]=0.0;
    subgrid->rhs[i]=0.0;
    subgrid->residual[i]=0.0;
    subgrid->Verr[i]=0.0;
    subgrid->Veff[i]=0.0;
    subgrid->Veffold[i]=0.0;
    subgrid->Heff[i]=0.0;
    subgrid->Heffold[i]=0.0;
    subgrid->Acceff[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acceffold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acveffold[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acveffold2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateSubgrid");
    subgrid->Acveff[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateSubgrid");
    for(k=0;k<grid->Nk[i];k++)
    {
      subgrid->Acceff[i][k]=0;
      subgrid->Acceffold[i][k]=0;
      subgrid->Acveff[i][k]=0;
      subgrid->Acveffold[i][k]=0; 
      subgrid->Acveffold2[i][k]=0;   
    }
    subgrid->Acveff[i][grid->Nk[i]]=0;
    subgrid->Acveffold[i][grid->Nk[i]]=0;
    subgrid->Acveffold2[i][grid->Nk[i]]=0;
  }

  // allocate flux variable to check flow condition
  subgrid->fluxp = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  subgrid->fluxn = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateSubgrid");
  
  for(i=0;i<grid->Nc;i++){
    subgrid->fluxp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateSubgrid");
    subgrid->fluxn[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateSubgrid");
  }

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
    {
      subgrid->fluxn[i][k]=0;
      subgrid->fluxp[i][k]=0;
    }
}

/* 
 * Function: OpenSubgridFiles
 * Usage: OpenFiles(myproc);
 * ------------------------------
 * Open all of the files used for i/o to store the file pointers.
 * include Veff and Aceff for each cell at different time steps
 * 
 */
void OpenSubgridFiles(int computeSediments, int mergeArrays, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  MPI_GetFile(filename,DATAFILE,"AeffFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  subgrid->AceffFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"VeffFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  subgrid->VeffFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  if(computeSediments){
    MPI_GetFile(filename,DATAFILE,"VsediFile","OpenFiles",myproc);
    if(mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    subgrid->VsediFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  
    MPI_GetFile(filename,DATAFILE,"AsediFile","OpenFiles",myproc);
    if(mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    subgrid->AsediFID = MPI_FOpen(str,"w","OpenFiles",myproc);

    MPI_GetFile(filename,DATAFILE,"subErosionFile","OpenFiles",myproc);
    if(mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    subgrid->subErosionFID = MPI_FOpen(str,"w","OpenFiles",myproc);

    MPI_GetFile(filename,DATAFILE,"subDepositionFile","OpenFiles",myproc);
    if(mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    subgrid->subDepositionFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  }
}

/*
 * Function: OutputSubgridVariables
 * Usage: OutputSubgridVariables(grid,phys,prop,myproc,numprocs,blowup,comm);
 * ---------------------------------------------------------------------------
 * Output the data every ntout steps as specified in suntans.dat
 * now include veff and aceff
 *
 */
 void OutputSubgridVariables(gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm)
 {
   int i, j, jptr, k, nwritten, arraySize, writeProc,nc1,nc2;
   char str[BUFFERLENGTH], filename[BUFFERLENGTH];

   if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart) {
     Write2DData(subgrid->Aceff,prop->mergeArrays,subgrid->AceffFID,"Error outputting subgrid Aceff data!\n",
       grid,numprocs,myproc,comm);
     Write2DData(subgrid->Veff,prop->mergeArrays,subgrid->VeffFID,"Error outputting subgrid Veff data!\n",
       grid,numprocs,myproc,comm);      
     
     if(prop->computeSediments)
     {
       Write2DData(subgrid->Vsedi,prop->mergeArrays,subgrid->VsediFID,"Error outputting subgrid Vsedi data!\n",
         grid,numprocs,myproc,comm);
       Write2DData(subgrid->Asedi,prop->mergeArrays,subgrid->AsediFID,"Error outputting subgrid Asedi data!\n",
         grid,numprocs,myproc,comm);
       for(i=0;i<grid->Nc;i++)
         subgrid->varNc[i]=sediments->Erosion[0][i][0];
       Write2DData(subgrid->varNc,prop->mergeArrays,subgrid->subErosionFID,"Error outputting subgrid erosion data!\n",
         grid,numprocs,myproc,comm);  
       for(i=0;i<grid->Nc;i++)
         subgrid->varNc[i]=sediments->Deposition[0][i];
       Write2DData(subgrid->varNc,prop->mergeArrays,subgrid->subDepositionFID,"Error outputting subgrid deposition data!\n",
         grid,numprocs,myproc,comm);
     }        

     if(prop->n==prop->nsteps+prop->nstart && myproc==0) {
       fclose(subgrid->AceffFID);
       fclose(subgrid->VeffFID);
       if(prop->computeSediments){
         fclose(subgrid->AsediFID);
         fclose(subgrid->VsediFID);
         fclose(subgrid->subErosionFID);
         fclose(subgrid->subDepositionFID);
       }
     }
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
      base=grid->cells[i*grid->maxfaces];
      dx1=(grid->xp[grid->cells[i*grid->maxfaces+1+j]]-grid->xp[base])/subgrid->segN;
      dx2=(grid->xp[grid->cells[i*grid->maxfaces+2+j]]-grid->xp[base])/subgrid->segN;
      dy1=(grid->yp[grid->cells[i*grid->maxfaces+1+j]]-grid->yp[base])/subgrid->segN;
      dy2=(grid->yp[grid->cells[i*grid->maxfaces+2+j]]-grid->yp[base])/subgrid->segN;        
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
 * Function:  CalculateCellSubgridXY
 * Usage: calculate all the points XY for each subcell
 * ----------------------------------------------------
 * use for grid.c to interpolate depth should input 
 * cell index n, segment number segN, x, y
 */
void CalculateCellSubgridXY(REAL *x, REAL *y, int n, int segN, gridT *grid, int myproc)
{
  int i,j,k,l,m,base;
  REAL dx1,dx2,dy1,dy2;

  // calculate x y for sub cell
  i=n;
  m=0;
  for(j=0;j<(grid->nfaces[i]-2);j++)
  {
    if(j==0)
    {
      base=grid->cells[i*grid->maxfaces];
      dx1=(grid->xp[grid->cells[i*grid->maxfaces+1]]-grid->xp[base])/segN;
      dx2=(grid->xp[grid->cells[i*grid->maxfaces+2]]-grid->xp[base])/segN;
      dy1=(grid->yp[grid->cells[i*grid->maxfaces+1]]-grid->yp[base])/segN;
      dy2=(grid->yp[grid->cells[i*grid->maxfaces+2]]-grid->yp[base])/segN;        
    } else {
      base=grid->cells[i*grid->maxfaces+2];
      dx1=(grid->xp[grid->cells[i*grid->maxfaces+3]]-grid->xp[base])/segN;
      dx2=(grid->xp[grid->cells[i*grid->maxfaces]]-grid->xp[base])/segN;
      dy1=(grid->yp[grid->cells[i*grid->maxfaces+3]]-grid->yp[base])/segN;
      dy2=(grid->yp[grid->cells[i*grid->maxfaces]]-grid->yp[base])/segN; 
    }
    for(k=0;k<=segN;k++)
    {
      for(l=k;l<=segN;l++)
      {
        x[m]=grid->xp[base]+dx1*(l-k)+dx2*(k);
        y[m]=grid->yp[base]+dy1*(l-k)+dy2*(k);
        m++;
      }  
    }    
  }

  // check m correct
  if(m!=((grid->nfaces[i]-2)*(segN+1)*(segN+2)/2)){
    printf("Error computing subgrid->xp/yp in function CalculateCellSubgridXY\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
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
 * Function: InterpolateSubgridHmarsh
 * Usage: calculate vegetation height for all subpoints
 * ----------------------------------------------------
 * input the same file as in grid.c to interpolate 
 * depth for all sub points
 *
 */
void InterpolateSubgridHmarsh(gridT *grid, physT *phys, int myproc, int numprocs)
{
  int i, j, n, nc,nc1,nc2,ne, Nd, scaledepth,ncount,base,k,intdd,ntotal,base1;
  REAL *xd, *yd, *d,intd,min;
  char str[BUFFERLENGTH],str1[BUFFERLENGTH];
  FILE *ifile;

  if(subgrid->hmarshint==1)
  {   
    MPI_GetFile(str,DATAFILE,"InputhmarshFile","InterpolateSubgridHmarsh",myproc);
    Nd = MPI_GetSize(str,"InterpolateSubgridHmarsh",myproc);
    xd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridHmarsh");
    yd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridHmarsh");
    d = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridHmarsh");

    ifile = MPI_FOpen(str,"r","InterpolateSubgridHmarsh",myproc);
    for(n=0;n<Nd;n++) {
      xd[n]=getfield(ifile,str);
      yd[n]=getfield(ifile,str);
      d[n]=getfield(ifile,str);//fabs(getfield(ifile,str));   
    }
    fclose(ifile);
    // interpolate subcell depth
    ncount=grid->Nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
    Interp(xd,yd,d,Nd,&(subgrid->xp[0]), &(subgrid->yp[0]),&(subgrid->hmarshp[0]),ncount,grid->maxfaces);

    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        if(subgrid->hmarshp[base+n]>=0.95*1)
          subgrid->hmarshp[base+n]=1;
        else
          subgrid->hmarshp[base+n]=0.0;
      }
    }

    // interpolate subedge depth
    ncount=grid->Ne*(subgrid->segN+1);
    Interp(xd,yd,d,Nd,&(subgrid->xpe[0]), &(subgrid->ype[0]),&(subgrid->hmarshpe[0]),ncount,grid->maxfaces);

    for(n=0;n<ncount;n++)
    {
      if(subgrid->hmarshpe[n]>=0.95*1)
        subgrid->hmarshpe[n]=1;
      else
        subgrid->hmarshpe[n]=0;
    }

  } else if(subgrid->hmarshint==0) {
    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        subgrid->hmarshp[base+n]=ReturnMarshHeight(subgrid->xp[base+n],subgrid->yp[base+n]);
        intd=subgrid->hmarshp[base+n];
        k=0;
        // adjust accuracy as the output grid.c
        if(intd>=1)
        {
          while(floor(intd)!=0)
          {
            intd=intd/10.0;
            k++;
          }
          intdd=floor(subgrid->hmarshp[base+n]*1e6/pow(10,k-1));
          if((subgrid->hmarshp[base+n]*1e6/pow(10,k-1)-intdd)>=0.5)
            intdd+=1;

          subgrid->hmarshp[base+n]=intdd/1e6*pow(10,k-1);
        } else if(intd<1 && intd!=0) {

          while(floor(intd)==0)
          {
            intd=intd*10.0;
            k++;
          }
          intdd=floor(subgrid->hmarshp[base+n]*1e6*pow(10,k));
          if((subgrid->hmarshp[base+n]*1e6*pow(10,k)-intdd)>=0.5)
            intdd+=1;
          subgrid->hmarshp[base+n]=intdd/1e6/pow(10,k);                
        }
      }
    }

    ncount=grid->Ne*(subgrid->segN+1);
    for(n=0;n<ncount;n++)
    {
      subgrid->hmarshpe[n]=ReturnMarshHeight(subgrid->xpe[n],subgrid->ype[n]);

      intd=subgrid->hmarshpe[n];
      k=0;
      if(intd>=1)
      {
        while(floor(intd)!=0)
        {
          intd=intd/10.0;
          k++;
        }
        intdd=floor(subgrid->hmarshpe[n]*1e6/pow(10,k-1));
        if((subgrid->hmarshpe[n]*1e6/pow(10,k-1)-intdd)>=0.5)
          intdd+=1;

        subgrid->hmarshpe[n]=intdd/1e6*pow(10,k-1);

      } else if(intd<1 && intd!=0){
        while(floor(intd)==0)
        {
          intd=intd*10.0;
          k++;
        }
        intdd=floor(subgrid->hmarshpe[n]*1e6*pow(10,k));
        if((subgrid->hmarshpe[n]*1e6*pow(10,k)-intdd)>=0.5)
          intdd+=1;
 
        subgrid->hmarshpe[n]=intdd/1e6/pow(10,k);  
      }      

    }

  } else {
    MPI_GetFile(str1,DATAFILE,"subhmarshinFile","InterpolateSubgridHmarsh",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ifile = MPI_FOpen(str,"r","InterpolateSubgridHmarsh",myproc);
    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        getfield(ifile,str);
        getfield(ifile,str);
        getfield(ifile,str);
        subgrid->hmarshp[base+n]=getfield(ifile,str);
      }
    }
    fclose(ifile);

    MPI_GetFile(str1,DATAFILE,"subhmarsheinFile","InterpolateSubgridHmarsh",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ifile = MPI_FOpen(str,"r","InterpolateSubgridHmarsh",myproc);
    ncount=grid->Ne*(subgrid->segN+1);
    for(n=0;n<ncount;n++)
    {
      getfield(ifile,str);
      getfield(ifile,str);
      getfield(ifile,str);
      subgrid->hmarshpe[n]=getfield(ifile,str);      
    }
    fclose(ifile);
  }

  // calculate mean vegetation height for sub cell and edge
  for(i=0;i<grid->Nc;i++)
  {
    ntotal=(grid->nfaces[i]-2)*subgrid->segN*subgrid->segN;
    base=i*(grid->maxfaces-2)*subgrid->segN*subgrid->segN;
    for(j=0;j<ntotal;j++)  
    {
      //subgrid->hmarshc[base+j]=(subgrid->hmarshp[subgrid->cellp[base*3+j*3]]\
      //  +subgrid->hmarshp[subgrid->cellp[base*3+j*3+1]]\
      //  +subgrid->hmarshp[subgrid->cellp[base*3+j*3+2]])/3;
      subgrid->hmarshc[base+j]=Min(subgrid->hmarshp[subgrid->cellp[base*3+j*3]],subgrid->hmarshp[subgrid->cellp[base*3+j*3+1]]);
      subgrid->hmarshc[base+j]=Min(subgrid->hmarshc[base+j],subgrid->hmarshp[subgrid->cellp[base*3+j*3+2]]);
    }
  }

  for(i=0;i<grid->Ne;i++)
  {
    base=i*subgrid->segN;
    base1=i*(subgrid->segN+1);
    for(j=0;j<subgrid->segN;j++)
       subgrid->hmarshe[base+j]=Min(subgrid->hmarshpe[base1+j],subgrid->hmarshpe[base1+j+1]);
  }
}

/*
 * Function: InterpolateSubgridCdV
 * Usage: calculate vegetation drag coefficient for all subpoints
 * ----------------------------------------------------
 * input the same file as in grid.c to interpolate 
 * depth for all sub points
 *
 */
void InterpolateSubgridCdV(gridT *grid, physT *phys, int myproc, int numprocs)
{
  int i, j, n, nc,nc1,nc2,ne, Nd, scaledepth,ncount,base,k,intdd,ntotal,base1;
  REAL *xd, *yd, *d,intd,min;
  char str[BUFFERLENGTH],str1[BUFFERLENGTH];
  FILE *ifile;

  if(subgrid->cdvint==1)
  {   
    MPI_GetFile(str,DATAFILE,"InputCdVFile","InterpolateSubgridCdV",myproc);
    Nd = MPI_GetSize(str,"InterpolateSubgridCdv",myproc);
    xd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridCdV");
    yd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridCdV");
    d = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridCdV");

    ifile = MPI_FOpen(str,"r","InterpolateSubgridCdv",myproc);
    for(n=0;n<Nd;n++) {
      xd[n]=getfield(ifile,str);
      yd[n]=getfield(ifile,str);
      d[n]=getfield(ifile,str);//fabs(getfield(ifile,str));   
    }
    fclose(ifile);
    // interpolate subcell depth
    ncount=grid->Nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
    Interp(xd,yd,d,Nd,&(subgrid->xp[0]), &(subgrid->yp[0]),&(subgrid->cdvp[0]),ncount,grid->maxfaces);

    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        if(subgrid->cdvp[base+n]>=0.95*1.3)
          subgrid->cdvp[base+n]=1.3;
        else
          subgrid->cdvp[base+n]=0;
      }
    }

    // interpolate subedge depth
    ncount=grid->Ne*(subgrid->segN+1);
    Interp(xd,yd,d,Nd,&(subgrid->xpe[0]), &(subgrid->ype[0]),&(subgrid->cdvpe[0]),ncount,grid->maxfaces);

    for(n=0;n<ncount;n++)
    {
      if(subgrid->cdvpe[n]>=0.95*1.3)
        subgrid->cdvpe[n]=1.3;
      else
        subgrid->cdvpe[n]=0;
    }

  } else if(subgrid->cdvint==0) {
    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        subgrid->cdvp[base+n]=ReturnMarshDragCoefficient(subgrid->xp[base+n],subgrid->yp[base+n]);
        intd=subgrid->cdvp[base+n];
        k=0;
        // adjust accuracy as the output grid.c
        if(intd>=1)
        {
          while(floor(intd)!=0)
          {
            intd=intd/10.0;
            k++;
          }
          intdd=floor(subgrid->cdvp[base+n]*1e6/pow(10,k-1));
          if((subgrid->cdvp[base+n]*1e6/pow(10,k-1)-intdd)>=0.5)
            intdd+=1;

          subgrid->cdvp[base+n]=intdd/1e6*pow(10,k-1);
        } else if(intd<1 && intd!=0) {

          while(floor(intd)==0)
          {
            intd=intd*10.0;
            k++;
          }
          intdd=floor(subgrid->cdvp[base+n]*1e6*pow(10,k));
          if((subgrid->cdvp[base+n]*1e6*pow(10,k)-intdd)>=0.5)
            intdd+=1;
          subgrid->cdvp[base+n]=intdd/1e6/pow(10,k);                
        }
      }
    }

    ncount=grid->Ne*(subgrid->segN+1);
    for(n=0;n<ncount;n++)
    {
      subgrid->cdvpe[n]=ReturnMarshDragCoefficient(subgrid->xpe[n],subgrid->ype[n]);

      intd=subgrid->cdvpe[n];
      k=0;
      if(intd>=1)
      {
        while(floor(intd)!=0)
        {
          intd=intd/10.0;
          k++;
        }
        intdd=floor(subgrid->cdvpe[n]*1e6/pow(10,k-1));
        if((subgrid->cdvpe[n]*1e6/pow(10,k-1)-intdd)>=0.5)
          intdd+=1;

        subgrid->cdvpe[n]=intdd/1e6*pow(10,k-1);

      } else if(intd<1 && intd!=0){
        while(floor(intd)==0)
        {
          intd=intd*10.0;
          k++;
        }
        intdd=floor(subgrid->cdvpe[n]*1e6*pow(10,k));
        if((subgrid->cdvpe[n]*1e6*pow(10,k)-intdd)>=0.5)
          intdd+=1;
 
        subgrid->cdvpe[n]=intdd/1e6/pow(10,k);  
      }      

    }

  } else {
    MPI_GetFile(str1,DATAFILE,"subcdvinFile","InterpolateSubgridHmarsh",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ifile = MPI_FOpen(str,"r","InterpolateSubgridCdV",myproc);
    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        getfield(ifile,str);
        getfield(ifile,str);
        getfield(ifile,str);
        subgrid->cdvp[base+n]=getfield(ifile,str);
      }
    }
    fclose(ifile);

    MPI_GetFile(str1,DATAFILE,"subcdveinFile","InterpolateSubgridCdV",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ifile = MPI_FOpen(str,"r","InterpolateSubgridCdV",myproc);
    ncount=grid->Ne*(subgrid->segN+1);
    for(n=0;n<ncount;n++)
    {
      getfield(ifile,str);
      getfield(ifile,str);
      getfield(ifile,str);
      subgrid->cdvpe[n]=getfield(ifile,str);      
    }
    fclose(ifile);
  }

  // calculate mean vegetation height for sub cell and edge
  for(i=0;i<grid->Nc;i++)
  {
    ntotal=(grid->nfaces[i]-2)*subgrid->segN*subgrid->segN;
    base=i*(grid->maxfaces-2)*subgrid->segN*subgrid->segN;
    for(j=0;j<ntotal;j++)  
      subgrid->cdvc[base+j]=(subgrid->cdvp[subgrid->cellp[base*3+j*3]]\
        +subgrid->cdvp[subgrid->cellp[base*3+j*3+1]]\
        +subgrid->cdvp[subgrid->cellp[base*3+j*3+2]])/3;
  }

  for(i=0;i<grid->Ne;i++)
  {
    base=i*subgrid->segN;
    base1=i*(subgrid->segN+1);
    for(j=0;j<subgrid->segN;j++)
       subgrid->cdve[base+j]=Min(subgrid->cdvpe[base1+j],subgrid->cdvpe[base1+j+1]);
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
void InterpolateSubgridDepth(gridT *grid, physT *phys, propT *prop, int myproc,int numprocs)
{
  int i, j, n, nc,nc1,nc2,ne, Nd, scaledepth,ncount,base,base1,k,intdd,ntotal;
  REAL *xd, *yd, *d, scaledepthfactor,depthelev,intd,min,d1,d2,d3,d0;
  char str[BUFFERLENGTH],str1[BUFFERLENGTH];
  FILE *ifile;

  scaledepth=(int)MPI_GetValue(DATAFILE,"scaledepth","InterpolateSubgridDepth",myproc);
  scaledepthfactor=MPI_GetValue(DATAFILE,"scaledepthfactor","InterpolateSubgridDepth",myproc);
  depthelev=MPI_GetValue(DATAFILE,"depthelev","InterpolateSubgridDepth",myproc); 

  if(subgrid->dpint==1)
  {
    MPI_GetFile(str,DATAFILE,"subgriddFile","InterpolateSubgridDepth",myproc);
    Nd = MPI_GetSize(str,"InterpolateSubgridDepth",myproc);
    xd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridDepth");
    yd = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridDepth");
    d = (REAL *)SunMalloc(Nd*sizeof(REAL),"InterpolateSubgridDepth");

    ifile = MPI_FOpen(str,"r","InterpolateSubgridDepth",myproc);
    for(n=0;n<Nd;n++) {
      xd[n]=getfield(ifile,str);
      yd[n]=getfield(ifile,str);
      d[n]=getfield(ifile,str);//fabs(getfield(ifile,str));
      if(scaledepth)
        d[n]*=scaledepthfactor;
      d[n]+=depthelev;    
    }
    fclose(ifile);
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


    // get subgrid->hmin subgrid->hmax
    for(i=0;i<grid->Nc;i++)
    {
      for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
      {
        if(j==0)
          subgrid->hmax[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];

        if(subgrid->hmin[i]>-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmin[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];

        if(subgrid->hmax[i]<-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmax[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
      }
      if(grid->dv[i]>-subgrid->hmin[i]){
        for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
          if(subgrid->hmin[i]==-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
            subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j]=grid->dv[i];        
        subgrid->hmin[i]=-grid->dv[i];      
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
  } else if(subgrid->dpint==0) {
    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        subgrid->dp[base+n]=ReturnSubgridPointDepth(subgrid->xp[base+n],subgrid->yp[base+n],grid->xv[nc],grid->yv[nc]);
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
        if(j==0)
          subgrid->hmax[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
        if(subgrid->hmax[i]<-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmax[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];        
        if(subgrid->hmin[i]>-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmin[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
      }
      if(grid->dv[i]>-subgrid->hmin[i]){
        for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
          if(subgrid->hmin[i]==-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
            subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j]=grid->dv[i];        
        subgrid->hmin[i]=-grid->dv[i];      
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

    for(nc=0;nc<grid->Ne;nc++){
      nc1=grid->grad[2*nc];
      nc2=grid->grad[2*nc+1];
      if(nc1==-1)
        nc1=nc2;
      if(nc2==-1)
        nc2=nc1;
    }

  } else {
    MPI_GetFile(str1,DATAFILE,"subpointinFile","InterpolateSubgridDepth",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ifile = MPI_FOpen(str,"r","InterpolateSubgridDepth",myproc);
    for(nc=0;nc<grid->Nc;nc++){
      ncount=(grid->nfaces[nc]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      base=nc*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2;
      for(n=0;n<ncount;n++)
      { 
        getfield(ifile,str);
        getfield(ifile,str);
        getfield(ifile,str);
        subgrid->dp[base+n]=getfield(ifile,str);
        if(scaledepth)
          subgrid->dp[base+n]*=scaledepthfactor;
        subgrid->dp[base+n]+=depthelev; 
        if(prop->culvertmodel)
          if(culvert->top[nc]!=INFTY){
            subgrid->dp[base+n]=grid->dv[nc]; 
          }
      }
    }
    fclose(ifile);

    // get subgrid->hmin
    for(i=0;i<grid->Nc;i++)
    {
      for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
      {
        if(j==0)
          subgrid->hmax[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
        if(subgrid->hmax[i]<-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmax[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
        if(subgrid->hmin[i]>-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
          subgrid->hmin[i]=-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j];
      }
      if(grid->dv[i]>-subgrid->hmin[i]){
        for(j=0;j<((grid->nfaces[i]-2)*(subgrid->segN+1)*(subgrid->segN+2)/2);j++)
          if(subgrid->hmin[i]==-subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j])
            subgrid->dp[i*(grid->maxfaces-2)*(subgrid->segN+1)*(subgrid->segN+2)/2+j]=grid->dv[i];        
        subgrid->hmin[i]=-grid->dv[i];      
      }
    }

    MPI_GetFile(str1,DATAFILE,"subpointeinFile","InterpolateSubgridDepth",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ifile = MPI_FOpen(str,"r","InterpolateSubgridDepth",myproc);
    ncount=grid->Ne*(subgrid->segN+1);
    for(n=0;n<ncount;n++)
    {
      getfield(ifile,str);
      getfield(ifile,str);
      getfield(ifile,str);
      subgrid->dpe[n]=getfield(ifile,str);   
      if(scaledepth)
        subgrid->dpe[n]*=scaledepthfactor;
      subgrid->dpe[n]+=depthelev;       
    }
    for(nc=0;nc<grid->Ne;nc++){
      nc1=grid->grad[2*nc];
      nc2=grid->grad[2*nc+1];
      if(nc1==-1)
        nc1=nc2;
      if(nc2==-1)
        nc2=nc1;
      if(prop->culvertmodel)
        if(culvert->top[nc1]!=INFTY || culvert->top[nc2]!=INFTY)
        {  
          if(culvert->top[nc1]!=INFTY && culvert->top[nc2]==INFTY)
            for(n=0;n<(subgrid->segN+1);n++)
              subgrid->dpe[nc*(subgrid->segN+1)+n]=grid->dv[nc1];
          if(culvert->top[nc2]!=INFTY && culvert->top[nc1]==INFTY)
            for(n=0;n<(subgrid->segN+1);n++)
              subgrid->dpe[nc*(subgrid->segN+1)+n]=grid->dv[nc2];
          if(culvert->top[nc2]!=INFTY && culvert->top[nc1]!=INFTY)
            for(n=0;n<(subgrid->segN+1);n++)
              subgrid->dpe[nc*(subgrid->segN+1)+n]=Min(grid->dv[nc2],grid->dv[nc1]);
        }
    }
    fclose(ifile);
  }

  // calculate mean depth for sub cell and edge
  for(i=0;i<grid->Nc;i++)
  {
    ntotal=(grid->nfaces[i]-2)*subgrid->segN*subgrid->segN;
    base=i*(grid->maxfaces-2)*subgrid->segN*subgrid->segN;
    for(j=0;j<ntotal;j++) { 
      d1=subgrid->dp[subgrid->cellp[base*3+j*3]];
      d2=subgrid->dp[subgrid->cellp[base*3+j*3+1]];
      d3=subgrid->dp[subgrid->cellp[base*3+j*3+2]];
      if(d1>d3)
      {
        d0=d1;
        d1=d3;
        d3=d0;
      }
      if(d2>d3)
      {
        d0=d2;
        d2=d3;
        d3=d0;
      }
      if(d2<d1)
      {
        d0=d2;
        d2=d1;
        d1=d0;
      }
      subgrid->d1[base+j]=d1;
      subgrid->d2[base+j]=d2;
      subgrid->d3[base+j]=d3;
      subgrid->dc[base+j]=(d1+d2+d3)/3;
    }
  }

  for(i=0;i<grid->Ne;i++)
  {
    base=i*subgrid->segN;
    base1=i*(subgrid->segN+1);
    for(j=0;j<subgrid->segN;j++)
       subgrid->de[base+j]=(subgrid->dpe[base1+j]+subgrid->dpe[base1+j+1])/2;
  }
}

/*
 * Function: CalculateAreaRatio
 * Usage: calculate the area ratio for each sub-triangle for a cell
 * since the cell is divided into several triangles and then to do 
 * subgrid discretization. This ratio is used for the calculation of wet 
 * area
 * ----------------------------------------------------
 *
 */
void CalculateAreaRatio(gridT *grid, int myproc)
{
  REAL xt[3],yt[3];
  int n,nf,i;
  for(n=0;n<grid->Nc;n++) {
    for(i=0;i<(grid->nfaces[n]-2);i++){        
      xt[0]=grid->xp[grid->cells[n*grid->maxfaces+i+1]];
      xt[1]=grid->xp[grid->cells[n*grid->maxfaces+i+2]];
      xt[2]=grid->xp[grid->cells[n*grid->maxfaces]];
      yt[0]=grid->yp[grid->cells[n*grid->maxfaces+i+1]];
      yt[1]=grid->yp[grid->cells[n*grid->maxfaces+i+2]];
      yt[2]=grid->yp[grid->cells[n*grid->maxfaces]];
      subgrid->Acratio[n*(grid->maxfaces-2)+i]=GetArea(xt,yt,3)/grid->Ac[n];
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
  //ac1=subgrid->Acbackup[nc]/ntotal;
  Ac=0;
  for(i=0;i<ntotal;i++)
  { 
    ac1=subgrid->Acbackup[nc]/N/N*subgrid->Acratio[nc*(maxfaces-2)+i/(N*N)];
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
  REAL Ac;
  
  // get the profile of Ac for each cell   
  for(nc=0;nc<grid->Nc;nc++)
  {
    base=nc*(subgrid->disN+1);
    subgrid->Acprof[base+subgrid->disN]=subgrid->Acbackup[nc];
    
    for(i=0;i<subgrid->disN;i++)
    {
      Ac=CalculateWetArea(nc, grid->nfaces[nc], grid->maxfaces, subgrid->hprof[base+i]);
      if(Ac>subgrid->Acbackup[nc])
        Ac=subgrid->Acbackup[nc];
      subgrid->Acprof[base+i]=Ac;
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
          r=sqrt(dis);
      } else {
        if(h<=-d2)
          r=0;
        else if(h>-d2 && h<=-d1)
        { 
          r1=(h+d2)/(d2-d1);
          r=r1*sqrt(dis); //pow((pow((d1-d2),2)+dis),0.5);
        } else {
          r=sqrt(dis);//pow((pow((d1-d2),2)+dis),0.5);
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
        r=sqrt(dis);//pow((pow((d2-d1),2)+dis),0.5);
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
    if(subgrid->Wetperiprof[base+subgrid->disN]<grid->df[ne])
      subgrid->Wetperiprof[base+subgrid->disN]=grid->df[ne];
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
  int i, k, j, base,n;
  REAL h;
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
      if(j!=subgrid->disN && subgrid->Acprof[base+j]>subgrid->Acprof[base+j+1])
      {
        printf("i %d cell j %dth Acprof is %f is bigger than j+1 %d th %f\n",i,j,subgrid->Acprof[base+j],(j+1),subgrid->Acprof[base+j+1]);
        printf("i %d nface %d Ac %e Acratio %e",i,grid->nfaces[i],grid->Ac[i],subgrid->Acratio[i*(grid->maxfaces-2)]);
      }
      if(j!=subgrid->disN && subgrid->Vprof[base+j]>subgrid->Vprof[base+j+1])
        printf("i %d cell j %dth Vprof is %f is bigger than j+1 %d th\n",i,j,subgrid->Vprof[base+j],(j+1));
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

  // check updatefreesurface function
  for(i=0;i<grid->Nc;i++){
    base=i*(subgrid->disN+1);
    for(j=0;j<=subgrid->disN;j++){
      h=UpdateFreeSurface(i,subgrid->Vprof[base+j]);
      if(fabs(subgrid->hprof[base+j]-h)>1e-5)
        printf("something wrong in update free surface in subgrid check\n");
    }
  } 

  // check subgrid dmax
  n=0;
  for(i=0;i<grid->Nc;i++)
  {
    if(fabs(grid->dv[i]+subgrid->hmin[i])>1e-5){
      n++;
      printf("dv is not consistent with dmax in subgrid bathymetry\n");
      printf("proc %d i %d xv %e yv %e dv %e hmin %e\n",myproc,i,grid->xv[i],grid->yv[i],grid->dv[i],subgrid->hmin[i]);
    }
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
  ReadSubgridProperties(prop, myproc);

  // allocate space for subgrid method
  AllocateandInitializeSubgrid(grid, prop, myproc);

  //open subgrid outputfiles
  OpenSubgridFiles(prop->computeSediments,prop->mergeArrays,myproc);

  // calculate all the x y for subcell and subedege
  CalculateSubgridXY(grid, myproc);

  // calculate the index for each subcell
  CalculateSubgridCellp(grid, myproc);

  // calculate all the for all the points in subgrid
  InterpolateSubgridDepth(grid, phys,prop,myproc,numprocs);

  // calculate the max and min depth to get a reasonable range for h
  CalculateHProfile(grid, myproc);

  // calculate the Area ratio between different sub triangles for each cell A_sub/Ac
  CalculateAreaRatio(grid, myproc);
  
  // calculate the wet area profile for future use
  CalculateAcProfile(grid, myproc);

  // calculate the cell volume profile for future use
  CalculateVolumeProfile(grid, myproc);

  // calculate the edge flux height profile for future use
  CalculateFluxHeightProfile(grid, myproc);
  
  // calculate wet perimeter profile for the last layer
  CalculateWetperimeterProfile(grid,myproc);

  // calculate subgrid marsh setup for all subcell and subedge
  if(prop->marshmodel && subgrid->dragpara)
  {
    InterpolateSubgridHmarsh(grid,phys,myproc,numprocs);
    InterpolateSubgridCdV(grid,phys,myproc,numprocs);
  }

  // check subgrid
  SubgridCheck(grid,myproc);
  
  // output subgrid for test
  OutputSubgrid(grid,phys,prop,myproc,numprocs,comm);
  
  if(myproc==0)
    printf("finish subgrid basic module preparation\n");
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
    //if(h+grid->dv[nc]<DRYCELLHEIGHT)
      //h=-grid->dv[nc]+DRYCELLHEIGHT;
    // no need to add the following lines by the new connection with subgrid
    // and culvert model
    //if(prop->culvertmodel)
     // if(h>culvert->top[nc])
      //  h=culvert->top[nc];
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
  REAL dh,V,dV,dac,ac;

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
    dac=subgrid->Acprof[base+i]-subgrid->Acprof[base+i-1]; 
    ac=subgrid->Acprof[base+i-1]+(h-subgrid->hprof[base+i-1])/dh*dac;
    V=subgrid->Vprof[base+i-1]+0.5*(subgrid->Acprof[base+i-1]+ac)*(h-subgrid->hprof[base+i-1]);
  }
  return V;
}

/*
 * Function: UpdateFreeSurface
 * Usage: basic function to calculate Free Surface based on water volume
 * ----------------------------------------------------
 * use hprof to get water free surface;
 *
 */
REAL UpdateFreeSurface(int nc, REAL V)
{
  int  i,base;
  REAL dh,h,dV,dac,ac,temp;

  base=nc*(subgrid->disN+1); 
  if(V==0)
    h=subgrid->hprof[base];
  else if(V>=subgrid->Vprof[base+subgrid->disN])
  {
    h=subgrid->hprof[base+subgrid->disN]+(V-subgrid->Vprof[base+subgrid->disN])/subgrid->Acbackup[nc];
  } else {
    i=1;
    while(i<=subgrid->disN)
    {
      if(V<subgrid->Vprof[base+i])
        break;
      i++;
    }
    dV=V-subgrid->Vprof[base+i-1];
    dac=subgrid->Acprof[base+i]-subgrid->Acprof[base+i-1]; 
    ac=subgrid->Acprof[base+i-1];
    if(dac!=0){
      dh=subgrid->hprof[base+i]-subgrid->hprof[base+i-1];
      temp=sqrt(ac*ac+2*dac/dh*dV);
      h=subgrid->hprof[base+i-1]+(-ac+temp)/dac*dh;
    }else{
      h=subgrid->hprof[base+i-1]+dV/ac;
    }
    if(h!=h)
      printf("nc %d dac %e dh %e i %d Acprof %e %e \n",nc,dac,dh,i,subgrid->Acprof[base+i],subgrid->Acprof[base+i-1]);
  }
  return h;
}

/*
 * Function: UpdateSubgridVeff
 * Usage: update phys->h based on subgrid->Veff
 * ----------------------------------------------------
 * use hprof to get free surface;
 *
 */
void UpdateSubgridFreeSurface(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int nc, i,base;
  REAL dh,h,V,dV;
  for(nc=0;nc<grid->Nc;nc++)
  {
    V=subgrid->Veff[nc];
    phys->h[nc]=UpdateFreeSurface(nc,V);
    if(prop->culvertmodel){
      // compare to culvert pressure to ensure culvert pressure is larger than h
      if(culvert->top[nc]!=INFTY)
      {
        if(culvert->pressure[nc]<phys->h[nc])
          phys->h[nc]=culvert->pressure[nc];
        else
          if(culvert->pressure[nc]<=culvert->top[nc])
          {
            culvert->pressure[nc]=phys->h[nc];
            culvert->pressure2[nc]=phys->h[nc];
          } else 
            phys->h[nc]=culvert->pressure[nc];
      }

      if(phys->h[nc]>culvert->top[nc])
        phys->h[nc]=culvert->top[nc];
      // for no culvert-> cell, also update culvert-> pressure to the free surface
      if(culvert->top[nc]==INFTY){
        culvert->pressure2[nc]=phys->h[nc];
        culvert->pressure[nc]=phys->h[nc];
      }
    }
    if(phys->h[nc]<=(-grid->dv[nc]+DRYCELLHEIGHT))
    {
      phys->h[nc]=-grid->dv[nc]+DRYCELLHEIGHT;
      phys->active[nc]=0;
    }
    h=phys->h[nc];
    V=UpdateVeff(nc,h);
    if(fabs(subgrid->Veff[nc]-V)>1e-3){
      if(!prop->culvertmodel)
        printf("something wrong to calculate free surface from volume myproc %d nc %d V %e Veff %e diff %e h %e\n",myproc,nc,V,subgrid->Veff[nc],V-subgrid->Veff[nc],h);
      else
        if(culvert->top[nc]==INFTY)
          printf("something wrong to calculate free surface from volume myproc %d nc %d V %e Veff %e diff %e h %e\n",myproc,nc,V,subgrid->Veff[nc],V-subgrid->Veff[nc],h);
    }
  }
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
    //subgrid->Aceffold[nc]=subgrid->Aceff[nc];
    // for wetting and drying and get rid of aceff=0
    //if(h+grid->dv[nc]<DRYCELLHEIGHT)
      //h=-grid->dv[nc]+DRYCELLHEIGHT;
    subgrid->Aceff[nc]=UpdateAceff(nc,h);
    if(prop->culvertmodel)
      if(h>culvert->top[nc])
        subgrid->Aceff[nc]=UpdateAceff(nc,culvert->top[nc]);
    if(prop->wetdry)
      if(subgrid->Aceff[nc]<=grid->Ac[nc]/1000)
        subgrid->Aceff[nc]=grid->Ac[nc]/1000;
  }
}

/*
 * Function: SubgridCulverttopArea
 * Usage: Calculate the value for Qcoef for outer loop of Culvert model
 * ----------------------------------------------------
 * Only calculate the wet area for where there is a culvert
 * otherwise 0
 *
 */
void SubgridCulverttopArea(gridT *grid, propT *prop, int myproc)
{
  int nc;
  for(nc=0;nc<grid->Nc;nc++)
    if(culvert->top[nc]!=INFTY)
      culvert->toparea[nc]=UpdateAceff(nc,culvert->top[nc]);
    else  
      culvert->toparea[nc]=0;
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
  REAL v,ac,ztop,zbot,dznew,Vc;
  int i,k,ktop;
   
  // prepare for Updatescalars
  if(option==0){
    for(i=0;i<grid->Nc;i++)
    {
      Vc=0;
      // store the old values
      for(k=0;k<grid->Nk[i];k++)
        subgrid->Acceffold[i][k]=subgrid->Acceff[i][k];
 
      for(k=0;k<=grid->Nk[i];k++){
        subgrid->Acveffold2[i][k]=subgrid->Acveffold[i][k];
        subgrid->Acveffold[i][k]=subgrid->Acveff[i][k];
      }

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
        //if(grid->dzz[i][k]>=0.9*DRYCELLHEIGHT)
          subgrid->Acceff[i][k]=(UpdateVeff(i,ztop)-UpdateVeff(i,zbot))/grid->dzz[i][k];
        //else
          //subgrid->Acceff[i][k]=0;
        Vc+=subgrid->Acceff[i][k]*grid->dzz[i][k];

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
        //zbot=-grid->dv[i];
        zbot=subgrid->hmin[i]; // use the subgrid hmin to get rid of 0 acceff
        ztop=zbot+dznew;
      }

      //if(ktop==(grid->Nk[i]-1))
        subgrid->Acceff[i][ktop]=(subgrid->Veff[i]-Vc)/dznew;//
      //else{
        //if(subgrid->Verr[i]<0)
      //subgrid->Acceff[i][ktop]=(subgrid->Veff[i]-Vc)/dznew;//
        //else
          //subgrid->Acceff[i][ktop]=(UpdateVeff(i,ztop)-UpdateVeff(i,zbot))/dznew;
      //}

      if(subgrid->Acceff[i][ktop]<0)
      {
        printf("something wrong in calculate subgrid acceff which induces negative acceff myproc %d n %d i %d Veff %e Veffold %e Vc %e ktop %d Nk %d Vlimit %e H %e h+dv %e dzz %e\n",
          myproc,prop->n,i,subgrid->Veff[i],subgrid->Veffold[i],Vc,ktop,grid->Nk[i],grid->Ac[i]*DRYCELLHEIGHT,phys->h[i],phys->h[i]+grid->dv[i],grid->dzz[i][grid->Nk[i]-1]);
        exit(1);
      }
      // combine all new cells to one cell at top
      //if(dznew>=0.9*DRYCELLHEIGHT)
      //else
        //subgrid->Acceff[i][ktop]=0;
      
      // update acveff to new ctop
      for(k=ktop;k>=grid->ctop[i];k--){
        ztop=zbot+grid->dzz[i][k];
        subgrid->Acveff[i][k]=UpdateAceff(i,ztop);
        zbot=ztop;
      }

    }
  } else {
    for(i=0;i<grid->Nc;i++)
    {
      Vc=0;
      if(grid->ctop[i]<grid->ctopold[i])
      {   
        // from bottom to top
        zbot=-grid->dv[i];
        ztop=zbot+grid->dzz[i][grid->Nk[i]-1];

        for(k=grid->Nk[i]-1;k>grid->ctop[i];k--)
        {  
          // update mean cell center wet area based on volume
          //if(grid->dzz[i][k]>=0.9*DRYCELLHEIGHT)
          subgrid->Acceff[i][k]=(UpdateVeff(i,ztop)-UpdateVeff(i,zbot))/grid->dzz[i][k];
          Vc+=subgrid->Acceff[i][k]*grid->dzz[i][k];
          //else
            //subgrid->Acceff[i][k]=0;

          //subgrid->Acveff[i][k]=UpdateAceff(i,ztop);
          zbot=ztop;
          ztop=zbot+grid->dzz[i][k-1];
        }
        subgrid->Acceff[i][grid->ctop[i]]=(subgrid->Veff[i]-Vc)/grid->dzz[i][grid->ctop[i]];

        if(subgrid->Acceff[i][grid->ctop[i]]<0)
          printf("something wrong in acceff\n");
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
  int ne, nc,nc1,nc2, nc3, i,base,k,j,nc_b;
  REAL dh,h,fh,fhr,dfhr,dfh,culverttop,u,hbot,hc,wetperi,Af,zc;
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
    {
      nc=nc1;
      nc_b=nc2;
    } else {
      nc=nc2;
      nc_b=nc1;
    }
    
    // h is the free surface height used to calculate flux height can be upwind or central differencing
    if(subgrid->dzfmeth==2)
      h=0.5*(phys->h[nc1]+phys->h[nc2]); 
    else
      h=UpWind(u,phys->h[nc1],phys->h[nc2]);

    
    if(prop->culvertmodel)
    {
      culverttop=UpWind(phys->u[ne][0],culvert->top[nc1],culvert->top[nc2]);
      if(h>culverttop)
        h=culverttop;
    }

    // for drag term we always use h=mean(hnc1,hnc2) check phys.c
    subgrid->he[ne]=h;
    
    for(k=0;k<grid->etop[ne];k++)
      grid->dzf[ne][k]=0;
    
    if(prop->vertcoord==1){
      hbot=hc;    
      for(k=grid->ctop[nc];k<grid->etop[ne];k++)
        hbot-=grid->dzz[nc][k];

      for(k=grid->etop[ne];k<grid->Nke[ne]-1;k++)
      {
        hbot-=grid->dzz[nc][k];//grid->dzf[ne][k]; // still old dzf from the original dzf
        grid->dzf[ne][k]=CalculateFluxHeight(ne,h)-CalculateFluxHeight(ne,hbot);//UpdateFluxHeight(ne,h)-UpdateFluxHeight(ne,hbot);
        h=hbot;
      }
      grid->dzf[ne][grid->Nke[ne]-1]=CalculateFluxHeight(ne,h);
    } else {
      for(k=grid->etop[ne];k<grid->Nke[ne];k++)
      {
        if(subgrid->dzfmeth==2)
          grid->dzf[ne][k]=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
        else{
          hbot=UpWind(phys->u[ne][k],vert->zc[nc1][k]-grid->dzz[nc1][k]/2,vert->zc[nc2][k]-grid->dzz[nc2][k]/2);
          h=UpWind(phys->u[ne][k],vert->zc[nc1][k]+grid->dzz[nc1][k]/2,vert->zc[nc2][k]+grid->dzz[nc2][k]/2);
          grid->dzf[ne][k]=CalculateFluxHeight(ne,h)-CalculateFluxHeight(ne,hbot);//UpdateFluxHeight(ne,h)-UpdateFluxHeight(ne,hbot);
          if(k==vert->Nkeb[ne])
            wetperi=CalculateWetperimeter(ne,h);
        }
      }      
    }

   // if(prop->culvertmodel==1)
     // if(culvert->top[nc1]!=INFTY && culvert->top[nc2]!=INFTY)
       // grid->dzf[ne][grid->Nke[ne]-1]=grid->dzf[ne][grid->Nke[ne]-1];  

    if(((grid->dv[nc]+phys->h[nc])<=2*DRYCELLHEIGHT))
    {
      grid->dzf[ne][grid->Nke[ne]-1]=0;
    }
    
    if(prop->vertcoord==1)
      wetperi=CalculateWetperimeter(ne,h);

    // added part
    k=grid->Nke[ne]-1;
    
    if(grid->etop[ne]==k && grid->mark[ne]==2 && grid->dzf[ne][k]<0.01)
    {  
       grid->dzf[ne][k]=0.01;
       wetperi=grid->df[ne];
    }

    for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
      if(grid->dzf[ne][k]<=2*DRYCELLHEIGHT)
        grid->dzf[ne][k]=0;
    
    k=grid->Nke[ne]-1;
    if(grid->etop[ne]==k && grid->mark[ne]!=2)
      if(subgrid->Veff[nc]<=2*grid->Ac[nc]*DRYCELLHEIGHT)
        grid->dzf[ne][k]=0;

    if(subgrid->Veff[nc1]>0 && subgrid->Veff[nc2]>0)
      nc=1;
    else 
      nc=0;
    
    if(prop->vertcoord==1)
      if(grid->dzf[ne][grid->Nke[ne]-1]==0)
        wetperi=grid->df[ne];

    subgrid->dzboteff[ne]=grid->dzf[ne][grid->Nke[ne]-1]*grid->df[ne]/wetperi;
    
    if(prop->vertcoord!=1){
      Af=0;
      for(k=vert->Nkeb[ne];k<grid->Nke[ne];k++)
        Af+=grid->dzf[ne][k]*grid->df[ne];
      if(Af==0)
        wetperi=grid->df[ne];
      subgrid->dzboteff[ne]=Af/wetperi;
    }
    
    if(subgrid->dzboteff[ne]<2*DRYCELLHEIGHT)
      subgrid->dzboteff[ne]=DRYCELLHEIGHT;
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
  char str[BUFFERLENGTH], str1[BUFFERLENGTH];
  FILE *ofile; 

  // output hprof for each cell
  MPI_GetFile(str1,DATAFILE,"HprofFile","OutputSubgrid",myproc);
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
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
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
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
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
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
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
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
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
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
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
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
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
  
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*(pow(subgrid->segN,2)*3)*(grid->maxfaces-2);
    for(j=0;j<((grid->nfaces[i]-2)*(pow(subgrid->segN,2)*3));j++)
      fprintf(ofile,"%d ",subgrid->cellp[base+j]);
    fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output xp yp dp for each cell
  MPI_GetFile(str1,DATAFILE,"subpointFile","OutputSubgrid",myproc);
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Nc;i++) 
  {
    base=i*((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->maxfaces-2);
    for(j=0;j<(((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->nfaces[i]-2));j++){
      fprintf(ofile,"%d %e %e %e\n",i,subgrid->xp[base+j],subgrid->yp[base+j],subgrid->dp[base+j]);
    }
    //fprintf(ofile,"\n");
  }
  fclose(ofile);

  // output xp yp hmarshp cdvp hmarshpe cdvpe for each cell
  if(prop->marshmodel)
  {
    MPI_GetFile(str1,DATAFILE,"subhmarshFile","OutputSubgrid",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
    for(i=0;i<grid->Nc;i++) 
    {
      base=i*((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->maxfaces-2);
      for(j=0;j<(((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->nfaces[i]-2));j++)
        fprintf(ofile,"%d %e %e %e\n",i,subgrid->xp[base+j],subgrid->yp[base+j],subgrid->hmarshp[base+j]);
    }
    fclose(ofile);

    MPI_GetFile(str1,DATAFILE,"subcdvFile","OutputSubgrid",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
    for(i=0;i<grid->Nc;i++) 
    {
      base=i*((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->maxfaces-2);
      for(j=0;j<(((subgrid->segN+1)*(subgrid->segN+2)/2)*(grid->nfaces[i]-2));j++)
        fprintf(ofile,"%d %e %e %e\n",i,subgrid->xp[base+j],subgrid->yp[base+j],subgrid->cdvp[base+j]);
    }
    fclose(ofile);

    // output xpe ype hmarshpe for each edge
    MPI_GetFile(str1,DATAFILE,"subhmarsheFile","OutputSubgrid",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
    for(i=0;i<grid->Ne;i++)
    { 
      base=i*(subgrid->segN+1);
      for(j=0;j<=subgrid->segN;j++)
        fprintf(ofile,"%d %e %e %e\n",i,subgrid->xpe[base+j],subgrid->ype[base+j],subgrid->hmarshpe[base+j]);    
    }
    fclose(ofile);  

    // output xpe ype cdvpe for each edge
    MPI_GetFile(str1,DATAFILE,"subcdveFile","OutputSubgrid",myproc);
    if(numprocs>1)
      sprintf(str,"%s.%d",str1,myproc);
    else
      sprintf(str,"%s",str1);
    ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
    for(i=0;i<grid->Ne;i++)
    { 
      base=i*(subgrid->segN+1);
      for(j=0;j<=subgrid->segN;j++)
        fprintf(ofile,"%d %e %e %e\n",i,subgrid->xpe[base+j],subgrid->ype[base+j],subgrid->cdvpe[base+j]);    
    }
    fclose(ofile);    
  }

  // output xpe ype dpe for each edge
  MPI_GetFile(str1,DATAFILE,"subpointeFile","OutputSubgrid",myproc);
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
  for(i=0;i<grid->Ne;i++)
  { 
    base=i*(subgrid->segN+1);
    for(j=0;j<=subgrid->segN;j++)
      fprintf(ofile,"%d %e %e %e\n",i,subgrid->xpe[base+j],subgrid->ype[base+j],subgrid->dpe[base+j]);    
    //fprintf(ofile,"\n");  
  }
  fclose(ofile);

  // output subgrid->hmin for each cell
  MPI_GetFile(str1,DATAFILE,"subdmaxFile","OutputSubgrid",myproc);
  if(numprocs>1)
    sprintf(str,"%s.%d",str1,myproc);
  else
    sprintf(str,"%s",str1);
  ofile = MPI_FOpen(str,"w","OutputSubgrid",myproc);
  
  for(i=0;i<grid->Nc;i++) 
    fprintf(ofile,"%e %e %e\n",grid->xv[i],grid->yv[i],(-subgrid->hmin[i]));
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


/*
 * Function: StoreSubgridAceffandVeffOld
 * Usage: Store the old value of Aceff and Veff
 * ----------------------------------------------------
 * Veffold used in Culvert model
 */
void StoreSubgridOldAceffandVeff(gridT *grid, int myproc)
{
  int i;
  for(i=0;i<grid->Nc;i++)
  {
    subgrid->Aceffold[i]=subgrid->Aceff[i];
    subgrid->Veffold[i]=subgrid->Veff[i];
  }
}


/*
 * Function: UpdateSubgridHeff
 * Usage: Calculate the average total depth of a cell 
 * ----------------------------------------------------
 */
void UpdateSubgridHeff(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int i;
  for(i=0;i<grid->Nc;i++)
  {
    subgrid->Heffold[i]=subgrid->Heff[i];
    subgrid->Heff[i]=subgrid->Veff[i]/subgrid->Aceff[i];
    if(subgrid->Aceff[i]/grid->Ac[i]<1e-5)
      subgrid->Heff[i]=DRYCELLHEIGHT;
  }
}

/*
 * Function: CalculateSubgridActiveErosion(latest version with tau_c_mean)
 * Usage: Calculate the active erosion for 1D model
 * assume a balance between freesurface gradient and bottom drag
 * use derived parameterization
 * ----------------------------------------------------
 * Can only be used for grid->Nkmax=1
 */
void CalculateSubgridActiveErosion(gridT *grid, physT *phys, propT *prop, int myproc)
{
  REAL utmp,ratio, em, alpha, erosionmax, nettmp; 
  REAL taubtmp2, taubtmp1, epslon, taub, taub0, max, mean, Hmean,taub1;
  int i,j,k,l,base,ntotal,n; 
  // coefficient for hard erosion
  em=1.0;
  // coefficient for soft erosion
  alpha=4.2;  // 4.2 to 25.6 according to MIKE21
  
  // calculate taub for each sub cell
  for(j=0;j<grid->Nc;j++)
  {
    Hmean=subgrid->Veff[j]/subgrid->Aceff[j];
    utmp=sqrt(pow(phys->uc[j][grid->Nk[j]-1], 2)+pow(phys->vc[j][grid->Nk[j]-1], 2));
    //taub0=(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow((log(30*Hmean/sediments->kb[i])+30*sediments->kb[i]/Hmean-1)/0.41,-2)*utmp*utmp;  
    base=j*(grid->maxfaces-2)*subgrid->segN*subgrid->segN;
    ntotal=(grid->nfaces[j]-2)*subgrid->segN*subgrid->segN;
    if(subgrid->erosionpara)
      taub0=(1+phys->rho[j][grid->Nk[j]-1])*RHO0*subgrid->Cdmean[j]*utmp*utmp;

    if(Hmean<BUFFERHEIGHT || !phys->active[j])
      taub0=0;
 
    if(Hmean<1e-3)
      taub0=0;

    if(prop->culvertmodel)
      if(culvert->top[j]!=INFTY)
        taub0=0;

    //if(taub0>0.1)
      //printf("n %d nc %d dv %e v %e taub %e cd %e u %e\n",prop->n,j,grid->dv[j]+phys->h[j],subgrid->Veff[j],taub0,subgrid->Cdmean[j],utmp);
          
    max=0.0;
    mean=0.0;
    n=0;

    for(k=0;k<ntotal;k++)
    {
      taub=taub0*subgrid->Cdratio[base+k];
      if(prop->wavemodel)
      {
        if(wprop->FetchModel)
        {
          if(taub==0){
            taub=0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*
              2*PI/wave->Twsig[j]),2);   
          } else {
            taubtmp1=sqrt(pow((taub*phys->uc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*
              RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*cos(wave->Winddir[j])),2)+
              pow((taub*phys->vc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*
              pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*sin(wave->Winddir[j])),2));
            taubtmp2=sqrt(pow((taub*phys->uc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*
              RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*(-cos(wave->Winddir[j]))),2)+
              pow((taub*phys->vc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*
              pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*(-sin(wave->Winddir[j]))),2));
            taub=Max(taubtmp1,taubtmp2);
          }
        } else {
          taubtmp1=0.5*wave->fw[j]*RHO0*pow(wave->ub[j],2);
          taubtmp2=taub;
          taub = sqrt(pow(taubtmp2, 2)+pow(taubtmp1, 2)+ 2*taubtmp1*taubtmp2*cos(wave->Cr[j]));
        }
      } 
      if(prop->culvertmodel)
        if(culvert->top[j]!=INFTY)
          taub=0;
      subgrid->taubsub[base+k]=taub;  
      
      if(max<taub)
        max=taub;
      
      if(taub>0)
      {
        n++;
        mean+=taub;
      }
    }
    if(n==0)
    {
      sediments->Seditbmax[j]=0.0;
      sediments->Seditb[j]=0.0;      
    } else {
      sediments->Seditbmax[j]=max;
      sediments->Seditb[j]=taub0;
    }
  } 

  for(i=0;i<sediments->Nsize;i++)
  {
    for(j=0;j<grid->Nc;j++)
    {
      base=j*(grid->maxfaces-2)*subgrid->segN*subgrid->segN;
      ntotal=(grid->nfaces[j]-2)*subgrid->segN*subgrid->segN;
      erosionmax=0;
      // calculate erosion for each layer
      for(k=0;k<sediments->Nlayer;k++)
      {
        erosionmax=sediments->Layerthickness[i][j][k]*sediments->Drydensity[i][k]*subgrid->Asedi[j]/grid->Ac[j]/prop->dt;
        sediments->Erosion_old[i][j][k]=sediments->Erosion[i][j][k];
        sediments->Erosion[i][j][k]=0.0;
        n=0;
        // only calculate when the mean bottom shear stress >0
        if(sediments->Seditbmax[j]>sediments->Taue[i][k])
          for(l=0;l<ntotal;l++)
          {
            if(subgrid->taubsub[base+l]>sediments->Taue[i][k])
            {
              n++;
              if(sediments->Softhard[k]==0){
                // soft erosion
                sediments->Erosion[i][j][k]+=subgrid->A_sub[base+l]*sediments->E0[i][k]*
                  exp(alpha*(subgrid->taubsub[base+l]-sediments->Taue[i][k]));
              } else {
                // hard erosion
                ratio=subgrid->taubsub[base+l]/sediments->Taue[i][k]-1;
                sediments->Erosion[i][j][k]+=subgrid->A_sub[base+l]*
                  sediments->E0[i][k]*pow(ratio,em);
              }               
            }
          }

        // check the total mass for each layer
        if(n==0)
          sediments->Erosion[i][j][k]=0.0;
        else
          sediments->Erosion[i][j][k]/=subgrid->Asedi[j];

        if(sediments->Erosion[i][j][k]>erosionmax)
          sediments->Erosion[i][j][k]=erosionmax; 

        if(sediments->bedComplex==1){
          nettmp=0;
          //assume last layer cannot be flushed away
          for(k=sediments->Nlayer-2;k>0;k--){
            nettmp=sediments->Layerthickness[i][j][k]*sediments->Drydensity[i][k]/prop->dt-sediments->Erosion[i][j][k]+sediments->Erosion[i][j][k+1]+sediments->Consolid[k-1]-sediments->Consolid[k]; //here can use theta method.
            if(nettmp<0)
              sediments->Erosion[i][j][k]+=nettmp;
          }

          // calculate erosion limit for the first layer
          nettmp=0;
          if(sediments->Nlayer>1)
            nettmp=sediments->Layerthickness[i][j][0]*sediments->Drydensity[i][0]/prop->dt-sediments->Erosion[i][j][0]+sediments->Erosion[i][j][1]+sediments->Deposition[i][j]-sediments->Consolid[0];
          if(nettmp<0)   
            sediments->Erosion[i][j][0]+=nettmp; 
        }  
      }
    }
  }
}

/*
 * Function: CalculateSubgridBottomDragCoef
 * Usage: Calculate the drag coefficient for all subedge and
 * then modify the original phys->CdB
 * ----------------------------------------------------
 * Can only be used for grid->Nkmax=1 and subgrid=1
 */
void CalculateSubgridDragCoef(gridT *grid, physT *phys, propT *prop)
{
   int ne,i,j,k,base,base1,nc1,nc2,n;
   REAL h, Cd, d1,d2,dtmp, Htmp, Cdbtmp, CdVtmp, Cdtmp, p, H_mean,Cd1,p1,p2,ac_sub,df_sub;
   
   for(ne=0;ne<grid->Ne;ne++)
   {
     p1=0;
     p2=0;
     nc1 = grid->grad[2*ne];
     nc2 = grid->grad[2*ne+1];
     if(nc1==-1)
       nc1=nc2;
     if(nc2==-1)
       nc2=nc1;
     n=0;
     h=subgrid->he[ne];
     df_sub=grid->df[ne]/subgrid->segN;
     base=ne*subgrid->segN;
     base1=ne*(subgrid->segN+1);
     //H_mean=0.5*(grid->dzz[nc1][grid->Nke[ne]-1]+grid->dzz[nc2][grid->Nke[ne]-1]);
     Cd=0;
     if(subgrid->dzboteff[ne]>=BUFFERHEIGHT && grid->etop[ne]==grid->Nke[ne]-1 && grid->mark[ne]!=2) //ignore inflow boundary
     {
       for(i=0;i<subgrid->segN;i++)
       {
         // store small depth at d1 and large at d2
         d1=subgrid->dpe[base1+i];
         d2=subgrid->dpe[base1+i+1];
         p=1;
         Cdbtmp=1;
         CdVtmp=0;

         if(d1>d2)
         { 
           dtmp=d2;
           d2=d1;
           d1=dtmp;
         }
         if(h<(-d2+BUFFERHEIGHT))
         {
           // dry subedge Htmp=0 so no effects on total Cdb
           Htmp=0;
         } else{
           n++;
           if(h>-d1)
           {
             Htmp=subgrid->de[base+i]+h;
             ac_sub=Htmp*df_sub;
           }
           else
           {
             Htmp=0.5*(h+d2);//*(h+d2)/(d2-d1);
             ac_sub=Htmp*(h+d2)/(d2-d1)*df_sub;
           }
           
           //when no bottom roughness specified use CdB directly
           if(prop->z0B==0)
             Cdbtmp=prop->CdB;
           else{
             Cdbtmp=(log(Htmp/phys->z0B[ne])+phys->z0B[ne]/Htmp-1)/0.41;
             Cdbtmp=1/Cdbtmp/Cdbtmp;
           }
           
           if(prop->marshmodel)
           {
             p=subgrid->hmarshe[base+i]/Htmp;
             if(p>1)
               p=1;
             p1+=Min(subgrid->hmarshe[base+i],Htmp)/Htmp*ac_sub;
             p2+=ac_sub;
             CdVtmp=marsh->Alphav*marsh->Rv*subgrid->cdve[base+i]*marsh->Na*p*p*Htmp\
               /(1-PI*pow(marsh->Rv,2)*marsh->Na*p)/pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
           }

         }
         // subgrid parameterization for total drag at each edge
         Cd+=sqrt(Htmp/(Cdbtmp+CdVtmp))*ac_sub;
       }

       Cd=pow(Cd,2);
       Cd=(grid->dzf[ne][0]*grid->df[ne])*(grid->dzf[ne][0]*grid->df[ne])*subgrid->dzboteff[ne]/Cd;
       Cdbtmp=(log(subgrid->dzboteff[ne]/phys->z0B[ne])+phys->z0B[ne]/subgrid->dzboteff[ne]-1)/0.41;
       Cdbtmp=1/Cdbtmp/Cdbtmp;
       if(prop->marshmodel)
       {
         p=p1/p2;
         CdVtmp=marsh->Alphav*marsh->Rv*1*marsh->Na*p*p*subgrid->dzboteff[ne]\
              /(1-PI*pow(marsh->Rv,2)*marsh->Na*p)/pow((1-2*marsh->Rv*sqrt(marsh->Na)),2); 
         Cd1=CdVtmp+Cdbtmp;
       } else {
         Cd1=Cdbtmp;
       }

       if(subgrid->dragpara==1)
         phys->CdB[ne]=Cd;
       else
         phys->CdB[ne]=Cd1;
     }
   }
}

/*
 * Function: CalculateSubgridCdmean
 * Usage: Calculate the total drag coefficient calculated
 * by mean depth and vegetation height
 * ----------------------------------------------------
 * Can only be used for erosion calculation when subgrid=1
 * This function is replaced by the new method to calculate
 * the cell center average drag forcing
 */
void CalculateSubgridCdmean(gridT *grid, physT *phys, propT *prop)
{
  int i, j, ntotal, base,n;
  REAL h,hleft,Htmp, Hmean, p, Cdb, Cdv,Cdsum,tmp,A,R,d1,d2,d3,d,Cdmax,Cdmin, Vc,Ac,ac1,r;//maxratio;
  for(i=0;i<grid->Nc;i++)
  {
    h=phys->h[i];
    Cdmax=0;
    ntotal=(grid->nfaces[i]-2)*subgrid->segN*subgrid->segN;
    base=i*(grid->maxfaces-2)*subgrid->segN*subgrid->segN;
    Vc=0;
    Ac=0;
    //Hmean=subgrid->Veff[i]/subgrid->Aceff[i];
    tmp=0;
    subgrid->Cdmean[i]=0;

    // calculate Cdratio CdB/Cdtotal*h(i,j)/Hmean
    for(j=0;j<ntotal;j++)
    { 
      d1=subgrid->d1[base+j];
      d2=subgrid->d2[base+j];
      d3=subgrid->d3[base+j];
      
      // calculate H_mean for each subcell
      if(h<=-d2)
        subgrid->H_sub[base+j]=(h+d3)/3;
      else if(h>-d2 && h<=-d1)
        subgrid->H_sub[base+j]=((d3-d1)*(d2-d1)*(d3+h)-(d2+h)*(d2+h)*(d2+h))/((d3-d1)*(d2-d1)-(d2+h)*(d2+h))/3;
      else
        subgrid->H_sub[base+j]=(d1+h)+(d2+d3-2*d1)/3;

      
      // calculate effective wet area for each subcell
      ac1=subgrid->Acbackup[i]/subgrid->segN/subgrid->segN*subgrid->Acratio[i*(grid->maxfaces-2)+j/(subgrid->segN*subgrid->segN)];
      
      // for different condition use different method to calculate ac
      if(d1==d2 && d2==d3 )//|| fabs(d3-d1)<=0.01)
      {
        if(h<-d3)
          subgrid->A_sub[base+j]=0;
        else
          subgrid->A_sub[base+j]=ac1;      
      } else if(d1==d2){
        if(h>=-d1)
          subgrid->A_sub[base+j]=ac1;
        else if(h<-d1 && h>=-d3){
          r=pow(((h+d3)/(d3-d1)),2);
          subgrid->A_sub[base+j]=ac1*r;      
        } else {
          subgrid->A_sub[base+j]=0;
        }
      } else if(d2==d3){
        if(h>=-d1)
          subgrid->A_sub[base+j]=ac1;
        else if(h<-d1 && h>=-d2){
          r=pow((-(h+d1)/(d3-d1)),2);
          subgrid->A_sub[base+j]=ac1*(1-r);      
        } else {
          subgrid->A_sub[base+j]=0;
        }
      } else {
        if(h>=-d1)
          subgrid->A_sub[base+j]=ac1;
        else if(h<-d1 && h>=-d2)
        {
          r=pow((-(h+d1)/(d2-d1)),2);
          subgrid->A_sub[base+j]=ac1*(d3-d2)/(d3-d1)+ac1*(d2-d1)/(d3-d1)*(1-r);
        } else if(h<-d2 && h>=-d3) {
          r=pow(((h+d3)/(d3-d2)),2);
          subgrid->A_sub[base+j]=ac1*(d3-d2)/(d3-d1)*r;
        } else { 
          subgrid->A_sub[base+j]=0;
        }
      }
      // calculate the total volume and wet area
      Ac+=subgrid->A_sub[base+j];
      Vc+=subgrid->A_sub[base+j]*subgrid->H_sub[base+j];

      subgrid->Cdratio[base+j]=0.0;
      //Htmp=phys->h[i]+subgrid->dc[base+j];
       
      Cdsum=0.0;
      R=0.0;
      if(subgrid->H_sub[base+j]>BUFFERHEIGHT) 
      {
        R=1;
        if(prop->z0B==0)
          Cdb=prop->CdB;
        else{
          A=(log(subgrid->H_sub[base+j]/prop->z0B)+prop->z0B/subgrid->H_sub[base+j]-1)/0.41;
          Cdb=1/A/A; // should add more for generalize case       
        }
        Cdsum=Cdb;
        
        if(prop->marshmodel)
        {
          if(subgrid->hmarshc[base+j]>0)
          {
            if(subgrid->hmarshc[base+j]<subgrid->H_sub[base+j])
              hleft=subgrid->hmarshc[base+j]; 
            else
              hleft=0;
            p=hleft/subgrid->H_sub[base+j];
            if(p>1)
              p=1;
            Cdv=marsh->Alphav*marsh->Rv*subgrid->cdvc[base+j]*marsh->Na*p*p*subgrid->H_sub[base+j]\
                 /(1-PI*pow(marsh->Rv,2)*marsh->Na*p)/pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
            Cdsum+=Cdv;
            A=Cdb/Cdsum;
            R*=A;
          }
        }
        subgrid->Cdratio[base+j]=R;
      }

      /*if(j==0){
        Cdmax=Cdsum;
        Cdmin=Cdsum;
      }

      if(Cdsum>Cdmax)
        Cdmax=Cdsum;
      if(Cdmin>Cdsum)
        Cdmin=Cdsum;*/

      if(Cdsum>0)
        tmp+=sqrt(subgrid->H_sub[base+j]/Cdsum)*subgrid->H_sub[base+j]*subgrid->A_sub[base+j];
    }
    
    for(j=0;j<ntotal;j++)
    { 
      if(subgrid->H_sub[base+j]>BUFFERHEIGHT)
       subgrid->Cdratio[base+j]*=subgrid->H_sub[base+j]/Vc*Ac;
    }
    subgrid->Asedi[i]=Ac;
    subgrid->Vsedi[i]=Vc;
    if(tmp>0)
      subgrid->Cdmean[i]=Vc*Vc*Vc/Ac/tmp/tmp;
    
    /*if(subgrid->Cdmean[i]<Cdmin)
      subgrid->Cdmean[i]=Cdmin;
    if(subgrid->Cdmean[i]>Cdmax)
      subgrid->Cdmean[i]=Cdmax;*/
  }
}

/*
 * Function: FreeSubgrid
 * Usage: free space for all the variables
 * ----------------------------------------------------
 * Basic sunfree function
 *
 */
void FreeSubgrid(gridT *grid, int myproc) {
  int i,j;
  free(subgrid->xp);
  free(subgrid->yp);
  free(subgrid->dp);
  free(subgrid->dc);
  free(subgrid->varNc);
  free(subgrid->rhs);
  free(subgrid->residual);
  free(subgrid->d1);
  free(subgrid->H_sub);
  free(subgrid->A_sub);
  free(subgrid->d2);
  free(subgrid->d3);
  free(subgrid->hmarshpe);
  free(subgrid->hmarshe);
  free(subgrid->de);
  free(subgrid->dpe);
  free(subgrid->hmarshp);
  free(subgrid->hmarshc);
  free(subgrid->he);
  free(subgrid->hprof);
  free(subgrid->hprofe);
  free(subgrid->hmin);
  free(subgrid->hmax);
  free(subgrid->Vprof);
  free(subgrid->Acprof);
  free(subgrid->Wetperiprof);
  free(subgrid->Acbackup);
  free(subgrid->Aceff);
  //free(subgrid->Cdmean);
  free(subgrid->Cdratio);
  free(subgrid->taubsub);
  //free(subgrid->delta);
  free(subgrid->Aceffold);
  free(subgrid->Heff);
  free(subgrid->Heffold);
  free(subgrid->Acwet);
  free(subgrid->Veff);
  //free(subgrid->Vsedi);
  //free(subgrid->Asedi);
  free(subgrid->Veffold);
  free(subgrid->dzboteff);
  free(subgrid->fluxhprof);
  free(subgrid->cellp);
  for(i=0;i<grid->Nc;i++)
  {
    free(subgrid->Acceff[i]);
    free(subgrid->fluxp[i]);
    free(subgrid->fluxn[i]);
    free(subgrid->Acceffold[i]);
    free(subgrid->Acveff[i]);
    free(subgrid->Acveffold[i]);
    free(subgrid->Acveffold2[i]);
  } 
  free(subgrid->fluxn);
  free(subgrid->fluxp);  
  free(subgrid->Acceffold);
  free(subgrid->Acceff);
  free(subgrid->Acveffold);
  free(subgrid->Acveffold2);
  free(subgrid->Acveff);
}


/*
 * Function: SubgridCellAverageTau
 * Usage: Finish the calculation for cell-center averaged 
 * drag forcing. In this function, only windstress and phys->uc_new/dt
 * are added. For coriolis, baroclinic, momentum advection and horizontal
 * diffusion is in Function HorizontalSource, For barotropic part is in 
 * Function Upredictor
 * ----------------------------------------------------
 * Function is turned on when subgrid sediment transport 
 * parameterization is turned on and grid->Nkmax==1
 */
/*void SubgridCellAverageTau(gridT *grid,physT *phys,propT *prop, metT *met, int myproc, int numprocs)
{
  int iptr, i;
  REAL tmp;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) 
  {
    i=grid->cellp[iptr];
    subgrid->tau_c_x[i]-=phys->uc[i][0]/prop->dt;
    subgrid->tau_c_y[i]-=phys->vc[i][0]/prop->dt;

    // wind stress works with met model 
    // need a new function for regular wind stress
    if(prop->metmodel){
      subgrid->tau_c_x[i]+=met->tau_x[i]*subgrid->Aceff[i]/subgrid->Veff[i];
      subgrid->tau_c_y[i]+=met->tau_y[i]*subgrid->Aceff[i]/subgrid->Veff[i];
    } else {
      tmp=ReturnCellWindStress(grid->xv[i], grid->yv[i], prop->n, prop->dt, 1);
      subgrid->tau_c_x[i]+=tmp*subgrid->Aceff[i]/subgrid->Veff[i];
      tmp=ReturnCellWindStress(grid->xv[i], grid->yv[i], prop->n, prop->dt, 0);
      subgrid->tau_c_y[i]+=tmp*subgrid->Aceff[i]/subgrid->Veff[i];
    }

    // before we calculate tau/H
    // here multifly back H_mean
    subgrid->tau_c_x[i]*=subgrid->Veff[i]/subgrid->Aceff[i];
    subgrid->tau_c_y[i]*=subgrid->Veff[i]/subgrid->Aceff[i];
    tmp=sqrt(pow(phys->uc[i][0],2)+pow(phys->vc[i][0],2));
    if(grid->dzz[i][0]<=BUFFERHEIGHT || tmp<0.01)
    {
      subgrid->tau_c_x[i]=0;
      subgrid->tau_c_y[i]=0;
    }
  }
}*/

/*
 * Function: SubgridFluxCheck
 * Usage: Calculate the inward and outward flux for each cell at each layer
 * to check the wetting drying condition for scalar transport to ensure the bounded
 * concentration 
 * ----------------------------------------------------
 * Function is turned on when subgrid sediment transport 
 * parameterization is turned on and grid->Nkmax==1
 */
void SubgridFluxCheck(gridT *grid, physT *phys, propT *prop,int myproc){
  int i,iptr,ne,ktop,nf,k;
  REAL sum0,normal,sum,fac1,fac2,fac3,flux;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    // get the effective top layer number
    if(grid->ctop[i]>=grid->ctopold[i]) {
      ktop=grid->ctop[i];
    } else {
      ktop=grid->ctopold[i];    
    }    
    
    // initialize flux variable
    for(k=0;k<grid->Nk[i];k++){
      subgrid->fluxn[i][k]=0;
      subgrid->fluxp[i][k]=0;
    }
    sum0=0;
    sum=0;

    // calculate horizontal flux
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      for(k=grid->etop[ne];k<grid->Nke[ne];k++){ 
        flux=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+
          fac3*phys->u_old2[ne][k])*grid->df[ne]*normal*grid->dzf[ne][k];
        if(flux>0){
          // outflux for each cell layer
          subgrid->fluxp[i][k]+=fabs(flux);
          // total flux
          sum0+=fabs(flux);
        }else{
          // influx for each cell layer
          subgrid->fluxn[i][k]+=fabs(flux);
          // total flux
          sum+=fabs(flux);
        }
      }
    }

    // if outward flux is bigger than volume, scalar transport is stopped to ensure bounded scalar concentration
    //if(fabs(sum0)>subgrid->Veff[i])
      //phys->active[i]=0;
    if(fabs(sum)>subgrid->Veff[i] || fabs(sum0)>subgrid->Veffold[i])
      phys->active[i]=0;
    
    // add back vertical flux for each layer of each cell
    if(ktop!=grid->Nk[i]-1)
    {
      flux=prop->dt*(fac1*phys->wnew[i][ktop+1]*subgrid->Acveff[i][ktop+1]+fac2*phys->w_old[i][ktop+1]*subgrid->Acveffold[i][ktop+1]+
        fac3*phys->w_old2[i][ktop+1]*subgrid->Acveffold2[i][ktop+1]);

      if(flux>0)
        subgrid->fluxp[i][ktop]+=fabs(flux);
      else
        subgrid->fluxn[i][ktop]+=fabs(flux);               
        
      for(k=ktop+1;k<grid->Nk[i];k++)
      {
         // bottom layer
        flux=prop->dt*(fac1*phys->wnew[i][k+1]*subgrid->Acveff[i][k+1]+
            fac2*phys->w_old[i][k+1]*subgrid->Acveffold[i][k+1]+
            fac3*phys->w_old2[i][k+1]*subgrid->Acveffold2[i][k+1]);

        if(flux<0)
          subgrid->fluxp[i][k]+=fabs(flux);
        else
          subgrid->fluxn[i][k]+=fabs(flux);                
          
        // top layer
        flux=prop->dt*(fac1*phys->wnew[i][k]*subgrid->Acveff[i][k]+
          fac2*phys->w_old[i][k]*subgrid->Acveffold[i][k]+
          fac3*phys->w_old2[i][k]*subgrid->Acveffold2[i][k]);
        if(flux>0)
          subgrid->fluxp[i][k]+=fabs(flux);
          else
          subgrid->fluxn[i][k]+=fabs(flux);               
      }
    }
  }
}


