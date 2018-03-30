/*
 * File: sediment.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file contains physically-based functions for cohesive sediment calculation.
 * Cohesive suspended sediments calculation is done in updatescalars when prop->sedi==1
 * Here sediment settling velocity, erosion rate and consolidation in different layers are calculated
 * 
 *
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
#include "culvert.h"
#include "timer.h"
 
void ReadSediProperties(int myproc);
void InitializeSediment(gridT *grid, physT *phys, propT *prop,int myproc);
void AllocateSediment(gridT *grid, int myproc);
void FreeSediment(gridT *grid, int myproc);
void BedChange(gridT *grid, physT *phys, propT *prop, int myproc);
void SedimentSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop,int Nosize, REAL theta);
void SedimentVerticalVelocity(gridT *grid, physT *phys,int Nosize,int symbol, int myproc);
void CalculateSedimentKb(gridT *grid, physT *phys, int myproc);
void RouseCurveAlpha(gridT *grid, physT *phys, propT *prop, int myproc);
void OpenSediFiles(propT *prop, int myproc);
void OutputSediment(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
void CalculateSediDiffusivity(gridT *grid, physT *phys,int Nosize,int myproc); 
void ISendRecvSediBedData3D(REAL **celldata, gridT *grid, int nlayer, int myproc, MPI_Comm comm);
void SettlingVelocity(gridT *grid, physT *phys, propT *prop, int myproc);
void CalculateErosion(gridT *grid, physT *phys, propT *prop,  int myproc);
void CalculateDeposition(gridT *grid, physT *phys, propT *prop, int myproc); //used by calculateerosion and Bedchange

/*
 * Function: ReadSediProperties
 * Usage: ReadSediProperties(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * Based on sedi.dat, load in the important parameters for 
 * cohesive sediment transport model
 *
 */
void ReadSediProperties(int myproc) {    
  int m,n;
  char str[BUFFERLENGTH];
  sediments->Nlayer = MPI_GetValue(DATAFILE,"Nlayer","ReadSediProperties",myproc);
  sediments->Nsize = MPI_GetValue(DATAFILE,"Nsize","ReadSediProperties",myproc);   
  sediments->WSconstant = MPI_GetValue(DATAFILE,"WSconstant","ReadSediProperties",myproc);
  sediments->sscvprof = MPI_GetValue(DATAFILE,"SSCvprof","ReadSediProperties",myproc);
  sediments->readSediment= MPI_GetValue(DATAFILE,"readSediment","ReadSediProperties",myproc);
  sediments->bedInterval= MPI_GetValue(DATAFILE,"bedInterval","ReadSediProperties",myproc);
  sediments->bedComplex= MPI_GetValue(DATAFILE,"bedComplex","ReadSediProperties",myproc);
  sediments->ParabolKappa= MPI_GetValue(DATAFILE,"ParabolKappa","ReadSediProperties",myproc);
  sediments->layerprop=MPI_GetValue(DATAFILE,"sedilayerprop","ReadSediProperties",myproc);
  if(sediments->bedComplex==1){
    printf("because bedComplex==1, so set bedInterval=1 automatically");
    sediments->bedInterval=1;
  }
  
  // condition check
  if(sediments->readSediment==1 && sediments->Nsize>3){
    printf("Nsize = %d>1, but readSediment==1 which means Nsize==1. You should set readSediment as 0 or Nsize==1.\n",sediments->Nsize);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  //allocate spaces for all sediments properties
  sediments->Ds=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Ws0=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Gsedi=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Anglerepos=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Prt=(REAL *)SunMalloc(sediments->Nsize*sizeof(REAL), "ReadSedimentProperties");
  sediments->Consolid=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  // fraction specified
  sediments->E0=(REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "ReadSedimentProperties");
  sediments->Taue=(REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "ReadSedimentProperties");
  sediments->Taud=(REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "ReadSedimentProperties");
  sediments->Drydensity=(REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "ReadSedimentProperties");
  for(m=0;m<sediments->Nsize;m++){
    sediments->E0[m]=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
    sediments->Taue[m]=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
    sediments->Taud[m]=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
    sediments->Drydensity[m]=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  }
  sediments->Thickness=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  sediments->Softhard=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "ReadSedimentProperties");
  for(m=1;m<=sediments->Nsize;m++) {
    sprintf(str,"Ds%d",m);
    sediments->Ds[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Ws0%d",m);
    if(sediments->WSconstant==1)
      sediments->Ws0[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    else
      sediments->Ws0[m-1]=0.0;
    sprintf(str,"Angle%d",m);
    sediments->Anglerepos[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc); 
    sprintf(str,"Gsedi%d",m);
    sediments->Gsedi[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc); 
    sprintf(str,"Prt%d",m);
    sediments->Prt[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
  }
  for(m=1;m<=sediments->Nlayer;m++) {
    sprintf(str,"Consolid%d",m);
    sediments->Consolid[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"Thickness%d",m);
    sediments->Thickness[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
    sprintf(str,"softhard%d",m);
    sediments->Softhard[m-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
  }
  if(!sediments->layerprop)
    for(m=1;m<=sediments->Nsize;m++){
      for(n=1;n<=sediments->Nlayer;n++){
        sprintf(str,"E%d%d",n,m);
        sediments->E0[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
        sprintf(str,"Taue%d%d",n,m);
        sediments->Taue[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
        sprintf(str,"Taud%d%d",n,m);
        sediments->Taud[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
        sprintf(str,"Drydensity%d%d",n,m);
        sediments->Drydensity[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
      }
    }
  else
    for(m=1;m<=sediments->Nsize;m++){
      for(n=1;n<=sediments->Nlayer;n++){
        sprintf(str,"E0%d",n);
        sediments->E0[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
        sprintf(str,"Taue%d",n);
        sediments->Taue[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
        sprintf(str,"Taud%d",n);
        sediments->Taud[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
        sprintf(str,"Drydensity%d",n);
        sediments->Drydensity[m-1][n-1]=MPI_GetValue(DATAFILE,str,"ReadSediProperties",myproc);
      }
    }    
  sediments->Kbed=MPI_GetValue(DATAFILE,"kbed","ReadSediProperties",myproc);
  sediments->Chind=MPI_GetValue(DATAFILE,"Chind","ReadSediProperties",myproc);
  sediments->Cfloc=MPI_GetValue(DATAFILE,"Cfloc","ReadSediProperties",myproc);
}

/*
 * Function: AllocateSediment
 * Usage: allocate space for sediment variables
 * ----------------------------------------------------
 * Based on the value from ReadSedimentProperties
 *
 */
void AllocateSediment(gridT *grid, int myproc) { 
  int i,j,jptr,k;

  // add bottom roughness kb z0r z0s z0b
  sediments->kb = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  sediments->z0b = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  sediments->z0r = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  sediments->z0s = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");

  sediments->SediC = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->SediC_old = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Erosion = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Erosion_old = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");   
  sediments->boundary_sediC = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Ws = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Deposition = (REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  sediments->alphaSSC = (REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  sediments->Deposition_old = (REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  sediments->Toplayerratio = (REAL **)SunMalloc(sediments->Nsize*sizeof(REAL *), "AllocateSediVariables");
  sediments->Layerthickness = (REAL ***)SunMalloc(sediments->Nsize*sizeof(REAL **), "AllocateSediVariables");
  sediments->Seditb = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  sediments->Totalthickness = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  sediments->Reposangle = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  sediments->Layertop = (int *)SunMalloc(grid->Nc*sizeof(int), "AllocateSediVariables");
  sediments->Seditbmax = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
  for(i=0;i<sediments->Nsize;i++){
    sediments->SediC[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    sediments->SediC_old[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    sediments->Deposition[i]=(REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    sediments->alphaSSC[i]=(REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    sediments->Deposition_old[i]=(REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    sediments->Toplayerratio[i]=(REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateSediVariables");
    // boundary part  
    sediments->boundary_sediC[i] = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocateSediment");
    sediments->Ws[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    sediments->Erosion[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    sediments->Erosion_old[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    sediments->Layerthickness[i] = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
    for(j=0;j<grid->Nc;j++){
      sediments->SediC[i][j] = (REAL *)SunMalloc(grid->Nk[j]*sizeof(REAL), "AllocateSediVariables");
      sediments->SediC_old[i][j] = (REAL *)SunMalloc(grid->Nk[j]*sizeof(REAL), "AllocateSediVariables");
      sediments->Erosion[i][j] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables"); 
      sediments->Erosion_old[i][j] = (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables"); 
      sediments->Layerthickness[i][j]= (REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
      // allocate boundary value
      for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
        k=grid->edgep[jptr];  
        sediments->boundary_sediC[i][jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[k]*sizeof(REAL), "AllocateSediVariables");
      }
      sediments->Ws[i][j] = (REAL *)SunMalloc((grid->Nk[j]+1)*sizeof(REAL), "AllocateSediVariables"); 
    }
  }
  sediments->Thicknesslayer = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->SediKappa_tv = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  sediments->Wnewsedi= (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "AllocateSediVariables");
  for(i=0;i<grid->Nc;i++){
    sediments->Thicknesslayer[i]=(REAL *)SunMalloc(sediments->Nlayer*sizeof(REAL), "AllocateSediVariables");
    sediments->Wnewsedi[i]= (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL), "AllocateSediVariables");
    sediments->SediKappa_tv[i]= (REAL *)SunMalloc((grid->Nk[i])*sizeof(REAL), "AllocateSediVariables");
  }
}

/*
 * Function: InitializeSediment
 * Usage: give initial value for sediment variables
 * ----------------------------------------------------
 * Based on the value from ReadSedimentProperties
 * here assume Erosion and Deposition initial value is 0, may have error
 */
void InitializeSediment(gridT *grid, physT *phys, propT *prop,  int myproc) { 
  int i,j,k,jptr,ne;
  REAL *stmp,z,alphassc;
  FILE *InitSedimentFID;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  // give 0 for SediC boundary_sediC Ws
  for(i=0;i<sediments->Nsize;i++) {
    for(j=0;j<grid->Nc;j++){
      for(k=0;k<grid->Nk[j];k++) {
        sediments->SediC[i][j][k]=0;   
        sediments->SediC_old[i][j][k]=0;             
        sediments->Ws[i][j][k]=0;
      }
      sediments->Ws[i][j][grid->Nk[j]]=0;
    }

    for(j=0;j<(grid->edgedist[5]-grid->edgedist[2]);j++)
    {
      ne=grid->edgep[j+grid->edgedist[2]]; 
      for(k=0;k<grid->Nke[ne];k++)
        sediments->boundary_sediC[i][j][k]=0; // boundary_sediC will be set using BoundaryScalars function  
    }
  }
  
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i]+1;k++){
      sediments->Wnewsedi[i][k]=0;
      if(k!=grid->Nk[i])
        sediments->SediKappa_tv[i][k]=0;
    }
  }

  // SSC vertical profile option
  if(grid->Nkmax==1)
  {
    if(sediments->sscvprof==2)
    {
      alphassc=MPI_GetValue(DATAFILE,"alphaSSC","InitializeSediment",myproc);
      for(i=0;i<sediments->Nsize;i++)
        for(j=0;j<grid->Nc;j++)
           sediments->alphaSSC[i][j]=alphassc;
    } else {
      RouseCurveAlpha(grid, phys,prop, myproc);
    }
  } else {
    for(i=0;i<sediments->Nsize;i++)
      for(j=0;j<grid->Nc;j++)
        sediments->alphaSSC[i][j]=1;
  }

  // give 0 for deposition erosion layerthickness
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++){
      sediments->Deposition[i][j]=0;
        
      sediments->Toplayerratio[i][j]=0;          
      sediments->Deposition_old[i][j]=0;
      for(k=0;k<sediments->Nlayer;k++){
        sediments->Erosion[i][j][k]=0;
        sediments->Erosion_old[i][j][k]=0;
        sediments->Layerthickness[i][j][k]=0;
      }
    }
  
  // set first tb
  for(i=0;i<grid->Nc;i++) {
    sediments->Seditb[i]=0;
    sediments->Seditbmax[i]=0;
    sediments->Layertop[i]=0;
    sediments->Totalthickness[i]=0;
    sediments->Reposangle[i]=0;
    for(j=0;j<sediments->Nlayer;j++)
      sediments->Thicknesslayer[i][j]=0;
  }

  // read or use function to set initial condition for SediC
  if(sediments->readSediment) { //it means Nsize==1
    MPI_GetFile(filename,DATAFILE,"InitSedimentFile","InitializeSediment",myproc);
    InitSedimentFID = MPI_FOpen(filename,"r","InitializeSediment",myproc);
    stmp = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializeSediment");
    if(fread(stmp,sizeof(REAL),grid->Nkmax,InitSedimentFID) != grid->Nkmax)
      printf("Error reading stmp second\n");
    fclose(InitSedimentFID);    
    
    for(i=0;i<grid->Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        sediments->SediC[0][i][k]=stmp[k];
    
    SunFree(stmp,grid->Nkmax,"InitializeSediment");
  } else {
    for(j=0;j<sediments->Nsize;j++){
      for(i=0;i<grid->Nc;i++) {
        z = 0;
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          z-=grid->dz[k]/2;
          sediments->SediC[j][i][k]=ReturnSediment(grid->xv[i],grid->yv[i],z,j); // add new functions in Initialization.c
          z-=grid->dz[k]/2;
        }
      }
    }
  }

  // set SediC_old value assume sediC^n-1=sediC^n for prop->n==1
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      for(k=0;k<grid->Nk[j];k++)
        sediments->SediC_old[i][j][k]=sediments->SediC[i][j][k];
  
  // give the initial value for thickness
  for(j=0;j<sediments->Nsize;j++)
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<sediments->Nlayer;k++){
        sediments->Layerthickness[j][i][k]=sediments->Thickness[k]*ReturnBedSedimentRatio(grid->xv[i],grid->yv[i],k,j,sediments->Nsize);
      }

  // calculate Thicknesslayer
  for(j=0;j<grid->Nc;j++)
    for(k=0;k<sediments->Nlayer;k++){
      sediments->Thicknesslayer[j][k]=sediments->Thickness[k];
      sediments->Totalthickness[j]+=sediments->Thickness[k];
    }

  //calculate Toplayerratio
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      sediments->Toplayerratio[i][j]=ReturnBedSedimentRatio(grid->xv[i],grid->yv[i],0,i,sediments->Nsize);

  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      sediments->Reposangle[j]+=sediments->Toplayerratio[i][j]*sediments->Anglerepos[i];


  // calculate initial settling velocity 
  if(sediments->WSconstant==1){
    for(i=0;i<sediments->Nsize;i++)
      for(j=0;j<grid->Nc;j++)
        for(k=0;k<grid->Nk[j]+1;k++)
          sediments->Ws[i][j][k]=sediments->Ws0[i];
  } else {
    SettlingVelocity(grid,phys,prop,myproc);
  }
  // initialize bottom roughness
  sediments->ds50=MPI_GetValue(DATAFILE,"Ds50","InitializeSediment",myproc);
  for(i=0;i<grid->Nc;i++){
    sediments->z0s[i] = 2.5*sediments->ds50/30;
    sediments->z0b[i] = 0;
    sediments->z0r[i] = 0;
    sediments->kb[i] = 30*sediments->z0s[i];
  }
    
  // calculate initial Deposition & deposition
  CalculateErosion(grid,phys,prop,myproc);
  CalculateDeposition(grid,phys,prop,myproc);
}

/*
 * Function: RouseCurveAlpha 
 * usage: calculate the ratio between bottom SSC and depth-averaged SSC
 * by assuming the existence of Rouse Curve (divide into 10 layers to calulcate)
 * integeral
 * ------------------------------------------------------
 * Use when Nkmax=1 and sscvprof==1
 *
 */
void RouseCurveAlpha(gridT *grid, physT *phys, propT *prop, int myproc){
  REAL dz,rn,ustar,integ,z_bot=2*sediments->ds50,H,z;
  int i,j,k,nk=10;
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
    { 
      H=phys->h[j]+grid->dv[j];
      dz=(H-z_bot)/nk;
      integ=0;
      if(prop->z0B==0)     
        ustar=sqrt(pow(phys->uc[j][grid->Nk[j]-1],2)+pow(phys->vc[j][grid->Nk[j]-1],2))*prop->CdB;
      else
        ustar=sqrt(pow(phys->uc[j][grid->Nk[j]-1],2)+pow(phys->vc[j][grid->Nk[j]-1],2))*pow((1/0.41*(log(H/prop->z0B)+prop->z0B/H-1)),-2);
      rn=0.41*ustar/sediments->Ws0[i];
      for (k = 0; k < nk; k++)
      {
        z=0.5*dz+(k-1)*dz;
        integ+=dz*pow((z_bot/(H-z_bot)*(H-z)/z),rn);
      }
      sediments->alphaSSC[i][j]=H/integ;
      if(H<DRYCELLHEIGHT)
        sediments->alphaSSC[i][j]=0;
    } 
}

/*
 * Function: FreeSediment
 * Usage: free space for all the variables
 * ----------------------------------------------------
 * Basic sunfree function
 *
 */
void FreeSediment(gridT *grid, int myproc) {
  int i,j;

  for(i=0;i<sediments->Nsize;i++){
    for(j=0;j<grid->Nc;j++){
      free(sediments->SediC[i][j]);
      free(sediments->Ws[i][j]);
      free(sediments->Layerthickness[i][j]);
      free(sediments->Erosion[i][j]);
      free(sediments->Erosion_old[i][j]);
    }
    free(sediments->SediC[i]);
    free(sediments->Ws[i]);
    free(sediments->Layerthickness[i]);
    free(sediments->Toplayerratio[i]);
    free(sediments->Erosion[i]);
    free(sediments->Erosion_old);
    free(sediments->Deposition[i]);
    free(sediments->alphaSSC[i]);
    free(sediments->Deposition_old[i]);
  }
  free(sediments->SediC);
  free(sediments->Erosion);
  free(sediments->Erosion_old);
  free(sediments->Ws);
  free(sediments->Deposition);
  free(sediments->alphaSSC);
  free(sediments->Deposition_old);
  for(i=0;i<grid->Nc;i++){
    free(sediments->Thicknesslayer[i]);
    free(sediments->Wnewsedi[i]);
    free(sediments->SediKappa_tv[i]);
  }
  free(sediments->Thicknesslayer);
  free(sediments->SediKappa_tv);
  free(sediments->Wnewsedi);
  free(sediments->Layerthickness);
  free(sediments->Toplayerratio);
  free(sediments->Reposangle);
  free(sediments->Ds);
  free(sediments->Ws0);
  free(sediments->Gsedi);
  free(sediments->Prt);
  free(sediments->Consolid);
  for(i=0;i<sediments->Nsize;i++){
    free(sediments->E0[i]);
    free(sediments->Taue[i]);
    free(sediments->Taud[i]);
    free(sediments->Drydensity[i]);
  }
  free(sediments->E0);
  free(sediments->Taue);
  free(sediments->Taud);
  free(sediments->Drydensity);
  free(sediments->Thickness);
  free(sediments->Softhard);
  free(sediments->Seditb);
  free(sediments->Seditbmax);
  free(sediments->Layertop);
  free(sediments->Totalthickness);
}

/*
 * Function: SettlingVelocity
 * this function still need to be tested
 * Usage: calculate flocculation and hindered settling velocity
 * ----------------------------------------------------
 * Richarson and Zaki(1954) equation
 * only work when WSconstant=0
 *
 */
void SettlingVelocity(gridT *grid, physT *phys, propT *prop,  int myproc) { 
  int i,j,k;
  REAL tmp,sum,Cgel=1800000,wsn=2.0;

  // roughly estimate by Stokes Law
  for(i=0;i<sediments->Nsize;i++)
    sediments->Ws0[i] = (sediments->Gsedi[i]-1)*prop->grav*pow(sediments->Ds[i], 2)/(18*prop->nu);
  // coefficient for Richarson and Zaki(1954) equation
  //Cgel = 1800000;
  //wsn = 2.0;

  for(j=0;j<grid->Nc;j++){  
    for(k=0;k<grid->Nk[j]+1;k++){   
      sum=0;
      if(k==0){
        for(i=0;i<sediments->Nsize;i++)
          sum+=sediments->SediC[i][j][k];
      } else if(k==grid->Nk[j]){
        for(i=0;i<sediments->Nsize;i++)
          sum+=sediments->SediC[i][j][k-1];
      } else {
        for(i=0;i<sediments->Nsize;i++)
          sum+=0.5*(sediments->SediC[i][j][k-1]+sediments->SediC[i][j][k+1]);
      }
      //for(i=0;i<sediments->Nsize;i++)
      //sum+=0.5*(sediments->SediC[i][j][k-1]+sediments->SediC[i][j][k+1]);
      for(i=0;i<sediments->Nsize;i++){
        if(sum < sediments->Cfloc)
          sediments->Ws[i][j][k] = sediments->Ws0[i];
        else if(sum>= sediments->Cfloc && sum < sediments->Chind)
          sediments->Ws[i][j][k] =sediments-> Ws0[i]*sum/sediments->Gsedi[i]/1000/1000;  ///??????????
        else if(sum>= sediments->Chind && sum < Cgel){
          tmp = Min(1,sum/Cgel);
          tmp = 1-tmp;
          sediments->Ws[i][j][k] = sediments->Ws0[i]*sediments->Chind/sediments->Gsedi[i]/1000/1000*pow(tmp,wsn);
        } else {
          sediments->Ws[i][j][k] = 0;
          //sediments->SediC[i][j][k] = Cgel;
        }
      }
    }
    /*for(i=0;i<Nsize;i++){
      sediments->Ws[i][j][grid->Nk[j]] = sediments->Ws[i][j][grid->Nk[i]-1];
      sediments->Ws[i][j][0]=Ws[i][j][1];
      }*/
  }
}

/*
 * Function: CalculateErosion
 * Usage: calculate Erosion Rate for suspended sediment transport
 *        calculate Deposition for n time step
 * ----------------------------------------------------
 * for two types, hard and soft
 * coefficients based on Mike21
 */
void CalculateErosion(gridT *grid, physT *phys, propT *prop, int myproc) { 
  int i,j,k;
  REAL taub,utmp,taubtmp1,taubtmp2,em,alpha,ratio, depotmp, nettmp, erosionmax,ds90,dzh;
  // assume ks=3*ds90
  //ds90=MPI_GetValue(DATAFILE,"Ds90","CalculationErosion",myproc);
  
  // use sediments->kb 
  // coefficient for hard erosion
  em=1.0;
  // coefficient for soft erosion
  alpha=4.2;  // 4.2 to 25.6 according to MIKE21

  for(j=0;j<grid->Nc;j++){
    // calculate taub
    //taub = prop->CdB*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*(pow(phys->uc[j][grid->Nk[j]-1], 2)+pow(phys->vc[j][grid->Nk[j]-1], 2));
    utmp=sqrt(pow(phys->uc[j][grid->Nk[j]-1], 2)+pow(phys->vc[j][grid->Nk[j]-1], 2));

    // need to add more specification for 2d model
    if(!prop->subgrid)
      if(grid->Nk[j]>1)    
        //taub=pow(log(30*0.5*grid->dzz[j][grid->Nk[j]-1]/sediments->kb[j])/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
        taub=pow(log(0.5*grid->dzz[j][grid->Nk[j]-1]/phys->z0B[j])/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
      else
        taub=pow((log(grid->dzz[j][grid->Nk[j]-1]/phys->z0B[j])+prop->z0B/grid->dzz[j][grid->Nk[j]-1]-1)/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
    else{
      if(grid->Nk[j]>1){
        dzh=subgrid->Acceff[j][grid->Nk[j]-1]*grid->dzz[j][grid->Nk[j]-1]/subgrid->Acveff[j][grid->Nk[j]-1];
        //taub=pow(log(30*subgrid->Heff[j]/sediments->kb[j])/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
        taub=pow(log(0.5*dzh/phys->z0B[j])/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
      } else
        //taub=pow((log(30*subgrid->Heff[j]/sediments->kb[j])+sediments->kb[j]/subgrid->Heff[j]/30-1)/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
        taub=pow((log(subgrid->Heff[j]/phys->z0B[j])+phys->z0B[j]/subgrid->Heff[j]-1)/0.41,-2)*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*pow(utmp,2);
    }

    // add wave effects here is fetch right now
    // here may amplify the wave effects and 
    if(prop->wavemodel){
      if(wprop->FetchModel)
      {
        if(taub==0){
          taub=0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2);   
        } else {
          taubtmp1=sqrt(pow((taub*phys->uc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*cos(wave->Winddir[j])),2)+pow((taub*phys->vc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*sin(wave->Winddir[j])),2));
          taubtmp2=sqrt(pow((taub*phys->uc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*(-cos(wave->Winddir[j]))),2)+pow((taub*phys->vc[j][grid->Nk[j]-1]/utmp+0.5*(1+phys->rho[j][grid->Nk[j]-1])*RHO0*wave->Fw[j]*pow((wave->Waveexcur[j]*2*PI/wave->Twsig[j]),2)*(-sin(wave->Winddir[j]))),2));
          taub=Max(taubtmp1,taubtmp2);
        }
      } else {
        taubtmp1=0.5*wave->fw[j]*RHO0*pow(wave->ub[j],2);
        taubtmp2=taub;
        taub = sqrt(pow(taubtmp2, 2)+pow(taubtmp1, 2)+ 2*taubtmp1*taubtmp2*cos(wave->Cr[j]));
      }
    }
    
    if((phys->h[j]+grid->dv[j])<BUFFERHEIGHT)
      taub=0;

    // can be eliminated by the setup of initialize vertical grid
    if(!prop->subgrid)
      if(grid->dzz[j][grid->Nk[j]-1]<BUFFERHEIGHT)
        taub=0;

    if(prop->subgrid){
      if(subgrid->Heff[j]<BUFFERHEIGHT || subgrid->Heff[j]<1e-3)
        taub=0;
      //if(subgrid->dzboteff[j]<1e-3)
        //taub=0;
    }
    
    if(!phys->active[j])
      taub=0;

    // inside culvert no erosion
    if(prop->culvertmodel)
      if(culvert->top[j]!=INFTY)
        taub=0;

    sediments->Seditb[j]=taub;

    if(sediments->Seditbmax[j]<=taub)
      sediments->Seditbmax[j]=taub;
  }

  for(i=0;i<sediments->Nsize;i++){
    erosionmax=0;
    for(j=0;j<grid->Nc;j++){
      taub=sediments->Seditb[j];
      for(k=0;k<sediments->Nlayer;k++){
        erosionmax=sediments->Layerthickness[i][j][k]*sediments->Drydensity[i][k]/prop->dt;
        sediments->Erosion_old[i][j][k]=sediments->Erosion[i][j][k];
        sediments->Erosion[i][j][k]=0;
        if(taub>sediments->Taue[i][k]){
          if(sediments->Softhard[k]==0){
	          // soft erosion
            sediments->Erosion[i][j][k]=sediments->E0[i][k]*exp(alpha*(taub-sediments->Taue[i][k]));
          } else {
            // hard erosion
            ratio=Max(taub/sediments->Taue[i][k]-1,0);
            sediments->Erosion[i][j][k]=sediments->E0[i][k]*pow(ratio,em);
          }
        }
        if(sediments->Erosion[i][j][k]>erosionmax)
          sediments->Erosion[i][j][k]=erosionmax;
      }     

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
/*
 * Function: CalculateDeposition
 * Usage: calculate deposition
 * ----------------------------------------------------
 * based on settling velocity and SediConcentration
 */
void CalculateDeposition(gridT *grid, physT *phys, propT *prop, int myproc) {
  int i,j;
  REAL alpha;
  for(i=0;i<sediments->Nsize;i++){
    for(j=0;j<grid->Nc;j++){
      sediments->Deposition_old[i][j]=sediments->Deposition[i][j];
      // here deposition is explicit n time steps while Erosion is n+1 time steps
      sediments->Deposition[i][j]=sediments->alphaSSC[i][j]*sediments->SediC[i][j][grid->Nk[j]-1]*sediments->Ws[i][j][grid->Nk[j]];
      alpha=1.0;
      if(prop->subgrid)
      {
        if(subgrid->erosionpara)
          alpha=subgrid->delta[j];
        if(subgrid->Acceff[j][grid->Nk[j]-1]==0)
          alpha=0;
      }
      sediments->Deposition[i][j]*=alpha;
      if(!phys->active[j])
        sediments->Deposition[i][j]=0;
    }
  }
}

/*
 * Function: BedChange
 * Usage: calculate Cbed and Thickness change
 * ----------------------------------------------------
 * update by intervals, calculate Cbed fraction to 
 * get erosion for each fraction
 * assume dry density is constant for each layer
 * assume the last layer cannot be flushed away 
 * when bedComplex=1, Bedchange will be calculated every time step
 * there may some error in extreme condition
 */
void BedChange(gridT *grid, physT *phys, propT *prop, int myproc) {
 
  int i,j,k,nf,nc1,nc2;
  REAL grad,angle,bedflux,ratio,alpha;
  // update thickness
  // assume Drydensity for each fraction is constant for each layer 
  // the only thing change is the thickness(mass) for each fraction
  // the erosion will be divided into different parts according to the ration of thickness for each layer
  // the gravity effects only happen for the layertop part
  for(i=0;i<sediments->Nsize;i++){
    for(j=0;j<grid->Nc;j++){
      alpha=1.0;
      if(prop->subgrid)
        alpha=subgrid->Aceff[j]/grid->Ac[j];
      // no gravity effects
      if(sediments->Nlayer>1){
        // here can use theta method
        sediments->Layerthickness[i][j][0]+=sediments->bedInterval*prop->dt/sediments->Drydensity[i][0]*(alpha*(-0.5*sediments->Erosion[i][j][0]-0.5*sediments->Erosion_old[i][j][0]+0.5*sediments->Erosion[i][j][1]+0.5*sediments->Erosion_old[i][j][1]+0.5*sediments->Deposition[i][j]+0.5*sediments->Deposition_old[i][j])-sediments->Consolid[0]);
      } else {
        // if there is one layer, no consolidation and erosion[j][1];
        // because here we use deposition n+1 step, so there may be some error to make Layerthickness<0
        sediments->Layerthickness[i][j][0]+=0.5*sediments->bedInterval*prop->dt/sediments->Drydensity[i][0]*(alpha*(-sediments->Erosion[i][j][0]-sediments->Erosion_old[i][j][0]+sediments->Deposition[i][j]+sediments->Deposition_old[i][j])-sediments->Consolid[0]); 
      }

      // update other layers without gravity and deposition
      for(k=1;k<sediments->Nlayer-1;k++){
        sediments->Layerthickness[i][j][k]+=sediments->bedInterval*prop->dt/sediments->Drydensity[i][k]*(alpha*(-sediments->Erosion[i][j][k]+sediments->Erosion[i][j][k+1])+sediments->Consolid[k-1]-sediments->Consolid[k]);
        if(sediments->Layerthickness[i][j][k]<0)
          sediments->Layerthickness[i][j][k]=0;
      }
    }
  }
 
  // add gravity effects
  for(i=0;i<grid->Ne;i++){
    nc1=grid->grad[2*i];
    nc2=grid->grad[2*i+1];
    if(nc1!=-1 && nc2!=-1){
      grad=(sediments->Totalthickness[nc1]-sediments->Totalthickness[nc2])/grid->dg[i];
      if(grad>=0){
        ratio=(atan(grad)-sediments->Reposangle[nc1])/sediments->Reposangle[nc1];
        if(ratio>0){
          bedflux=sediments->Kbed*ratio*grad*grid->df[i];
          for(j=0;j<sediments->Nsize;j++){
            sediments->Layerthickness[j][nc1][0]-=bedflux/grid->Ac[nc1]*sediments->Toplayerratio[j][nc1];
            sediments->Layerthickness[j][nc2][0]+=bedflux/grid->Ac[nc2]*sediments->Toplayerratio[j][nc1];
          }
        }
      } else {
        ratio=(atan(-grad)-sediments->Reposangle[nc2])/sediments->Reposangle[nc2];
        if(ratio>0){
          bedflux=sediments->Kbed*ratio*grad*grid->df[i];
          for(j=0;j<sediments->Nsize;j++){
            sediments->Layerthickness[j][nc1][0]-=bedflux/grid->Ac[nc1]*sediments->Toplayerratio[j][nc2];
            sediments->Layerthickness[j][nc2][0]+=bedflux/grid->Ac[nc2]*sediments->Toplayerratio[j][nc2];
            if(sediments->Layerthickness[j][nc1][0]<0)
              sediments->Layerthickness[j][nc1][0]=0;
            if(sediments->Layerthickness[j][nc2][0]<0)
              sediments->Layerthickness[j][nc2][0]=0;                         
          }   
        }    
      }
    }
  }

  // update Layertop, Thicknesslayer and Totalthickness
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++){
      if(i==0)
        sediments->Totalthickness[j]=0;
      for(k=0;k<sediments->Nlayer;k++){
        if(i==0)
          sediments->Thicknesslayer[j][k]=0;
        sediments->Thicknesslayer[j][k]+=sediments->Layerthickness[i][j][k];
        if(i==sediments->Nsize-1)
          sediments->Totalthickness[j]+=sediments->Thicknesslayer[j][k];
      }
      if(i==sediments->Nsize-1){
        sediments->Layertop[j]=0;
        while(sediments->Thicknesslayer[j][sediments->Layertop[j]]<=0 && sediments->Layertop[j]<sediments->Nlayer){
          sediments->Layertop[j]++;
        }
      }
    }

  // update Toplayerratio
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      sediments->Toplayerratio[i][j]=sediments->Layerthickness[i][j][sediments->Layertop[j]]/sediments->Thicknesslayer[j][sediments->Layertop[j]];

  // update Reposangle
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      sediments->Reposangle[j]+=sediments->Toplayerratio[i][j]*sediments->Anglerepos[i];

}

/*
 * Function: SedimentSource
 * Usage: dT/dt + u dot grad T = d/dz ( kappaT dT/dz) + A + B*T
 *--------------------------------------------------------------
 * same use as Heatsource, add new parameter Nosize to represent the No. of Nsize for fraction
 *
 */
void SedimentSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop,int Nosize, REAL theta) {
  int i, k;
  REAL alpha;
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i]-1;k++) {
      B[i][k]=0;
      A[i][k]=0;
    }

    //advection
    alpha=1.0;
    if(prop->subgrid)
    {
      alpha=subgrid->Aceff[i]/subgrid->Acceff[i][grid->Nk[i]-1];
      if(subgrid->erosionpara)
        alpha*=subgrid->delta[i];

      if(subgrid->Acceff[i][grid->Nk[i]-1]==0)
        alpha=0;
    }
    
    //if(grid->Nkmax==1)
    // alpha is included in deposition np matter whether there is only one layer
    alpha*=sediments->alphaSSC[Nosize][i];

    B[i][grid->Nk[i]-1]=alpha*sediments->Ws[Nosize][i][grid->Nk[i]]/grid->dzz[i][k]; 

    //erosion
    alpha=1.0;
    if(prop->subgrid)
    {
      alpha=subgrid->Aceff[i]/subgrid->Acceff[i][grid->Nk[i]-1];
      if(subgrid->Acceff[i][grid->Nk[i]-1]==0)
        alpha=0;
    }
    A[i][grid->Nk[i]-1]=alpha*((1-theta)*sediments->Erosion_old[Nosize][i][0]+theta*sediments->Erosion[Nosize][i][0])/grid->dzz[i][grid->Nk[i]-1];
  }
}

/*
 * Function: SedimentVerticalVelocity
 * Usage: wnewsedi=phys->wnew-ws, woldsedi=phys->w_old-ws, ws is explicit
 *--------------------------------------------------------------
 * provide the vertical velocity field for updatescalar function
 *
 */
void SedimentVerticalVelocity(gridT *grid, physT *phys,int Nosize,int symbol, int myproc) {
  int i,k;
  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i]+1;k++)
      sediments->Wnewsedi[i][k]=phys->w_im[i][k]-sediments->Ws[Nosize][i][k];
}

/*
 * Function: OpenSediFiles
 * read sediment output file names and open sediment files
 *--------------------------------------------------------
 *
 */
void OpenSediFiles(propT *prop, int myproc) {
  int i;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  sediments->SedimentFID=(FILE **)SunMalloc(sediments->Nsize*sizeof(FILE *), "OpenSediFiles"); 
  for(i=0;i<sediments->Nsize;i++){
    sprintf(str,"Sediment%dFile",i+1);
    MPI_GetFile(filename,DATAFILE,str,"OpenSediFiles",myproc);
    if(prop->mergeArrays)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    sediments->SedimentFID[i] = MPI_FOpen(str,"w","OpenFiles",myproc); 
  }
  MPI_GetFile(filename,DATAFILE,"LayerFile","OpenSediFiles",myproc);
  if(prop->mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  sediments->LayerthickFID = MPI_FOpen(str,"w","OpenFiles",myproc);
  
  MPI_GetFile(filename,DATAFILE,"tbFile","OpenSediFiles",myproc);
  if(prop->mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  sediments->SeditbFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"tbmaxFile","OpenSediFiles",myproc);
  if(prop->mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  sediments->SeditbmaxFID = MPI_FOpen(str,"w","OpenFiles",myproc);
}

/*
 * Function: OutputSediment
 * output sediment data for each time step
 *--------------------------------------------------------------
 *
 */
void OutputSediment(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, int blowup, MPI_Comm comm) {
  int i, j, jptr, k, nwritten, nosize, nolayer;
  char str[BUFFERLENGTH];
  FILE *ofile;
  REAL thicktmp;

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {
    
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      sprintf(str,"Error outputting SSC data for size class %d of %d.\n",nosize+1,sediments->Nsize);
      Write3DData(sediments->SediC[nosize],phys->htmp,prop->mergeArrays,sediments->SedimentFID[nosize],
		  str,grid,numprocs,myproc,comm);
    }
    
    Write2DData(sediments->Totalthickness,prop->mergeArrays,sediments->LayerthickFID,"Error outputting bed thickness data!\n",
    		grid,numprocs,myproc,comm);
    
    Write2DData(sediments->Seditb,prop->mergeArrays,sediments->SeditbFID,"Error outputting bed shear stress data!\n",
		  grid,numprocs,myproc,comm);
  }
  
  if(prop->n==prop->nsteps+prop->nstart) {
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      fclose(sediments->SedimentFID[nosize]);
    }
    fclose(sediments->LayerthickFID);
   
    fclose(sediments->SeditbFID);

    Write2DData(sediments->Seditbmax,prop->mergeArrays,sediments->SeditbmaxFID,"Error outputting max bed shear stress data!\n",
		  grid,numprocs,myproc,comm);    
    fclose(sediments->SeditbmaxFID);
  }
}

/*
 * Function: CalculateSediDiffusivity
 * calculate sediment diffusivity 
 *--------------------------------------------------------------
 * may use phys->kappa_tv directly from my25 function or 
 * use parabolic diffusivity
 */
void CalculateSediDiffusivity(gridT *grid, physT *phys,int Nosize,int myproc) {  
  int ii, kk;
  REAL z;
  
  if(sediments->ParabolKappa==1){
    for(ii=0;ii<grid->Nc;ii++){
      z=grid->dv[ii]+phys->h[ii]-0.5*grid->dzz[ii][0];
      for(kk=0;kk<grid->Nk[ii];kk++){
	if(phys->CdB[ii]!=-1)
	  sediments->SediKappa_tv[ii][kk]=z*(1-z/(grid->dv[ii]+phys->h[ii]))*phys->uc[ii][grid->Nk[ii]-1]*sqrt(phys->CdB[ii])*0.41/sediments->Prt[Nosize];
	if(kk!=grid->Nk[ii]-1)
	  z-=0.5*(grid->dzz[ii][kk]+grid->dzz[ii][kk+1]);
      }
    }
  } else {
    for(ii=0;ii<grid->Nc;ii++)
      for(kk=0;kk<grid->Nk[ii];kk++)
	sediments->SediKappa_tv[ii][kk]=phys->kappa_tv[ii][kk]/sediments->Prt[Nosize];
  }
}

/*
 * Function: ISendRecvSediBedData3D
 * Usage: ISendRecvSediBedData3D(SediCbed or Thickness,grid,myproc,comm);
 * ----------------------------------------------------
 * This function will transfer the 3D cell data for sediment bed back and forth between
 * processors using nonblocking sends/recvs.
 *
 */
void ISendRecvSediBedData3D(REAL **celldata, gridT *grid, int nlayer,int myproc,MPI_Comm comm) {
  int k, n, nstart, neigh, neighproc;
  REAL t0=Timer();
  
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    
    nstart=0;
    for(n=0;n<grid->num_cells_send[neigh];n++) {
      for(k=0;k<nlayer;k++) 
        grid->send[neigh][nstart+k]=celldata[grid->cell_send[neigh][n]][k];
      nstart+=nlayer;
    }
    
    MPI_Isend((void *)(grid->send[neigh]),grid->total_cells_send[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[neigh])); 
  }
  
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    neighproc = grid->myneighs[neigh];
    MPI_Irecv((void *)(grid->recv[neigh]),grid->total_cells_recv[neigh],MPI_DOUBLE,neighproc,1,
	      comm,&(grid->request[grid->Nneighs+neigh]));
  }
  MPI_Waitall(2*grid->Nneighs,grid->request,grid->status);
  
  for(neigh=0;neigh<grid->Nneighs;neigh++) {
    nstart=0;
    for(n=0;n<grid->num_cells_recv[neigh];n++) {
      for(k=0;k<nlayer;k++) 
        celldata[grid->cell_recv[neigh][n]][k]=grid->recv[neigh][nstart+k];
      nstart+=nlayer;
    }
  }
  // t_comm+=Timer()-t0;
}

/*
 * Function: CalculateSedimentKb
 * usage: CalculateSedimentKb(gridT *grid, physT *phys, int myproc);
 * --------------------------------------------------------------------
 * This function calculate/update the bottom roughness for sediment 
 * transport and Yi-Ju's Wave model
 * 
 */
void CalculateSedimentKb(gridT *grid, physT *phys, int myproc)
{
  int i,n;
  REAL alpha, tauw, lambda_ano, eta_ano, d0, lambda_orb, lamda;
  REAL A1, A2, A3, B1, B2, B3, F, eta, ar, a1, a2, Tstar, tmp;

  alpha = 0.056;
  lambda_ano = 535*sediments->ds50;
  eta_ano = 0.171*lambda_ano;
 
  a1 = 0.068;
  a2 = 0.0204*pow(log(sediments->ds50),2)+0.022*log(sediments->ds50)+0.0709;
  ar = 0.267;
  
  A1 = 0.095;
  A2 = 0.442;
  A3 = 2.28;
  
  B1 = 1/A1;
  B2 = 0.5*(1+A2)*B1;
  B3 = pow(B2, 2)-A3*B1;

  for(i=0;i<grid->Nc;i++)
  {
    // ripple roughness
    //Obtain the ripple field
    /*
    if(prop->wavemodel && wprop->FetchModel==0){
      tauw=0.5*wave->fw[i]*1000*pow(wave->ub[i],2);
      if (tauw > sediments->Taue[i][sediments->Layertop[i]]){
        d0 = 2*wave->ab[i];
        lambda_orb = 0.62*d0;
        if (d0/eta_ano < 20)
          lambda = 0.62*d0;
        else if (d0/eta_ano <= 100)
          lambda = lambda_ano*exp(-log(lambda_orb/lambda_ano)*log(0.01*d0/eta_ano)/log(5));
        else
          lambda = lambda_ano;
      
        F = B3-B1*log(d0/lambda);
      
        if (F < 0)
          eta = 0;
        else
          eta = d0/(exp(B2-sqrt(F)));
        sedi->z0r[i] = ar*pow(eta, 2)/lambda;
      }
    }
    */
    // roughness due to shear stress
    Tstar=0;
    for(n=0;n<sediments->Nsize;n++)
      Tstar += sediments->Seditb[i]/sediments->Taue[n][sediments->Layertop[i]]/sediments->Nsize;

    if (Tstar > 1.0)
      sediments->z0b[i] = alpha*sediments->ds50*a1*Tstar/(1+a2*Tstar);
    
    //update kb
    sediments->kb[i] = 30*Max(sediments->z0s[i], sediments->z0b[i]+sediments->z0r[i]);
  }

}


/*
 * Function: ComputeSedimentsRestart
 * Usage: ComputeSedimentsRestart(grid,phys,prop, myproc, numproc, blowup, comm);
 * ----------------------------------------------------
 * This function is the function to prepare for the restart run of SUNTANS
 * It allocates and initializes all variables of sediment calculation
 *
*/
void ComputeSedimentsRestart(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int k,i;
  // when computeSediments=2 means the suntans use restart run
  sediments=(sedimentsT *)SunMalloc(sizeof(sedimentsT),"ComputeSedimentsRestart");
  // allocate and initialize all the sediment variables
  ReadSediProperties(myproc); 
  OpenSediFiles(prop,myproc); 
  AllocateSediment(grid,myproc);  
  InitializeSediment(grid,phys,prop,myproc);
  BoundarySediment(grid,phys,prop);   
}

/*
 * Function: ComputeSedimentsBedRestart
 * Usage: ComputeSedimentsBedRestart(grid,phys,prop, myproc, numproc, blowup, comm);
 * ----------------------------------------------------
 * This function is the function to prepare for the restart run of SUNTANS
 * It calculates the preparation for bed model
 *
*/
void ComputeSedimentsBedRestart(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int k,i,j;

  // set Sedic_old
  for(j=0;j<sediments->Nsize;j++)
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
        sediments->SediC_old[j][i][k]=sediments->SediC[j][i][k];

  // compute bed layer model
  // update Layertop, Thicknesslayer and Totalthickness
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++){
      if(i==0)
        sediments->Totalthickness[j]=0;
      for(k=0;k<sediments->Nlayer;k++){
          if(i==0)
          sediments->Thicknesslayer[j][k]=0;
        sediments->Thicknesslayer[j][k]+=sediments->Layerthickness[i][j][k];
        if(i==sediments->Nsize-1)
          sediments->Totalthickness[j]+=sediments->Thicknesslayer[j][k];
      }
      if(i==sediments->Nsize-1){
        sediments->Layertop[j]=0;
        while(sediments->Thicknesslayer[j][sediments->Layertop[j]]<=0 && sediments->Layertop[j]<sediments->Nlayer){
          sediments->Layertop[j]++;
        }
      }
    }

  // update Toplayerratio
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      sediments->Toplayerratio[i][j]=sediments->Layerthickness[i][j][sediments->Layertop[j]]/sediments->Thicknesslayer[j][sediments->Layertop[j]];

  // update Reposangle
  for(i=0;i<sediments->Nsize;i++)
    for(j=0;j<grid->Nc;j++)
      sediments->Reposangle[j]+=sediments->Toplayerratio[i][j]*sediments->Anglerepos[i];          

  // calculate settling velocity
  if(sediments->WSconstant==0)
    SettlingVelocity(grid,phys,prop,myproc);
  // calculate initial erosion and deposition
  CalculateErosion(grid,phys,prop,myproc);    
  CalculateDeposition(grid,phys,prop,myproc);
}


/*
 * Function: ComputeSediments
 * Usage: ComputeSediments(grid,phys,prop, myproc, numproc, blowup, comm);
 * ----------------------------------------------------
 * This function is the main function for sediment model part
 * include all the calculation for sediment transport
 * and called by phys.c
 *
*/
void ComputeSediments(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int k,i;

  // when computeSediments=2 means the suntans use restart run
  if(prop->n==1+prop->nstart && prop->computeSediments!=2){

    sediments=(sedimentsT *)SunMalloc(sizeof(sedimentsT),"ComputeSediments");
    // allocate and initialize all the sediment variables
    ReadSediProperties(myproc); 
    OpenSediFiles(prop,myproc); 
    AllocateSediment(grid,myproc);  
    InitializeSediment(grid,phys,prop,myproc);
    BoundarySediment(grid,phys,prop);   
  }
  
  // update wave->Cr 
  if(prop->wavemodel && wprop->FetchModel==0)
  {
    UpdateWaveCr(grid,phys,myproc);
    UpdateWaveFw(grid,myproc);
  }

  // calculate n+theta Erosion for boundary 
  if(prop->subgrid && grid->Nkmax==1 && subgrid->erosionpara)
  {
    CalculateSubgridCdmean(grid, phys, prop);
    CalculateSubgridActiveErosion(grid, phys, prop, myproc);
  } else {
    CalculateErosion(grid,phys,prop,myproc);
  }

  // update kb
  CalculateSedimentKb(grid,phys,myproc);

  // update SSCalpha
  if(grid->Nkmax==1 && sediments->sscvprof==1)
    RouseCurveAlpha(grid, phys,prop, myproc);

  // calculate n+1 Sediment concentration field
  for(k=0;k<sediments->Nsize;k++){
    SedimentSource(phys->wtmp,phys->uold,grid,phys,prop,k,prop->theta);
    SedimentVerticalVelocity(grid,phys,k,1,myproc);
    CalculateSediDiffusivity(grid,phys,k,myproc);
    UpdateScalars(grid,phys,prop,sediments->Wnewsedi,sediments->SediC[k],sediments->SediC_old[k],
      sediments->boundary_sediC[k],phys->Cn_T,
      0,0,sediments->SediKappa_tv,prop->theta,
      phys->uold,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,prop->TVDtemp);
    //SedimentVerticalVelocity(grid,phys,k,-1,myproc);
    ISendRecvCellData3D(sediments->SediC[k],grid,myproc,comm);
  }          
  if(sediments->WSconstant==0)
    SettlingVelocity(grid,phys,prop,myproc);
  // update Deposition
  CalculateDeposition(grid,phys,prop,myproc);
  // update bed change
  if(prop->n%sediments->bedInterval==0 && sediments->bedInterval>0)
    BedChange(grid,phys,prop,myproc);   
  // get the boundary value for the next time step
  BoundarySediment(grid,phys,prop);
  // output sediment results
  OutputSediment(grid,phys,prop,myproc,numprocs,blowup,comm);
  // free space
  //if(prop->n==prop->nstart+prop->nsteps)
  //FreeSediment(grid,sediments,myproc);
}

