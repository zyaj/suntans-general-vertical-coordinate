/*
 * File: marsh.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file contains physically-based functions for marsh model parameterization.
 * The marsh model effects will be included when prop->marsh==1
 * The marsh data is set in marsh.dat
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
#include "marsh.h"
#include "physio.h"
#include "subgrid.h"

void InterpMarsh(gridT *grid, physT *phys, propT *prop, int myproc,int numprocs);
void ReadMarshProperties(int myproc);
void AllocateMarsh(gridT *grid, int myproc);
void OutputMarsh(gridT *grid,physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm);

/*
 * Function: ReadMarshProperties
 * Usage: ReadMarshProperties(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * Based on marsh.dat, load in the important parameters for 
 * marsh model
 *
 */
void ReadMarshProperties(int myproc)
{ 
  marsh->IntCdV = MPI_GetValue(DATAFILE,"IntCdV","ReadMarshProperties",myproc);
  marsh->Inthmarsh = MPI_GetValue(DATAFILE,"Inthmarsh","ReadMarshProperties",myproc);
  marsh->cdV = MPI_GetValue(DATAFILE,"CdV","ReadMarshProperties",myproc);
  marsh->Alphav = MPI_GetValue(DATAFILE,"alphav","ReadMarshProperties",myproc);
  marsh->Hmarsh = MPI_GetValue(DATAFILE,"hmarsh","ReadMarshProperties",myproc);
  marsh->Na = MPI_GetValue(DATAFILE,"Na","ReadMarshProperties",myproc);
  marsh->Rv = MPI_GetValue(DATAFILE,"rv","ReadMarshProperties",myproc);
}

/*
 * Function: AllocateMarsh
 * Allocate space for marsh variables
 * ----------------------------------------------------
 *
 */
void AllocateMarsh(gridT *grid, int myproc)
{
  marsh->hmarshcenter = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateMarsh");
  marsh->CdVcenter = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateMarsh");
  marsh->marshtop = (int *)SunMalloc(grid->Ne*sizeof(int),"AllocateMarsh");
  marsh->hmarsh = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"AllocateMarsh");
  marsh->hmarshleft = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"AllocateMarsh");
  marsh->CdV = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"AllocateMarsh");
}

/*
 * Function: FreeMarsh
 * Usage: free space for all the variables
 * ----------------------------------------------------
 * Basic sunfree function
 *
 */
void FreeMarsh(gridT *grid, int myproc)
{
  free(marsh->CdV);
  free(marsh->marshtop);
  free(marsh->hmarsh);
  free(marsh->hmarshleft);
  free(marsh->hmarshcenter);
  free(marsh->CdVcenter);
}

/*
 * Function: OutputMarsh
 * Usage: output hmarsh and cdv file for sunplot
 * ----------------------------------------------------
 * 
 */
void OutputMarsh(gridT *grid,physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm)
{
  int n,nf,ne;
  REAL *CdVtmp, *hmarshtmp;
  char str[BUFFERLENGTH], str1[BUFFERLENGTH], str2[BUFFERLENGTH];
  FILE *ofile;
  // for CdV

  MPI_GetFile(str1,DATAFILE,"CdVFile","OutputMarsh",myproc);
  if(prop->mergeArrays)
    strcpy(str,str1);
  else
    sprintf(str,"%s.%d",str1,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str); 
  ofile = MPI_FOpen(str,"w","OutputMarsh",myproc);
  for(n=0;n<grid->Nc;n++) {
    marsh->CdVcenter[n]=0;
    for(nf=0;nf<grid->nfaces[n];nf++) {
      ne = grid->face[n*grid->maxfaces+nf];
      marsh->CdVcenter[n]+=marsh->CdV[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
    }
    marsh->CdVcenter[n]/=2*grid->Ac[n];  
  }

  Write2DData(marsh->CdVcenter,prop->mergeArrays,ofile,"Error outputting marsh drag coefficient data!\n",
		  grid,numprocs,myproc,comm);
  fclose(ofile);

  // for hmarsh
  MPI_GetFile(str2,DATAFILE,"hmarshFile","OutputMarsh",myproc);
  if(prop->mergeArrays)
    strcpy(str,str2);
  else
    sprintf(str,"%s.%d",str2,myproc);
  if(VERBOSE>2) printf("Outputting %s...\n",str); 
  ofile = MPI_FOpen(str,"w","OutputMarsh",myproc);
  for(n=0;n<grid->Nc;n++) {
    marsh->hmarshcenter[n]=0;
    for(nf=0;nf<grid->nfaces[n];nf++) {
      ne = grid->face[n*grid->maxfaces+nf];
      marsh->hmarshcenter[n]+=marsh->hmarsh[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
    }
    marsh->hmarshcenter[n]/=2*grid->Ac[n];
  }
  Write2DData(marsh->hmarshcenter,prop->mergeArrays,ofile,"Error outputting marsh height data!\n",
		  grid,numprocs,myproc,comm);
  fclose(ofile);
  
 // for hmarsh
 ofile = MPI_FOpen("hmarshedge.dat","w","OutputData",myproc);

 for(n=0;n<grid->Ne;n++)
   fprintf(ofile,"%e %e %e\n",grid->xe[n],grid->ye[n],marsh->hmarsh[n]);
 fclose(ofile);
}

/*
* Function: InterpMarsh(grid,phys,prop,myproc)
* usage: interpolate the value for hmarsh and CdV
* --------------------------------------------     
* Author: Yun Zhang   
* Interpolate hmarsh and CdV according to CdVint.dat and hmarshint.dat
* IntCdV==1, interpolate value
* IntCdV==2, read center point data directly, same for Inthmarsh
* only works when marshmodel==1
*/
void InterpMarsh(gridT *grid, physT *phys, propT *prop,int myproc, int numprocs)
{  
   int n,Nv,Nh;
   REAL maxcdv,maxhmarsh,*xv,*yv,*CdVtmp,*xh,*yh,*hmarshtmp;
   char str[BUFFERLENGTH],str1[BUFFERLENGTH];
   FILE *fid;
   // for CdV
   MPI_GetFile(str1,DATAFILE,"InputCdVFile","ReadMarshProperties",myproc);

   if(marsh->IntCdV==1){
     maxcdv=0;
     Nv = MPI_GetSize(str1,"InterpMarsh",myproc); // change to z0TFILE
     xv = (REAL *)SunMalloc(Nv*sizeof(REAL),"InterpMarsh");
     yv = (REAL *)SunMalloc(Nv*sizeof(REAL),"InterpMarsh");
     CdVtmp = (REAL *)SunMalloc(Nv*sizeof(REAL),"InterpMarsh");
     fid = MPI_FOpen(str1,"r","InterpMarsh",myproc);

     for(n=0;n<Nv;n++) {
       xv[n]=getfield(fid,str);
       yv[n]=getfield(fid,str);
       CdVtmp[n]=getfield(fid,str); // but cd should be positive already
       if(maxcdv<CdVtmp[n])
         maxcdv=CdVtmp[n];
     }
     fclose(fid);
     Interp(xv,yv,CdVtmp,Nv,&(grid->xe[0]), &(grid->ye[0]),&(marsh->CdV[0]),grid->Ne,grid->maxfaces);

     for(n=0;n<grid->Ne;n++)
	     if(marsh->CdV[n]>=0.95*1.3)
	       marsh->CdV[n]=1.3;
       else
	       marsh->CdV[n]=0;
     free(xv);
     free(yv);
     free(CdVtmp);
   } else if(marsh->IntCdV==2){
     if(numprocs==1)
       sprintf(str,"%s-edge",str1);
     else
       sprintf(str,"%s-edge.%d",str1,myproc);
     fid = MPI_FOpen(str,"r","InterpDrag",myproc);
     for(n=0;n<grid->Ne;n++) {
       getfield(fid,str);
       getfield(fid,str);
       marsh->CdV[n]=getfield(fid,str);
     }
     fclose(fid);
   } else if(marsh->IntCdV==0){
     for(n=0;n<grid->Ne;n++) {
       marsh->CdV[n]=ReturnMarshDragCoefficient(grid->xe[n],grid->ye[n]);
     }
   } else if(marsh->IntCdV==3){
     for(n=0;n<grid->Ne;n++) {
       marsh->CdV[n]=marsh->cdV;
     }
   } else{
     printf("IntCdV=%d, IntCdV can only be 0, 1, 2 and 3\n",marsh->IntCdV);
     MPI_Finalize();
     exit(EXIT_FAILURE);
   }


   // for hamrsh
   MPI_GetFile(str1,DATAFILE,"InputhmarshFile","ReadMarshProperties",myproc);
   if(marsh->Inthmarsh==1){
     maxhmarsh=0;
     Nh = MPI_GetSize(str1,"InterpMarsh",myproc); 
     xh = (REAL *)SunMalloc(Nh*sizeof(REAL),"InterpMarsh");
     yh = (REAL *)SunMalloc(Nh*sizeof(REAL),"InterpMarsh");
     hmarshtmp = (REAL *)SunMalloc(Nh*sizeof(REAL),"InterpMarsh");
     fid = MPI_FOpen(str1,"r","InterpMarsh",myproc);
     for(n=0;n<Nh;n++) {
       xh[n]=getfield(fid,str);
       yh[n]=getfield(fid,str);
       hmarshtmp[n]=getfield(fid,str); // but cd should be positive already
       if(maxhmarsh<hmarshtmp[n])
         maxhmarsh=hmarshtmp[n];
     }
     fclose(fid);
     Interp(xh,yh,hmarshtmp,Nh,&(grid->xe[0]), &(grid->ye[0]),&(marsh->hmarsh[0]),grid->Ne,grid->maxfaces);
     for(n=0;n<grid->Ne;n++)
       if(marsh->hmarsh[n]>=0.95*1)
	       marsh->hmarsh[n]=1;
       else
	       marsh->CdV[n]=0;
     free(xh);
     free(yh);
     free(hmarshtmp);

   } else if(marsh->Inthmarsh==2){
     if(numprocs==1)
       sprintf(str,"%s-edge",str1);
     else
       sprintf(str,"%s-edge.%d",str1,myproc);
     fid = MPI_FOpen(str,"r","InterpMarsh",myproc);
     for(n=0;n<grid->Ne;n++) {
       getfield(fid,str);
       getfield(fid,str);
       marsh->hmarsh[n]=getfield(fid,str);
     }
     fclose(fid);
   } else if(marsh->Inthmarsh==0){
     for(n=0;n<grid->Ne;n++) {
       marsh->hmarsh[n]=ReturnMarshHeight(grid->xe[n],grid->ye[n]);
     }
   } else if(marsh->Inthmarsh==3){
     for(n=0;n<grid->Ne;n++) {
       marsh->hmarsh[n]=marsh->Hmarsh;
     }
   } else {
     printf("Inthmarsh=%d, Inthmarsh can only be 0, 1 and 2\n",marsh->Inthmarsh);
     MPI_Finalize();
     exit(EXIT_FAILURE);
   }
}

/* 
* Function: SetMarshTop(grid,phys,myproc)
* usage: find the top layer which have marsh part
* ------------------------------------------------
* Author: Yun Zhang
* dzf here = 0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k])
* store the layer number at phys->marshtop
* only calculate the top layer for computational edge
*/
void SetMarshTop(gridT *grid, physT *phys, int myproc)
{
  int i,jptr,j,nc1,nc2,n;
  REAL d;
  // initialize
  for(i=0;i<grid->Ne;i++){
    marsh->marshtop[i]=grid->Nke[i]-1;
    marsh->hmarshleft[i]=0;
  } 
  // for computational edges 
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    n=grid->Nke[j]-1;
    d=marsh->hmarsh[j];
    // for edges with multiple layers
    if(grid->Nke[j]-grid->etop[j]>1){
      while(d>0.5*(grid->dzz[nc1][n]+grid->dzz[nc2][n]) && n>grid->etop[j]){
        d-=0.5*(grid->dzz[nc1][n]+grid->dzz[nc2][n]);
        n--;  
      }
      if(n==grid->etop[j]){
        marsh->hmarshleft[j]=Min(0.5*(grid->dzz[nc1][n]+grid->dzz[nc2][n]),d);
      } else {    
        marsh->hmarshleft[j]=d;
      }
      marsh->marshtop[j]=n;
    } else {
      marsh->hmarshleft[j]=Min(0.5*(grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]]),d);
      marsh->marshtop[j]=grid->etop[j];
    }
  }
}

/* 
* Function: Marshexplicitterm(grid,phys,prop,edgep,theta,dt, myproc)
* usage: calculate marsh source explicit term for 1 edge including all layers
* ------------------------------------------------
* Author: Yun Zhang
* called by Upredictor function to calculate marsh effest
*/
void MarshExplicitTerm(gridT *grid, physT *phys, propT *prop,int edgep, REAL theta, REAL dt, int myproc)
{ 
  int iv,j,nc1,nc2,base,base1,i;
  REAL p,p1,d1,d2,h;
  j=edgep;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1)
    nc1=nc2;
  if(nc2==-1)
    nc2=nc1;
  // FOR INNER LAYERS
  for(iv=marsh->marshtop[j]+1;iv<grid->Nke[j];iv++){
    phys->utmp[j][iv]-=(1-theta)*dt*marsh->CdV[j]*fabs(phys->u[j][iv])*phys->u[j][iv]*marsh->Alphav\
      *2*marsh->Rv*marsh->Na/(1-pow(marsh->Rv,2)*PI*marsh->Na)/pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
  }
  // FOR TOP LAYERS
  p=marsh->hmarshleft[j]/(0.5*(grid->dzz[nc1][marsh->marshtop[j]]+grid->dzz[nc2][marsh->marshtop[j]]));
  if(prop->subgrid && (grid->Nke[j]-grid->etop[j])==1)
  { 
    p=0;
    p1=0;
    h=subgrid->he[j];
    base=j*subgrid->segN;
    base1=j*(subgrid->segN+1);
    for(i=0;i<subgrid->segN;i++)
    {
       d1=subgrid->dpe[base1+i];
       d2=subgrid->dpe[base1+i+1];
       if(((d1+d2)/2+h)>BUFFERHEIGHT)
       {
         p+=Min((d1+d2)/2+h,subgrid->hmarshe[base+i]);
         p1+=(d1+d2)/2+h;        
       }
    } 
    p=p/p1;
    if(p1==0)
      p=0; 
  }

  phys->utmp[j][marsh->marshtop[j]]-=(1-theta)*dt*marsh->CdV[j]*fabs(phys->u[j][marsh->marshtop[j]])\
    *phys->u[j][marsh->marshtop[j]]*marsh->Alphav*marsh->Rv*marsh->Na*pow(p,2)\
    /(1-PI*pow(marsh->Rv,2)*marsh->Na*p)/pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
}

/* 
* Function: Marshimplicitterm(grid,phys,prop,edgep,layer,theta,dt,myproc)
* usage: calculate marsh source implicit effect and add on the trisolve coefficient
* ------------------------------------------------
* Author: Yun Zhang
* called by Upredictor function to calculate marsh effest
*/
REAL MarshImplicitTerm(gridT *grid, physT *phys, propT *prop,int edgep, int layer, REAL theta, REAL dt, int myproc)
{ 
  int jv,j,nc1,nc2,i,base,base1;
  REAL term,p,p1,h,d1,d2;
  j=edgep;
  jv=layer;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1)
    nc1=nc2;
  if(nc2==-1)
    nc2=nc1;
  p=marsh->hmarshleft[j]/(0.5*(grid->dzz[nc1][marsh->marshtop[j]]+grid->dzz[nc2][marsh->marshtop[j]]));

  // add subgrid calculation
  if(prop->subgrid && (grid->Nke[j]-grid->etop[j])==1)
  { 
    p=0;
    p1=0;
    h=subgrid->he[j];
    base=j*subgrid->segN;
    base1=j*(subgrid->segN+1);
    for(i=0;i<subgrid->segN;i++)
    {
       d1=subgrid->dpe[base1+i];
       d2=subgrid->dpe[base1+i+1];
       if(((d1+d2)/2+h)>BUFFERHEIGHT)
       {
         p+=Min((d1+d2)/2+h,subgrid->hmarshe[base+i]);
         p1+=(d1+d2)/2+h;        
       }
    } 

    p=p/p1;
    if(p1==0)
      p=0; 
  }

  if(grid->Nke[j]-grid->etop[j]>1){
     if(jv>marsh->marshtop[j]) 
       // inner layers
       term=theta*dt*marsh->CdV[j]*fabs(phys->u[j][jv])*marsh->Rv*marsh->Na*marsh->Alphav\
         /(1-pow(marsh->Rv,2)*PI*marsh->Na)/pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
     else     
       // top layers
       term=theta*dt*marsh->CdV[j]*fabs(phys->u[j][marsh->marshtop[j]])*marsh->Rv*marsh->Na\
         *marsh->Alphav*pow(p,2)/(1-PI*pow(marsh->Rv,2)*marsh->Na*p)\
         /pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
  } else {
    term=theta*dt*marsh->CdV[j]*fabs(phys->u[j][marsh->marshtop[j]])*marsh->Rv*marsh->Na\
      *marsh->Alphav*pow(p,2)/(1-PI*pow(marsh->Rv,2)*marsh->Na*p)\
      /pow((1-2*marsh->Rv*sqrt(marsh->Na)),2);
  }
  return term;
}

/* 
* Function: SetupMarshmodel
* usage: the setup function for marsh model
* ------------------------------------------------
* Author: Yun Zhang
*/
void SetupMarshmodel(gridT *grid,physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm)
{
  if(myproc==0)
    printf("\n\nmarsh model has been started\n\n");
  marsh=(marshT *)SunMalloc(sizeof(marshT),"SetupMarshmodel");
  ReadMarshProperties(myproc);
  AllocateMarsh(grid,myproc);
  InterpMarsh(grid,phys,prop,myproc,numprocs);
  OutputMarsh(grid,phys,prop,myproc,numprocs,comm);
}