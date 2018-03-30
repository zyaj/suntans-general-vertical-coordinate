/*
 * File: culvert.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This contains all the functions needed when there are culverts
 * in calculation domain
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
#include "culvert.h"
#include "subgrid.h"

void ReadCulvertProperties(int myproc);
void AllocateCulvert(gridT *grid,int myproc);
void InitializeCulvert(gridT *grid, physT *phys,propT *prop, int myproc);

/*
 * Function: SetupCulvertmodel
 * Usage: first function to setup all the properties of culvert model
 * ----------------------------------------------------
 *
 */
void SetupCulvertmodel(gridT *grid, physT *phys, propT *prop, int myproc){
  culvert=(culvertT *)SunMalloc(sizeof(culvertT),"SetupCulvertmodel");
  ReadCulvertProperties(myproc);      
  AllocateCulvert(grid,myproc);      
  InitializeCulvert(grid,phys,prop,myproc);
}

/*
 * Function: ReadCulvertProperties
 * Usage: ReadCulvertProperties(myproc);
 * ----------------------------------------------------
 * Based on culvert.dat, load in the important parameters for 
 * wave model
 *
 */
void ReadCulvertProperties(int myproc)
{
  culvert->CdC=MPI_GetValue(DATAFILE,"cdculvert","ReadCulvertProperties",myproc);
  culvert->constant= MPI_GetValue(DATAFILE,"constantculvert","ReadCulvertProperties",myproc);
  culvert->eps= MPI_GetValue(DATAFILE,"culverteps","ReadCulvertProperties",myproc);
}

/*
 * Function: AllocateCulvert
 * Usage: allocate space for culvert variables
 * ----------------------------------------------------
 * Based on the value from ReadCulvertProperties
 *
 */
void AllocateCulvert(gridT *grid,int myproc)
{
  culvert->top = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->pressure = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->pressure2 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->pressure3 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->pressure4 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->toparea = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->Qcoef = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->condition = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->condition2 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->source1 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->htmp = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->source2 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->source3 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  culvert->source4 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
}

/*
 * Function: FreeCulvert
 * Usage: free space for all the variables
 * ----------------------------------------------------
 * Basic sunfree function
 *
 */
void FreeCulvert(int myproc)
{
  free(culvert->top);
  free(culvert->condition);
  free(culvert->condition2);
  free(culvert->pressure);
  free(culvert->pressure2);
  free(culvert->pressure3);
  free(culvert->pressure4);
  free(culvert->Qcoef);
  free(culvert->source1);
  free(culvert->source2);
  free(culvert->source3);
  free(culvert->source4);
  free(culvert->htmp);
}

/*
 * Function: InitializeCulvert
 * usage: initial culvert variables
 * -----------------------------------------------------
 * May use Returnculvert->top function in Initialize.c
 * 
 */
void InitializeCulvert(gridT *grid, physT *phys, propT *prop, int myproc)
{ 
  int i;
  REAL culvertheight;
  culvertheight=MPI_GetValue(DATAFILE,"culvertheight","InitializeCulvert",myproc);
  for(i=0;i<grid->Nc;i++){
    if(culvert->constant==1)
      culvert->top[i]=culvertheight-grid->dv[i];
    else
      culvert->top[i]=ReturnCulvertTop(grid->xv[i],grid->yv[i],grid->dv[i]);
    culvert->condition[i]=0; //assume free surface first
    culvert->condition2[i]=0; //assume free surface first
    culvert->htmp[i]=0;
    culvert->source1[i]=0;
    culvert->source2[i]=0;
    culvert->source3[i]=0;
    culvert->source4[i]=0;
    culvert->pressure[i]=phys->h[i];
    culvert->pressure2[i]=phys->h[i];
    culvert->pressure3[i]=phys->h[i];
    culvert->pressure4[i]=phys->h[i];
    // change  10/01/2014
    if(phys->h[i]>culvert->top[i])
      phys->h[i]=culvert->top[i];
    if(culvert->top[i]!=INFTY)
      culvert->toparea[i]=grid->Ac[i];
    else 
      culvert->toparea[i]=0;

    culvert->Qcoef[i]=0;
  }
  culvert->sum=0;
  UpdateDZ(grid,phys,prop,1);
}

/*
 * Function: CheckCulvertCondition
 * usage: check whether the culvert assumption (free surface, 0 or pressurized flow, 1)
 * whether meets results
 * ---------------------------------------------------------------------------------
 * calculate the residual of the newton iteration method 
 *
 */ 
void CheckCulvertCondition(gridT *grid,physT *phys, propT *prop, int myproc)
{  
   int i,iptr,j;
   REAL sum;
   for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    if(prop->subgrid)
      culvert->condition[i]+=subgrid->Veff[i]-subgrid->Aceff[i]*phys->h[i];
    else
    {
      if(phys->h[i]>=-grid->dv[i] && culvert->htmp[i]>=-grid->dv[i])
        culvert->condition[i]+=grid->dv[i]*grid->Ac[i];
      else if(phys->h[i]>=-grid->dv[i] && culvert->htmp[i]<-grid->dv[i])
        culvert->condition[i]+=grid->Ac[i]*(grid->dv[i]+phys->h[i]);
      else if(phys->h[i]<-grid->dv[i] && culvert->htmp[i]>=-grid->dv[i])
        culvert->condition[i]-=grid->Ac[i]*phys->h[i];
      else if(phys->h[i]<-grid->dv[i] && culvert->htmp[i]<-grid->dv[i])
        culvert->condition[i]=0;
      //else if(phys->h[i]>Culverttop[i] && culvert->htmp[i]>=-grid->dv[i])
      //  culvert->htmp[i]+=grid->Ac[i]*(Culverttop[i]+grid->dv[i]-phys->h[i]);
      //else if(phys->h[i]>Culverttop[i] && culvert->htmp[i]<-grid->dv[i])
      //  culvert->htmp[i]+=grid->Ac[i]*(Culverttop[i]+grid->dv[i]);
      else
        printf(" %dth cell has something wrong with the residual\n",i);
    }
   }
}

/*
 * Function: SetCulvertDragCoefficient 
 * usage: SetCulvertDragCoefficient(gridT grid, physT phys, int myproc)
 * ----------------------------------------------------------------------
 * set cdb and cdt FOR culvert edge, if CulvertCodition(grad(edge))=1 set 
 * cdT and cdB as culvert->CdC, otherwise only cdB will be changed to culvert->CdC
 */
void SetCulvertDragCoefficient(gridT *grid, physT *phys, int myproc)
{
  int i,ne1,ne2;
  for(i=0;i<grid->Ne;i++){
    ne1 = grid->grad[2*i];
    ne2 = grid->grad[2*i+1];
    if(ne1!=-1 && ne2!=-1) // not for bc edge, because the velocity at bc is not calculated but set
      if(culvert->top[ne1]!=INFTY || culvert->top[ne2]!=INFTY){
        phys->CdB[i]=culvert->CdC;
        if(phys->h[ne1]>culvert->top[ne1] || phys->h[ne2]>culvert->top[ne2])
          phys->CdT[i]=culvert->CdC;   
      } 

    if(grid->dzf[i][grid->Nke[i]-1]<BUFFERHEIGHT && grid->etop[i]==grid->Nke[i]-1){
      phys->CdB[i]=100;
      //printf("Making CdB a large value due to small cell!\n");
    }
  }
}

/*
 * Function: CulvertHCoefficient
 * usage: change the coefficients of free surface eqn. 
 * -----------------------------------------------------------
 * Called by CGSolve function, change Hcoefficient with new continuity eqn. without time gradient
 * of free surface. it also changes htmp, to delete time gradient part
 *
 */
void CulvertHCoefficients(REAL *coef, REAL *fcoef, gridT *grid, physT *phys, propT *prop, int myproc)
{
  int i,iptr,ne, jptr, nf,k,check;
  REAL sum,normal, tmp,fac;
   
  fac=prop->imfac1;

  tmp=prop->grav*pow(fac*prop->dt,2);
  // calculate fcoef and coef
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    check=1;
    if(!prop->subgrid){   

      if(phys->h[i]<-grid->dv[i] && culvert->pressure[i]<culvert->top[i])
        coef[i] = 0-culvert->Qcoef[i];
      else
        coef[i] = grid->Ac[i]-culvert->Qcoef[i];

      if(culvert->pressure[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
      //if(culvert->htmp[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
        coef[i]=0;

    } else {
      coef[i]=subgrid->Aceff[i]-culvert->Qcoef[i];

     if(culvert->pressure[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
        coef[i]=0;

    }

   for(nf=0;nf<grid->nfaces[i];nf++) 
     if(grid->neigh[i*grid->maxfaces+nf]!=-1) {
       ne = grid->face[i*grid->maxfaces+nf];
       fcoef[i*grid->maxfaces+nf]=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne];
       coef[i]+=fcoef[i*grid->maxfaces+nf];

      if(fcoef[i*grid->maxfaces+nf]>0)
        check=0;
     }
    
   // when 0=0 exists make sure hnew=hold
   if(check)
   {
     coef[i]=1.0;
     phys->htmp[i]=phys->h[i];
   }
  }
}

/*
 * Function: Storeculvert->pressure
 * usage: store phys->h results with culvert pressure, change phys->h to culvertop for 
 * culvert->condition=1 cell 
 * -----------------------------------------------------------------------------------
 * no=1 to store pressure, change phys->h, no=0 change phys->h to be pressure field
 */
void StoreCulvertPressure(REAL *h, int Nc, int no, int myproc)
{ 
  int i;
  if(no==1){
    for(i=0;i<Nc;i++){
      culvert->pressure[i]=h[i]; //store new pressure field (for culvert part maybe not freesurface)
      culvert->pressure2[i]=h[i]; //store new pressure field (for culvert part maybe not freesurface)
      if(h[i]>culvert->top[i] && culvert->top[i]!=INFTY)
        h[i]=culvert->top[i]; // transfer it into Culvertop for fully discharge part
    }
  } else {
    for(i=0;i<Nc;i++)
      h[i]=culvert->pressure[i];
  }  
}

/*
 * Function: CulvertIterationSource
 * usage: calculate source term for each iteration -----------------------------------------------------------------------------------
 * based on Newton-nested method 
 */
void CulvertIterationSource(gridT *grid, physT *phys,  propT *prop, REAL theta,REAL dt, int myproc)
{ 
  int iptr,i,nf,ne,k,check=0;
  REAL sum,normal;
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    phys->htmp[i]=culvert->source1[i];
    if(!prop->subgrid){

      if(phys->h[i]<-grid->dv[i] && culvert->pressure[i]<culvert->top[i])
        phys->htmp[i]+=grid->Ac[i]*grid->dv[i];
   
      if(culvert->pressure2[i]>culvert->top[i])
        phys->htmp[i]-=culvert->toparea[i]*culvert->top[i];

      if(culvert->pressure[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
      //if(culvert->htmp[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
      {
        phys->htmp[i]=culvert->source1[i]-culvert->toparea[i]*culvert->top[i]; 
        sum = 0;
        // add back explicit part
        for(nf=0;nf<grid->nfaces[i];nf++) {
          ne = grid->face[i*grid->maxfaces+nf];
          normal = grid->normal[i*grid->maxfaces+nf];
          for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
            sum+=(1-theta)*phys->u[ne][k]*grid->df[ne]*normal*grid->dzf[ne][k];
        } 
        phys->htmp[i]+=dt*sum;
      }
    } else {
      phys->htmp[i]-=(subgrid->Veff[i]-phys->h[i]*subgrid->Aceff[i]);
 
      if(culvert->pressure2[i]>culvert->top[i])
        phys->htmp[i]-=culvert->toparea[i]*culvert->top[i]; 

      if(culvert->pressure[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
      //if(culvert->htmp[i]>culvert->top[i] && phys->h[i]>culvert->top[i])
      {
        phys->htmp[i]=culvert->source1[i]-subgrid->Veffold[i];//+culvert->toparea[i]*(culvert->pressure[i]-culvert->top[i]); 
        sum = 0;
        // add back explicit part
        for(nf=0;nf<grid->nfaces[i];nf++) {
          ne = grid->face[i*grid->maxfaces+nf];
          normal = grid->normal[i*grid->maxfaces+nf];
          for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
            sum+=(1-theta)*phys->u[ne][k]*grid->df[ne]*normal*grid->dzf[ne][k];
        } 
        phys->htmp[i]+=dt*sum;
      }
    }     
  }
}

/*
 * Function: CulvertInitIteration
 * usage: (restore htmp2 htmp3 hold and store h guess) or (store htmp2 htmp3 hold) 
 -----------------------------------------------------------------------------------
 * no=-1 store no=1 restore 
 */
void CulvertInitIteration(gridT *grid, physT *phys, propT *prop, int no, int myproc)
{ 
  int iptr,i;

  if(no==-1){
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      if(prop->subgrid)
        culvert->condition[i]=-subgrid->Veff[i]+subgrid->Aceff[i]*phys->h[i];
      else
      {
        if(phys->h[i]>=-grid->dv[i])
          culvert->condition[i]=-grid->Ac[i]*grid->dv[i];
        else
          culvert->condition[i]=0;
      }
      // store formal iteration value
      culvert->htmp[i]=phys->h[i];
    }

    // store each guess
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      phys->hold[i]=culvert->source2[i];
      phys->htmp2[i]=culvert->source3[i];
      phys->htmp3[i]=culvert->source4[i];
    }
  } else {
    // store each guess
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      if(phys->h[i]>culvert->top[i]){ 
        if(!prop->subgrid)
          phys->htmp[i]=phys->htmp[i]-culvert->toparea[i]*phys->h[i]+culvert->toparea[i]*culvert->top[i];
        else 
          phys->htmp[i]=phys->htmp[i];//-culvert->toparea[i]*(phys->h[i]-culvert->top[i]);
      }
      culvert->source1[i]=phys->htmp[i];
      culvert->source2[i]=phys->hold[i];
      culvert->source3[i]=phys->htmp2[i];
      culvert->source4[i]=phys->htmp3[i];
    }
  }
}

/*
 * Function: UpdateculvertQcoef
 * usage: Update Qcoef for outer loop based on culvert->pressure2
 -----------------------------------------------------------------------------------
 * 
 */
void UpdateCulvertQcoef(gridT *grid, propT *prop,int no, int myproc)
{
   int nc;
   if(no==0)
     for(nc=0;nc<grid->Nc;nc++)
       if(culvert->pressure2[nc]>culvert->top[nc])
         culvert->Qcoef[nc]=culvert->toparea[nc];
       else
         culvert->Qcoef[nc]=0;
   else 
     for(nc=0;nc<grid->Nc;nc++)
       culvert->Qcoef[nc]=0;   
}