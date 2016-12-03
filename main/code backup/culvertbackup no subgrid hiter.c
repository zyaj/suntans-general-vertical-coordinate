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
  ConstantCulvert= MPI_GetValue(DATAFILE,"constantculvert","ReadCulvertProperties",myproc);
  Culverteps= MPI_GetValue(DATAFILE,"culverteps","ReadCulvertProperties",myproc);
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
  Culverttop = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  Culvertpressure = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  Culvertcondition = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  Culvertsource1 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  Culvertsource2 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  Culvertsource3 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
  Culvertsource4 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateCulvert");
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
  free(Culverttop);
  free(Culvertcondition);
  free(Culvertpressure);
  free(Culvertsource1);
  free(Culvertsource2);
  free(Culvertsource3);
  free(Culvertsource4);
}

/*
 * Function: InitializeCulvert
 * usage: initial culvert variables
 * -----------------------------------------------------
 * May use ReturnCulvertTop function in Initialize.c
 * 
 */
void InitializeCulvert(gridT *grid,physT *phys, int myproc)
{ 
  int i;
  REAL culvertheight;
  culvertheight=MPI_GetValue(DATAFILE,"culvertheight","InitializeCulvert",myproc);
  for(i=0;i<grid->Nc;i++){
    if(ConstantCulvert==1)
      Culverttop[i]=culvertheight-grid->dv[i];
    else
      Culverttop[i]=ReturnCulvertTop(grid->xv[i],grid->yv[i],grid->dv[i]);
    Culvertcondition[i]=0; //assume free surface first
    Culvertsource1[i]=0;
    Culvertsource2[i]=0;
    Culvertsource3[i]=0;
    Culvertsource4[i]=0;
    Culvertpressure[i]=phys->h[i];
  }
  Culvertsum=0;
}

/*
 * Function: CheckCulvertCondition
 * usage: check whether the culvert assumption (free surface, 0 or pressurized flow, 1)
 * whether meets results
 * ---------------------------------------------------------------------------------
 * if yes, CulvertCheck will become 1 to stop try, and keep CulvertCondition based on 
 * results
 *
 */ 
void CheckCulvertCondition(REAL *h, gridT *grid, int myproc)
{  
   int i,iptr,j;
   REAL sum;
   sum=0;
   for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
     i = grid->cellp[iptr];
     if(Culvertcondition[i]<=Culverttop[i]){
       if(h[i]>Culverttop[i])
         Culvertcondition[i]=-(h[i]-Culverttop[i])/(grid->dv[i]+Culverttop[i]);
       else
         Culvertcondition[i]=0;
     } else {
       if(h[i]>Culverttop[i])
         Culvertcondition[i]=0;
       else
         Culvertcondition[i]=-(h[i]-Culverttop[i])/(grid->dv[i]+Culverttop[i]);
     }
   }
}

/*
 * Function: SetCulvertDragCoefficient 
 * usage: SetCulvertDragCoefficient(gridT grid, physT phys, int myproc)
 * ----------------------------------------------------------------------
 * set cdb and cdt FOR culvert edge, if CulvertCodition(grad(edge))=1 set 
 * cdT and cdB as Cdculvert, otherwise only cdB will be changed to cdculvert
 */
void SetCulvertDragCoefficient(gridT *grid, physT *phys, int myproc)
{
  int i,ne1,ne2;
  REAL cdculvert;
  cdculvert=MPI_GetValue(DATAFILE,"cdculvert","SetCulvertDragCoefficient",myproc);
  for(i=0;i<grid->Ne;i++){
    ne1 = grid->grad[2*i];
    ne2 = grid->grad[2*i+1];
    if(ne1!=-1 && ne2!=-1) // not for bc edge, because the velocity at bc is not calculated but set
      if(Culverttop[ne1]!=INFTY || Culverttop[ne2]!=INFTY){
        phys->CdB[i]=cdculvert;
        if(phys->h[ne1]>=Culverttop[ne1] || phys->h[ne2]>=Culverttop[ne2]!=0)
          phys->CdT[i]=cdculvert;      
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
  int i,iptr,ne, jptr, nf,k;
  REAL sum,normal, tmp = prop->grav*pow(prop->theta*prop->dt,2);
  // calculate fcoef and coef
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];    
    if(phys->h[i]<=Culverttop[i])
      coef[i] = grid->Ac[i];
    else
      coef[i] = 0;

   for(nf=0;nf<grid->nfaces[i];nf++) 
     if(grid->neigh[i*grid->maxfaces+nf]!=-1) {
       ne = grid->face[i*grid->maxfaces+nf];
       fcoef[i*grid->maxfaces+nf]=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne];
       coef[i]+=fcoef[i*grid->maxfaces+nf];
     }
  }
}

/*
 * Function: StoreCulvertPressure
 * usage: store phys->h results with culvert pressure, change phys->h to culvertop for 
 * culvertcondition=1 cell 
 * -----------------------------------------------------------------------------------
 * no=1 to store pressure, change phys->h, no=0 change phys->h to be pressure field
 */
void StoreCulvertPressure(REAL *h, int Nc, int no, int myproc)
{ 
  int i;
  if(no==1){
    for(i=0;i<Nc;i++){
      Culvertpressure[i]=h[i]; //store new pressure field (for culvert part maybe not freesurface)
      if(h[i]>Culverttop[i] && Culverttop[i]!=INFTY)
        h[i]=Culverttop[i]; // transfer it into Culvertop for fully discharge part
    }
  } else {
    for(i=0;i<Nc;i++)
      h[i]=Culvertpressure[i];
  }  
}

/*
 * Function: CulvertIterationSource
 * usage: calculate source term for each iteration -----------------------------------------------------------------------------------
 * based on Newton-nested method 
 */
void CulvertIterationSource(gridT *grid, physT *phys, REAL theta,REAL dt, int myproc)
{ 
  int iptr,i,nf,ne,k;
  REAL sum,normal;
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    phys->htmp[i]=Culvertsource1[i];
    /*if(phys->h[i]>Culverttop[i] && Culvertpressure[i]>Culverttop[i]){
      phys->htmp[i]=Culvertsource1[i]-grid->Ac[i]*Culverttop[i]; 
    } 
    if(phys->h[i]>Culverttop[i] && Culvertpressure[i]<=Culverttop[i]){
      phys->htmp[i]=Culvertsource1[i]-grid->Ac[i]*Culvertpressure[i]; 
      sum = 0;
      // add back explicit part
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        normal = grid->normal[i*grid->maxfaces+nf];
        for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
          sum+=(1-theta)*phys->u[ne][k]*grid->df[ne]*normal*grid->dzf[ne][k];
      } 
      phys->htmp[i]+=dt*sum;
    }*/

    if(phys->h[i]>Culverttop[i] && Culvertpressure[i]>Culverttop[i]){
      phys->htmp[i]=Culvertsource1[i]-grid->Ac[i]*Culverttop[i]; 
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
    if(phys->h[i]>Culverttop[i] && Culvertpressure[i]<=Culverttop[i]){
      phys->htmp[i]-=grid->Ac[i]*Culverttop[i];
    }   
  }
}

/*
 * Function: CulvertInitIteration
 * usage: (restore htmp2 htmp3 hold and store h guess) or (store htmp2 htmp3 hold) 
 -----------------------------------------------------------------------------------
 * no=-1 store no=1 restore 
 */
void CulvertInitIteration(gridT *grid, physT *phys, int no, int myproc)
{ 
  int iptr,i;
  if(no==-1){
    // store each guess
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      Culvertcondition[i]=phys->h[i];
      phys->hold[i]=Culvertsource2[i];
      phys->htmp2[i]=Culvertsource3[i];
      phys->htmp3[i]=Culvertsource4[i];
    }
  } else {
    // store each guess
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      if(phys->h[i]>Culverttop[i]){
        phys->htmp[i]=phys->htmp[i]-grid->Ac[i]*phys->h[i]+grid->Ac[i]*Culverttop[i];
      }    
      Culvertsource1[i]=phys->htmp[i];
      Culvertsource2[i]=phys->hold[i];
      Culvertsource3[i]=phys->htmp2[i];
      Culvertsource4[i]=phys->htmp3[i];
    }
  }
}
