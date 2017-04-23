#include "math.h"
#include "fileio.h"
#include "suntans.h"
#include "initialization.h"
#include "sediments.h"
#include "wave.h"
#include "culvert.h"
#define sech 1/cosh
/*
 * Function: GetDZ
 * Usage: GetDZ(dz,depth,Nkmax,myproc);
 * ------------------------------------
 * Returns the vertical grid spacing in the array dz.
 *
 */
int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc) {
  int k, status,Nk_noiso=10;
  REAL z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status),sum=0;
  REAL a0=120;
  if(dz!=NULL) {
    if(r==1){
      for(k=0;k<Nkmax;k++)
        dz[k]=depth/Nkmax;
        
    }else if(r>1 && r<=1.1) {    
      dz[0] = depth*(r-1)/(pow(r,Nkmax)-1);
      if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[0]);
      for(k=1;k<Nkmax;k++) 
  dz[k]=r*dz[k-1];
    } else if(r>-1.1 && r<-1) {    
      r=fabs(r);
      dz[Nkmax-1] = depth*(r-1)/(pow(r,Nkmax)-1);
      if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[Nkmax-1]);
      for(k=Nkmax-2;k>=0;k--) 
  dz[k]=r*dz[k+1];
    } else {
      printf("Error in GetDZ when trying to create vertical grid:\n");
      printf("Absolute value of stretching parameter rstretch must  be in the range (1,1.1).\n");
      exit(1);
    }
  } else {
    r=fabs(r);
    if(r!=1)
      dz0 = depth*(r-1)/(pow(r,Nkmax)-1);
    else
      dz0 = depth/Nkmax;
    z = dz0;
    for(k=1;k<Nkmax;k++) {
      dz0*=r;
      z+=dz0;
      if(z>=localdepth) {
  return k;
      }
    }
  }
}
  
/*
 * Function: ReturnDepth
 * Usage: grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
 * --------------------------------------------------------
 * Helper function to create a bottom bathymetry.  Used in
 * grid.c in the GetDepth function when IntDepth is 0.
 *
 */
REAL ReturnDepth(REAL x, REAL y) {
  REAL d;
  REAL dshallow=400,ddeep=2000,x_start=400000, Lp=12000;
  d=3000-2600*exp(-(x-x_start)*(x-x_start)/2/Lp/Lp);
  return d;
}
 
 
/*
  * Function: ReturnFreeSurface
  * Usage: grid->h[n]=ReturnFreeSurface(grid->xv[n],grid->yv[n]);
  * -------------------------------------------------------------
  * Helper function to create an initial free-surface. Used
  * in phys.c in the InitializePhysicalVariables function.
  *
  */
REAL ReturnFreeSurface(REAL x, REAL y, REAL d) {
  REAL h;
  h=0;
  return h;
}

/*
 * Function: ReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnSalinity(REAL x, REAL y, REAL z) {
  REAL drho=0.00522,h0=224,p=1.15,s1=-0.00018,s2=-3.64e-8,beta=1e-3,s;
  s=drho/beta*(1-pow(exp(z/h0),p)+s1*z+s2*z*z);
  return s;
}

/*
 * Function: ReturnTemperature
 * Usage: grid->T[n]=ReturnTemperaturegrid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnTemperature(REAL x, REAL y, REAL z, REAL depth) {
  REAL h;
  return 0;
}

/*
 * Function: IsoReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field under iso 
 * pycnal coordinate. Used in phys.c in the 
 * InitializePhysicalVariables function.
 *
 */
REAL IsoReturnSalinity(REAL x, REAL y, REAL z, int i,int k) {
  REAL drho=0.00522,h0=224,p=1.15,s1=-0.00018,s2=-3.64e-8,beta=1e-3,s,zp;
  // here rho0=1024.75
    if(k==0)
      s=3.171505;
    if(k==1)
      s=5.31515;
  return s;
}

/*
 * Function: IsoReturnTemperature
 * Usage: grid->T[n]=ReturnTemperaturegrid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field under iso 
 * pycnal coordinate. Used in phys.c in the 
 * InitializePhysicalVariables function.
 *
 */
REAL IsoReturnTemperature(REAL x, REAL y, REAL z, REAL depth, int i, int k) {
  REAL h;
  return 0;
}

/*
 * Function: ReturnHorizontalVelocity
 * Usage: grid->u[n]=ReturnHorizontalVelocity(grid->xv[n],grid->yv[n],
 *                                            grid->n1[n],grid->n2[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial velocity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z) {
  REAL U=0,D=3000;
  return U*D/ReturnDepth(x,y)*n1;
}

/*
 * Function: ReturnSediment
 * Usage: SediC[Nsize][n][Nk]=ReturnSediment(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial sediment concentration field.  Used
 * in sediment.c IntitalizeSediment function
 *
 */
REAL ReturnSediment(REAL x, REAL y, REAL z, int sizeno) {
  return 0;
}

/*
 * Function: ReturnBedSedimentRatio
 * Usage: SediC[Nsize][n][Nk]=ReturnBedSedimentRatio(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial bed sediment concentration field.  Used
 * in sediment.c IntitalizeSediment function
 * the sum of ratio should be 1
 */
REAL ReturnBedSedimentRatio(REAL x, REAL y, int layer, int sizeno, int nsize) {
  return 1.00/nsize;
}

/*
 * Function: ReturnWindSpeed 
 * Usage: Uwind[n]=ReturnWindSpeed(grid->xv[n],grid->yv[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial wind velocity field.  Used
 * in wave.c InitializeWave
 *
 */
REAL ReturnWindSpeed(REAL x, REAL y) {
  return 0;
}

/*
 * Function: ReturnWindDirection
 * Usage: Winddir[n]=ReturnWindDirection(grid->xv[n],grid->yv[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial wind velocity field.  Used
 * in wave.c InitializeWave
 *
 */
REAL ReturnWindDirection(REAL x, REAL y) {
  return 0;
}

/*
 * Function: ReturnCulvertTop
 * usage: Culvertheight[i]=ReturnCulvertTop(grid->xv[i],grid->yv[i],myproc);
 * -------------------------------------------------------------------------------
 * provide the culvert size for culvert cell, for non culvert cell, assume Culvertheight[i]=EMPTY
 *
 */
REAL ReturnCulvertTop(REAL x, REAL y, REAL d)
{ 
  REAL sum;
  if(x>350 && x<650)
    return -1.2;
  else
    return INFTY;
}

/*
 * Function: ReturnMarshHeight
 * Usage: hmarsh[n]=ReturnMarshHeight(grid->xe[n],grid->ye[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial hmarsh field.  Used
 * in marsh.c Interpmarsh
 *
 */
REAL ReturnMarshHeight(REAL x, REAL y)
{
  return 0;
}
/*
 * Function: ReturnMarshHeight
 * Usage: CdV[n]=ReturnMarshDragCoefficient(grid->xe[n],grid->ye[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial CdV field.  Used
 * in marsh.c Interpmarsh
 *
 */
REAL ReturnMarshDragCoefficient(REAL x, REAL y)
{
  return 0.0;
}
/*
 * Function: ReturnSubgridPointDepth
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * Helper function to give the depth of each point for subgrid method
 *
 */
REAL ReturnSubgridPointDepth(REAL x, REAL y, REAL xv, REAL yv)
{
  return 1;
}

/*
 * Function: ReturnSubgridPointDepth
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * Helper function to give the depth of each point for subgrid method
 *
 */
REAL ReturnSubgridPointeDepth(REAL x, REAL y)
{
  return 1;
}



/*
 * Function: ReturnSubCellArea
 * Usage: calculate area for each subcell for each cell 
 * ------------------------------------------------------------
 * Helper function to give the area of each subcell for subgrid method
 *
 */
REAL ReturnSubCellArea(REAL x1, REAL y1, REAL x2, REAL y2, REAL x3, REAL y3, REAL h)
{
  return 0;
}

/*
 * Function: ReturnFluxHeight
 * Usage: calculate flux height for each subedge 
 * ------------------------------------------------------------
 * Helper function to give flux height of each subedge for subgrid method
 *
 */
REAL ReturnFluxHeight(REAL x1,REAL y1, REAL x2, REAL y2, REAL h)
{
  return 0;
}

/*
 * Function: ReturnSubgridErosionParameterizationEpslon
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * give the value of epslon when erosion parameterization is on
 *
 */
REAL ReturnSubgridErosionParameterizationEpslon(REAL x, REAL y)
{
  return 0;
}
