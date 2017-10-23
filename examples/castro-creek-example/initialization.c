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
  int k,Nkk, status;
  REAL dztop1,dztop,z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status);
  
  // set average dz distribution to which layer 
  Nkk=31;
  dztop=0.105;
  dztop1=3.5;
  if(dz!=NULL) {
    if(Nkmax>Nkk){
      dz[0]=dztop1;
      for(k=1;k<Nkk;k++)
        dz[k]=dztop;
      dz[Nkk] =(depth-(Nkk-1)*dztop-dztop1)/(Nkmax-Nkk);
      if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[0]);
      for(k=Nkk+1;k<Nkmax;k++) 
	dz[k]=dz[Nkk];
    }else{
      for(k=0;k<Nkmax;k++)
        dz[k]=depth/Nkmax;
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
  return 10;
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
  return -5;
}

/*
 * Function: ReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnSalinity(REAL x, REAL y, REAL z) 
{
  REAL x0,y0,x1,y1,x2,y2,x3,y3,dis,sx1,sx2,sx3,sx4,sx,s;
  int winter;
  winter=1; 
  //summer
  if(!winter)
  {
    sx1=1/212;
    sx2=5.0/1600;
    sx3=9.0/600;
    sx4=9.0/1000;
    x0=2026.66;
    y0=2051.04;
    x1=3254.1;
    y1=1087.2;
    y2=1842;
    x2=2159.2;
    y3=1407;
    x3=2892.4;
    // wildcat creek
    if(y<=2900 && y>1405 && x>2890){
      dis=sqrt(pow((x-x3),2)+pow((y-y3),2));
      s=23-sx4*dis;
    } else {
        // inner part of castro creek
      if(y<=y1){
        dis=sqrt(pow((x-x1),2)+pow((y-y1),2));
        s=21-dis*sx3;
      // interior point of castro cove
      } else if (y>y2 & x>x0 & x<x2 && y<2900) { 
        dis=sqrt(pow((x-x0),2)+pow((y-y0),2));
        s=27-dis*sx1;
      // interior point to site c
      } else if (y<y2 & x>x2 & y>y1){
        dis=sqrt(pow((x-x2),2)+pow((y-y2),2));
        s=26-dis*sx2;
      // outer point
      } else { 
        dis=sqrt(pow((x-x0),2)+pow((y-y0),2));
        s=27+dis*sx1/100;
      }
    };
    if(s<0)
      s=0;
    if(s>29)
      s=29;
  }else{
    // winter
    sx1=5/212;
    sx2=5.0/1600;
    sx3=9.0/600;
    sx4=9.0/1000;
    x0=2026.66;
    y0=2051.04;
    x1=3254.1;
    y1=1087.2;
    y2=1842;
    x2=2159.2;
    y3=1407;
    x3=2892.4;
    // wildcat creek
    if(y<=2900 && y>1405 && x>2890){
        dis=sqrt(pow((x-x3),2)+pow((y-y3),2));
        s=18-sx4*dis;
    } else {
      // inner part of castro creek
      if(y<=y1){
        dis=sqrt(pow((x-x1),2)+pow((y-y1),2));
        s=16-dis*sx3;
      // interior point of castro cove
      } else if (y>y2 & x>x0 & x<x2 && y<2900) { 
        dis=sqrt(pow((x-x0),2)+pow((y-y0),2));
        s=26-dis*sx1;
      // interior point to site c
      } else if (y<y2 & x>x2 & y>y1){
        dis=sqrt(pow((x-x2),2)+pow((y-y2),2));
        s=21-dis*sx2;
      // outer point
      } else { 
        dis=sqrt(pow((x-x0),2)+pow((y-y0),2));
        s=26+dis*sx1/100;
      }
    };
    if(s<0)
      s=0;
    if(s>28)
      s=28;
  }
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
  if(y<=550 && y>=500){
     return 0;
  } else {
     return 0;
  }
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
  return 0;
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
  //double a;
  //a=1/3;
  //printf("%f %d\n",a,Nsize);
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
  if(x<=3214.9 && x>=3196.9 && y>=1299.8 && y<=1305.7){
     sum=-4.2; 
     return sum;
  } else if(x<=3253.7 && x>=3237.1 && y>=1299.3 && y<=1304) {
     sum=-4.2;
     return sum;
  } else {
    return INFTY; // for non culvert cell return INFTY
  }
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
  if(x<=200 && x>=100)
    return 0;
  else
    return 1;
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
  if(x<=200 && x>=100)
    return 0;
  else
    return 1.6;
}

/*
 * Function: ReturnSubgridPointDepth
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * Helper function to give the depth of each point at subcell for subgrid method
 *
 */
REAL ReturnSubgridPointDepth(REAL x, REAL y,REAL xv,REAL yv)
{
  return 0;
}

/*
 * Function: ReturnSubgridPointeDepth
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * Helper function to give the depth of each point at subedge for subgrid method
 *
 */
REAL ReturnSubgridPointeDepth(REAL x, REAL y)
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
 * Function: IsoReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field under iso 
 * pycnal coordinate. Used in phys.c in the 
 * InitializePhysicalVariables function.
 *
 */
REAL IsoReturnSalinity(REAL x, REAL y, REAL z, int i,int k) {
  return ReturnSalinity(x,y,z);
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
