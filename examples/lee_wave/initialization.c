#include "math.h"
#include "fileio.h"
#include "suntans.h"
#include "initialization.h"
#include "util.h"

#define sech 1/cosh
/*
 * Function: GetDZ
 * Usage: GetDZ(dz,depth,Nkmax,myproc);
 * ------------------------------------
 * Returns the vertical grid spacing in the array dz.
 *
 */
int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc) {
  int k, status;
  REAL z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status);

  if(dz!=NULL) {
    if(r==1) 
      for(k=0;k<Nkmax;k++)
	dz[k]=depth/Nkmax;
    else if(r>1 && r<=1.1) {    
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
  REAL Lhill = 1000, h=65, D=1885+0*3142+0*5027, res = 10, lambda = 1000, k = 2*PI/lambda;
  REAL xc =  1000/2;//10935;//+0*5e4+res/2;
  
  // Gaussian
  // return D - h*exp(-pow((x-xc)/Lhill,2));
   // return D - h/(1+pow((x-xc)/Lhill,2));

  //sin wave with period of 7290m... used with 12-19_sponge where res = 270m
  return Max(D - h*cos(k*(x)), Max(D - h*cos(k*(x-res/2)), D - h*cos(k*(x+res/2))));
  // return D - h*cos(k*(x));

  // Witch of Agnesi (being sure to return the maximum depth for any cell) with tanh to smooth the edges (current settings on tanh are not smnooth!)
  
  // here's one made smooth by adding in the witches next door (it is periodic, afterall!)
  // REAL distant_witches = - (h/(1+pow((x-xc - 2*xc)/Lhill,2)))- (h/(1+pow((x-xc + 2*xc)/Lhill,2)));
  // if(x-xc < 0)
  //   return D - h/(1+pow(((x-res/2)-xc)/Lhill,2)) + distant_witches;
  // else 
  //   return D - h/(1+pow(((x+res/2)-xc)/Lhill,2)) + distant_witches;

  // and one made smooth using tanh functions near the boundary. This might violate the linear solution more.
  // if(x-xc < 0)
  //   return D - h/(1+pow(((x-res/2)-xc)/Lhill,2)) * 1/2*(1+tanh((x-2*res)/res)) * 1/2*(1+tanh((2*xc-2*res-x)/res));
  // else if (x-xc > 0)
  //   return D - h/(1+pow(((x+res/2)-xc)/Lhill,2)) * 1/2*(1+tanh((x-2*res)/res)) * 1/2*(1+tanh((2*xc-2*res-x)/res));
  // else
  //   return D - (h/(1+pow((x-xc)/Lhill,2))) * 1/2*(1+tanh((x-2*res)/res)) * 1/2*(1+tanh((2*xc-2*res-x)/res));



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
  REAL L = 2.0, a = 0.0, k;
  k = 2*PI/L;
  return a/2*cos(k*x);
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
  return 0;
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
  REAL alpha=1e-4, g = 9.81, dT, c0=0.25, Fr;
    REAL U=0.2, N=0.002, D=5027, Lx = 10000, h0=1, xc=40000/2, lambda = 10000, k = 2*PI/lambda;

  dT = N*N*D/alpha/g;

  // uniform stratification
  return 20+z*dT/D;

  //gill's sin solution added to uniform stratification
  // return 20+z*dT/D + U*U*h0/g/alpha*(N*N/U/U - k*k)*cos(k*x + sqrt(N*N/U/U - k*k)*(z+D));

  // queney's solution added to uniform stratification
  // return 20+z*dT/D - N*N*h0*Lx/alpha/g / (Lx*Lx + (x-xc)*(x-xc)) *(Lx*cos(N*(z+D+h0)/U) - (x-xc)*sin(N*(z+D+h0)/U));

  // no stratification
  // return 1;
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
  REAL U=0.2, N=0.002, D=6911, Lx = 10000, h0=1, xc=40000/2, lambda = 10000, k = 2*PI/lambda;;

  // depth averaged uniform flow
  // return U*D/ReturnDepth(x,y)*n1;

  //gill's sin solution added to depth averaged uniform flow
  // return (U*D/ReturnDepth(x,y) + U*h0*sqrt(N*N/U/U - k*k)*sin(k*x + sqrt(N*N/U/U - k*k)*(z+D)))*n1;

  // queney's solution
  // return (U*D/ReturnDepth(x,y) + N*h0*Lx/(Lx*Lx + (x-xc)*(x-xc))*(Lx*sin(N*(z+D+h0)/U) + (x-xc)*cos(N*(z+D+h0)/U)))*n1;

  // return U*n2; // periodic forcing via 3D coriolis
  return 0;
}

 
REAL ReturnSediment(REAL x, REAL y, REAL z, int sizeno) {}
REAL ReturnBedSedimentRatio(REAL x, REAL y, int layer, int sizeno,int nsize) {}

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
  return 0;
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
  // REAL xc = 2e5, Lx = 3159.5, h = 20, D = 7000, res = 250;
  // // Gaussian
  // // return D - h*exp(-pow((x-xc)/Lx,2));

  // // Witch of Agnesi (being sure to return the maximum depth for any cell)
  // if(x-xc < 0)
  //   return D - h/(1+pow(((x-res/2)-xc)/Lx,2));
  // else if (x-xc > 0)
  //   return D - h/(1+pow(((x+res/2)-xc)/Lx,2));
  // else
  //   return D - h/(1+pow((x-xc)/Lx,2));
  return ReturnDepth(x,y);
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
  // REAL xc = 2e5, Lx = 3159.5, h = 20, D = 7000, res = 250;
  // // Gaussian
  // // return D - h*exp(-pow((x-xc)/Lx,2));

  // // Witch of Agnesi (being sure to return the maximum depth for any cell)
  // if(x-xc < 0)
  //   return D - h/(1+pow(((x-res/2)-xc)/Lx,2));
  // else if (x-xc > 0)
  //   return D - h/(1+pow(((x+res/2)-xc)/Lx,2));
  // else
  //   return D - h/(1+pow((x-xc)/Lx,2));
  return ReturnDepth(x,y);
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
