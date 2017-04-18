/*
 * File: initialization.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for initialization.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _initialization_h
#define _initialization_h

int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc);
REAL ReturnDepth(REAL x, REAL y);
REAL ReturnFreeSurface(REAL x, REAL y, REAL d);
REAL ReturnSalinity(REAL x, REAL y, REAL z);
REAL ReturnTemperature(REAL x, REAL y, REAL z, REAL depth);
REAL IsoReturnSalinity(REAL x, REAL y, REAL z,int i,int k);
REAL IsoReturnTemperature(REAL x, REAL y, REAL z, REAL depth,int i,int k);
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z);
REAL ReturnSediment(REAL x, REAL y, REAL z, int sizeno);
REAL ReturnBedSedimentRatio(REAL x, REAL y, int layer, int sizeno,int nsize);
REAL ReturnMarshHeight(REAL x, REAL y);
REAL ReturnMarshDragCoefficient(REAL x, REAL y);
REAL ReturnCulvertTop(REAL x, REAL y, REAL d);
REAL ReturnWindSpeed(REAL x, REAL y);
REAL ReturnWindDirection(REAL x, REAL y);
REAL ReturnSubgridPointDepth(REAL x, REAL y, REAL xv, REAL yv);
REAL ReturnSubgridPointeDepth(REAL x, REAL y);
REAL ReturnSubCellArea(REAL x1, REAL y1, REAL x2, REAL y2, REAL x3, REAL y3, REAL h);
REAL ReturnFluxHeight(REAL x1,REAL y1, REAL x2, REAL y2, REAL h);   
REAL ReturnSubgridErosionParameterizationEpslon(REAL x, REAL y);
#endif
