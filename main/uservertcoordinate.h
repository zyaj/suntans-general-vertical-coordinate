/*
 * File: vertcoordinate.h
 * Author : Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for vertcoordinate.c
 * including all the variables for the general vertical coordinate
 *
 */
#ifndef _uservertcoordinate_h
#define _uservertcoordinate_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

void UserDefinedVerticalCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
void InitializeVerticalCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
void InitializeIsopycnalCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
void InitializeSigmaCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
void InitializeVariationalCoordinate(gridT *grid, propT *prop, physT *phys, int myproc);
void MonitorFunctionForVariationalMethod(gridT *grid, propT *prop, physT *phys, int myproc, int numprocs, MPI_Comm comm);
void MonitorFunctionForAverageMethod(gridT *grid, propT *prop, physT *phys, int myproc);
#endif