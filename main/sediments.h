/*
 * File: sediments.h
 * Author : Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for sediment.c
 * including all the variables for sediments
 *
 */
#ifndef _sediments_h
#define _sediments_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _sedimentsT {
REAL ***SediC,// sediment concentration for time step n [fraction][cell][Nkmax]
     ***SediC_old,// sediment concentration for time step n-1 [fraction][cell][Nkmax]
     //***SediCbed, // sediment concentration in bed [fraction][cell][Nlayer]
     ***boundary_sediC, // sediment transport boundary condition [fraction][cell][Nkmax]
     ***Ws, // settling velocity [fraction][cell][Nkmax]
     *Ws0, // constant settling velocity [fraction] -> given in sedi.dat
     *Ds,  // sediment diameter [fraction] -> given in sedi.dat
     *kb, // bottom roughnesss [cell], kb=3*ds90
     *z0b,
     *z0r,
     *z0s, // z0x to calculate kb
     *Anglerepos, // sediment angle of repose [fraction] -> given in sedi.dat
     *Reposangle, // average angle of repos for each cell at top layer [cell]
     *Gsedi, // density of sediment/ density of water [fraction] -> given in sedi.dat
     *Seditb, // store the tb for each cell [cell]
     *Seditbmax, // store the tb max
     ***Erosion, // erosion term in first layer is for suspended sediment transport [fraction][cell][Nlayer]
     ***Erosion_old, //added
     //**Woldsedi, //vertical velocity for sediment particles [cell][Nkmax+1]
     **Wnewsedi, // calculated by SedimentVerticalVelocity
     **SediKappa_tv, // tubulent sediment diffusivity [cell][Nkmax]
     //**Erosiontotal, // total erosion for each cell each layer [cell][Nlayer]
     //**Erosiontotal_old, // store former step [cell][Nlayer]
     //**Neterosion,  // store the net erosion in Bedinterval steps [cell][Nlayer]
     **Deposition, // deposition term in first layer is for suspended sediment transport [fraction][cell]
     **alphaSSC, // represent the ration between bottom SSC and depth-averaged SSC [fraction][cell]
     **Deposition_old, //added [fraction][cell]
     **Toplayerratio, // the ratio for each fraction at top layer [fraction][cell]
     *Consolid, // consolidation rate for each layer [Nlayer]-> given in sedi.dat
     //*Bedmudratio, // the ratio of mud for each layer [Nlayer]-> given in sedi.dat
     **E0, // constant/basic erosion rate [Nsize][Nlayer]-> given in sedi.dat
     **Taue, // Erosion critical stress [Nsize][Nlayer]-> given in sedi.dat
     **Taud, // deposition critical stress [Nsize][Nlayer]-> given in sedi.dat
     *Prt,  // Prandtl number [fraction]-> now just assume 1
     //**Layermass, // Layer total mass for each layer [cell][Nlayer]
     ***Layerthickness, // Layer thickness [Nsize][cell][Nlayer] //added
     **Thicknesslayer,// the thickness for each layer [cell][Nlayer] sum of Layerthickness[i][const][const]
     // use this one in order to increase speed, may cause memory waste
     *Totalthickness,  // total bed thickness [cell]
     *Thickness, // Layerthickness [Nlayer]-> given in sedi.dat
     *Softhard,  // the layer is soft or hard-> given in sedi.dat [Layer] 
     **Drydensity, // dry density for each layer, constant [fraction][Nlayer]
     ds50,
     Kbed, // diffusion coefficient for bed layerthickness
     Chind, // concentration for hindered settling velocity
     Cfloc; // concentration for flocuation for settling velocity

int *Layertop,  // store the first layer number for each cell [cell]
    Nlayer, // number of bed layer -> given in sedi.dat
    Nsize,  // number of fractions -> given in sedi.dat
    sscvprof, // the switch to decide the vertical profile of SSC for 2d OR 1d model, 1 rouse curve 2 constant
    WSconstant, // if 1,just use constant settling velocity, 0 otherwise -> given in sedi.dat
    ParabolKappa, // whether to use parabolic tubulent diffusivity
    bedInterval, // the interval steps to update bed change
    bedComplex,  // whether consider the possibility to flush away a whole layer
    restart, // whether use restart
    layerprop, // whether the properies (E0,taue,taud,drydensity) of each layer is constant for all sediment classes
    readSediment; // if 1, we will read sediment file as the IC for sediment Concentration, Now just support Nsizemax=1
FILE *LayerthickFID, **SedimentFID, *SeditbFID, *SeditbmaxFID;
} sedimentsT;

// Globally allocate the pointer to the sediments structure
sedimentsT *sediments;

void ComputeSediments(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
void ComputeSedimentsRestart(gridT *grid, physT *phys, propT *prop, int myproc);
void ComputeSedimentsBedRestart(gridT *grid, physT *phys, propT *prop, int myproc);
#endif


