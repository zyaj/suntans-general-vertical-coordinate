/*
 * File: subgrid.h
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for subgrid.c.
 *
 */
#ifndef _subgrid_h
#define _subgrid_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _subgridT {
REAL *xp, // the x coordinate for subgrid points [Nc*(N+1)(N+2)/2*(maxfaces-2)] in each cell
     *yp, // the y coordinate for subgrid points [Nc*(N+1)(N+2)/2*(maxfaces-2)]
     *dp, // depth for each sub point [Nc*(N+1)(N+2)/2*(maxfaces-2)]
     *xpe, // the x coordinate for subgrid points at each edge [Ne*(N+1)]
     *dpe, // the y coordinate for subgrid points at each edge [Ne*(N+1)]
     *ype, // depth for subgrid points at each edge [Ne*(N+1)]
     *hprof, // the discrete h value for effective area/Volume profile [Nc*(disN+1)]
     *hprofe, // the discrete h value for effective flux height profile [Ne*(disN+1)]
     *hmin, // store the minimum h value for subgrid cell [Nc];
     *Vprof, // the total volume profile for each cell [Nc*(disN+1)]
     *Acprof, // the effective wet area profile for each cell [Nc*(disN+1)]
     *Wetperiprof, // the wetted parameter for the last layer [Ne*(disN+1)]
     *Acbackup, // store the original cell Area [Nc]
     *Aceff,  // store the top wet area at h=phys->h[nc] to use when use iteration method 
     *Aceffold,
     *Heff,   // average total depth for a cell (veff/aceff)
     *Heffold,   // average total depth for a cell     
     **Acceff, // the effective wet area for each layer center of each cell [Nc][Nk]
     **Acceffold, // the formal time step effective wet area for each layer center of each cell [Nc][Nk]
     **Acveff, // the formal time step effective wet area for each layer top and bottom center of each cell [Nc][Nk+1]
     **Acveffold, // the effective wet area for each layer center of each cell [Nc][Nk] 
     *Acwet, // the wet area for specified h for each cell [Nc]
     *Veff, // effective volume for each cell [Nc]
     *Veffold, // store the effective volume for the last time step
     *dzboteff,     // effective height for bottom layer to calulcate drag coefficient [Ne] 
     *fluxhprof; // the effective flux height profile for each cell [Ne*disN]

int *cellp, // index of sub points for each subgrid [(Nc*N^2+Nquad*N^2)*3] (Nquad number of quad grid)
    segN, // the number of segments of each edge to generate subgrid 
    disN, // the number of discrete points to calculate effective Ac, cell Volume  
    meth, // 1 for regular subgrid method with edge segment, 2 for simple subgrid model with sub cell center depth 3 for user defined function to define V, Ac and fluxh profile
    dpint, // 1 for interpolation, 2 for user defined function
    dzfmeth; // the method to calculate h for flux height calculation*/
} subgridT;

// Globally allocate the pointer to the subgrid structure
subgridT *subgrid;

void SubgridBasic(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void UpdateSubgridVeff(gridT *grid, physT *phys, propT *prop, int myproc);
void StoreSubgridVeff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridAceff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridVerticalAceff(gridT *grid, physT *phys, propT *prop,int option, int myproc);
void UpdateSubgridFluxHeight(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridHeff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridDZ(gridT *grid, physT *phys, propT *prop, int myproc);
void StoreSubgridOldAceffandVeff(gridT *grid, int myproc);
void SubgridCulverttopArea(gridT *grid, propT *prop, int myproc);
void CalculateSubgridActiveErosion(gridT *grid, physT *phys, propT *prop, int myproc);
#endif

