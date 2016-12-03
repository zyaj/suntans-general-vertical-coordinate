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
#include "met.h"
typedef struct _subgridT {
REAL *varNc, // variable [Nc] to store any thing you want to output
     *rhs, // store the right hand side for subgrid free surface solver
     *residual, // store the residual of the free surface solver
     *Verr, // store the error of the continuity equation when using free surface solver
     *xp, // the x coordinate for subgrid points [Nc*(N+1)(N+2)/2*(maxfaces-2)] in each cell
     *yp, // the y coordinate for subgrid points [Nc*(N+1)(N+2)/2*(maxfaces-2)]
     *dp, // depth for each sub point [Nc*(N+1)(N+2)/2*(maxfaces-2)]
     *dc, // depth for sub each cell [(Nc*N^2+2*Nquad*N^2)]
     *d1, // minimim depth for sub each cell [(Nc*N^2+2*Nquad*N^2)]
     *d2, // intermediate depth for sub each cell [(Nc*N^2+2*Nquad*N^2)]
     *d3, // maximum depth for sub each cell [(Nc*N^2+2*Nquad*N^2)]     
     *H_sub, // the total depth for each subcell= V_sub/A_sub [(Nc*N^2+2*Nquad*N^2)] 
     *A_sub, // the wet area for each subcell [(Nc*N^2+2*Nquad*N^2)] 
     *Asedi, // the wet area calculated for Cdmean [Nc]
     *Vsedi, // the volume calculated for Cdmean [Nc]
     //*Asedi_nv, // the wet area calculated for Cdmean [Nc] for nonvegetated part
     //*Vsedi_nv, // the volume calculated for Cdmean [Nc] for nonvegetated part
     //*Asedi_v, // the wet area calculated for Cdmean [Nc] for vegetated part
     //*Vsedi_v, // the volume calculated for Cdmean [Nc] for vegetated part
     *xpe, // the x coordinate for subgrid points at each edge [Ne*(N+1)]
     *ype, // the y coordinate for subgrid points at each edge [Ne*(N+1)]
     *hmarshpe, // vegetation height for each subedge point [Ne*(N+1)]
     *hmarshe, // vegetation height for each sub edge center [Ne*N]
     *hmarshp, // vegetation height for each subcell point [Nc*(N+1)(N+2)/2*(maxfaces-2)]
     *hmarshc, // vegetation height for each subcell [(Nc*N^2+2*Nquad*N^2)]
     *cdvpe, // vegetation drag coefficient for each subedge point [Ne*(N+1)]
     *cdve, // vegetation drag coefficient for each sub edge center [Ne*N]
     *cdvp, // vegetation drag coefficient for each subcell point [Nc*(N+1)(N+2)/2*(maxfaces-2)]
     *cdvc, // vegetation drag coefficient for each subcell [(Nc*N^2+2*Nquad*N^2)]
     *de, // the depth for each sub edge [Ne*N+1]
     *dpe, // depth for subgrid points at each edge [Ne*(N+1)]
     *he, // free surface height at each edge (can be upwind or center differencing)
     *hprof, // the discrete h value for effective area/Volume profile [Nc*(disN+1)]
     *hprofe, // the discrete h value for effective flux height profile [Ne*(disN+1)]
     *hmin, // store the minimum h value for subgrid cell [Nc];
     *hmax, // store the max h (-dmin) value for subgrid cell [Nc];
     *Vprof, // the total volume profile for each cell [Nc*(disN+1)]
     *Acprof, // the effective wet area profile for each cell [Nc*(disN+1)]
     *Wetperiprof, // the wetted parameter for the last layer [Ne*(disN+1)]
     *Acbackup, // store the original cell Area [Nc]
     *Aceff,  // store the top wet area at h=phys->h[nc] to use when use iteration method 
     *Acratio, // store the ratio of Ac_subpolygon/Ac_total when nf>3 for unstructured grid [Nc*(grid->maxfaces-2)]
     *Cdmean, // store the total drag coefficient calculated by mean depth and vegetation height [Ncell]
     //*Cdmean_nv, // store the total drag coefficient calculated by mean depth and vegetation height [Ncell] for nonvegetated part
     //*Cdmean_v, // store the total drag coefficient calculated by mean depth and vegetation height [Ncell] for vegetated part
     *taubsub,
     *Cdratio, // store the ratio between exact bottom drag coefficient at sub cell and Cdmean [(maxfaces-2)*segN*segN]
     //*epslon, // store eplson value for each cell for sediment transport parameterization [cell]
     //*stmp1_x, // store the explicit (time step n) term in the momentum conservation for each cell [Nc] including horizontal and vertical momentum advection, baroclinic pressure, surface pressure, coriolis 
     //*stmp2_x, // store the explicit term (time step n-1)
     //*stmp1_y, // store the explicit (time step n) term in the momentum conservation for each cell [Nc] including horizontal and vertical momentum advection, baroclinic pressure, surface pressure, coriolis 
     //*stmp2_y, // store the explicit term (time step n-1)
     //*tau_c_x, // store cell volume averaged tau_c include vegetation and bottom drag [Nc]
     //*tau_c_y, // store cell volume averaged tau_c include vegetation and bottom drag [Nc]
     **fluxp, // sum of flow flux at each cell layer positive, outflow
     **fluxn, // sum of flow flux at each cell layer negative, inflow
     *delta, // the ratio between arithmetic mean and volume average SSC in one cell [Nc]
     *Aceffold,
     *Heff,   // average total depth for a cell (veff/aceff)
     *Heffold,   // average total depth for a cell     
     **Acceff, // the effective wet area for each layer center of each cell [Nc][Nk]
     **Acceffold, // the formal time step effective wet area for each layer center of each cell [Nc][Nk]
     **Acveff, // the formal time step effective wet area for each layer top and bottom center of each cell [Nc][Nk+1] n+1
     **Acveffold, // the effective wet area for each layer center of each cell [Nc][Nk] n 
     **Acveffold2,// the effective wet area for each layer center of each cell [Nc][Nk] n-1
     *Acwet, // the wet area for specified h for each cell [Nc]
     *Veff, // effective volume for each cell [Nc]
     *Veffold, // store the effective volume for the last time step
     *dzboteff,     // effective height for bottom layer to calulcate drag coefficient [Ne] 
     *hiter, // store the iterative h for each iteration
     *hiter_min, // store the iterative h for minimum residual    
     *fluxhprof; // the effective flux height profile for each cell [Ne*disN]

REAL eps; // the convergence requirement for maximum residual

int *cellp, // index of sub points for each subgrid [(Nc*N^2+2*Nquad*N^2)*3] (Nquad number of quad grid)
    segN, // the number of segments of each edge to generate subgrid 
    disN, // the number of discrete points to calculate effective Ac, cell Volume  
    meth, // 1 for regular subgrid method with edge segment, 2 for simple subgrid model with sub cell center depth 3 for user defined function to define V, Ac and fluxh profile
    dpint, // 1 for interpolation, 2 for user defined function
    dzfmeth, // the method to calculate h for flux height calculation*/
    hmarshint, // 1 for interpolation, 2 for user defined function
    cdvint, // 1 for interpolation, 2 for user defined function
    erosionpara, // whether to use the new method to calculate sediment erosion (1) or not(0)
    dragpara; // whether to use the new method to calculate bottom shear stress

FILE *VeffFID, *AceffFID,*AsediFID,*VsediFID, *subDepositionFID,*subErosionFID;
} subgridT;

// Globally allocate the pointer to the subgrid structure
subgridT *subgrid;

void SubgridBasic(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void UpdateSubgridVeff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridFreeSurface(gridT *grid, physT *phys, propT *prop, int myproc);
void StoreSubgridVeff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridAceff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridVerticalAceff(gridT *grid, physT *phys, propT *prop,int option, int myproc);
void UpdateSubgridFluxHeight(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridHeff(gridT *grid, physT *phys, propT *prop, int myproc);
void UpdateSubgridDZ(gridT *grid, physT *phys, propT *prop, int myproc);
void CalculateSubgridDragCoef(gridT *grid, physT *phys, propT *prop);
void StoreSubgridOldAceffandVeff(gridT *grid, int myproc);
void SubgridCulverttopArea(gridT *grid, propT *prop, int myproc);
void CalculateSubgridActiveErosion(gridT *grid, physT *phys, propT *prop, int myproc);
void CalculateSubgridCdmean(gridT *grid, physT *phys, propT *prop);
void FreeSubgrid(gridT *grid, int myproc);
void CalculateCellSubgridXY(REAL *x, REAL *y, int n, int segN, gridT *grid, int myproc);
//void SubgridCellAverageTau(gridT *grid,physT *phys,propT *prop, metT *met, int myproc, int numprocs);
void OpenSubgridFiles(int computeSediments, int mergeArrays,int myproc);
void OutputSubgridVariables(gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm);
void SubgridFluxCheck(gridT *grid, physT *phys, propT *prop,int myproc);
#endif


