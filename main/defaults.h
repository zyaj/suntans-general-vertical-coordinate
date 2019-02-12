/*
 * File: defaults.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains default values of variables that may not be defined in suntans.dat.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "phys.h"


/* prettyplot:
 * uses quadratic interpolation for output values if 1, otherwise uses whichever interp method specified
 * use this by default to get better approximations for values on skewed grids
*/
const int prettyplot_DEFAULT=0;

/* linearFS:
   default value is for a nonlinear free surface, linearFS=0
*/
const int linearFS_DEFAULT=0;

/* interp:
   default value uses classic SUNTANS Perot method (type 0) vs Quadratic (type 1)
*/
const int interp_DEFAULT=0;

/* kinter
   use o.klepsova interpolation method (type 1) vs original perot interpolation(type 0)
*/
const int kinterp_DEFAULT=0;


/* gravity:
   default value is that of Earth's gravitational constant (SI units)
*/
const REAL grav_DEFAULT=9.81;

/* minumum_depth:
    0: Do nothing
    Positive value: Will be the minimum allowable depth.  
    Negative value: Sets the minimum depth to the depth of the upper layer.
*/
const REAL minimum_depth_DEFAULT=0;   

/* fixdzz:
   0: Do not adjust bottom cell height.
   1: Do the adjustment but assume bottom cell height must be greater than dzsmall*dz[Nkmax-1];
  -1: Do the adjustment but assume bottom cell height must be greater than dzsmall (i.e. not relative).
*/
const int fixdzz_DEFAULT=1;   

/* TVDsalt, TVDturb, TVDtemp:
   0: No TVD scheme
   1: First-order upwind (Psi(r)=0)
   2: Lax-Wendroff (Psi(r)=1)
   3: Superbee
   4: Van Leer

   Defaults for salt and temperature are Van Leer, for turbulence model use first-order upwind.
*/
const int TVDsalt_DEFAULT=4;
const int TVDtemp_DEFAULT=4;
const int TVDturb_DEFAULT=0;

/* laxWendroff:
   0: Nothing
   1: Set eddy-viscosity values when using nonlinear=2 to those dictated by the lax-wendroff
   scheme.
*/
const int laxWendroff_DEFAULT = 0;
   
/* laxWendroff_Vertical: 
   0: Do not employ Lax-Wendroff coefficient for vertical advection.
   1: Employ it.
*/
const REAL laxWendroff_Vertical_DEFAULT = 0;

/* hprecond:
   0: No preconditioner for free-surface solver
   1: Jacobi preconditioner
*/
const int hprecond_DEFAULT = 1;

/* ntoutStore:
   How often to save restart data.  If 0 then just save at the last time step.
*/
const int ntoutStore_DEFAULT = 0;

/* TVDmomentum
   TVD for advection of momentum, default is vanleer
*/
const int TVDmomentum_DEFAULT = 3;

/* conserveMomentum
   Use conservative momentum advection scheme by default.
*/
const int conserveMomentum_DEFAULT = 1;

/* thetaM
   Implicit vertical advection of horizontal momentum when thetaM>0.5.
   A value of -1 implies that the original conservative scheme is used in UPredictor().
*/
const REAL thetaM_DEFAULT = -1;

/* wetdry
   Don't do wetting and drying by default
*/
const int wetdry_DEFAULT = 0;

/* smoothbot:
   Treatment in SetFluxHeight and ComputeVelocityVector for smooth
   bottom flow when partial stepping is used
*/
const REAL smoothbot_DEFAULT = 0.0;

/* mergeArrays
   If mergeArrays=1 then merge output data into one file.  Otherwise output into separate
   files on each processor with suffix file.processor_number
*/
const int mergeArrays_DEFAULT = 1;

/* computeSediments
   Whether or not to compute sediments.  Off by default.
*/
const int computeSediments_DEFAULT = 0;

/* 
 *  Heat flux model and meteorological IO netcdf Parameters
 */
// Latitude - required by solar radiation function
const int latitude_DEFAULT = 29.0;

// 0 - no meteorological input; 1 - COARE3.0, short and longwave radiation calculated
const int metmodel_DEFAULT = 0; 

// Time offset parameter in days
const REAL toffSet_DEFAULT = 0.0;

// GMT offset parameter in hours used to correct solar radiation term
const REAL gmtoffset_DEFAULT = 0.0;

// Interpolation model. 0 - inverse distance weighting; 1 - kriging with spherical variogram
const int varmodel_DEFAULT = 1;

// variogram nugget parameter. Covariance = 1 - nugget @ distance = 0
const REAL nugget_DEFAULT = 0.1;

// variogram sill parameter. Covariance = 1 - sill @ distance = range
const REAL sill_DEFAULT = 0.9;

// variogram range parameter. Decorrelation length scale.
const REAL range_DEFAULT = 1e5;

//Output data to netcdf format (0 - binary, 1 - netcdf)
const int outputNetcdf_DEFAULT = 0;

// Number of steps to write to each netcdf file (mergeArray=1 only)
const int nstepsperncfile_DEFAULT = 999999;

// File number of first netcdf output file (mergeArray=1 only)
const int ncfilectr_DEFAULT = 0;

//Light extinction depth [m]
const REAL Lsw_DEFAULT = 2.0;

//Drag and heat flux coefficients
const REAL Cda_DEFAULT = 1.1e3;
const REAL Ch_DEFAULT = 1.4e3;
const REAL Ce_DEFAULT = 1.4e3;

//Start and base time string
const char *starttime_DEFAULT = "19900101.000000";
const char *basetime_DEFAULT =  "19900101.000000";

// NetCDF boundary conidtion default
const int netcdfBdy_DEFAULT = 0;

// Read initial condition netcdf
const int readinitialnc_DEFAULT = 0;

// Calculate Age variables
const int calcage_DEFAULT = 0;

// Age calculation method: 1 use river boundaries, 2 - internal source
const int agemethod_DEFAULT = 1;

// Calculate average quantities
const int calcaverage_DEFAULT = 0;

// Maximum number of faces 
const int maxFaces_DEFAULT = DEFAULT_NFACES;

/* Intz0T and Intz0B
   whether to read and interpolate Cd from z0tint.dat and z0bint.dat in Rundata
*/
const int Intz0T_DEFAULT = 0;

const int Intz0B_DEFAULT = 0;

/* IntCdV
   whether to read and interpolate CdV due to marsh from CdVint.dat in Rundata
*/
const int IntCdV_DEFAULT = 0;

/* Inthmarsh
   whether to read and interpolate hmarsh due to marsh from hmarshint.dat in Rundata
*/
const int Inthmarsh_DEFAULT = 0;

/* outvwgt
   whether to output vwgt for each cell 
*/
const int outvwgt_DEFAULT = 0;

/* wavemodel
   whether to consider wavemodel 
*/
const int wavemodel_DEFAULT = 0;

/* culvertmodel
   whether to consider culvertmodel 
*/
const int culvertmodel_DEFAULT = 0;

/* marshmodel
   whether to consider marsh model
*/
const int marshmodel_DEFAULT = 0;

/* subgrid
   whether to use subgrid method to calculate free surface
   casulli subgrid method
*/
const int subgrid_DEFAULT = 0;

const REAL subgrideps_DEFAULT = 1e-3;
/* im
   which implicit method to use for momentum equation 
   0 as theta method, 1 as AM2, 2 as AI2
*/
const int im_DEFAULT = 0;

/* ex:
   1: AX2 2:AB2 3:AB3
*/
const int ex_DEFAULT = 2;


/* udrag
  udrag=u^n-1/2dt*dh/dx*g
  whether to use udrag to calculate tau_b=cdb*|udrag|u^n+1
*/
const int udrag_DEFAULT = 0;

/* vertcoord
  default 1 which means z-level coordinate
*/
const int vertcoord_DEFAULT = 1;

/* modifydzf
  default 0 whether use u^im to modify dzf 
  now only works for the general vertical coordinate
*/
const int modifydzf_DEFAULT = 0;

/* dJdtmeth
  default 0 to use implicit method to resolve u/JdJdt term
  otherwise use explicit method
*/
const int dJdtmeth_DEFAULT = 0;

/* thetaT
  the time averaging coefficient for variational mesh
*/
const REAL thetaT_DEFAULT = 1;

/* output_user_var
  whether output the user defined variable
  0: no output 1: nc 2: nc_nk
*/
const int output_user_var_DEFAULT=0;

/* vertdzmin
  the minimum layerthickness when the generalized
  vertical coordinate is applied
*/
const REAL vertdzmin_DEFAULT=1e-3;

/* dzfmeth
  how to calculate flux height 1 upwind, 2 laxwendroff, 3 superbee, 4 van deer, 5 central differencing
  only use when vertcoord!=1
  default upwind
*/
const int dzfmeth_DEFAULT=1;

/* depthelev
  add constand depth for all grid cell
  default value is set as 0
*/
const REAL depthelev_DEFAULT=0.0;

/* SSCvprof
  the method to estimate the ratio between bottom SSC and depth-average SSC 
  to calculate deposition
  1 Rouse profile, 2 user defined alpha
*/
const REAL SSCvprof_DEFAULT=0.0;

/* sedilayerprop
  whether the sediment layer properties are constant for all sediment classes
  0 not constant, requires input values for each layer each sediment class
  1 only need to provide values for each layer, the properties are constant for all other sediment classes
*/
const REAL sedilayerprop_DEFAULT=1;


/* periodicbc
  whether there is periodic boundary condition in the simulation
  1: there is otherwise 0. 
  When periodicbc==1, there should be periodic_point.dat to store the pair between 
*/
const int periodicbc_DEFAULT=0;
