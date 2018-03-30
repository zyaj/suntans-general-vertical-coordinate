/*
 * File: wave.c
 * Author: Yi-Ju Chou & Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file contains physically-based functions for wave energy density model like swan.
 * The wave model is only turned on when prop->wavemodel=1
 * There are two models, one is the simple fetchmodel with the equations from Navy manual
 * another model is from Yi-Ju's codes, same as swan model which can be coupled with 
 * hydrodynamic mdoel
 *
 */

#include "grid.h"
#include "wave.h"
#include "util.h"
#include "tvd.h"
#include "phys.h"
#include "physio.h"
#include "mympi.h"
#include "memory.h"
#include "sediments.h"
#include "suntans.h"
#include "timer.h"
#include "boundaries.h"
#include "initialization.h"

static REAL InnerProduct2(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm);
static REAL InterpCgToFace(int m, int n, int j, gridT *grid);
static REAL semivariogram(REAL Cov0, REAL Dmax, REAL D);
static void OperatorN(int m, int n, REAL *x, REAL *y, gridT *grid, propT *prop);
static void ObtainWaveVelocity(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
static void ObtainEdgeUwField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
static void CenterRadiationStress(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
static void EdgeRadiationStressToFlow(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
static REAL InterpWavePropToFace(int j, REAL *value, gridT *grid);
static void ClampK(gridT *grid, int nc1, int nc2, int j, int k, int *k_nc1, int *k_nc2, int *k_j);
static void UpdateActionDensitySinkSource(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs);
static void UpdateActionDensitySpectral(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);
static void ImplicitUpdateGeographic(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs);
static void ExplicitUpdateGeographic(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs);
static REAL InterpWindShearToFace(int j, physT *phys, gridT *grid);
static void BiCGSolveN(int m0, int n0, gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void ReadWaveVariables(gridT *grid, propT *prop, int myproc, MPI_Comm comm);
static void OpenWaveFiles(int merge, int myproc);
static void OutputWaveData(gridT *grid, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
static void ObtainKrigingCoef(gridT *grid, int myproc, int numprocs);
static void WindField(propT *prop, gridT *grid, physT *phys, int myproc, int numprocs);
static void WindSurfaceShear(gridT *grid, physT *phys, propT *prop,MPI_Comm comm, int myproc);
static void RadiationStressToFlow(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);

// Fetch model
void FetchAllocateWave(gridT *grid, int myproc);
void FetchInitializeWave(gridT *grid, physT *phys, int myproc);
void FetchFreeWave(gridT *grid, int myproc);
void CalculateFetch(gridT *grid, int myproc);
void FetchCalculateHwsigTwsig(gridT *grid,physT *phys, int myproc);
void FetchCalculateWaveexcur(gridT *grid,physT *phys, int myproc);
void FetchCalculateFw(gridT *grid,int myproc);
void FetchOutputWave(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
void FetchOpenWaveFiles(int merge, int myproc);


void InitializeWaveProperties(wpropT **wprop, propT *prop, int myproc)
{
  int i, Ns;
  char xstr[BUFFERLENGTH], ystr[BUFFERLENGTH];
  char xw[5], yw[5];
  char fd[BUFFERLENGTH], ifile[BUFFERLENGTH];
  REAL m;

  *wprop = (wpropT *)SunMalloc(sizeof(wpropT), "InitializeWaveProperties");

  (*wprop)->FetchModel = MPI_GetValue(DATAFILE,"fetchmodel","InitializeWaveProperties",myproc); 

  (*wprop)->Mw = MPI_GetValue(DATAFILE,"Mw","InitializeWaveProperties",myproc); 
  (*wprop)->Nw = MPI_GetValue(DATAFILE,"Nw","InitializeWaveProperties",myproc);
  (*wprop)->wnstep = MPI_GetValue(DATAFILE,"wnstep","InitializeWaveProperties",myproc);
  (*wprop)->sgmin = MPI_GetValue(DATAFILE,"sgmin","InitializeWaveProperties",myproc);
  (*wprop)->sgmax = MPI_GetValue(DATAFILE,"sgmax","InitializeWaveProperties",myproc);
  (*wprop)->wind_dt = MPI_GetValue(DATAFILE,"wind_dt","InitializeWaveProperties",myproc);
  (*wprop)->implicit_whitecap = MPI_GetValue(DATAFILE,"implicit_whitecap","InitializeWaveProperties",myproc);
  (*wprop)->implicit_advection = MPI_GetValue(DATAFILE,"implicit_advection","InitializeWaveProperties",myproc);
  (*wprop)->wind_forcing = MPI_GetValue(DATAFILE,"wind_forcing","InitializeWaveProperties",myproc);
  (*wprop)->nstation = MPI_GetValue(DATAFILE,"nstation","InitializeWaveProperties",myproc);
  (*wprop)->tail_opt = MPI_GetValue(DATAFILE,"tail_opt","InitializeWaveProperties",myproc);
  (*wprop)->tail_pow = MPI_GetValue(DATAFILE,"tail_pow","InitializeWaveProperties",myproc);
  (*wprop)->wind_shear = MPI_GetValue(DATAFILE,"wind_shear","InitializeWaveProperties",myproc);
  (*wprop)->rad_stress = MPI_GetValue(DATAFILE,"rad_stress","InitializeWaveProperties",myproc);
  (*wprop)->form_drag = MPI_GetValue(DATAFILE,"form_drag","InitializeWaveProperties",myproc);
  (*wprop)->btm_mud = MPI_GetValue(DATAFILE,"btm_mud","InitializeWaveProperties",myproc);
  (*wprop)->btm_sedi_erosion = MPI_GetValue(DATAFILE,"btm_sedi_erosion","InitializeWaveProperties",myproc);
  (*wprop)->depth_fw_cutoff = MPI_GetValue(DATAFILE,"depth_fw_cutoff","InitializeWaveProperties",myproc);
  (*wprop)->fw_drag = MPI_GetValue(DATAFILE,"fw_drag","InitializeWaveProperties",myproc);
  (*wprop)->btm_conc = MPI_GetValue(DATAFILE,"btm_conc","InitializeWaveProperties",myproc);
  (*wprop)->btm_vis = MPI_GetValue(DATAFILE,"btm_vis","InitializeWaveProperties",myproc);
  (*wprop)->btm_mud_thickness = MPI_GetValue(DATAFILE,"btm_mud_thickness","InitializeWaveProperties",myproc);
  (*wprop)->depth_brk_cutoff = MPI_GetValue(DATAFILE,"depth_brk_cutoff","InitializeWaveProperties",myproc);
  (*wprop)->depth_brk_indx = MPI_GetValue(DATAFILE,"depth_brk_indx","InitializeWaveProperties",myproc);
  (*wprop)->NLtriad = MPI_GetValue(DATAFILE,"NLtriad","InitializeWaveProperties",myproc);
  (*wprop)->NLquad = MPI_GetValue(DATAFILE,"NLquad","InitializeWaveProperties",myproc);
  (*wprop)->BRKdepth = MPI_GetValue(DATAFILE,"BRKdepth","InitializeWaveProperties",myproc);


//The tail after sg_max; (based on SWAN)
  //N_0.01 = N_end (sg99/sgmax)^(-m)
  //sg99=sgmax*(0.01)^(-1/m)
  (*wprop)->sg99 = (*wprop)->sgmax*pow(0.01, -1/(*wprop)->tail_pow);
  //sg_tail is obtained by dividing the integral of N*sg dsg by integral of N dsg. The integral is taken
  //from sg_max to sg99.
  m = (*wprop)->tail_pow;
  (*wprop)->sgtail = (1.0-m)/(2.0-m)*(pow((*wprop)->sg99, 2.0-m)-pow((*wprop)->sgmax, 2.0-m))
    /(pow((*wprop)->sg99, 1.0-m)-pow((*wprop)->sgmax, 1.0-m));

  Ns = (*wprop)->nstation;

  (*wprop)->Nwind = floor((prop->nsteps*prop->dt)/(*wprop)->wind_dt)+1;
  (*wprop)->nwind = floor((*wprop)->wind_dt/prop->dt);
  
  (*wprop)->xw = (REAL *)SunMalloc(Ns*sizeof(REAL), "InitializeWaveProperties");
  (*wprop)->yw = (REAL *)SunMalloc(Ns*sizeof(REAL), "InitializeWaveProperties");

  for (i=1; i<= Ns; i++){
    sprintf(xstr, "xw.%d", i);
    (*wprop)->xw[i-1] = MPI_GetValue(DATAFILE,xstr,"InitializeWaveProperties",myproc);
    sprintf(ystr, "yw.%d", i);
    (*wprop)->yw[i-1] = MPI_GetValue(DATAFILE,ystr,"InitializeWaveProperties",myproc);
    
  }
  
     
}


void AllocateWaveVariables(gridT *grid, waveT **wave, propT *prop, wpropT *wprop)
{
  int flag=0, i, j, m, n; 
  int Nc=grid->Nc, Ne=grid->Ne, Mw=wprop->Mw, Nw=wprop->Nw,
      MN = Max(Mw, Nw), nstation = wprop->nstation, Nwind=wprop->Nwind;

  *wave = (waveT *)SunMalloc(sizeof(waveT), "AllocateWaveVariables");
  
  //intrinsic wave parameters
  
  (*wave)->sg = (REAL *)SunMalloc(Mw*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->dsg = (REAL *)SunMalloc(Mw*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->thtaw = (REAL *)SunMalloc(Nw*sizeof(REAL), "AllocateWaveVariables");
  //significant wave height and orbital velocity
  (*wave)->Hw =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->ub =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->fz =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->fphi =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->uscx =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->uscy =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->ab =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->wind_spfx =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->wind_spfy =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->wind_dgf =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->wind_spf =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->Etot = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->Etmp = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->kmean = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->sgmean = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->T0 = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables"); 
  (*wave)->Hs = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->T01 = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->ktail = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->sg_PM = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->Cr = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->thtamean = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->tmp = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->fw = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->Etail = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->Ntail = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->klambda = (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  //wave variables in cell centers
  (*wave)->kw = (REAL **)SunMalloc(Mw*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->cgx = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->cgy = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->cph = (REAL **)SunMalloc(Mw*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->cs = (REAL ***)SunMalloc((Mw+1)*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->ct = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->N = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->Ntmp = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->Nold = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->ssrc = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->tsrc = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  (*wave)->src = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  
  //wave velocity stored at cell centers
  (*wave)->ux =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->uy =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->uz =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->divScx =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->divScy =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");

  //wave properties stored at cell centers
  (*wave)->Uw =  (REAL **)SunMalloc(Ne*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->divSe =  (REAL **)SunMalloc(Ne*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->ab_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->sgmean_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->thtamean_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->kmean_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->kw_edge =  (REAL **)SunMalloc(Ne*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->use =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->Uwind =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  //wind variables
  (*wave)->wind_sp = (REAL **)SunMalloc(nstation*sizeof(REAL *), "AllocateWaveVariables");
  (*wave)->wind_dg = (REAL **)SunMalloc(nstation*sizeof(REAL *), "AllocateWaveVariables");

  for(i=0; i< Nc; i++){
    (*wave)->klambda[i] = (REAL *)SunMalloc(nstation*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->fz[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->uz[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->ux[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->uy[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->divScx[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->divScy[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
  }

  for(j=0; j< Ne; j++){    
    (*wave)->Uw[j] =  (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->divSe[j] =  (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->kw_edge[j] =  (REAL *)SunMalloc(Mw*sizeof(REAL), "AllocateWaveVariables");
  }

  for(i=0; i< nstation; i++){
    (*wave)->wind_sp[i] = (REAL *)SunMalloc(Nwind*sizeof(REAL), "AllocateWaveVariables"); 
    (*wave)->wind_dg[i] = (REAL *)SunMalloc(Nwind*sizeof(REAL), "AllocateWaveVariables"); 
  }

  //coordinates of wind stations


  //matrix operators
  (*wave)->a = (REAL *)SunMalloc(MN*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->b = (REAL *)SunMalloc(MN*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->c = (REAL *)SunMalloc(MN*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->sp = (REAL *)SunMalloc((Mw+1)*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->sm = (REAL *)SunMalloc((Mw+1)*sizeof(REAL), "AllocateWaveVariables");

  (*wave)->tp = (REAL *)SunMalloc(Nw*sizeof(REAL), "AllocateWaveVariables");
  (*wave)->tm = (REAL *)SunMalloc(Nw*sizeof(REAL), "AllocateWaveVariables");

  for (m=0; m<Mw; m++){
    (*wave)->kw[m] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
    (*wave)->cgx[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->cgy[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");    
    (*wave)->cs[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->ct[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->N[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->Ntmp[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->Nold[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->ssrc[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->tsrc[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->src[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    (*wave)->cph[m] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");

    for (n=0; n<Nw; n++){
      (*wave)->cgx[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->cgy[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->cs[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->ct[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->N[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->Ntmp[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->Nold[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->ssrc[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->tsrc[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      (*wave)->src[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
    }
  }
  (*wave)->cs[Mw] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
  for (n=0; n<Nw; n++){
    (*wave)->cs[Mw][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");

  }

  //wave variables in cell edges
  (*wave)->cg = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  for (m=0; m<Mw; m++){
    (*wave)->cg[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");   
    for (n=0; n<Nw; n++){
      (*wave)->cg[m][n] = (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
    }
  }
}


void FreeWaveVariables(gridT *grid, waveT *wave, propT *prop, wpropT *wprop)
{
  int i, j, m, n, Nc=grid->Nc, Ne=grid->Ne, Mw=wprop->Mw, Nw=wprop->Nw,
    nstation=wprop->nstation, Nwind=wprop->Nwind;
    
  
  for(i=0; i< nstation; i++){
    free(wave->wind_sp[i]);
    free(wave->wind_dg[i]);
  }
  for(i = 0; i< Nc; i++){
    free(wave->klambda[i]);
    free(wave->fz[i]);
    free(wave->uz[i]);
    free(wave->ux[i]);
    free(wave->uy[i]);
    free(wave->divScx[i]);
    free(wave->divScy[i]);
  }
  for(j = 0; j< Ne; j++){
    free(wave->Uw[j]);
    free(wave->divSe[j]);
    free(wave->kw_edge[j]);
  }
  free(wave->Uw);
  free(wave->divSe);
  free(wave->klambda);
  free(wave->wind_sp);
  free(wave->wind_dg);
  free(wave->Hw);
  free(wave->wind_spfx);
  free(wave->wind_spfy);
  free(wave->wind_spf);
  free(wave->wind_dgf);
  free(wave->ub);
  free(wave->fz);
  free(wave->fphi);
  free(wave->ab);
  free(wave->fw);
  free(wave->ab_edge);
  free(wave->kmean_edge);
  free(wave->sgmean_edge);
  free(wave->thtamean_edge);
  free(wave->kw_edge);

  for(m=0;m<Mw;m++){
    for(n=0;n<Nw;n++){
      free(wave->cgx[m][n]);
      free(wave->cgy[m][n]);
      free(wave->cs[m][n]);
        free(wave->ct[m][n]);
	free(wave->N[m][n]);
	free(wave->Ntmp[m][n]);
	free(wave->Nold[m][n]);
	free(wave->ssrc[m][n]);
	free(wave->tsrc[m][n]);
	free(wave->src[m][n]); 
      }
      free(wave->kw[m]);
      free(wave->cgx[m]);
      free(wave->cgy[m]);
      free(wave->cph[m]);
      free(wave->cs[m]);
      free(wave->ct[m]);
      free(wave->N[m]);
      free(wave->Ntmp[m]);
      free(wave->Nold[m]);
      free(wave->tsrc[m]);
      free(wave->ssrc[m]);
      free(wave->src[m]);

  }
  for(n=0;n<Nw;n++){
    free(wave->cs[Mw][n]);
  }
  free(wave->cs[Mw]);
  free(wave->kw);
  free(wave->cgx);
  free(wave->cgy);
  free(wave->cph);
  free(wave->cs);
  free(wave->ct);
  free(wave->N);
  free(wave->Ntmp);
  free(wave->Nold);

  for(m=0;m<Mw;m++){
    for(n=0;n<Nw;n++){
      free(wave->cg[m][n]);
    }
    free(wave->cg[m]);
  }
  free(wave->cg);
  free(wave->thtaw);
  free(wave->sg);
  free(wave->dsg);
  free(wave->a);
  free(wave->b);
  free(wave->c);
  free(wave->sp);
  free(wave->sm);
  free(wave->ssrc);
  free(wave->tp);
  free(wave->tm);
  free(wave->tsrc);
  free(wave->src);
  free(wave->Etot);
  free(wave->Etmp);
  free(wave->T0);
  free(wave->Hs);
  free(wave->T01);
  free(wave->ktail);
  free(wave->sg_PM);
  free(wave->Cr);
  free(wave->tmp);
  free(wave->kmean);
  free(wave->thtamean);
  free(wave->Etail);
  free(wave->Ntail);
  free(wave->sgmean);
  free(wave->ux);
  free(wave->uy);
  free(wave->uz);
  free(wave->divScx);
  free(wave->divScy);
  free(wave->uscx);
  free(wave->uscy);
  free(wave->use);
  free(wave->Uwind);
  free(wave);

}

void InitializeWaveVariables(gridT *grid, propT *prop, int myproc, MPI_Comm comm)
{
  int i, j, k, m, n, Nc=grid->Nc, Ne=grid->Ne, Mw=wprop->Mw, Nw=wprop->Nw, MN=Max(Mw, Nw);
  int is, nw, nstation = wprop->nstation, Nwind=wprop->Nwind;
  REAL dsw, dtw, logsg_min, logsg_max, sg_face;
  
  for(i=0;i<Nc;i++){
    for (j=0; j<nstation; j++){
      wave->klambda[i][j] = 0;
    }
    for (k=0; k<grid->Nk[i]; k++){
      wave->fz[i][k] = 0;
      wave->uz[i][k] = 0;
      wave->uy[i][k] = 0;
      wave->ux[i][k] = 0;
      wave->divScx[i][k] = 0;
      wave->divScy[i][k] = 0;
    }
    wave->fphi[i] = 0;
    wave->uscx[i] = 0;
    wave->uscy[i] = 0;

  } 

  for(j=0;j<Ne;j++){
    for (k = 0; k < grid->Nkc[j]; k++){
      wave->Uw[j][k] = 0;
      wave->divSe[j][k] = 0;
    }
    for(m=0; m<Mw; m++){
      wave->kw_edge[j][m] = 0;
    }    
    wave->ab_edge[j] = 0;
    wave->sgmean_edge[j] = 0;
    wave->thtamean_edge[j] = 0;
    wave->kmean_edge[j] = 0;
    wave->use[j] = 0;
    wave->Uwind[j] = 0;
  }

  for(m=0; m<Mw; m++){
    for(n=0; n<Nw; n++){
      for(i=0;i<Nc;i++){
	wave->cgx[m][n][i]=0;
        wave->cgy[m][n][i]=0;
        wave->cs[m][n][i]=0;
        wave->ct[m][n][i]=0;	
	wave->N[m][n][i]=0;
	wave->Ntmp[m][n][i]=0;
	wave->Nold[m][n][i]=0;
        wave->ssrc[m][n][i]=0;
        wave->tsrc[m][n][i]=0;
        wave->src[m][n][i]=0;
      }
    }
  }

  for(m=0; m<Mw; m++){
    for(i=0;i<Nc;i++){
      wave->kw[m][i]=100.0;//1.0/grid->dv[i];  //initial guess for the wave number field
      wave->cph[m][i]=0;
    }
  }


  for(n=0; n<Nw; n++){
    for(i=0;i<Nc;i++){
      wave->cs[Mw][n][i]=0;
    }
  }


  for(m=0; m<Mw; m++){
    for(n=0; n<Nw; n++){
      for(i=0;i<Ne;i++){
	wave->cg[m][n][i]=0;
      }
    }
    wave->sp[m] = 0;
    wave->sm[m] = 0;
  }
  wave->sp[Mw] = 0;
  wave->sm[Mw] = 0;
 
  for (n=0; n<Nw; n++){
    wave->tp[n] = 0;
    wave->tm[n] = 0;

  }  

  for(m=0; m<MN; m++){
    wave->a[m] = 0;
    wave->b[m] = 0;
    wave->c[m] = 0;    
  }
  
  logsg_min = log(wprop->sgmin);
  logsg_max = log(wprop->sgmax);
  dsw = (logsg_max-logsg_min)/(double)Mw;
  sg_face = wprop->sgmin;
  
  for(m=0; m<Mw; m++){
    wave->sg[m]=exp(logsg_min + 0.5*dsw + (double)m*dsw);
    wave->dsg[m] = exp(logsg_min + (double)(m+1)*dsw)
      - exp(logsg_min + (double)(m)*dsw);
    
  }

  //dsw = (wprop->sgmax-wprop->sgmin)/(double)Mw;
  //  sg_face = wprop->sgmin;
  
  //wave->sg[0] = wprop->sgmin+0.5*dsw;
  //wave->dsg[0] = dsw;
  //for(m=1; m<Mw; m++){
  //  wave->sg[m] = wave->sg[m-1]+dsw;
  //  wave->dsg[m] = dsw;    
  //}



  dtw = 2*PI/Nw;
  wave->thtaw[0] = dtw;
  for(n=1; n<Nw; n++){
    wave->thtaw[n] = wave->thtaw[n-1]+dtw;
  }

  for(i=0; i<Nc; i++){
    wave->wind_spfx[i] = 0;
    wave->wind_spfy[i] = 0;
    wave->wind_dgf[i] = 0;
    wave->wind_spf[i] = 0;
    wave->Hw[i] = 0;
    wave->ub[i] = 0;
    wave->ab[i] = 0;
    wave->T0[i] = 0;
    wave->T01[i] = 0;
    wave->Hs[i] = 0;
    wave->ktail[i] = 1000;
    wave->sg_PM[i] = 0;
    wave->Cr[i] = 0;
    wave->thtamean[i] = 0;
    wave->Etail[i] = 0;
    wave->Ntail[i] = 0;
    wave->Etot[i] = 0;
    wave->Etmp[i] = 0;
    wave->tmp[i] = 0;
    wave->kmean[i] = 0;
    wave->sgmean[i] = 0;
    wave->fw[i] = 0;
  }

  for(is=0; is < nstation; is++){
    for(nw=0; nw<Nwind; nw++){
      wave->wind_sp[is][nw] = 0;
      wave->wind_dg[is][nw] = 0;
    }
  } 
}

void InputWind(int station, propT *prop, int myproc, int numprocs)
{
  char str[BUFFERLENGTH], c, str2[BUFFERLENGTH], str0[BUFFERLENGTH];
  char tmp[BUFFERLENGTH], tmp2[BUFFERLENGTH];
  char  fd[BUFFERLENGTH], filename[BUFFERLENGTH];
  FILE *ifile,*ofile;
  int Nwind = wprop->Nwind, i, j, n;
  int Noffset;

  Noffset = floor((prop->nstart*prop->dt)/wprop->wind_dt);
  sprintf(fd, "WindFile.%d", station+1);
  MPI_GetFile(filename, DATAFILE, fd, "InputWind",myproc);
  ifile = MPI_FOpen(filename, "r", "InputWind", myproc);

  str0[0] = ' ';

  /*while (str0[0] != '-'){
    getline(&str0, "...");
  }

  c = fgetc(ifile);

  if (c == EOF){
    return;
  }

  for(n=0; n<Noffset; n++)
    getline(ifile, str0, "");
  
  c = fgetc(ifile);

  if (c == EOF){
    return;
  }*/


  for(n=0; n<Nwind; n++){
    i = 0;
    while(c != '\n' & c != EOF & c != '\t'){
      c = fgetc(ifile);
      i++;
      if (i >= 24 & i <= 29)
	str[i-24] = c;
      if (i >= 31 & i <= 37)
	str2[i-31] = c;      
    }
    if (str[3] == ' '){
      for (j = 0; j < 29-24+1; j++)
	str[j] = tmp[j];
    }

    if (str2[3] == ' '){
      for (j = 0; j < 37-31+1; j++)
	str2[j] = tmp2[j];
    }
    for (j=0; j < 29-24+1; j++)
      tmp[j] = str[j];
    for (j = 0; j < 37-31+1; j++)
      tmp2[j] = str2[j];

    wave->wind_sp[station][n] = strtod(str, (char **)NULL);
    wave->wind_dg[station][n] = strtod(str2, (char **)NULL);
    wave->wind_dg[station][n]*=PI/180;
    c = fgetc(ifile);
 
  }
  fclose(ifile);

   // output input data to check
  /*sprintf(str,"windfile.%d",station);
  ofile = MPI_FOpen(str,"w","OutputMarsh",myproc);
  for(n=0;n<Nwind;n++) {
    fprintf(ofile,"%e %e\n",wave->wind_sp[station][n],wave->wind_dg[station][n]);
  }
  fclose(ofile);*/
}
REAL NewtonToGetK(REAL init, REAL omg, REAL depth)
{
  int indx = 1;
  REAL rsd = 1, g, gp, k, k1;

  k = init;
  while(fabs(rsd) > SMALL){
    g = pow(omg, 2)-GRAV*k*tanh(k*depth);
    gp = -GRAV*tanh(k*depth)-depth*GRAV*k*pow(cosh(k*depth), -2);
    k1 = k - g/gp;
    rsd = k1-k;
    indx = indx + 1;
  }
  return k1;
}

void SourceByWind(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc;
  REAL dsw, Cd, ustar, H, A, B, coss, tmp, spd;

  dsw = (wprop->sgmax-wprop->sgmin)/(double)Mw;
  

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	ustar = pow(wave->uscx[i],2)+pow(wave->uscy[i],2);
	ustar = sqrt(ustar);
	if(ustar <= SMALL){
	  wave->sg_PM[i] = 10000.0;
	  H = 0;
	}else{
	  wave->sg_PM[i] = 0.13*GRAV/28/ustar*2*PI;
	  tmp = wave->sg[m]/wave->sg_PM[i];
	  tmp = pow(tmp, -4);
	  H = exp(-tmp);

	}
	H = 1.0;
        coss = cos(wave->thtaw[n]-wave->wind_dgf[i]);
        A = 80/(2*PI)*pow(RHOair/RHO0/GRAV*wave->cph[m][i], 2)/wave->sg[m]*pow(ustar*Max(0, coss), 4);
	B = Max(0, 5*RHOair/RHO0*(wave->wind_spf[i]/wave->cph[m][i]*coss-0.9))*wave->sg[m];
        wave->ssrc[m][n][i] = A/wave->sg[m];
	wave->src[m][n][i] = B;
      
      }
    }
  }
}


void SinkByWhitecapping(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc;
  REAL dsw, dthta, tmp, Cds, delta, p, Gamma, s_PM, s, alpha, alpha_PM;

  dsw = (wprop->sgmax-wprop->sgmin)/(double)Mw;
  dthta = 2*PI/(double)Nw;
  //  Cds = 0.000033;
  Cds = 0.0000236;
  delta = 0.0;
  p = 4;
  alpha_PM = 0.00457;
  s_PM = sqrt(0.00302);


  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	s = wave->kmean[i]*sqrt(wave->Etot[i]);
	if (wave->kmean[i] != 0){
	  //Gamma = Cds*((1-delta)+delta*wave->kw[m][i]/wave->kmean[i])*pow(s/s_PM, p);
	  //Gamma = Cds*pow(s/s_PM, p);
	  //wave->src[m][n][i] -= Gamma*wave->sgmean[i]*wave->kw[m][i]/wave->kmean[i];//*wave->N[m][n][i];
	  alpha = wave->Etot[i]*pow(wave->sgmean[i], 4)/pow(GRAV, 2);
	  Gamma = Cds*pow(wave->sg[m]/wave->sgmean[i], 2)*pow(alpha/alpha_PM, 2);
	  //wave->src[m][n][i] -= Gamma*wave->sgmean[i]*wave->kw[m][i]/wave->kmean[i];
	  wave->src[m][n][i] -= Gamma*wave->sgmean[i];
	}
      }
    }
  }

}

void SinkByWhitecapping_implicit(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, mm, nn, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc, p;
  REAL dthta, trcDM, tmp2, J, dt = prop->dt*wprop->wnstep, Cwh, alpha_PM, Nmn, **RHS, DMmn;
  REAL alpha, Gamma, s, tmp, tmpE;

  dthta = 2*PI/(double)Nw;
  Cwh = 0.000033;
  alpha_PM = 0.00457;

  RHS = (REAL **)SunMalloc(Mw*sizeof(REAL *), "SinkByWhitecapping_implicit");
  for (m = 0; m < Mw; m++)
    RHS[m] = (REAL *)SunMalloc(Nw*sizeof(REAL), "SinkByWhitecapping_implicit");

  //First update all the sink/source except the whitecapping sink  
  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      for (i= 0; i < Nc; i++){        
	wave->Ntmp[m][n][i] = wave->N[m][n][i] + wave->ssrc[m][n][i]*dt*wprop->wnstep;
	wave->Ntmp[m][n][i] = wave->Ntmp[m][n][i]*exp(wave->src[m][n][i]*dt*wprop->wnstep);
      }
    }
  }
  //Second, update the intermediate part of the whitecapping sink.
  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	wave->Nold[m][n][i] = wave->Ntmp[m][n][i];
	s = wave->kmean[i]*sqrt(wave->Etot[i]);
	if (wave->kmean[i] != 0){
	  alpha = wave->Etot[i]*pow(wave->sgmean[i], 4)/pow(GRAV, 2);
	  Gamma = Cwh*pow(wave->sg[m]/wave->sgmean[i], 2)*pow(alpha/alpha_PM, 2);
	  wave->src[m][n][i] = -0.5*Gamma*wave->sgmean[i];
	  wave->Ntmp[m][n][i] = wave->Ntmp[m][n][i]*exp(wave->src[m][n][i]*dt*wprop->wnstep);
	}	
      }
    }
  }

  //Calculate the intermediate total energy
  for(i=0; i< Nc; i++){
    tmpE = 0;
    for(m = 0; m < Mw;m++){
      tmp = 0;
      for(n=0; n< Nw; n++){
      	if(wave->Ntmp[m][n][i] < 0) wave->Ntmp[m][n][i] = 0;
	//E(sigma) = INTEGRAL E(sigma, theta)*dtheta = SUM N*sigma*dtheta
      	tmp += wave->Ntmp[m][n][i]*dthta;
      }
      //E = INTEGRAL E(sigma)*dsigma
      tmpE+=tmp*wave->sg[m]*wave->dsg[m];      
    }
    wave->Etmp[i] = tmpE;
  }

  //Obtain the inverse of the matrix at LHS:
  //The semi-implicit scheme for whitecapping becomes (I + A)*N = RHS, then N = (I + A)^-1*RHS, where A is    
  //the matrix for the manipulation that involves Etot = SUM(N*sg*dsg), i.e. A*N = J*DIA(Nstar)*SUM(N*sg*dsg), 
  //in which N here is an M*N vector and DIA(N) is an M*N x M*N diagonal matrix that has Nstar as each of its
  //diagonal element. Due to the structure of the singular matrix A, the inverse of I+A can be obtained by
  // (I + A)^-1 = I - 1/[trace(A)+1]*A. 
  for(i = 0; i<Nc; i++){
    trcDM = 0;
    //J is the parameterization for whitecapping that is a function of global wave
    //variables, sgmean and Etot (Koman et al., 1984).
    J = Cwh*pow(wave->sgmean[i], 7)*pow(GRAV, -4)*pow(alpha_PM, -2)*wave->Etmp[i];
    //First, find the trace of A, stored in "tmp", and the RHS
    for (m=0; m < Mw; m++){
      for (n=0; n < Nw; n++){
	Nmn = wave->Nold[m][n][i];
	//Nmn = wave->Ntmp[m][n][i];
	trcDM += pow(wave->sg[m], 3)*Nmn*wave->dsg[m]*dthta;
	RHS[m][n] = Nmn - 0.25*dt*J*pow(wave->sg[m], 2)*wave->Etot[i]*Nmn;
      }
    }    

    //Second, update N = (I+A)^-1*RHS = RHS + 1/[trace(A)+1]*(-A)*RHS
    for (m=0; m < Mw; m++){
      for (n=0; n < Nw; n++){
	tmp2 = 0.25*dt*J/(1.0 + 0.25*dt*J*trcDM); //Each element in -1/[trace(A0)+1]*A0
	Nmn = RHS[m][n]; //This is equivalent to multiplication of an identity matrix with RHS, stored in a temporary variable Nmn.  
	for (mm = 0; mm < Mw; mm++)
	  for (nn = 0; nn < Nw; nn++){
	    DMmn = wave->Nold[m][n][i]*pow(wave->sg[m], 2)*(dthta*wave->sg[mm]*wave->dsg[mm]);
	    Nmn -= tmp2*DMmn*RHS[mm][nn];
	  }

	   
	if (Nmn >= 0.0)	  
	  wave->Ntmp[m][n][i] = Nmn;
	//else
	  //printf("Warning!! Skip negative N at i = %d, m = %d, n = %d, N = %f, RHS = %f\n", i, m, n, Nmn, RHS[m][n]);
      }
    }
  }
  for (m=0; m < Mw; m++)
    for (n=0; n < Nw; n++)
      ISendRecvCellData2D(wave->Ntmp[m][n], grid, myproc, comm);      
  
  for (m = 0; m < Mw; m++)
    SunFree(RHS[m], Nw*sizeof(REAL), "SinkByWhitecapping_implicit");
  SunFree(RHS, Mw*sizeof(REAL), "SinkByWhitecapping_implicit");
}


void SinkByMud(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, l, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc, Nl;
  REAL dpth, kd, St_thickness, dm, gamma, cg, Dw, conc, fac, Cb, kwd, fw;

  dm = wprop->btm_mud_thickness;  
  if (prop->computeSediments)
    for(i=0;i<sediments->Nsize;i++)
      conc = sediments->Drydensity[i][0]/sediments->Gsedi[i]/sediments->Nsize*0.000001;
  else
    conc = wprop->btm_conc;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	dpth = grid->dv[i]+phys->h[i];
	kd = wave->kw[m][i]*dpth;
	if(wave->sg[m] > 0){
	  St_thickness = sqrt(2*wprop->btm_vis/wave->sg[m]);

	  if (prop->computeSediments)
	    if (wprop->btm_sedi_erosion)
	      dm = sediments->Thicknesslayer[i][0];
	  
	  if(dm < 0) dm = 0;
	  dm /= St_thickness;
	  if(dm < 4)
	    fac = (sinh(dm)*cosh(dm)-sin(dm)*cos(dm))/
	       (pow(cosh(dm), 2)*pow(cos(dm), 2)+pow(sinh(dm), 2)*pow(sin(dm), 2));
	  else
	    fac = 1.0;
	}else{
	  St_thickness = 0.0;
	  dm = 0;
	}	
	if (wprop->btm_sedi_erosion)
	  gamma = pow(1.0+conc*(sediments->ds50-1.0), -1);
	else
	  gamma = pow(1.0+conc*(2.65-1.0), -1);
	  
	if(kd > SMALL){
	  Dw = St_thickness*gamma*pow(wave->kw[m][i], 2)/(sinh(2*kd)+2*kd)*fac;
          
	  cg = sqrt(pow(wave->cgx[m][n][i], 2) + pow(wave->cgy[m][n][i], 2));
	  wave->src[m][n][i] -= 2*Dw*cg;
	}
      }
    }
  }
  
}

void SinkByDrag(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc;
  REAL Cb, fw, dpth, Cf, kwd; 

  
  Cf = 0.15;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){       
	dpth = grid->dv[i]+phys->h[i];	
	if(dpth <= wprop->depth_fw_cutoff)
	  fw = wprop->fw_drag;
	else
	  fw = wave->fw[i];	
	Cb = fw*GRAV/sqrt(2)*wave->ub[i];	
	kwd = wave->kw[m][i]*(dpth);
	if(kwd > SMALL){
	  wave->src[m][n][i]-=Cb*pow(wave->sg[m], 2)*pow(GRAV*sinh(kwd), -2);
	}else{
	  wave->src[m][n][i] = 0.0;
	}
      }
    }
  }
  
  
}

void SinkByBreaking(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, Nc=grid->Nc, Nw=wprop->Nw, Mw=wprop->Mw;
  REAL B, Qb, Hrms, Dtot, gamma, nb, dpth, beta, Q0, tmp, alphaBJ, Hmax;
  
  B = 0.5;
  gamma = 0.75;
  nb = 1;
  alphaBJ = 1;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	dpth = grid->dv[i]+phys->h[i];
	if (dpth > 0){
	  Hrms = wave->Hs[i]/4.0;
	  Hmax = dpth*gamma;
	  beta = Hrms/Hmax;
	  if (beta > 0.2) {
	    if (beta < 1) {
	      if (beta < 0.5) Q0 = 0;
	      else Q0 = pow(2*beta-1, 2);
	      tmp = exp(Q0-1)/pow(beta, 2);
	      Qb = -(Q0 - pow(beta, 2)*(Q0-tmp)/(pow(beta, 2)-tmp));  
	      if (Qb > 1) Qb = 1;
	    }else Qb = 1;
	    Dtot = alphaBJ*Qb*wave->sgmean[i]*pow(Hmax, 2)/(8*PI);
	    wave->src[m][n][i]-=alphaBJ*Qb*wave->sgmean[i]/pow(beta, 2)/(8*PI);
	  }	  
	}
	  
      }
    }
  }
}


void SourceByTriad(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc;
  int indx2, indx3, indx1;
  REAL sg1, sg2, Ur, dpth, alpha_EB, J1, J2, dsw, beta;
  REAL sgmax = wprop->sgmax, sgmin =wprop->sgmin;
  REAL Splus, Sminus, cg, tmp, tmp0;
  REAL E, E05, E2, G, k2;



  dsw = (log(sgmax)-log(sgmin))/(double)Mw;
  alpha_EB = 0.2;

  for (m=0; m < Mw; m++){
    
    indx1 = m + floor(log(0.5)/dsw);
    indx2 = m + floor(log(2)/dsw) + 1;
        
    sg1 = wave->sg[m]/2.0;
    sg2 = wave->sg[m]*2.0;

    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	dpth = grid->dv[i]+phys->h[i];
	if (dpth > 0.1){
	  Ur = GRAV/(8*sqrt(2)*pow(PI, 2))*wave->Hs[i]*pow(wave->T01[i], 2)/pow(dpth, 2);
	  if (Ur >=0 && Ur <= 1){
	    if (Ur < 0.02)
	      beta = 0;
	    else
	      beta = -PI/2 + PI/2*tanh(0.2/Ur);
	    
	    E = wave->N[m][n][i]*wave->sg[m];
	    if (indx1 < 0){
	      E05 = 0;
	      J1 = 0;
	    }else{
	      E05 = wave->N[indx1][n][i]*sg1;
	      J1 = pow(wave->kw[indx1][i], 2)*(GRAV*dpth+2*pow(wave->cph[indx1][i], 2))
	      /(wave->kw[m][i]*dpth*(GRAV*dpth + 2/15*GRAV*pow(dpth, 3)*pow(wave->kw[m][i], 2) - 2/5*pow(wave->sg[m]*dpth, 2)));
	    }
	    if (indx2 >= Mw){
	      E2 = wave->N[Mw-1][n][i]*pow(sg2/sgmax, -wprop->tail_pow)*sg2;
	      G = pow(sg2, 2)*dpth/GRAV;
	      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
	      tmp = G + 1/tmp;
	      k2 = sg2*sqrt(tmp/(GRAV*dpth));
	      J2 = pow(wave->kw[m][i], 2)*(GRAV*dpth+2*pow(wave->cph[m][i], 2))
	      	/(k2*dpth*(GRAV*dpth + 2/15*GRAV*pow(dpth, 3)*pow(k2, 2) - 2/5*pow(sg2*dpth, 2)));
	    }else{
	      E2 = wave->N[indx2][n][i]*sg2;
	      k2 = wave->kw[indx2][i];
	      J2 = pow(wave->kw[m][i], 2)*(GRAV*dpth+2*pow(wave->cph[m][i], 2))
		/(k2*dpth*(GRAV*dpth + 2/15*GRAV*pow(dpth, 3)*pow(k2, 2) - 2/5*pow(sg2*dpth, 2)));
	    }
	    
	    
	    cg = sqrt(pow(wave->cgx[m][n][i], 2) + pow(wave->cgy[m][n][i], 2));
	    tmp0 = alpha_EB*2*PI*wave->cph[m][i]*cg*pow(J1, 2)*fabs(sin(beta));
	    tmp = tmp0*(pow(E05, 2)-2*E*E05); 
	    if (tmp > 0) Splus = tmp;
	    else Splus = 0.0;
	    tmp0 = alpha_EB*2*PI*wave->cph[m][i]*cg*pow(J2, 2)*fabs(sin(beta));
	    tmp = tmp0*(pow(E, 2)-2*E*E2);
	    if (tmp > 0) Sminus = -2*tmp;
	    else Sminus = 0.0; 
	    wave->ssrc[m][n][i] += (Sminus + Splus)/wave->sg[m];
	  }
	    
	}
	
      }
    }
  }
}
void SourceByQuad(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc;
  int nf, nb;
  REAL dsw, lambda = 0.25, sg_pls, Cnl4, sgmax = wprop->sgmax, sgmin =wprop->sgmin;
  REAL E1, E1pls, E1mns, E2, E2pls, E2mns, E3, E3pls, E3mns;
  REAL R, Csh1 = 5.5, Csh2 = 6/7, Csh3 = -1.25, kp, dpth;
  REAL ang1, ang2, ang3, ang4, dEdsg, dEdthta, dtw, coef;
  REAL S1, S2, S3, Snl4_1, Snl4_2, Snl4;  

  Cnl4 = 3*10000000;
  dsw = (log(sgmax)-log(sgmin))/(double)Mw;    
  dtw = 2*PI/Nw;
  ang1 = -11.5/180*PI;
  ang2 =  33.6/180*PI;
  ang3 =  11.5/180*PI;
  ang4 = -33.6/180*PI;


  for (m=0; m < Mw; m++){
    sg_pls = wave->sg[m]*(1.0+lambda);
    for (n=0; n < Nw; n++){      
      for(i = 0; i<Nc; i++){
	E1 = wave->N[m][n][i]*wave->sg[m];
	if (m == Mw-1){
	  sgmax = wave->sg[m];
	  E2 = wave->N[m][n][i]*pow(sg_pls/sgmax, -wprop->tail_pow)*sg_pls;
	  dEdsg = (wave->N[m][n][i]*wave->sg[m] - wave->N[m-1][n][i]*wave->sg[m-1])
	    /(wave->sg[m]-wave->sg[m-1]);
	  E3 = wave->N[m][n][i]*wave->sg[m] - dEdsg*lambda*wave->sg[m];
	}else if (m == 0){
          dEdsg = (wave->N[m+1][n][i]*wave->sg[m+1] - wave->N[m][n][i]*wave->sg[m])
	    /(wave->sg[m+1]-wave->sg[m]);
	  E2 = wave->N[m][n][i]*wave->sg[m] + dEdsg*lambda*wave->sg[m];
	  E3 = 0;
	}else{
	  dEdsg = (wave->N[m+1][n][i]*wave->sg[m+1] - wave->N[m-1][n][i]*wave->sg[m-1])
	    /(wave->sg[m+1]-wave->sg[m-1]);
	  E2 = wave->N[m][n][i]*wave->sg[m] + dEdsg*lambda*wave->sg[m];
	  E3 = wave->N[m][n][i]*wave->sg[m] - dEdsg*lambda*wave->sg[m];
	}
	if (n == Nw-1){
	  nf = 0;
	  nb = n-1;
	}else if (n == 0){
	  nf = n+1;
	  nb = Nw-1;
	}else{
	  nf = n+1;
	  nb = n-1;
	}
	dEdthta = (wave->N[m][nf][i] - wave->N[m][nb][i])/(2.0*dtw);

	E1pls = E2 + dEdthta*ang1;
	E1mns = E3 + dEdthta*ang2;

	E2pls = E1 + dEdsg*(2.0*lambda + pow(lambda, 2.0))*wave->sg[m]
	  + dEdthta*ang1;
	E2mns = E1 + dEdsg*(-pow(lambda, 2.0))*wave->sg[m] + dEdthta*ang2;

	E3pls = E1 + dEdsg*(-pow(lambda, 2.0))*wave->sg[m] + dEdthta*ang1;
	E3mns = E1 + dEdsg*(-2.0*lambda + pow(lambda, 2.0))*wave->sg[m]
	  + dEdthta*ang2;
	
	coef = Cnl4*pow(2.0*PI, 2)*pow(GRAV, -4)*pow(wave->sg[m]/(2*PI), 11);
	S1 = coef*(pow(E1, 2.0)*(E1pls*pow(1.0+lambda, -4) + E1mns*pow(1.0-lambda, -4))
		   -2*E1*E1pls*E1mns*pow(1.0-lambda*lambda, -4));
	S2 = coef*(pow(E2, 2.0)*(E2pls*pow(1.0+lambda, -4) + E2mns*pow(1.0-lambda, -4))
		   -2*E2*E2pls*E2mns*pow(1.0-lambda*lambda, -4));
	S3 = coef*(pow(E3, 2.0)*(E3pls*pow(1.0+lambda, -4) + E3mns*pow(1.0-lambda, -4))
		   -2*E3*E3pls*E3mns*pow(1.0-lambda*lambda, -4));
	Snl4_1 = 2.0*S1 - S2 - S3;

	E1pls = E2 + dEdthta*ang3;
	E1mns = E3 + dEdthta*ang4;

	E2pls = E1 + dEdsg*(2.0*lambda + pow(lambda, 2.0))*wave->sg[m]
	  + dEdthta*ang3;
	E2mns = E1 + dEdsg*(-pow(lambda, 2))*wave->sg[m] + dEdthta*ang4;


	E3pls = E1 + dEdsg*(-pow(lambda, 2.0))*wave->sg[m] + dEdthta*ang3;
	E3mns = E1 + dEdsg*(-2.0*lambda + pow(lambda, 2.0))*wave->sg[m]
	  + dEdthta*ang4;

	S1 = coef*(pow(E1, 2.0)*(E1pls*pow(1.0+lambda, -4) + E1mns*pow(1.0-lambda, -4))
		   -2*E1*E1pls*E1mns*pow(1.0-lambda*lambda, -4));
	S2 = coef*(pow(E2, 2.0)*(E2pls*pow(1.0+lambda, -4) + E2mns*pow(1.0-lambda, -4))
		   -2*E2*E2pls*E2mns*pow(1.0-lambda*lambda, -4));
	S3 = coef*(pow(E3, 2.0)*(E3pls*pow(1.0+lambda, -4) + E3mns*pow(1.0-lambda, -4))
		   -2*E3*E3pls*E3mns*pow(1.0-lambda*lambda, -4));

	Snl4_2 = 2.0*S1 - S2 - S3;

	Snl4 = Snl4_1 + Snl4_2;

	kp = 0.75*wave->kmean[i];
        if(kp <= 0.5) kp = 0.5;
	dpth = phys->h[i] + grid->dv[i];
	if (dpth <= 0)
	  R = 0;
	else	  
	  R = 1+Csh1/(kp*dpth)*(1-Csh2*kp*dpth)*exp(Csh3*kp*dpth);

	wave->ssrc[m][n][i]+=R*Snl4/wave->sg[m];

      }
    }
  }

}


void ObtainTotalEnergy(gridT *grid, physT *phys, propT *prop,MPI_Comm comm, int myproc){
  
  int i, m, n, Nc=grid->Nc, Nw=wprop->Nw, Mw=wprop->Mw;
  REAL dsw, dthta, tmp, Nend, mm, dpth, max_wh, bf;
  REAL tmpE, tmpsg, tmpEsg;

  dsw = (wprop->sgmax-wprop->sgmin)/(double)Mw;
  dthta = 2*PI/(double)Nw;
  mm = wprop->tail_pow;

  for(i=0; i< Nc; i++){
    tmpE = 0;
    tmpEsg = 0;
    Nend = 0;
    for(m = 0; m < Mw;m++){
      tmp = 0;
      tmpsg = 0;
      for(n=0; n< Nw; n++){
      	if(wave->N[m][n][i] < 0) wave->N[m][n][i] = 0;
	//E(sigma) = INTEGRAL E(sigma, theta)*dtheta = SUM N*sigma*dtheta
      	tmp += wave->N[m][n][i]*dthta;
	//Esigma(sigma) = INTEGRAL E(sigma, theta)*sigma*dtheta = SUM N*sigma*sigma*dtheta
	tmpsg += wave->N[m][n][i]*dthta;
      }
      tmpE+=tmp*wave->sg[m]*wave->dsg[m];
      //E = INTEGRAL Esigma(sigma)*dsigma
      tmpEsg +=tmp*wave->sg[m]*wave->sg[m]*wave->dsg[m];
      if (m == Mw-1)
	Nend = tmp;    
      //E = INTEGRAL E(sigma)*dsigma  
    }
    wave->Etot[i] = tmpE;
    wave->T01[i] = tmpE/tmpEsg*2*PI; //The tail is not considered while obtaining T01    

    if(wprop->tail_opt){
      wave->Etail[i] = Nend*pow(wprop->sgmax,mm)*1.0/(2.0-mm)
	*(pow(wprop->sg99, 2.0-mm)-pow(wprop->sgmax, 2.0-mm));
      wave->Ntail[i] = wave->Etail[i]/wprop->sgtail;
      wave->Etot[i]+=wave->Etail[i];
    }

    wave->Hs[i] = 4*sqrt(wave->Etot[i]);

    if (wprop->depth_brk_cutoff){
      dpth = grid->dv[i]+phys->h[i];
      max_wh = wprop->depth_brk_indx*dpth;
      bf = tmpE/pow(max_wh, 2.0);
      if (bf > 1.0){
	wave->Etot[i] = pow(max_wh, 2.0);
	wave->Hs[i] = 4.0*max_wh;
	for (m = 0; m < Mw; m++){
	  for (n = 0; n < Nw; n++){
	    wave->N[m][n][i] /= bf;
	  }
	}
      }
    }
  }  
}


void ObtainMeanKSG(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Nc=grid->Nc, Nw=wprop->Nw, Mw=wprop->Mw;
  REAL dthta, tmp, invE, tmpsinh, dpth, m0, tmpthta, tail, kd;

  dthta = 2*PI/(double)Nw;

  for(i=0; i< Nc; i++){
    wave->kmean[i] = 0;
    wave->sgmean[i] = 0;
    wave->ub[i] = 0;
    wave->ab[i] = 0;
    wave->thtamean[i] = 0;
    wave->T0[i] = 0;
    dpth = grid->dv[i]+phys->h[i];  
    if (wave->Etot[i] > SMALL){
      invE = 1/wave->Etot[i];
      m0 = 0;
      for(m = 0; m < Mw;m++){
	tmp = 0;
	for(n=0; n< Nw; n++){
	  //E(sigma) = INTEGRAL E(sigma, theta) dthta
	  //= INTERGRAL N(sigma, theta)*sigma dthta
	  tmp += wave->N[m][n][i]*dthta;
	}
	wave->sgmean[i]+=tmp*wave->dsg[m];

	tmp *= wave->sg[m];
	wave->kmean[i] += tmp*wave->dsg[m]*1/sqrt(wave->kw[m][i]);
	m0 += tmp*wave->dsg[m];	

	if (dpth > 0){
	  kd = wave->kw[m][i]*(dpth);
	  tmpsinh = pow(sinh(kd), 2);
	  tmpsinh = 1/tmpsinh;
          wave->ab[i] += tmpsinh*tmp*wave->dsg[m];
	  wave->ub[i] += tmpsinh*tmp*pow(wave->sg[m], 2)*wave->dsg[m];
	}

      }

      if(wprop->tail_opt){       
	wave->sgmean[i]+=wave->Ntail[i];
	wave->kmean[i]+=wave->Etail[i]*1/sqrt(wave->ktail[i]);
	m0 += wave->Etail[i]*(wprop->sg99-wprop->sgmax);
	if (dpth > 0){
	  kd = wave->ktail[i]*(dpth);
	  tmpsinh = pow(sinh(kd), 2);
	  tmpsinh = 1/tmpsinh;
          wave->ab[i] += tmpsinh*wave->Etail[i];
	  wave->ub[i] += tmpsinh*wave->Etail[i]*pow(wprop->sgtail, 2.0);
	}	
      }
            
      wave->sgmean[i]*= invE;
      wave->sgmean[i] = 1/wave->sgmean[i];      
      
      wave->kmean[i]*= invE;
      wave->kmean[i] = pow(wave->kmean[i], -2);
      
      wave->T0[i] = m0;
      
      wave->ab[i] = wave->ab[i]*2;
      wave->ab[i] = sqrt(wave->ab[i]);

      //wave->ub[i] = wave->ub[i]*2;
      //wave->ub[i] = sqrt(wave->ub[i]);
      wave->ub[i] = wave->ab[i]*wave->sgmean[i];
      

      for(n=0; n< Nw; n++){
	tmp = 0;
	for(m = 0; m < Mw;m++){
	  tmp+=wave->N[m][n][i]*wave->sg[m]*wave->dsg[m];
	}
	wave->thtamean[i]+=tmp*wave->thtaw[n]*dthta;
      }
      wave->thtamean[i]*=invE;
    }
  }  

}


//Obtain the wave number field corresponding to each intrinsic frequency using the
//Taylor series expansion, as well as find the phase speed by k/omega. 
void ObtainKField(gridT *grid, physT *phys, propT *prop) 
{
  int i,m;
  REAL dpth, G, tmp, cph;


  for (m=0; m < wprop->Mw; m++){
    for(i = 0; i < grid->Nc; i++){
      if(phys->h[i] < -grid->dv[i])
	dpth = 0.1;
      else
	dpth = grid->dv[i]+phys->h[i];
      G = pow(wave->sg[m], 2)*dpth/GRAV;
      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
      tmp = G + 1/tmp;
      wave->cph[m][i] = 1.0/sqrt(tmp/(GRAV*dpth));
      wave->kw[m][i] = wave->sg[m]/wave->cph[m][i];
    }
  }
  if (wprop->tail_opt)
    for(i = 0; i < grid->Nc; i++){
      if(phys->h[i] < -grid->dv[i])
	dpth = 0.1;
      else
	dpth = grid->dv[i]+phys->h[i];
      G = pow(1.5*wprop->sgtail, 2)*dpth/GRAV;
      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
      tmp = G + 1/tmp;
      cph= 1.0/sqrt(tmp/(GRAV*dpth));
      wave->ktail[i] = wprop->sgtail/cph;    
    }
}
void ObtainEdgeKField(gridT *grid, physT *phys, propT *prop) 
{
  int k, j, m, Ne = grid->Ne, Mw = wprop->Mw;
  REAL dpth, G, tmp, cph;

  for(j = 0; j < Ne; j++){
    dpth=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      dpth+=grid->dzf[j][k];
    if (dpth < 0.1) dpth = 0.1;    
    for (m=0; m < Mw; m++){      
      G = pow(wave->sg[m], 2)*dpth/GRAV;
      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
      tmp = G + 1/tmp;
      cph = 1.0/sqrt(tmp/(GRAV*dpth));
      wave->kw_edge[j][m] = wave->sg[m]/cph;
    }
  } 
}

void ObtainCenterCgField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int i, m, n, Mw=wprop->Mw,Nw=wprop->Nw, Nc=grid->Nc;
  REAL tmp, depth, kw, kd;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i < Nc; i++){
	if(phys->h[i] < -grid->dv[i])
	  depth = 0.1;
	else
	  depth = phys->h[i] + grid->dv[i];
	
	kw = wave->kw[m][i];
	kd = kw*depth;
	
	if(kd > 0.1)
	  tmp = 0.5*(1+2*kd/sinh(2*kd))*wave->sg[m]/kw;
	else
	  tmp = sqrt(GRAV*depth);

	wave->cgx[m][n][i] = tmp*cos(wave->thtaw[n]) + phys->uc[i][grid->ctop[i]];
	wave->cgy[m][n][i] = tmp*sin(wave->thtaw[n]) + phys->vc[i][grid->ctop[i]];

      }
      ISendRecvCellData2D(wave->cgx[m][n], grid, myproc, comm);
      ISendRecvCellData2D(wave->cgy[m][n], grid, myproc, comm);
    }
  }
    
}

void ObtainEdgeCgField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int j, m, n, nc1, nc2, Mw=wprop->Mw,Nw=wprop->Nw, Ne=grid->Ne;
  REAL tmp, depth, kw;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(j = 0; j<Ne; j++){
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(nc1 == -1){
	  wave->cg[m][n][j] = wave->cgx[m][n][nc2]*grid->n1[j]+wave->cgy[m][n][nc2]*grid->n2[j];	  
	}else if(nc2 == -1){
	  wave->cg[m][n][j] = wave->cgx[m][n][nc1]*grid->n1[j]+wave->cgy[m][n][nc1]*grid->n2[j];
	}else{
	  wave->cg[m][n][j] = InterpCgToFace(m, n, j, grid);
	}
      }
      ISendRecvEdgeData2D(wave->cg[m][n], grid, myproc, comm);
     
    }
  }
    
}

static REAL InterpCgToFace(int m, int n, int j, gridT *grid){
  int nc1, nc2;
  REAL def1, def2, Dj, Cgnc1, Cgnc2;

  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  Dj = grid->dg[j];
  def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
  def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
  Cgnc1 = wave->cgx[m][n][nc1]*grid->n1[j] + wave->cgy[m][n][nc1]*grid->n2[j];
  Cgnc2 = wave->cgx[m][n][nc2]*grid->n1[j] + wave->cgy[m][n][nc2]*grid->n2[j];

  //  if (def1 == 0 || def2 == 0)
  // return UpWind(wave->cg[m][n][j], Cgnc1, Cgnc2);
  //else
  return (Cgnc1*def2+Cgnc2*def1)/(def1+def2);

}

void ObtainEdgeWaveProp(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int j, nc1, nc2, Ne=grid->Ne;
  REAL tmp, depth, kw;

  for(j = 0; j<Ne; j++){
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
	
    if(nc1 == -1){
      wave->ab_edge[j] = wave->ab[nc2];
      wave->kmean_edge[j] = wave->kmean[nc2];
      wave->sgmean_edge[j] = wave->sgmean[nc2];
      wave->thtamean_edge[j] = wave->thtamean[nc2];      
    }else if(nc2 == -1){
      wave->ab_edge[j] = wave->ab[nc1];
      wave->kmean_edge[j] = wave->kmean[nc1];
      wave->sgmean_edge[j] = wave->sgmean[nc1];
      wave->thtamean_edge[j] = wave->thtamean[nc1];
    }else{
      wave->ab_edge[j] = InterpWavePropToFace(j, wave->ab, grid);
      wave->kmean_edge[j] = InterpWavePropToFace(j, wave->kmean, grid);
      wave->sgmean_edge[j] = InterpWavePropToFace(j, wave->sgmean, grid);
      wave->thtamean_edge[j] = InterpWavePropToFace(j, wave->thtamean, grid);
    }
  }
  ISendRecvEdgeData2D(wave->ab_edge, grid, myproc, comm);
  ISendRecvEdgeData2D(wave->kmean_edge, grid, myproc, comm);
  ISendRecvEdgeData2D(wave->sgmean_edge, grid, myproc, comm);
  ISendRecvEdgeData2D(wave->thtamean_edge, grid, myproc, comm);

}

static REAL InterpWavePropToFace(int j, REAL *value, gridT *grid){
  int nc1, nc2;
  REAL def1, def2, Cgnc1, Cgnc2;

  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
  def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
  //  if (def1 == 0 || def2 == 0)
  // return UpWind(wave->cg[m][n][j], Cgnc1, Cgnc2);
  //else
  return (value[nc1]*def2+value[nc2]*def1)/(def1+def2);

}


void ObtainCenterCsCtField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{

  int i, m, n, nf, ne, normal, k, nc1, nc2;
  REAL Ac, df, kw, tmp0, tmp1, tmp00, tmp11, cost, sint, cost2, sint2, tanhkd;
  REAL dudx, dudy, dvdx, dvdy, dddx, dddy, depth;
  int Nc=grid->Nc, Mw=wprop->Mw, Nw=wprop->Nw;
  REAL dthta = 2*PI/(double)Nw;

  for(i = 0; i < Nc; i++){
    Ac = grid->Ac[i];
    dudx = 0;
    dudy = 0;
    dvdx = 0;
    dvdy = 0;
    dddx = 0;
    dddy = 0;
    k = grid->ctop[i];
    //loop through each face of the the cell to obtain the velocity and depth gradients at the cell center

    for(nf=0; nf<grid->nfaces[i]; nf++){
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      df = grid->df[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if (nc1 == -1) nc1 = nc2;
      if (nc2 == -1) nc2 = nc1;
      dudx += 1/Ac*0.5*(phys->uc[nc1][grid->ctop[nc1]]
		      + phys->uc[nc2][grid->ctop[nc2]])*grid->n1[ne]*normal*df;
      dudy += 1/Ac*0.5*(phys->uc[nc1][grid->ctop[nc1]]
		      + phys->uc[nc2][grid->ctop[nc2]])*grid->n2[ne]*normal*df;

      dvdx += 1/Ac*0.5*(phys->vc[nc1][grid->ctop[nc1]]
		      + phys->vc[nc2][grid->ctop[nc2]])*grid->n1[ne]*normal*df;
      dvdy += 1/Ac*0.5*(phys->vc[nc1][grid->ctop[nc1]]
		      + phys->vc[nc2][grid->ctop[nc2]])*grid->n2[ne]*normal*df;

      dddx += 1/Ac*0.5*(phys->h[nc1]+grid->dv[nc1] + 
                        phys->h[nc2]+grid->dv[nc2])*grid->n1[ne]*normal*df;
      dddy += 1/Ac*0.5*(phys->h[nc1]+grid->dv[nc1] +  
                        phys->h[nc2]+grid->dv[nc2])*grid->n2[ne]*normal*df;
            
    }
    depth = phys->h[i]+grid->dv[i];
    
    if(depth > 0.1){
      for (n = 0; n < Nw; n++){
	//Center value used for frequency velocity
	cost = cos(wave->thtaw[n]);
	sint = sin(wave->thtaw[n]);
	for (m = 0; m < Mw; m++){
	  kw = wave->kw[m][i];
	  tanhkd = tanh(kw*depth);
	  tmp1 = 0.5*sqrt(GRAV)*pow(kw, 1.5)*(1-pow(tanhkd, 2))/sqrt(tanhkd);	
	  if (m >= 1){
	    
	    tmp11 = tmp1*(phys->dhdt[i]+phys->uc[i][grid->ctop[i]]*dddx+
			  +phys->vc[i][grid->ctop[i]]*dddy)  
	    //	    tmp11 = tmp1*phys->dhdt[i]
	      - sqrt(pow(wave->cgx[m][n][i], 2)+pow(wave->cgy[m][n][i], 2))
	      * kw*(cost*(dudx*cost+dudy*sint)
		    + sint*(dvdx*cost+dvdy*sint));
	    
	    tmp00 = tmp0*(phys->dhdt[i]+phys->uc[i][grid->ctop[i]]*dddx+
			  +phys->vc[i][grid->ctop[i]]*dddy)
	      //tmp00 = tmp0*phys->dhdt[i]
		       - sqrt(pow(wave->cgx[m-1][n][i], 2)+pow(wave->cgy[m-1][n][i], 2))
	      * kw*(cost*(dudx*cost+dudy*sint)
		    + sint*(dvdx*cost+dvdy*sint));
	    wave->cs[m][n][i] = 0.5*(tmp11 + tmp00);
	  }else{
	    	    wave->cs[m][n][i] = tmp1*(phys->dhdt[i]+phys->uc[i][grid->ctop[i]]*dddx+
	    		      +phys->vc[i][grid->ctop[i]]*dddy)
	    // wave->cs[m][n][i] = tmp1*phys->dhdt[i]
	      - sqrt(pow(wave->cgx[m][n][i], 2)+pow(wave->cgy[m][n][i], 2))
	      * kw*(cost*(dudx*cost+dudy*sint)
		    + sint*(dvdx*cost+dvdy*sint));
	  }
	  tmp0 = tmp1;
	}
      }
      for (n = 0; n < Nw; n++){
	wave->cs[0][n][i] = 2*wave->cs[1][n][i] - wave->cs[2][n][i];
	wave->cs[Mw][n][i] = 2*wave->cs[Mw-1][n][i] - wave->cs[Mw-2][n][i];
      }
    }

    if(depth > 0.1){
      for (m = 0; m < Mw; m++){
	kw = wave->kw[m][i];
	tanhkd = tanh(kw*depth);
	tmp1 = 0.5*sqrt(GRAV)*pow(kw, 1.5)*(1-pow(tanhkd, 2))/sqrt(tanhkd);
	for (n = 0; n < Nw; n++){
	  //Edge value used for theta velocity
	  cost2 = cos(wave->thtaw[n]+0.5*dthta);
	  sint2 = sin(wave->thtaw[n]+0.5*dthta);
	  
	  wave->ct[m][n][i] = -1/kw*(tmp1*(-sint2*dddx + cost2*dddy) 
				     + kw*(cost2*(-dudx*sint2+dudy*cost2)
					   + sint2*(-dvdx*sint2+dvdy*cost2)));
	}
      }
    }
  }
    
}

void UpdateWave(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int blowup, int myproc, int numprocs){
  
  int i,fetchmodel,n,m;
  REAL t0,t1;
  if(prop->n == 1+prop->nstart){
    fetchmodel=MPI_GetValue(DATAFILE,"fetchmodel","UpdateWave",myproc);
    if(fetchmodel){
      wave=(waveT *)SunMalloc(sizeof(waveT),"UpdateWave");
      wprop=(wpropT *)SunMalloc(sizeof(wpropT),"UpdateWave");
      FetchOpenWaveFiles(prop->mergeArrays,myproc); 
      FetchAllocateWave(grid,myproc);
      FetchInitializeWave(grid,phys,myproc);

    } else {
      InitializeWaveProperties(&wprop, prop, myproc);
      AllocateWaveVariables(grid, &wave, prop, wprop);
      OpenWaveFiles(prop->mergeArrays,myproc);
      InitializeWaveVariables(grid, prop, myproc, comm);

      if (RESTART)
        ReadWaveVariables(grid, prop, myproc, comm);
      
      if(wprop->wind_forcing){
        ObtainKrigingCoef(grid, myproc, numprocs);
        for(i = 0; i < wprop->nstation; i++)
	  InputWind(i, prop, myproc, numprocs);
      }
    }
  }

  if(wprop->FetchModel){
    if(wprop->ConstantWind==0){
      FetchWindSpeedandDirection(grid,prop,myproc); // add new functions in Boundaries.c
      CalculateFetch(grid,myproc);
      FetchCalculateHwsigTwsig(grid,phys,myproc);
      FetchCalculateWaveexcur(grid,phys,myproc);
      FetchCalculateFw(grid,myproc);
    }

    FetchOutputWave(grid,phys,prop,myproc,numprocs,blowup,comm); 

    //if(prop->n==prop->nstart+prop->nsteps)
     //FetchFreeWave(grid,myproc);
  } else {
    if(wprop->wind_forcing){
      WindField(prop, grid, phys, myproc, numprocs);
      if(wprop->wind_shear)
        WindSurfaceShear(grid, phys, prop, comm, myproc);
    }
  
    if ((prop->n-1) % wprop->wnstep == 0){
      ObtainKField(grid, phys, prop);
      ObtainEdgeKField(grid, phys, prop); 
      ObtainCenterCgField(grid, phys, prop, comm, myproc);
      ObtainEdgeCgField(grid, phys, prop, comm, myproc);
      ObtainCenterCsCtField(grid, phys, prop, comm, myproc);
      if (wprop->wind_forcing)
        SourceByWind(grid, phys, prop, comm, myproc);
      if (wprop->form_drag)
        SinkByDrag(grid, phys, prop, comm, myproc);
      if (wprop->btm_mud && prop->n!=prop->nstart+1)
        SinkByMud(grid, phys, prop, comm, myproc);
      if (wprop->NLquad)
        SourceByQuad(grid, phys, prop, comm, myproc);
      if (wprop->NLtriad)
        SourceByTriad(grid, phys, prop, comm, myproc);
      if (wprop->BRKdepth)
        SinkByBreaking(grid, phys, prop, comm, myproc);  
      if (wprop->implicit_whitecap == 1)
        SinkByWhitecapping_implicit(grid, phys, prop, comm, myproc);
      else{
        SinkByWhitecapping(grid, phys, prop, comm, myproc);
        UpdateActionDensitySinkSource(grid, phys, prop, comm, myproc, numprocs);
      }
      if (wprop->implicit_advection == 1)
        ImplicitUpdateGeographic(grid, phys, prop, comm, myproc, numprocs);
      else
        ExplicitUpdateGeographic(grid, phys, prop, comm, myproc, numprocs); 
     
      UpdateActionDensitySpectral(grid, phys, prop, comm, myproc);
      ObtainTotalEnergy(grid, phys, prop,comm, myproc);
      ObtainMeanKSG(grid, phys, prop, comm, myproc);
      ObtainWaveVelocity(grid, phys, prop, comm, myproc);
      if(wprop->rad_stress){
        ObtainEdgeWaveProp(grid, phys, prop, comm, myproc);    
        ObtainEdgeUwField(grid, phys, prop, comm, myproc);      
        CenterRadiationStress(grid, phys, prop, comm, myproc);
        EdgeRadiationStressToFlow(grid, phys, prop, comm, myproc);
      }
    }

    //Output data
    OutputWaveData(grid, prop, myproc, numprocs, 0, comm);
  }  
}

static void UpdateActionDensitySpectral(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, Nc=grid->Nc, Mw=wprop->Mw, Nw=wprop->Nw, MN=Max(wprop->Mw, wprop->Nw);
  REAL dsw, fct, dang;
  REAL x[MN], a[MN], b[MN], c[MN], d[MN];

  dsw = (wprop->sgmax-wprop->sgmin)/(double)Mw;



  //First update accident density in the frequency space
  for(i=0; i < Nc; i++){
    for(n=0; n< Nw; n++){

      //Compute the source term for the frequency discretization
      //based on the Crank-Nicolson scheme
      for(m=0; m < Mw; m++){
	fct = prop->dt*wprop->wnstep/(2*wave->dsg[m]);
	a[m] = fct*wave->sp[m];
	b[m] = 1+(fct*wave->sm[m]-fct*wave->sp[m+1]);
	c[m] = -fct*wave->sm[m+1];
	if(m == 0)
	  d[m] = b[m]*wave->N[m][n][i]+c[m]*wave->N[m+1][n][i];
	else if(m == Mw-1){
	  d[m] = a[m]*wave->N[m-1][n][i]+b[m]*wave->N[m][n][i];
	  d[Mw-1] += c[Mw-1]*wave->N[m][n][i]*pow(1+0.5*wave->dsg[m]/wprop->sgmax, -wprop->tail_pow);
	}else
	  d[m] = a[m]*wave->N[m-1][n][i]+b[m]*wave->N[m][n][i]
                              + c[m]*wave->N[m+1][n][i];
      }

      //Update upwind velocity operators for the current time step
      for(m=0; m < Mw+1; m++){
	wave->sp[m] = 0.5*(wave->cs[m][n][i]+fabs(wave->cs[m][n][i]));
	wave->sm[m] = 0.5*(wave->cs[m][n][i]-fabs(wave->cs[m][n][i]));	  
      }

      //Generate matrix operator for the matrix solver
      for(m=0; m < Mw; m++){
	fct = prop->dt*wprop->wnstep/(2*wave->dsg[m]);
	a[m] = -fct*wave->sp[m];
	b[m] = 1-(fct*wave->sm[m]-fct*wave->sp[m+1]);
	c[m] = fct*wave->sm[m+1];
	x[m] = 0;
      }
      a[0]= 0;
      d[Mw-1] -= c[Mw-1]*wave->N[Mw-1][n][i]*pow(1+0.5*wave->dsg[m]/wprop->sgmax, -wprop->tail_pow);
      c[Mw-1] = 0;
      
      TriSolve(a, b, c, d, (&x[0]), Mw);

      //Copy to the action density field
      for(m=0; m < Mw; m++){
      	wave->N[m][n][i] = x[m];
      }
    }
  }
  

 
  dang = 360/Nw;
  fct = prop->dt*wprop->wnstep/(2*dang);

  for(i=0; i < Nc; i++){
    for(m=0; m< Mw; m++){
      for(n=0; n < Nw; n++){
	if (n == 0){
	  a[n] = fct*wave->tp[Nw-1];
	  b[n] = 1+(fct*wave->tm[Nw-1]-fct*wave->tp[n]);
	}else{
	  a[n] = fct*wave->tp[n-1];
	  b[n] = 1+(fct*wave->tm[n-1]-fct*wave->tp[n]);
	}
	c[n] = -fct*wave->tm[n];
	if(n == 0){
	  d[n] = b[n]*wave->N[m][n][i]+c[m]*wave->N[m][n+1][i];
	  d[n] += a[n]*wave->N[m][Nw-1][i]; //Source due to periodicity
  	}else if(n == Nw-1){
  	  d[n] = a[n]*wave->N[m][n-1][i]+b[m]*wave->N[m][n][i];
	  d[n] += c[n]*wave->N[m][1][i];   //Source due to periodicity
	}else
  	  d[n] = a[n]*wave->N[m][n-1][i]+b[m]*wave->N[m][n][i]
	    + c[n]*wave->N[m][n+1][i];
      }

      for(n=0; n < Nw; n++){
	wave->tp[n] = 0.5*(wave->ct[m][n][i]+fabs(wave->ct[m][n][i]));
	wave->tm[n] = 0.5*(wave->ct[m][n][i]-fabs(wave->ct[m][n][i]));	  
      }

      for(n=0; n < Nw; n++){
	if (n == 0){
	  a[n] = -fct*wave->tp[Nw-1];
	  b[n] = 1-(fct*wave->tm[Nw-1]-fct*wave->tp[n]);
	  d[n] -= a[n]*wave->N[m][Nw-1][i]; //Source due to periodicity
	}else{
	  a[n] = -fct*wave->tp[n-1];
	  b[n] = 1-(fct*wave->tm[n-1]-fct*wave->tp[n]);
	  d[n] -= c[n]*wave->N[m][1][i];   //Source due to periodicity
	}
	c[n] = fct*wave->tm[n];
	x[n] = 0;
	
      }
      TriSolve(a, b, c, d, (&x[0]), Nw);
      //Copy to the action density field
      for(n=0; n < Nw; n++){
        wave->N[m][n][i] = x[n];
      }
    }    
  }

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      ISendRecvCellData2D(wave->N[m][n], grid, myproc, comm);
    }
  }
  
  
}

void UpdateActionDensitySinkSource(gridT *grid, physT *phys, propT *prop,MPI_Comm comm, int myproc, int numprocs)
{
  int i, m, n, Nc = grid->Nc, Mw=wprop->Mw, Nw=wprop->Nw;
  REAL dt = prop->dt;

  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      //      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //	i = grid->cellp[iptr];
      for (i= 0; i <  Nc; i++){
	wave->Ntmp[m][n][i] = wave->N[m][n][i] + wave->ssrc[m][n][i]*dt*wprop->wnstep;
	wave->Ntmp[m][n][i] = wave->Ntmp[m][n][i]*exp(wave->src[m][n][i]*dt*wprop->wnstep);
      }
      ISendRecvCellData2D(wave->Ntmp[m][n], grid, myproc, comm);      
    }
  }
}


void ExplicitUpdateGeographic(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs)
{
  int i, iptr, m, nf, ne, n, nc1, nc2, Nc = grid->Nc, Mw=wprop->Mw, Nw=wprop->Nw, normal, flag = 0;
  REAL dt = prop->dt, source, Ac, df, dg;
  
  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //i = grid->cellp[iptr];
      for(i = 0; i < Nc; i++){
	Ac = grid->Ac[i];

	source = 0;
	flag = 0;
	for(nf=0;nf<grid->nfaces[i];nf++) {
	  ne = grid->face[i*grid->maxfaces+nf];
	  normal = grid->normal[i*grid->maxfaces+nf];
	  df = grid->df[ne];
	  dg = grid->dg[ne];
	  nc1 = grid->grad[2*ne];
	  nc2 = grid->grad[2*ne+1];
	  if(nc1==-1){
	    source -= 0.5*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]+fabs(wave->cg[m][n][ne]))*
	    				   wave->Ntmp[m][n][nc2]);
	    flag = 1;
	  }else if(nc2==-1){
	    source -= 0.5*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]-fabs(wave->cg[m][n][ne]))*
	    				   wave->Ntmp[m][n][nc1]);
	    flag = 1;
	  }else{
	    source -= 0.5*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]+fabs(wave->cg[m][n][ne]))*
					   wave->Ntmp[m][n][nc2]+
					   (wave->cg[m][n][ne]-fabs(wave->cg[m][n][ne]))*
					   wave->Ntmp[m][n][nc1]);
	  }  
	}
	
	wave->N[m][n][i] = wave->Ntmp[m][n][i] + source; 
       
	
      }
      ISendRecvCellData2D(wave->N[m][n], grid, myproc, comm);      
    }
  }
}

void ImplicitUpdateGeographic(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs)
{  
  int i, iptr, m, nf, ne, n, nc1, nc2, Nc = grid->Nc, Mw=wprop->Mw, Nw=wprop->Nw, normal, flag = 0;
  REAL dt = prop->dt, source, Ac, df, dg;
  
  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //i = grid->cellp[iptr];
      for(i = 0; i < Nc; i++){
	Ac = grid->Ac[i];

	source = 0;
	flag = 0;
	for(nf=0;nf<grid->nfaces[i];nf++) {
	  ne = grid->face[i*grid->maxfaces+nf];
	  normal = grid->normal[i*grid->maxfaces+nf];
	  df = grid->df[ne];
	  dg = grid->dg[ne];
	  nc1 = grid->grad[2*ne];
	  nc2 = grid->grad[2*ne+1];
          //if(i==559)
          //  if(wave->Ntmp[m][n][i]!=wave->Ntmp[m][n][i])
          ////    printf("%d %d %d %d %f %d %f \n",ne, i,m,n,wave->Ntmp[m][n][i],normal,df);
	  if(nc1==-1){
	    source -= 0.25*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]+fabs(wave->cg[m][n][ne]))*
	    				   wave->Ntmp[m][n][nc2]);
	    flag = 1;
	  }else if(nc2==-1){
	    source -= 0.25*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]-fabs(wave->cg[m][n][ne]))*
	    				   wave->Ntmp[m][n][nc1]);
	    flag = 1;
	  }else{
	    source -= 0.25*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]+fabs(wave->cg[m][n][ne]))*
					   wave->Ntmp[m][n][nc2]+
					   (wave->cg[m][n][ne]-fabs(wave->cg[m][n][ne]))*
					   wave->Ntmp[m][n][nc1]);
	  }  
	}
	wave->N[m][n][i] = wave->Ntmp[m][n][i] + source; 
	wave->Nold[m][n][i] = wave->N[m][n][i];
      }
      ISendRecvCellData2D(wave->N[m][n], grid, myproc, comm);      
    }
  }
  
  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){	
      BiCGSolveN(m, n, grid, prop, myproc, numprocs, comm);
    }
  }
}

// This function reconstructs the wind field using the Kriging spatial interpolation from wind stations
// and the linear interpolation in time. It also calculates the shear stress at the cell center for wind
// forcing input as well as at the edge if the surface shear induce by winds is considered. Since the
// measured wind data (usually, NOAA) have already been read and stored, this function is not responsible
// for any data reading. It is called every wave time step before solving the wind wave field. If forcing
// by a constant wind speed and direction is desired, one can just enforce constant values here. -YJ  
static void WindField(propT *prop, gridT *grid, physT *phys, int myproc, int numprocs)
{
  int Nc=grid->Nc, Ns = wprop->nstation, Ne = grid->Ne;
  int i, j;
  int N;
  REAL r1, r2, sx1, sx2,sy1,sy2, dg1, dg2, rtime_rel, Cd, dg1xy, dg2xy;

  rtime_rel = prop->rtime-(double)prop->nstart*prop->dt;
  N =  floor(rtime_rel/wprop->wind_dt);
  r1 = (double)((int)rtime_rel % (int)wprop->wind_dt);
  r1/=wprop->wind_dt;
  r2 = 1-r1;
  
  if (N+1 < wprop->Nwind){
    for (i=0; i< Nc; i++){
      wave->wind_spfx[i] = 0;
      wave->wind_spfy[i] = 0;
      wave->wind_dgf[i] = 0;
      sx1 = 0;
      sx2 = 0;
      sy1 = 0;
      sy2 = 0;
      dg1 = 0;
      dg2 = 0;
      dg1xy = 0;
      dg2xy = 0;
      for (j=0; j< Ns; j++){
	//Transform the wind coordinate to geographic (Cartesian) coordinate
        dg1xy = 0.5*PI - wave->wind_dg[j][N]+PI;
	dg2xy = 0.5*PI - wave->wind_dg[j][N+1]+PI;
  
	sx1 += wave->klambda[i][j]*wave->wind_sp[j][N]*cos(dg1xy);
	sx2 += wave->klambda[i][j]*wave->wind_sp[j][N+1]*cos(dg2xy);
	sy1 += wave->klambda[i][j]*wave->wind_sp[j][N]*sin(dg1xy);
	sy2 += wave->klambda[i][j]*wave->wind_sp[j][N+1]*sin(dg2xy);
	dg1 += wave->klambda[i][j]*dg1xy;
	dg2 += wave->klambda[i][j]*dg2xy;
      }
      wave->wind_spfx[i] = sx1*r2 + sx2*r1;
      wave->wind_spfy[i] = sy1*r2 + sy2*r1;
      wave->wind_spf[i] = sqrt(pow(wave->wind_spfx[i], 2) + pow(wave->wind_spfy[i], 2));

      wave->wind_dgf[i] = asin(wave->wind_spfy[i]/wave->wind_spf[i]);
      
      if(wave->wind_spfx[i] < 0)
	wave->wind_dgf[i] = PI - wave->wind_dgf[i];

      if(wave->wind_spf[i] < 7.5)
       	Cd = 1.2875*0.001;
      else
       	Cd = (0.8 + 0.065*wave->wind_spf[i])*0.001;	

      wave->uscx[i] = sqrt(Cd)*wave->wind_spf[i]*cos(wave->wind_dgf[i]);
      wave->uscy[i] = sqrt(Cd)*wave->wind_spf[i]*sin(wave->wind_dgf[i]);

    }
  }else if(N+1 == wprop->Nwind) {
      for (i=0; i< Nc; i++){
	wave->wind_spfx[i] = 0;
	wave->wind_spfy[i] = 0;
	wave->wind_dgf[i] = 0;
	for (j=0; j< Ns; j++){
	  dg1xy = 0.5*PI - wave->wind_dg[j][N]+PI;
	  wave->wind_spfx[i] += wave->klambda[i][j]*wave->wind_sp[j][N]*cos(dg1xy);
	  wave->wind_spfy[i] += wave->klambda[i][j]*wave->wind_sp[j][N]*sin(dg1xy);
	  wave->wind_dgf[i] += wave->klambda[i][j]*dg1xy;
	}
	wave->wind_spf[i] = sqrt(pow(wave->wind_spfx[i], 2) + pow(wave->wind_spfy[i], 2));
    
	wave->wind_dgf[i] = asin(wave->wind_spfy[i]/wave->wind_spf[i]);

	if(wave->wind_spfx[i] < 0)
	  wave->wind_dgf[i] = PI - wave->wind_dgf[i];

	if(wave->wind_spf[i] < 7.5)
	  Cd = 1.2875*0.001;
	else
	  Cd = (0.8 + 0.065*wave->wind_spf[i])*0.001;
	wave->uscx[i] = sqrt(Cd)*wave->wind_spf[i]*cos(wave->wind_dgf[i]);
	wave->uscy[i] = sqrt(Cd)*wave->wind_spf[i]*sin(wave->wind_dgf[i]);

      }
    }
  if (wprop->wind_shear == 1){
    for (j=0; j< Ne; j++)            
      phys->tau_T[j] = InterpWindShearToFace(j, phys, grid);
  }
}

void WaveTurbMixing(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int Nc = grid->Nc;
  int i,k;
  REAL lz, lz0, z, depth, fzprime;

  for(i = 0; i < Nc; i++){
    wave->fphi[i] = 1.22+0.22*cos(2*wave->Cr[i]);
    if (wave->ab[i] > 0.02){
      depth = grid->dv[i]+phys->h[i];
      lz0 = log10(sediments->kb[i]/30*wave->sgmean[i]/wave->ub[i]);
      z = 0;
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	z -= 0.5*grid->dzz[i][k];
	if (z < 0.0){
	  lz = log(fabs(z)*wave->sgmean[i]/fabs(wave->ub[i]));
	  wave->fz[i][k] = -0.0488 + 0.02917*lz + 0.01703*pow(lz, 2)
	    + (1.125*(lz0+5)+0.125*pow((lz0+5), 4))*(-0.0102-0.00253*lz+0.00273*pow(lz, 2));
	  fzprime = 0.02917 + 2.0*0.01703*lz
	    + (1.125*(lz0+5)+0.125*pow((lz0+5), 4))*(-0.00253+2.0*0.00273*lz);	  
	  if(wave->fz[i][k] < 0) wave->fz[i][k] = 0; 
	  if(fzprime <= 0) wave->fz[i][k] = 0;
	}else{
	  wave->fz[i][k] = 0;
	}
	z -= 0.5*grid->dzz[i][k];
      }
    }else{
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	wave->fz[i][k] = 0; 
      }
    }
      
  }

   
}


void ObtainKrigingCoef(gridT *grid, int myproc, int numprocs)
{
  int Nc = grid->Nc, Ns=wprop->nstation;
  int i, j, k, jj;

  REAL covMax = 0.78, Dmax = 160000, D, r;
  REAL **A = (REAL **)SunMalloc(Ns*sizeof(REAL *), "ObtainKrigingCoef");
  REAL **H = (REAL **)SunMalloc(Ns*sizeof(REAL *), "ObtainKrigingCoef");
  REAL *b = (REAL *)SunMalloc(Ns*sizeof(REAL), "ObtainKrigingCoef");

  for (k = 0; k< Ns; k++){
    A[k] = (REAL *)SunMalloc(Ns*sizeof(REAL), "ObtainKrigingCoef");
    H[k] = (REAL *)SunMalloc(Ns*sizeof(REAL), "ObtainKrigingCoef");
  }

  for(i = 0; i< Ns; i++){
    A[i][i] = covMax;
    H[i][i] = 1;
    for(j = i+1; j < Ns; j++){
      D = sqrt(pow(wprop->xw[i]-wprop->xw[j], 2)+
	       pow(wprop->yw[i]-wprop->yw[j], 2));
      A[i][j] = semivariogram(covMax, Dmax, D);
      A[j][i] = A[i][j];
      H[i][j] = 0;
      H[j][i] = 0;
    }
  }

  for(k = 0; k < Ns-1; k++){
    for(i = k+1; i < Ns; i++){
      r = A[i][k]/A[k][k];
      A[i][k] = 0;
      for (j = k+1; j < Ns; j++){
	A[i][j] = A[i][j]-A[k][j]*r;
      }

      for (jj = 0; jj < k+1; jj++){
	H[i][jj] = H[i][jj]-r*H[k][jj];
      }
    }
  }

  for(k = Ns-1; k >0; k--){
    for(i = k-1; i >= 0; i--){
      r = A[i][k]/A[k][k];
      A[i][k] = 0;
      for(j=0; j< Ns; j++){
	H[i][j]-=r*H[k][j];
      }
    }
  }

  
  for (i = 0; i < Ns; i++ ){
    for (j = 0; j < Ns; j++ ){
      H[i][j] /= A[i][i];   
    }
  }

  for (i = 0; i < Nc; i++){
    for (j = 0; j < Ns; j++ ){
      D = sqrt(pow(grid->xv[i]-wprop->xw[j], 2)+
	       pow(grid->yv[i]-wprop->yw[j], 2));
      b[j] = semivariogram(covMax, Dmax, D);   
      //if(b[j]==0)
        //printf("%d %f %f %f \n",i, covMax,Dmax,D);   
    }
    for (j = 0; j < Ns; j++ ){
      wave->klambda[i][j] = 0;
      for (k = 0; k < Ns; k++ ){
	wave->klambda[i][j]+=H[j][k]*b[k];
      }
    }
  }


  for (k = 0; k < Ns; k++){
    SunFree(A[k], Ns, "ObtainKrigingCoef");
    SunFree(H[k], Ns, "ObtainKrigingCoef");
  }
  SunFree(A, Ns, "ObtainKrigingCoef");
  SunFree(H, Ns, "ObtainKrigingCoef");
  SunFree(b, Ns, "ObtainKrigingCoef");
}

static REAL semivariogram(REAL Cov0, REAL Dmax, REAL D)
{
  REAL tmp;

  if (D > Dmax)
    return 0.0;
  else{
    tmp = D/Dmax;
    return (Cov0*(1-1.5*tmp+0.5*pow(tmp, 3)));
  }
} 

static void BiCGSolveN(int m0, int n0, gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, n, niters;
  REAL *x, *r, *p, *z, *r0, *v, *t, *s, eps, eps0, epsW = 10E-10, sum;
  REAL rho0 = 1, rho, alpha = 1, omg = 1, beta, tmp;

  z = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  v = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  s = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  t = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  r0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  x = wave->N[m0][n0];
  r = wave->Ntmp[m0][n0];
  p = wave->Nold[m0][n0];
  



  //  niters = prop->maxiters;
  niters = 1000;

  for(i = 0; i< grid->Nc; i++){
    z[i] = 0;
    v[i] = 0;
    s[i] = 0;
    t[i] = 0;
    r0[i] = 0;

  }

  ISendRecvCellData2D(x, grid, myproc, comm);
  OperatorN(m0, n0, x, z, grid, prop);

  for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
    i = grid->cellp[iptr];
    p[i] = p[i]-z[i]; //Now initial P stores the initial residual.
    r[i] = p[i];  //Initial residual
    r0[i] = r[i]; //Another initial residual for BiCG
  }

  for(iptr=grid->celldist[1]; iptr<grid->celldist[2]; iptr++){
    i = grid->cellp[iptr];
    p[i] = 0;
  }
  eps0=eps=InnerProduct2(r,r,grid,myproc,numprocs,comm);
  
  if (eps > epsW){
    if(!prop->resnorm) eps0 = 1;

    for(n = 0; n < niters && eps!=0; n++){ 
      rho = InnerProduct2(r0,r,grid,myproc,numprocs,comm);
      beta = rho/rho0*alpha/omg; //When n = 0, rho0, alpha, omg = 1.
 
     
      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
	i = grid->cellp[iptr];
	p[i] = r[i]+beta*(p[i]-omg*v[i]); 
      }


      ISendRecvCellData2D(p, grid, myproc, comm);
      OperatorN(m0, n0, p, v, grid, prop);
      tmp = InnerProduct2(r0,v,grid,myproc,numprocs,comm);
      tmp = 1/tmp;   
      alpha = rho*tmp;

    
     
      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
      	i = grid->cellp[iptr];
      	s[i] = r[i]-alpha*v[i];
      }
      
      eps = InnerProduct2(s, s, grid, myproc, numprocs, comm);
      
      if(VERBOSE>3) printf("BiCGSolve Iteration (1): %d, resid=%e, proc=%d\n",n,eps,myproc);
      if(eps<SMALL){
	x[i] += alpha*p[i];
	break;
      }
    
      ISendRecvCellData2D(s, grid, myproc, comm);
      OperatorN(m0, n0, s, t, grid, prop);
      
      tmp = InnerProduct2(t, t, grid, myproc, numprocs, comm);
      tmp = 1/tmp;
      omg = InnerProduct2(t, s, grid, myproc, numprocs, comm)*tmp;
    
      
      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
      	i = grid->cellp[iptr];
      	x[i] += alpha*p[i] + omg*s[i];
	r[i] = s[i]-omg*t[i];
      }
      rho0 = rho;
      
      eps = InnerProduct2(r, r, grid, myproc, numprocs, comm);
      
      if(VERBOSE>3) printf("BiCGSolve Iteration (2): %d, resid=%e, proc=%d\n",n,eps,myproc);
      if(eps < epsW)
	break;
    }
  }

  if(myproc==0 && VERBOSE>3)
  {
    if(eps0 < epsW)
      printf("Step %d, norm of action density (%d, %d) source at = %e is already small\n", prop->n, m0, n0, eps0);
    else
      if(n==niters) printf("Warning... Step %d, action density (%d, %d) iteration not converging after %d steps! RES=%e > %.2e\n",
			   prop->n, m0, n0, n, eps, SMALL);
      else printf("Step %d, BiCGSolve action density (%d, %d) converged after %d iterations, rsdl = %e < %e\n",
		  prop->n, m0, n0, n, eps, epsW);
  }  

  ISendRecvCellData2D(x, grid, myproc, comm);

  SunFree(z, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(v, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(s, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(t, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(r0, grid->Nc*sizeof(REAL), "BiCGSolveN");
}


static void OperatorN(int m, int n, REAL *x, REAL *y, gridT *grid, propT *prop){

  int i, nf, ne, nc1, nc2, iptr, Nc = grid->Nc, normal;
  REAL dt = prop->dt, source, Ac, df, dg;

  for (i = 0; i < Nc; i++){   
    Ac = grid->Ac[i];
    y[i] = x[i];
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1)
	y[i] += 0.25*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]+fabs(wave->cg[m][n][ne]))*x[nc2]);
      if(nc2==-1)
	y[i] += 0.25*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]-fabs(wave->cg[m][n][ne]))*x[nc1]);
      else
	y[i] += 0.25*dt*wprop->wnstep*df*normal/Ac*((wave->cg[m][n][ne]+fabs(wave->cg[m][n][ne]))*x[nc2]+
				     (wave->cg[m][n][ne]-fabs(wave->cg[m][n][ne]))*x[nc1]);         
	       
    }
  }


}

static REAL InnerProduct2(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm){
  int i, iptr;
  REAL sum, mysum = 0;

  for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
    i = grid->cellp[iptr];
      mysum+=x[i]*y[i];
  }

  MPI_Reduce(&mysum, &(sum), 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Bcast(&sum, 1, MPI_DOUBLE,0,comm);

  return sum;
}


static void CenterRadiationStress(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int i, k;
  REAL df, dg, Ac, dt=prop->dt, normal, volume, dz, sumUw;
  int Nc=grid->Nc, nf, ne, Nk_shallow, nc_deep;
  int nc1, nc2;
  int k_nc1, k_nc2, k_ne;

  for(i = 0; i < Nc; i++){
    Ac = grid->Ac[i];
    for(k = 0; k < grid->Nk[i]; k++){
      wave->divScx[i][k] = 0;
      wave->divScy[i][k] = 0;
      volume = 0;
      sumUw = 0;
      for(nf=0; nf<grid->nfaces[i]; nf++){
	ne = grid->face[i*grid->maxfaces+nf];
	normal = grid->normal[i*grid->maxfaces+nf];
	df = grid->df[ne];
	dg = grid->dg[ne];

	nc1 = grid->grad[2*ne];
	nc2 = grid->grad[2*ne+1];
	if (nc1 == -1){
	  //0.5 is due to phase average, i.e. phase average of cos(phi)^2 

	  wave->divScx[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	    *wave->ux[nc2][k]*grid->dzz[nc2][k];
	  wave->divScy[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	    *wave->uy[nc2][k]*grid->dzz[nc2][k];
	  volume += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))*grid->dzz[nc2][k];
	  sumUw += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]));

	}else if (nc2 == -1){

	  wave->divScx[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]-fabs(wave->Uw[ne][k]))
	    *wave->ux[nc1][k]*grid->dzz[nc1][k];
	  wave->divScy[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]-fabs(wave->Uw[ne][k]))
	    *wave->uy[nc1][k]*grid->dzz[nc1][k];
	  volume += (wave->Uw[ne][k]-fabs(wave->Uw[ne][k]))*grid->dzz[nc1][k];
	  sumUw += (wave->Uw[ne][k]-fabs(wave->Uw[ne][k]));

	}else{
	  Nk_shallow = grid->Nk[nc1];
	  nc_deep = nc2;
	  if (grid->Nk[nc2] < Nk_shallow){
	    Nk_shallow = grid->Nk[nc2];
	    nc_deep = nc1;
	  }
	  //if (k < Nk_shallow){

	  wave->divScx[i][k] -= (0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	      *wave->ux[nc2][k]*grid->dzz[nc2][k]
	      +0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]-fabs(wave->Uw[ne][k]))
				 *wave->ux[nc1][k]*grid->dzz[nc1][k]);

	  wave->divScy[i][k] -= (0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	      *wave->uy[nc2][k]*grid->dzz[nc2][k]
	      +0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]-fabs(wave->Uw[ne][k]))
				 *wave->uy[nc1][k]*grid->dzz[nc1][k]);

	    volume += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))*grid->dzz[nc2][k]
	      +(wave->Uw[ne][k]-fabs(wave->Uw[ne][k]))*grid->dzz[nc1][k];

	    sumUw += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))+(wave->Uw[ne][k]-fabs(wave->Uw[ne][k]));
	    /*}else if (nc_deep == nc2){
	    wave->divScx[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	      *wave->ux[nc2][k]*grid->dzz[nc2][k];
	    wave->divScy[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	      *wave->uy[nc2][k]*grid->dzz[nc2][k];
	    volume += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))*grid->dzz[nc2][k];
	    sumUw += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]));
	  }else{
	    wave->divScx[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	      *wave->ux[nc1][k]*grid->dzz[nc1][k];
	    wave->divScy[i][k] -= 0.5*0.5*dt*wprop->wnstep*df*normal/Ac*(wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))
	      *wave->uy[nc1][k]*grid->dzz[nc1][k];
	    volume += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]))*grid->dzz[nc2][k];
	    sumUw += (wave->Uw[ne][k]+fabs(wave->Uw[ne][k]));
	    }*/

	}

      }
     
      if (volume > 0.00001){
	wave->divScx[i][k] = wave->divScx[i][k]*sumUw/volume;
	wave->divScy[i][k] = wave->divScy[i][k]*sumUw/volume;
      }else{
	wave->divScx[i][k] = 0.0;
	wave->divScy[i][k] = 0.0;
      }
      if (wave->divScx[i][k] != wave->divScx[i][k]){
	printf("Bad Scx proc=%d Scx[i=%d][k=%d]=%f\n", myproc, i, k, wave->divScx[i][k]);
      }
      if (wave->divScy[i][k] != wave->divScy[i][k]){
	printf("Bad Scx proc=%d Scx[i=%d][k=%d]=%f\n", myproc, i, k, wave->divScy[i][k]);
      }
    }
  }
  // So far, radiation stress due to vertically oscillating motion is not included yet.
  // It will be included when interpolating to the cell edge.
  ISendRecvCellData3D(wave->divScx, grid, myproc, comm);    
  ISendRecvCellData3D(wave->divScy, grid, myproc, comm);        
}

static REAL InterpWindShearToFace(int j, physT *phys, gridT *grid){
  int nc1, nc2;
  REAL def1, def2, Dj, Senc1, Senc2, U, Ustar, usex, usey, cd;

  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  Dj = grid->dg[j];
  cd = sqrt(0.00128);
  
  if (nc1 == -1){
    U = sqrt(pow(wave->wind_spfx[nc2]-phys->uc[nc2][grid->ctop[nc2]], 2.0)
	    +pow(wave->wind_spfy[nc2]-phys->vc[nc2][grid->ctop[nc2]], 2.0));
    if (U <= 7.5)
      cd = 0.00128;
    else
      cd = (0.8+0.065*U)*0.001;

    return 0.00125*cd*U*((wave->wind_spfx[nc2]-phys->uc[nc2][grid->ctop[nc2]])*grid->n1[j] 
			 +(wave->wind_spfy[nc2]-phys->vc[nc2][grid->ctop[nc2]])*grid->n2[j]);
  }else if (nc2 == -1){
    U = sqrt(pow(wave->wind_spfx[nc1]-phys->uc[nc1][grid->ctop[nc1]], 2.0)
	    +pow(wave->wind_spfy[nc1]-phys->uc[nc1][grid->ctop[nc1]], 2.0));
    if (U <= 7.5)
      cd = 0.00128;
    else
      cd = (0.8+0.065*U)*0.001;
    return 0.00125*cd*U*((wave->wind_spfx[nc1]-phys->uc[nc1][grid->ctop[nc1]])*grid->n1[j]
			+(wave->wind_spfy[nc1]-phys->vc[nc1][grid->ctop[nc1]])*grid->n2[j]);

  }else{
    
    def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
    def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
  
    usex = ((wave->wind_spfx[nc2]-phys->uc[nc2][grid->ctop[nc2]])*def1+
	    (wave->wind_spfx[nc1]-phys->uc[nc1][grid->ctop[nc1]])*def2)/(def1+def2);
    usey = ((wave->wind_spfy[nc2]-phys->vc[nc2][grid->ctop[nc2]])*def1+
	    (wave->wind_spfy[nc1]-phys->vc[nc1][grid->ctop[nc1]])*def2)/(def1+def2);
    
    U = sqrt(pow(usex, 2.0) + pow(usey, 2.0));



    if (U <= 7.5)
      cd = 0.00128;
    else
      cd = (0.8+0.065*U)*0.001;

    return 0.00125*cd*U*(usex*grid->n1[j] + usey*grid->n2[j]);
  }
}

void EdgeRadiationStressToFlow(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int j, k, flag, jptr;
  int Ne=grid->Ne, nc1, nc2, nc;
  REAL weight1, dz, fab;  
  if(prop->n == 1+prop->nstart)
    fab = 1.0;
  else
    fab = 1.5;


  for(jptr = grid->edgedist[0]; jptr < grid->edgedist[1]; jptr++){
    j = grid->edgep[jptr];
    for(k=grid->etop[j];k<grid->Nkc[j];k++){
      phys->utmp[j][k]+=(1.0-fab)*wave->divSe[j][k];
      
    }
  }


  for(jptr = grid->edgedist[0]; jptr < grid->edgedist[1]; jptr++){
    j = grid->edgep[jptr];
    //  for(j = 0; j<Ne; j++){
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    dz = 0;

    for (k = 0; k < grid->etop[j]; k++) wave->divSe[j][k] = 0.0;
    flag = 0;
    
    for(k = grid->etop[j]; k < grid->Nke[j]; k++){
      if(nc1 == -1){
	wave->divSe[j][k] = wave->divScx[nc2][k]*grid->n1[j]+wave->divScy[nc2][k]*grid->n2[j];
	//Add radiation stress due to vertically oscillating motion
	//0.5 is due to phase average
	wave->divSe[j][k] += prop->dt * wprop->wnstep*0.5*(0 - pow(wave->uz[nc2][k], 2))/grid->dg[j];
	if (k == grid->etop[j]){
	  dz = grid->dzz[nc2][k];
	  if (dz > 0.01){
	    wave->divSe[j][k] -= prop->dt * wprop->wnstep * 1/16/dz*(0 - wave->Etot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}
	if(k == grid->etop[j] + 1){
	  if (flag == 0){
	    dz += grid->dzz[nc2][k];
	    wave->divSe[j][k] -= prop->dt * wprop->wnstep * 1/16/dz*(0 - wave->Etot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}
      }else if(nc2 == -1){
	wave->divSe[j][k] = wave->divScx[nc1][k]*grid->n1[j]+wave->divScy[nc1][k]*grid->n2[j];
	wave->divSe[j][k] += prop->dt * wprop->wnstep*0.5*(pow(wave->uz[nc1][k], 2) - 0)/grid->dg[j];
	if (k == grid->etop[j]){
	  dz = grid->dzz[nc1][k];
	  if (dz > 0.01){
	    wave->divSe[j][k] -= prop->dt * wprop->wnstep * 1/16/dz*(wave->Etot[nc1]*GRAV - 0)/grid->dg[j];
	    flag = 1;
	  }
	}
	if(k == grid->etop[j] + 1){
	  if (flag == 0){
	    dz += grid->dzz[nc1][k];
	    wave->divSe[j][k] -= prop->dt * wprop->wnstep * 1/16/dz*(wave->Etot[nc1]*GRAV - 0)/grid->dg[j];
	    flag = 1;
	  }
	}
      }else{
	weight1 = 1 - (grid->def[nc1*grid->maxfaces + grid->gradf[2*j]]/grid->dg[j]);

	//weight1 = 0.5;
	wave->divSe[j][k] = 
	  (1.0-weight1)*(wave->divScx[nc1][k] * grid->n1[j] + wave->divScy[nc1][k] * grid->n2[j]) + 
	  weight1*(wave->divScx[nc2][k] * grid->n1[j] + wave->divScy[nc2][k] * grid->n2[j]); 
	wave->divSe[j][k] += prop->dt * wprop->wnstep * 0.5 * (pow(wave->uz[nc1][k], 2) - pow(wave->uz[nc2][k], 2))/grid->dg[j];

	if (k == grid->etop[j]){
	  dz = 0.5*(grid->dzz[nc2][k] + grid->dzz[nc1][k]);
	  if (dz > 0.01){
	    wave->divSe[j][k] -= prop->dt * wprop->wnstep * 1/16/dz*(wave->Etot[nc1]*GRAV - wave->Etot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}  
       //if(wave->divSe[j][k] != wave->divSe[j][k] ){
	 // printf("dasdadasdas-proc=%d, divSe[j=%d][k=%d] = %f\n", myproc, j, k, wave->divSe[j][k]);
         // printf("dz%f nc1 %d %f nc2%d %f \n",dz,nc1, wave->Etot[nc1],nc2, wave->Etot[nc2]);
        //}
	if (k == grid->etop[j] + 1){
	  if (flag == 0){
	    dz = 0.5*(grid->dzz[nc2][k] + grid->dzz[nc1][k]);
	    wave->divSe[j][k] -= prop->dt * wprop->wnstep * 1/16/dz*(wave->Etot[nc1]*GRAV - wave->Etot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}      
      }
     // if(wave->divSe[j][k] != wave->divSe[j][k] )
	//  printf("----------------proc=%d, divSe[j=%d][k=%d] = %f\n", myproc, j, k, wave->divSe[j][k]);
    }
    for (k = grid->Nke[j]; k < grid->Nkc[j]; k++){
      wave->divSe[j][k] = wave->divSe[j][grid->Nke[j]-1];
    }

    j = grid->edgep[jptr];
    for(k=grid->etop[j];k<grid->Nkc[j];k++){
      phys->utmp[j][k]+=fab*wave->divSe[j][k];
      
    }

  }
  
  ISendRecvEdgeData3D(wave->divSe, grid, myproc, comm);        
}
static void WindSurfaceShear(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int jptr, j, nc1, nc2; 
  REAL h0;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    if(grid->etop[j]<grid->Nke[j]){
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
	
      if(nc1==-1)
	nc1=nc2;
      if(nc2==-1)
	nc2=nc1;
      h0 = 0.5*(grid->dzz[nc1][grid->etop[j]]+grid->dzz[nc2][grid->etop[j]]);
      // Add the shear stress to the top cell
      if (h0 < 0.01)
	h0 = 0.01;
      
      phys->Cn_U[j][grid->etop[j]] += phys->tau_T[j]/h0*prop->dt;
    }
  }
}

static void RadiationStressToFlow(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int jptr, j, k; 
  REAL fab;

  if (prop->n == 1)
    fab = 1;
  else
    fab = 1.5;
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    for(k=grid->etop[j];k<grid->Nkc[j];k++){
      phys->utmp[j][k]+=wave->divSe[j][k];
      
    }
  }      


}

void ObtainWaveVelocity(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int i, k;
  REAL z, depth, kD, kz, kw, kx, ky;
  int Nc=grid->Nc;
  REAL tmp, invsinh;

  for(i = 0; i < Nc; i++){
    z = 0;
    depth = grid->dv[i]+phys->h[i];
    kD = wave->kmean[i]*depth;
    kw  = wave->kmean[i];

    if (kD > 0.1){
      invsinh = 1/sinh(kD);
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	z -= 0.5*grid->dzz[i][k];
	kz = wave->kmean[i]*(z+depth);
	if (kD < 100)
	  tmp = wave->Hs[i]/8.0*wave->sgmean[i]*cosh(kz)*invsinh;
	else
	  tmp = wave->Hs[i]/8.0*wave->sgmean[i]*exp(kz-kD);
	
	wave->uz[i][k] = tmp*tanh(kz);
	wave->ux[i][k] = fabs(tmp*cos(wave->thtamean[i]));
	wave->uy[i][k] = fabs(tmp*sin(wave->thtamean[i]));
	z -= 0.5*grid->dzz[i][k];
	if (wave->ux[i][k] != wave->ux[i][k])
	  printf("Bad ux: p = %d i=%d k=%d ux=%f\n", myproc, i, k, wave->ux[i][k]);
	
	if (wave->uz[i][k] != wave->uz[i][k])
	  printf("Bad uy: p = %d i=%d k=%d uy=%f\n", myproc, i, k, wave->uy[i][k]);
	
	if (wave->uz[i][k] != wave->uz[i][k])
	  printf("Bad uz: p = %d i=%d k=%d uz=%f\n", myproc, i, k, wave->uz[i][k]);

      }
    }else{
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	wave->ux[i][k] = 0;
	wave->uy[i][k] = 0;
	wave->uz[i][k] = 0;
      }
    }
      
  }
  ISendRecvCellData3D(wave->ux, grid, myproc, comm);    
  ISendRecvCellData3D(wave->uy, grid, myproc, comm);    
  ISendRecvCellData3D(wave->uz, grid, myproc, comm);    
}

void ObtainEdgeUwField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int j, k, Ne=grid->Ne, nc1, nc2, jptr;
  REAL weight1;

  for(jptr = grid->edgedist[0]; jptr < grid->edgedist[1]; jptr++){
    j = grid->edgep[jptr];
    //  for(j = 0; j<Ne; j++){
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j + 1];
    
    weight1 = 0.5;
    for (k = 0; k < grid->etop[j]; k++) wave->Uw[j][k] = 0.0;
    for (k = grid->etop[j]; k < grid->Nke[j]; k++){
      if (nc1 == -1)
    	wave->Uw[j][k] = wave->ux[nc2][k]*grid->n1[j] + wave->uy[nc2][k]*grid->n2[j];
      else if (nc2 == -1)
    	wave->Uw[j][k] = wave->ux[nc1][k]*grid->n1[j] + wave->uy[nc1][k]*grid->n2[j];
      else
	//weight1 = 1 - (grid->def[nc1*grid->maxfaces + grid->gradf[2*j]]/grid->dg[j]);
    	wave->Uw[j][k] = weight1 * (wave->ux[nc1][k]*grid->n1[j] + wave->uy[nc1][k]*grid->n2[j])
    	  + (1 - weight1) * (wave->ux[nc2][k]*grid->n1[j] + wave->uy[nc2][k]*grid->n2[j]);
    }
    for (k = grid->Nke[j]; k < grid->Nkc[j]; k++)
      wave->Uw[j][k] = wave->Uw[j][grid->Nke[j] - 1];

    
  }
  ISendRecvEdgeData3D(wave->Uw, grid, myproc, comm);    
}


static void ReadWaveVariables(gridT *grid, propT *prop, int myproc, MPI_Comm comm) {

  int i, j;
  int m, n, s, l;

  //fread(&(prop->nstart), sizeof(int),1,wprop->StartWaveFID);

  for(m=0; m<wprop->Mw; m++)
    for(n=0; n<wprop->Nw; n++){
      fread(wave->N[m][n],sizeof(REAL),grid->Nc,wprop->StartWaveFID);
      //      ISendRecvCellData2D(wave->N[m][n], grid, myproc, comm);    
    }
  


  fclose(wprop->StartWaveFID);

}

static void OpenWaveFiles(int merge, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];


  MPI_GetFile(filename,DATAFILE,"WaveHeightFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->WaveHeightFID = MPI_FOpen(str,"w","OpenFile",myproc);

  MPI_GetFile(filename,DATAFILE,"WindSpeedFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->WindSpeedFID = MPI_FOpen(str,"w","OpenFile",myproc);

  MPI_GetFile(filename,DATAFILE,"WindDirectionFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->WindDirectionFID = MPI_FOpen(str,"w","OpenFile",myproc);

  MPI_GetFile(filename,DATAFILE,"WaveVelocityFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->WaveVelocityFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"StoreWaveFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->StoreWaveFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  if(RESTART) {
    MPI_GetFile(filename,DATAFILE,"StartWaveFile","OpenWaveFiles",myproc);
    if(merge)
      strcpy(str,filename);
    else
      sprintf(str,"%s.%d",filename,myproc);
    wprop->StartWaveFID = MPI_FOpen(str,"r","OpenWaveFiles",myproc);
  }

}


static void OutputWaveData(gridT *grid, propT *prop,
		int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int i, j, k, index, nwritten;
  int m, n, l, s; 
  REAL z;
  REAL *tmp = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"OutputWaveData");

  // change prop->n==prop->nstart+1 to prop->nstart to output initial condition 
  if(!(prop->n%prop->ntout) || prop->n==prop->nstart || blowup|| prop->n==prop->nsteps+prop->nstart) {

    Write2DData(wave->Hs,prop->mergeArrays,wprop->WaveHeightFID,"Error outputting wave height data!\n", grid,numprocs,myproc,comm);
  
    Write2DData(wave->wind_spf,prop->mergeArrays,wprop->WindSpeedFID,"Error outputting wind speed data!\n", grid,numprocs,myproc,comm);
   
    Write2DData(wave->wind_dgf,prop->mergeArrays,wprop->WindDirectionFID,"Error outputting wind direction data!\n", grid,numprocs,myproc,comm);

    Write2DData(wave->ub,prop->mergeArrays,wprop->WaveVelocityFID,"Error outputting wind velocity data!\n", grid,numprocs,myproc,comm);

  }

  //if(prop->n==prop->nsteps+prop->nstart || blowup) {
  if(!(prop->n%(prop->ntout)) || prop->n==prop->nsteps+prop->nstart || blowup){  
 
    if(VERBOSE>=1 && myproc==0) printf("Writing to wave rstore...\n");
    fseek( wprop->StoreWaveFID, 0, SEEK_SET );
    //nwritten=fwrite(&(prop->n),sizeof(int),1,wprop->StoreWaveFID);

    if (prop->wavemodel)
      for(m=0; m<wprop->Mw; m++)
	for(n=0; n<wprop->Nw; n++)
	  fwrite(wave->N[m][n],sizeof(REAL),grid->Nc,wprop->StoreWaveFID);
    fflush(wprop->StoreWaveFID);

  }
  
  if(prop->n==prop->nsteps+prop->nstart) {
    fclose(wprop->WaveHeightFID);
    fclose(wprop->WaveVelocityFID);
    fclose(wprop->WindSpeedFID);
    fclose(wprop->WindDirectionFID);
    fclose(wprop->StoreWaveFID);
  }
}

/*
 * Function: UpdateWaveCr
 * usage: calculate the angle between wave direction and 
 * bottom current direction
 * -----------------------------------------------------
 * include in ComputeSediments in sediments.c
 *
 */
void UpdateWaveCr(gridT *grid, physT *phys, int myproc)
{
   int i;
   REAL angle;
   for(i=0;i<grid->Nc;i++){
     if(phys->uc[i][grid->Nk[i]-1] != 0) 
       angle = atan(phys->vc[i][grid->Nk[i]-1]/phys->uc[i][grid->Nk[i]-1]);
     else
       if(phys->vc[i][grid->Nk[i]-1] > 0)
	 angle = PI/2;
       else
	 angle = 3*PI/2;
     wave->Cr[i] = angle-wave->thtamean[i];
   }
}

/*
 * Function: UpdateWaveFw
 * usage: calculate the bottom drag coefficient from wave
 * -----------------------------------------------------
 * include in ComputeSediments in sediments.c
 *
 */
void UpdateWaveFw(gridT *grid, int myproc)
{
  int i;
  for(i=0;i<grid->Nc;i++){
    if(wave->ab[i] <= 0.2*sediments->kb[i])
      wave->fw[i] = 0.3;
    else if(wave->ab[i] <= 100*sediments->kb[i])
      wave->fw[i] = exp(-8.82+7.02*pow((wave->ab[i]/sediments->kb[i]), -0.140));
    else
      wave->fw[i] = exp(-7.30+5.61*pow((wave->ab[i]/sediments->kb[i]), -0.209));
    if(wave->fw[i] > 0.3) wave->fw[i] = 0.3;
  }
}



/*----------------------------------------------------------------------------*/
//fetch model

/*
 * Function: AllocateWave
 * Usage: allocate space for wave variables
 * ----------------------------------------------------
 * Based on the value from ReadWaveProperties
 *
 */
void FetchAllocateWave(gridT *grid, int myproc){   
  char filename[BUFFERLENGTH];
  FILE *ifile;

  wprop->ConstantWind = MPI_GetValue(DATAFILE,"constantwind","ReadWaveProperties",myproc);
  wprop->FetchModel = MPI_GetValue(DATAFILE,"fetchmodel","ReadWaveProperties",myproc);
  wave->Hwsig = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  wave->Twsig = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  wave->Fw = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  MPI_GetFile(filename,DATAFILE,"InputCoarseDomain","AllocateWave",myproc); 
  wprop->Numdomain = MPI_GetSize(filename,"AllocateWave",myproc);
  wave->Waveexcur = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  wave->Fetch = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  wave->Coarsedomain = (REAL *)SunMalloc(wprop->Numdomain*2*sizeof(REAL), "AllocateWave");
 
  wave->Uwind = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  wave->Winddir = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");

}

/*
 * Function: InitializeWave
 * Usage: give initial value for Wave variables
 * ----------------------------------------------------
 * Based on the value from ReadWaveProperties
 */
void FetchInitializeWave(gridT *grid, physT *phys, int myproc)
{
  int i;
  REAL uwind, winddir;
  char str[BUFFERLENGTH],filename[BUFFERLENGTH];
  FILE *ifile;

  // give initial 0 for Hwsig Twsig Fw Fetch Uwind Winddir
  for(i=0;i<grid->Nc;i++){
    wave->Hwsig[i]=0;
    wave->Twsig[i]=0;
    wave->Fw[i]=0;
    wave->Fetch[i]=0;

   
    wave->Uwind[i]=0;
    wave->Winddir[i]=0;
    wave->Waveexcur[i]=0;
  } 
  for(i=0;i<wprop->Numdomain;i++)
  {
    wave->Coarsedomain[i*2]=0;
    wave->Coarsedomain[i*2+1]=0;
  }
  MPI_GetFile(filename,DATAFILE,"InputCoarseDomain","CalculateFetch",myproc);     
  ifile=MPI_FOpen(filename,"r","CalculateFetch",myproc);
  for(i=0;i<wprop->Numdomain;i++){
    wave->Coarsedomain[2*i]=getfield(ifile,filename);
    wave->Coarsedomain[2*i+1]=getfield(ifile,filename);
    // make sure coarsedomain include all the voro points!!!!
  }
  fclose(ifile);
 
  uwind= MPI_GetValue(DATAFILE,"uwind","InitializeWave",myproc);
  winddir= MPI_GetValue(DATAFILE,"winddir","InitializeWave",myproc);
  
  // give initial value for Uwind and winddir
  if(wprop->ConstantWind==0 || wprop->ConstantWind==2){  
    for(i=0;i<grid->Nc;i++){
      wave->Uwind[i]=ReturnWindSpeed(grid->xv[i],grid->yv[i]);
      wave->Winddir[i]=ReturnWindDirection(grid->xv[i],grid->yv[i]);
    }
  } else if(wprop->ConstantWind==1){
    for(i=0;i<grid->Nc;i++){
      wave->Uwind[i]=uwind;
      wave->Winddir[i]=winddir;
    }
  } else {
    MPI_GetFile(str,DATAFILE,"InitWindFile","InitializeWind",myproc);     
    ifile=MPI_FOpen(str,"r","InitializeWind",myproc);
    for(i=0;i<grid->Nc;i++){ 
      getfield(ifile,str);
      getfield(ifile,str);
      wave->Uwind[i]=(REAL)getfield(ifile,str);
      wave->Winddir[i]=(REAL)getfield(ifile,str);
    }
    fclose(ifile);
  }
  
  CalculateFetch(grid,myproc);

  FetchCalculateHwsigTwsig(grid,phys,myproc);

  FetchCalculateWaveexcur(grid,phys,myproc);

  FetchCalculateFw(grid,myproc);
}

/*
 * Function: FreeWave
 * Usage: free space for all the variables
 * ----------------------------------------------------
 * Basic sunfree function
 *
 */
void FetchFreeWave(gridT *grid, int myproc)
{
  free(wave->Hwsig);
  free(wave->Twsig);
  free(wave->Uwind);
  free(wave->Winddir);
  free(wave->Fw);
  free(wave->Fetch);
  free(wave->Coarsedomain);
  free(wave->Waveexcur);
}

/*
 * Function: CalculateFetch
 * usage: calculate fetch based on a coarse domain
 * --------------------------------------------------------
 * calculate intersection point and distance
 *
 */
void CalculateFetch(gridT *grid, int myproc)
{
  int i,j;
  REAL tmp,xinter,yinter,distmp,dismin;
  for(i=0;i<grid->Nc;i++){
    wave->Fetch[i]=10000000; // assign a big number
    dismin=10000000; // remember min distance
    for(j=0;j<wprop->Numdomain-1;j++){
      tmp=(wave->Coarsedomain[2*j]-wave->Coarsedomain[2*j+2])*sin(wave->Winddir[i])-(wave->Coarsedomain[2*j+1]-wave->Coarsedomain[2*j+3])*cos(wave->Winddir[i]);
      // if tmp!=0 two lines are not parallel and have an intersection point
      if(tmp!=0){
        xinter=wave->Coarsedomain[2*j]*wave->Coarsedomain[2*j+3]-wave->Coarsedomain[2*j+1]*wave->Coarsedomain[2*j+2];
        xinter=xinter*cos(wave->Winddir[i]);
        xinter=xinter+(-grid->yv[i]*cos(wave->Winddir[i])+sin(wave->Winddir[i])*grid->xv[i])*(wave->Coarsedomain[2*j]-wave->Coarsedomain[2*j+2]);
        xinter=xinter/tmp;
        // test whether is in the middle part
        if(xinter>=Min(wave->Coarsedomain[2*j],wave->Coarsedomain[2*j+2]) && xinter<=Max(wave->Coarsedomain[2*j],wave->Coarsedomain[2*j+2])){
          if((wave->Coarsedomain[2*j]-wave->Coarsedomain[2*j+2])!=0){
            yinter=((wave->Coarsedomain[2*j+1]-wave->Coarsedomain[2*j+3])*xinter+(wave->Coarsedomain[2*j]*wave->Coarsedomain[2*j+3]-wave->Coarsedomain[2*j+1]*wave->Coarsedomain[2*j+2]))/(wave->Coarsedomain[2*j]-wave->Coarsedomain[2*j+2]);
          } else {
            yinter=sin(wave->Winddir[i])/cos(wave->Winddir[i])*(xinter-grid->xv[i])+grid->yv[i];
          }
          distmp=sqrt(pow((grid->xv[i]-xinter),2)+pow((grid->yv[i]-yinter),2));
          if(distmp<dismin)
            dismin=distmp;
          if(distmp<wave->Fetch[i] && (xinter-grid->xv[i])*cos(wave->Winddir[i])>=0 && (yinter-grid->yv[i])*sin(wave->Winddir[i])>=0)
            wave->Fetch[i]=distmp;
        }
      }    
    }
    // output error
    if(wave->Fetch[i]==10000000){
      wave->Fetch[i]=dismin;
      printf("\n there is something wrong at cell No.%d when calculating Fetch, the reason may be the voro point is out of coarse domain\n",i);
    }
  }
}

 
/*
 * Function: CalculateHwsigTwsig
 * usage: Calculate significant wave height and wave period by Corps of Engineers
 * shore protection Manual Eqns. 3-39 3-40
 * --------------------------------------------------------------------------------
 * Using results of Fetch and Uwind
 *
 */
void FetchCalculateHwsigTwsig(gridT *grid, physT *phys, int myproc)
{
   int i;
   REAL tmp,tmp1,tmp2,grav,ua;
   grav=9.81; 
   for(i=0;i<grid->Nc;i++){
     ua=0.71*pow(wave->Uwind[i],1.23);
     tmp1=grav*wave->Fetch[i]/ua/ua;
     tmp2=grav*(phys->h[i]+grid->dv[i])/ua/ua;
     // for hwsig
     tmp=0;
     tmp=0.00565*pow(tmp1,0.5);
     tmp=tmp/tanh(0.530*pow(tmp2,0.75));
     tmp=tanh(tmp)*tanh(0.530*pow(tmp2,0.75))*0.283;
     wave->Hwsig[i]=tmp*ua*ua/grav;
     // for twsig
     tmp=0;
     tmp=0.0379*pow(tmp1,0.33333333333);
     tmp=tmp/tanh(0.833*pow(tmp2,0.375));
     tmp=tanh(tmp)*tanh(0.833*pow(tmp2,0.75))*7.54;
     wave->Twsig[i]=tmp*ua/grav; 
   }
}

/*
 * Function CalculateWaveexcur
 * usage: get wave semi-excursion (z=-D) for by solving dispersion relation
 * ------------------------------------------------------------
 * using Hwsig and Twsig
 *
 */
void FetchCalculateWaveexcur(gridT *grid,physT *phys, int myproc)
{  
  int i;
  REAL k1,k2,kguess,k,err,grav,omega,dk;
  grav=9.81;
  for(i=0;i<grid->Nc;i++){     
    omega=2*PI/wave->Twsig[i];
    k1=pow(omega,2)/grav;
    k2=pow((pow(omega,2)/grav/(phys->h[i]+grid->dv[i])),0.5);
    kguess=Min(k1,k2);  
    // initial guess
    err=fabs(pow(omega,2)/kguess-grav*tanh(kguess*(phys->h[i]+grid->dv[i])));      
    k=kguess;
    dk=fabs(k1-k2)/20;    
    // solve dispersion relation
    while(kguess<Max(k1,k2)){
      if(err>=fabs(pow(omega,2)/kguess-grav*tanh(kguess*(phys->h[i]+grid->dv[i])))){
        k=kguess;
        err=fabs(pow(omega,2)/kguess-grav*tanh(kguess*(phys->h[i]+grid->dv[i])));
      }
      kguess=kguess+dk; // you can set accuracy
    }
    wave->Waveexcur[i]=wave->Hwsig[i]/2/sinh(k*(phys->h[i]+grid->dv[i]));
  }
}

/*
 * Function CalculateFw
 * usage: calculate wave friction coefficient for bottom stress
 * ------------------------------------------------------------
 * based on the value of Wave renolds number
 * may change when using other case
 */
void FetchCalculateFw(gridT *grid,int myproc)
{  
  int i;
  REAL nu,rew;
  nu=1e-6;
  for(i=0;i<grid->Nc;i++){
    rew=wave->Waveexcur[i]*wave->Waveexcur[i]*2*PI/wave->Twsig[i]/nu;
    if(rew!=0){
      if(rew<3e5){
        wave->Fw[i]=2*pow(rew,-0.5);
      } else if(rew>=3e5 && rew<=6e5) {
        wave->Fw[i]=0.024*pow(rew,-0.123);
      } else {
        wave->Fw[i]=0.005;
      }
    } 
    if(rew==0) // you can set limit
      wave->Fw[i]=0;
  }
}

/*
 * Function: OutputWave
 * usage: output Hwsig, Twsig and Fetch for Sunplot
 * -----------------------------------------------
 * same as layerthick.dat
 * right now only for fetch model may change with different wave model
 */
void FetchOutputWave(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm)
{
  int i, j, jptr, k, nwritten,nosize;
  char str[BUFFERLENGTH];
  FILE *ofile;
  
  if(wprop->ConstantWind==0){
    if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) { 
       Write2DData(wave->Fetch,prop->mergeArrays,wprop->FetchFID,"Error outputting fetch data!\n", grid,numprocs,myproc,comm);
       Write2DData(wave->Hwsig,prop->mergeArrays,wprop->HwsigFID,"Error outputting significant wave height data!\n", grid,numprocs,myproc,comm);
       Write2DData(wave->Twsig,prop->mergeArrays,wprop->TwsigFID,"Error outputting significant wave period data!\n", grid,numprocs,myproc,comm);
    }
    if(prop->n==prop->nsteps+prop->nstart) {
      fclose(wprop->FetchFID);
      fclose(wprop->HwsigFID);
      fclose(wprop->TwsigFID);
    }
  } else {
    if(prop->n==1+prop->nstart){
      Write2DData(wave->Fetch,prop->mergeArrays,wprop->FetchFID,"Error outputting fetch data!\n", grid,numprocs,myproc,comm);
      Write2DData(wave->Hwsig,prop->mergeArrays,wprop->HwsigFID,"Error outputting significant wave height data!\n", grid,numprocs,myproc,comm);
      Write2DData(wave->Twsig,prop->mergeArrays,wprop->TwsigFID,"Error outputting significant wave period data!\n", grid,numprocs,myproc,comm);
      fclose(wprop->FetchFID);
      fclose(wprop->HwsigFID);
      fclose(wprop->TwsigFID);
    }
  }
}


/* 
 * Function: OpenWaveFile
 * usage: open wave output files for Hwsig, Twsig, and Fetch
 * ----------------------------------------------------------
 * same as OpenSediFile
 *
 */
void FetchOpenWaveFiles(int merge, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  int i;
  MPI_GetFile(filename,DATAFILE,"HwsigFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->HwsigFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"TwsigFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->TwsigFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);
  
  MPI_GetFile(filename,DATAFILE,"FetchFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  wprop->FetchFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);
}
