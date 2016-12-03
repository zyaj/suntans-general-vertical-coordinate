/*
 * File: wave.c
 * Author : Yun Zhang & Yi-Ju Chou
 * Institution: Stanford University
 * --------------------------------
 * main function: UpdateWave()
 * include Fetchmodel and Yi-Ju's model
 */

#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "initialization.h"
#include "boundaries.h"
#include "sediments.h"
#include "util.h"
#include "tvd.h"
#include "mympi.h"
#include "scalars.h"
#include "marsh.h"
#include "wave.h"
#include "physio.h"
#include "memory.h"

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
static void OpenWaveFiles(int myproc);
static void OutputWaveData(gridT *grid, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
static void ObtainKrigingCoef(gridT *grid, int myproc, int numprocs);
static void WindField(propT *prop, gridT *grid, physT *phys, int myproc, int numprocs);
static void WindSurfaceShear(gridT *grid, physT *phys, propT *prop,MPI_Comm comm, int myproc);
static void RadiationStressToFlow(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc);

// Fetch model
void FetchReadWaveProperties(int myproc);
void FetchAllocateWave(gridT *grid, int myproc);
void FetchInitializeWave(gridT *grid, physT *phys, int myproc);
void FetchFreeWave(gridT *grid, int myproc);
void FetchCalculateFetch(gridT *grid, int myproc);
void FetchCalculateHwsigTwsig(gridT *grid,physT *phys, int myproc);
void FetchCalculateWaveexcur(gridT *grid,physT *phys, int myproc);
void FetchCalculateFw(gridT *grid,int myproc);
void FetchOutputWave(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
void FetchOpenWaveFiles(int merge, int myproc);

void InitializeWaveProperties(propT *prop, int myproc)
{
  int i, Ns;
  char xstr[BUFFERLENGTH], ystr[BUFFERLENGTH];
  char xw[5], yw[5];
  char fd[BUFFERLENGTH], ifile[BUFFERLENGTH];
  REAL m;
  
  WaveMw = MPI_GetValue(DATAFILE,"Mw","InitializeWaveProperties",myproc); 
  WaveNw = MPI_GetValue(DATAFILE,"Nw","InitializeWaveProperties",myproc);
  Wavewnstep = MPI_GetValue(DATAFILE,"wnstep","InitializeWaveProperties",myproc);
  Wavesgmin = MPI_GetValue(DATAFILE,"sgmin","InitializeWaveProperties",myproc);
  Wavesgmax = MPI_GetValue(DATAFILE,"sgmax","InitializeWaveProperties",myproc);
  Wavewind_dt = MPI_GetValue(DATAFILE,"wind_dt","InitializeWaveProperties",myproc);
  Waveimplicit_whitecap = MPI_GetValue(DATAFILE,"implicit_whitecap","InitializeWaveProperties",myproc);
  Waveimplicit_advection = MPI_GetValue(DATAFILE,"implicit_advection","InitializeWaveProperties",myproc);
  Wavewind_forcing = MPI_GetValue(DATAFILE,"wind_forcing","InitializeWaveProperties",myproc);
  Wavenstation = MPI_GetValue(DATAFILE,"nstation","InitializeWaveProperties",myproc);
  Wavetail_opt = MPI_GetValue(DATAFILE,"tail_opt","InitializeWaveProperties",myproc);
  Wavetail_pow = MPI_GetValue(DATAFILE,"tail_pow","InitializeWaveProperties",myproc);
  Wavewind_shear = MPI_GetValue(DATAFILE,"wind_shear","InitializeWaveProperties",myproc);
  Waverad_stress = MPI_GetValue(DATAFILE,"rad_stress","InitializeWaveProperties",myproc);
  Waveform_drag = MPI_GetValue(DATAFILE,"form_drag","InitializeWaveProperties",myproc);
  Wavebtm_mud = MPI_GetValue(DATAFILE,"btm_mud","InitializeWaveProperties",myproc);
  Wavebtm_sedi_erosion = MPI_GetValue(DATAFILE,"btm_sedi_erosion","InitializeWaveProperties",myproc);
  Wavedepth_fw_cutoff = MPI_GetValue(DATAFILE,"depth_fw_cutoff","InitializeWaveProperties",myproc);
  Wavefw_drag = MPI_GetValue(DATAFILE,"fw_drag","InitializeWaveProperties",myproc);
  Wavebtm_conc = MPI_GetValue(DATAFILE,"btm_conc","InitializeWaveProperties",myproc);
  Wavebtm_vis = MPI_GetValue(DATAFILE,"btm_vis","InitializeWaveProperties",myproc);
  Wavebtm_mud_thickness = MPI_GetValue(DATAFILE,"btm_mud_thickness","InitializeWaveProperties",myproc);
  Wavedepth_brk_cutoff = MPI_GetValue(DATAFILE,"depth_brk_cutoff","InitializeWaveProperties",myproc);
  Wavedepth_brk_indx = MPI_GetValue(DATAFILE,"depth_brk_indx","InitializeWaveProperties",myproc);
  WaveNLtriad = MPI_GetValue(DATAFILE,"NLtriad","InitializeWaveProperties",myproc);
  WaveNLquad = MPI_GetValue(DATAFILE,"NLquad","InitializeWaveProperties",myproc);
  WaveBRKdepth = MPI_GetValue(DATAFILE,"BRKdepth","InitializeWaveProperties",myproc);


//The tail after sg_max; (based on SWAN)
  //N_0.01 = N_end (sg99/sgmax)^(-m)
  //sg99=sgmax*(0.01)^(-1/m)
  Wavesg99 = Wavesgmax*pow(0.01, -1/Wavetail_pow);
  //sg_tail is obtained by dividing the integral of N*sg dsg by integral of N dsg. The integral is taken
  //from sg_max to sg99.
  m = Wavetail_pow;
  Wavesgtail = (1.0-m)/(2.0-m)*(pow(Wavesg99, 2.0-m)-pow(Wavesgmax, 2.0-m))
    /(pow(Wavesg99, 1.0-m)-pow(Wavesgmax, 1.0-m));

  Ns = Wavenstation;

  WaveNwind = floor((prop->nsteps*prop->dt)/Wavewind_dt)+1;
  Wavenwind = floor(Wavewind_dt/prop->dt);
  
  Wavexw = (REAL *)SunMalloc(Ns*sizeof(REAL), "InitializeWaveProperties");
  Waveyw = (REAL *)SunMalloc(Ns*sizeof(REAL), "InitializeWaveProperties");

  for (i=1; i<= Ns; i++){
    sprintf(xstr, "xw.%d", i);
    Wavexw[i-1] = MPI_GetValue(DATAFILE,xstr,"InitializeWaveProperties",myproc);
    sprintf(ystr, "yw.%d", i);
    Waveyw[i-1] = MPI_GetValue(DATAFILE,ystr,"InitializeWaveProperties",myproc);
    
  }
  
     
}


void AllocateWaveVariables(gridT *grid, propT *prop)
{
  int flag=0, i, j, m, n; 
  int Nc=grid->Nc, Ne=grid->Ne, Mw=WaveMw, Nw=WaveNw,
      MN = Max(Mw, Nw), nstation = Wavenstation, Nwind=WaveNwind;
  
  //intrinsic wave parameters
  
  Wavesg = (REAL *)SunMalloc(Mw*sizeof(REAL), "AllocateWaveVariables");
  Wavedsg = (REAL *)SunMalloc(Mw*sizeof(REAL), "AllocateWaveVariables");
  Wavethtaw = (REAL *)SunMalloc(Nw*sizeof(REAL), "AllocateWaveVariables");
  //significant wave height and orbital velocity
  WaveHw =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Waveub =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavefz =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  Wavefphi =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Waveuscx =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Waveuscy =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Waveab =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavewind_spfx =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavewind_spfy =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavewind_dgf =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavewind_spf =  (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveEtot = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveEtmp = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavekmean = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavesgmean = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveT0 = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables"); 
  WaveHs = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveT01 = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavektail = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavesg_PM = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveCr = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavethtamean = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavetmp = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Wavefw = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveEtail = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  WaveNtail = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
  Waveklambda = (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  //wave variables in cell centers
  Wavekw = (REAL **)SunMalloc(Mw*sizeof(REAL *), "AllocateWaveVariables");
  Wavecgx = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  Wavecgy = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  Wavecph = (REAL **)SunMalloc(Mw*sizeof(REAL *), "AllocateWaveVariables");
  Wavecs = (REAL ***)SunMalloc((Mw+1)*sizeof(REAL **), "AllocateWaveVariables");
  Wavect = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  WaveN = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  WaveNtmp = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  WaveNold = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  Wavessrc = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  Wavetsrc = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  Wavesrc = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  
  //wave velocity stored at cell centers
  Waveux =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  Waveuy =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  Waveuz =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  WavedivScx =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");
  WavedivScy =  (REAL **)SunMalloc(Nc*sizeof(REAL *), "AllocateWaveVariables");

  //wave properties stored at cell centers
  WaveUw =  (REAL **)SunMalloc(Ne*sizeof(REAL *), "AllocateWaveVariables");
  WavedivSe =  (REAL **)SunMalloc(Ne*sizeof(REAL *), "AllocateWaveVariables");
  Waveab_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  Wavesgmean_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  Wavethtamean_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  Wavekmean_edge =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  Wavekw_edge =  (REAL **)SunMalloc(Ne*sizeof(REAL *), "AllocateWaveVariables");
  Waveuse =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  WaveUwind =  (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
  //wind variables
  Wavewind_sp = (REAL **)SunMalloc(nstation*sizeof(REAL *), "AllocateWaveVariables");
  Wavewind_dg = (REAL **)SunMalloc(nstation*sizeof(REAL *), "AllocateWaveVariables");

  for(i=0; i< Nc; i++){
    Waveklambda[i] = (REAL *)SunMalloc(nstation*sizeof(REAL), "AllocateWaveVariables");
    Wavefz[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    Waveuz[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    Waveux[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    Waveuy[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    WavedivScx[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
    WavedivScy[i] =  (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL), "AllocateWaveVariables");
  }

  for(j=0; j< Ne; j++){    
    WaveUw[j] =  (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL), "AllocateWaveVariables");
    WavedivSe[j] =  (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL), "AllocateWaveVariables");
    Wavekw_edge[j] =  (REAL *)SunMalloc(Mw*sizeof(REAL), "AllocateWaveVariables");
  }

  for(i=0; i< nstation; i++){
    Wavewind_sp[i] = (REAL *)SunMalloc(Nwind*sizeof(REAL), "AllocateWaveVariables"); 
    Wavewind_dg[i] = (REAL *)SunMalloc(Nwind*sizeof(REAL), "AllocateWaveVariables"); 
  }

  //coordinates of wind stations


  //matrix operators
  Wavea = (REAL *)SunMalloc(MN*sizeof(REAL), "AllocateWaveVariables");
  Waveb = (REAL *)SunMalloc(MN*sizeof(REAL), "AllocateWaveVariables");
  Wavec = (REAL *)SunMalloc(MN*sizeof(REAL), "AllocateWaveVariables");
  Wavesp = (REAL *)SunMalloc((Mw+1)*sizeof(REAL), "AllocateWaveVariables");
  Wavesm = (REAL *)SunMalloc((Mw+1)*sizeof(REAL), "AllocateWaveVariables");

  Wavetp = (REAL *)SunMalloc(Nw*sizeof(REAL), "AllocateWaveVariables");
  Wavetm = (REAL *)SunMalloc(Nw*sizeof(REAL), "AllocateWaveVariables");

  for (m=0; m<Mw; m++){
    Wavekw[m] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
    Wavecgx[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    Wavecgy[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");    
    Wavecs[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    Wavect[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    WaveN[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    WaveNtmp[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    WaveNold[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    Wavessrc[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    Wavetsrc[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    Wavesrc[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
    Wavecph[m] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");

    for (n=0; n<Nw; n++){
      Wavecgx[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      Wavecgy[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      Wavecs[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      Wavect[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      WaveN[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      WaveNtmp[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      WaveNold[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      Wavessrc[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      Wavetsrc[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
      Wavesrc[m][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");
    }
  }
  Wavecs[Mw] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");
  for (n=0; n<Nw; n++){
    Wavecs[Mw][n] = (REAL *)SunMalloc(Nc*sizeof(REAL), "AllocateWaveVariables");

  }

  //wave variables in cell edges
  Wavecg = (REAL ***)SunMalloc(Mw*sizeof(REAL **), "AllocateWaveVariables");
  for (m=0; m<Mw; m++){
    Wavecg[m] = (REAL **)SunMalloc(Nw*sizeof(REAL *), "AllocateWaveVariables");   
    for (n=0; n<Nw; n++){
      Wavecg[m][n] = (REAL *)SunMalloc(Ne*sizeof(REAL), "AllocateWaveVariables");
    }
  }
}


void FreeWaveVariables(gridT *grid, propT *prop)
{
  int i, j, m, n, Nc=grid->Nc, Ne=grid->Ne, Mw=WaveMw, Nw=WaveNw,
    nstation=Wavenstation, Nwind=WaveNwind;
    
  
  for(i=0; i< nstation; i++){
    free(Wavewind_sp[i]);
    free(Wavewind_dg[i]);
  }
  for(i = 0; i< Nc; i++){
    free(Waveklambda[i]);
    free(Wavefz[i]);
    free(Waveuz[i]);
    free(Waveux[i]);
    free(Waveuy[i]);
    free(WavedivScx[i]);
    free(WavedivScy[i]);
  }
  for(j = 0; j< Ne; j++){
    free(WaveUw[j]);
    free(WavedivSe[j]);
    free(Wavekw_edge[j]);
  }
  free(WaveUw);
  free(WavedivSe);
  free(Waveklambda);
  free(Wavewind_sp);
  free(Wavewind_dg);
  free(WaveHw);
  free(Wavewind_spfx);
  free(Wavewind_spfy);
  free(Wavewind_spf);
  free(Wavewind_dgf);
  free(Waveub);
  free(Wavefz);
  free(Wavefphi);
  free(Waveab);
  free(Wavefw);
  free(Waveab_edge);
  free(Wavekmean_edge);
  free(Wavesgmean_edge);
  free(Wavethtamean_edge);
  free(Wavekw_edge);

  for(m=0;m<Mw;m++){
    for(n=0;n<Nw;n++){
      free(Wavecgx[m][n]);
      free(Wavecgy[m][n]);
      free(Wavecs[m][n]);
        free(Wavect[m][n]);
	free(WaveN[m][n]);
	free(WaveNtmp[m][n]);
	free(WaveNold[m][n]);
	free(Wavessrc[m][n]);
	free(Wavetsrc[m][n]);
	free(Wavesrc[m][n]); 
      }
      free(Wavekw[m]);
      free(Wavecgx[m]);
      free(Wavecgy[m]);
      free(Wavecph[m]);
      free(Wavecs[m]);
      free(Wavect[m]);
      free(WaveN[m]);
      free(WaveNtmp[m]);
      free(WaveNold[m]);
      free(Wavetsrc[m]);
      free(Wavessrc[m]);
      free(Wavesrc[m]);

  }
  for(n=0;n<Nw;n++){
    free(Wavecs[Mw][n]);
  }
  free(Wavecs[Mw]);
  free(Wavekw);
  free(Wavecgx);
  free(Wavecgy);
  free(Wavecph);
  free(Wavecs);
  free(Wavect);
  free(WaveN);
  free(WaveNtmp);
  free(WaveNold);

  for(m=0;m<Mw;m++){
    for(n=0;n<Nw;n++){
      free(Wavecg[m][n]);
    }
    free(Wavecg[m]);
  }
  free(Wavecg);
  free(Wavethtaw);
  free(Wavesg);
  free(Wavedsg);
  free(Wavea);
  free(Waveb);
  free(Wavec);
  free(Wavesp);
  free(Wavesm);
  free(Wavessrc);
  free(Wavetp);
  free(Wavetm);
  free(Wavetsrc);
  free(Wavesrc);
  free(WaveEtot);
  free(WaveEtmp);
  free(WaveT0);
  free(WaveHs);
  free(WaveT01);
  free(Wavektail);
  free(Wavesg_PM);
  free(WaveCr);
  free(Wavetmp);
  free(Wavekmean);
  free(Wavethtamean);
  free(WaveEtail);
  free(WaveNtail);
  free(Wavesgmean);
  free(Waveux);
  free(Waveuy);
  free(Waveuz);
  free(WavedivScx);
  free(WavedivScy);
  free(Waveuscx);
  free(Waveuscy);
  free(Waveuse);
  free(WaveUwind);

}

void InitializeWaveVariables(gridT *grid, propT *prop, int myproc, MPI_Comm comm)
{
  int i, j, k, m, n, Nc=grid->Nc, Ne=grid->Ne, Mw=WaveMw, Nw=WaveNw, MN=Max(Mw, Nw);
  int is, nw, nstation = Wavenstation, Nwind=WaveNwind;
  REAL dsw, dtw, logsg_min, logsg_max, sg_face;
  
  for(i=0;i<Nc;i++){
    for (j=0; j<nstation; j++){
      Waveklambda[i][j] = 0;
    }
    for (k=0; k<grid->Nk[i]; k++){
      Wavefz[i][k] = 0;
      Waveuz[i][k] = 0;
      Waveuy[i][k] = 0;
      Waveux[i][k] = 0;
      WavedivScx[i][k] = 0;
      WavedivScy[i][k] = 0;
    }
    Wavefphi[i] = 0;
    Waveuscx[i] = 0;
    Waveuscy[i] = 0;

  } 

  for(j=0;j<Ne;j++){
    for (k = 0; k < grid->Nkc[j]; k++){
      WaveUw[j][k] = 0;
      WavedivSe[j][k] = 0;
    }
    for(m=0; m<Mw; m++){
      Wavekw_edge[j][m] = 0;
    }    
    Waveab_edge[j] = 0;
    Wavesgmean_edge[j] = 0;
    Wavethtamean_edge[j] = 0;
    Wavekmean_edge[j] = 0;
    Waveuse[j] = 0;
    WaveUwind[j] = 0;
  }

  for(m=0; m<Mw; m++){
    for(n=0; n<Nw; n++){
      for(i=0;i<Nc;i++){
	Wavecgx[m][n][i]=0;
        Wavecgy[m][n][i]=0;
        Wavecs[m][n][i]=0;
        Wavect[m][n][i]=0;	
	WaveN[m][n][i]=0;
	WaveNtmp[m][n][i]=0;
	WaveNold[m][n][i]=0;
        Wavessrc[m][n][i]=0;
        Wavetsrc[m][n][i]=0;
        Wavesrc[m][n][i]=0;
      }
    }
  }

  for(m=0; m<Mw; m++){
    for(i=0;i<Nc;i++){
      Wavekw[m][i]=100.0;//1.0/grid->dv[i];  //initial guess for the wave number field
      Wavecph[m][i]=0;
    }
  }


  for(n=0; n<Nw; n++){
    for(i=0;i<Nc;i++){
      Wavecs[Mw][n][i]=0;
    }
  }


  for(m=0; m<Mw; m++){
    for(n=0; n<Nw; n++){
      for(i=0;i<Ne;i++){
	Wavecg[m][n][i]=0;
      }
    }
    Wavesp[m] = 0;
    Wavesm[m] = 0;
  }
  Wavesp[Mw] = 0;
  Wavesm[Mw] = 0;
 
  for (n=0; n<Nw; n++){
    Wavetp[n] = 0;
    Wavetm[n] = 0;

  }  

  for(m=0; m<MN; m++){
    Wavea[m] = 0;
    Waveb[m] = 0;
    Wavec[m] = 0;    
  }
  
  logsg_min = log(Wavesgmin);
  logsg_max = log(Wavesgmax);
  dsw = (logsg_max-logsg_min)/(double)Mw;
  sg_face = Wavesgmin;
  
  for(m=0; m<Mw; m++){
    Wavesg[m]=exp(logsg_min + 0.5*dsw + (double)m*dsw);
    Wavedsg[m] = exp(logsg_min + (double)(m+1)*dsw)
      - exp(logsg_min + (double)(m)*dsw);
    
  }

  //dsw = (Wavesgmax-Wavesgmin)/(double)Mw;
  //  sg_face = Wavesgmin;
  
  //Wavesg[0] = Wavesgmin+0.5*dsw;
  //Wavedsg[0] = dsw;
  //for(m=1; m<Mw; m++){
  //  Wavesg[m] = Wavesg[m-1]+dsw;
  //  Wavedsg[m] = dsw;    
  //}



  dtw = 2*PI/Nw;
  Wavethtaw[0] = dtw;
  for(n=1; n<Nw; n++){
    Wavethtaw[n] = Wavethtaw[n-1]+dtw;
  }

  for(i=0; i<Nc; i++){
    Wavewind_spfx[i] = 0;
    Wavewind_spfy[i] = 0;
    Wavewind_dgf[i] = 0;
    Wavewind_spf[i] = 0;
    WaveHw[i] = 0;
    Waveub[i] = 0;
    Waveab[i] = 0;
    WaveT0[i] = 0;
    WaveT01[i] = 0;
    WaveHs[i] = 0;
    Wavektail[i] = 1000;
    Wavesg_PM[i] = 0;
    WaveCr[i] = 0;
    Wavethtamean[i] = 0;
    WaveEtail[i] = 0;
    WaveNtail[i] = 0;
    WaveEtot[i] = 0;
    WaveEtmp[i] = 0;
    Wavetmp[i] = 0;
    Wavekmean[i] = 0;
    Wavesgmean[i] = 0;
    Wavefw[i] = 0;
  }

  for(is=0; is < nstation; is++){
    for(nw=0; nw<Nwind; nw++){
      Wavewind_sp[is][nw] = 0;
      Wavewind_dg[is][nw] = 0;
    }
  } 
}

void InputWind(int station, propT *prop, int myproc, int numprocs)
{
  char str[BUFFERLENGTH], c, str2[BUFFERLENGTH], str0[BUFFERLENGTH];
  char tmp[BUFFERLENGTH], tmp2[BUFFERLENGTH];
  char  fd[BUFFERLENGTH], filename[BUFFERLENGTH];
  FILE *ifile;
  int Nwind = WaveNwind, i, j, n;
  int Noffset;

  Noffset = floor((prop->nstart*prop->dt)/Wavewind_dt);
  sprintf(fd, "WindFile.%d", station+1);
  MPI_GetFile(filename, DATAFILE, fd, "InputWind",myproc);
  ifile = MPI_FOpen(filename, "r", "InputWind", myproc);
    
  str0[0] = ' ';
  while (str0[0] != '-')

    getline(ifile, str0, "...");

  c = fgetc(ifile);

  if (c == EOF){
    return;
  }

  for(n=0; n<Noffset; n++)
    getline(ifile, str0, "");
  
  c = fgetc(ifile);

  if (c == EOF){
    return;
  }

  
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
    Wavewind_sp[station][n] = strtod(str, (char **)NULL);
    Wavewind_dg[station][n] = strtod(str2, (char **)NULL);
    Wavewind_dg[station][n]*=PI/180;
    c = fgetc(ifile);
    

  }
  fclose(ifile);
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

  int i, m, n, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc;
  REAL dsw, Cd, ustar, H, A, B, coss, tmp, spd;

  dsw = (Wavesgmax-Wavesgmin)/(double)Mw;
  

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	ustar = pow(Waveuscx[i],2)+pow(Waveuscy[i],2);
	ustar = sqrt(ustar);
	if(ustar <= SMALL){
	  Wavesg_PM[i] = 10000.0;
	  H = 0;
	}else{
	  Wavesg_PM[i] = 0.13*GRAV/28/ustar*2*PI;
	  tmp = Wavesg[m]/Wavesg_PM[i];
	  tmp = pow(tmp, -4);
	  H = exp(-tmp);

	}
	H = 1.0;
       coss = cos(Wavethtaw[n]-Wavewind_dgf[i]);
       A = 80/(2*PI)*pow(RHOair/RHO0/GRAV*Wavecph[m][i], 2)/Wavesg[m]*pow(ustar*Max(0, coss), 4);
	B = Max(0, 5*RHOair/RHO0*(Wavewind_spf[i]/Wavecph[m][i]*coss-0.9))*Wavesg[m];
        Wavessrc[m][n][i] = A/Wavesg[m];
	Wavesrc[m][n][i] = B;
      
      }
    }
  }
}


void SinkByWhitecapping(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc;
  REAL dsw, dthta, tmp, Cds, delta, p, Gamma, s_PM, s, alpha, alpha_PM;

  dsw = (Wavesgmax-Wavesgmin)/(double)Mw;
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
	s = Wavekmean[i]*sqrt(WaveEtot[i]);
	if (Wavekmean[i] != 0){
	  //Gamma = Cds*((1-delta)+delta*Wavekw[m][i]/Wavekmean[i])*pow(s/s_PM, p);
	  //Gamma = Cds*pow(s/s_PM, p);
	  //Wavesrc[m][n][i] -= Gamma*Wavesgmean[i]*Wavekw[m][i]/Wavekmean[i];//*WaveN[m][n][i];
	  alpha = WaveEtot[i]*pow(Wavesgmean[i], 4)/pow(GRAV, 2);
	  Gamma = Cds*pow(Wavesg[m]/Wavesgmean[i], 2)*pow(alpha/alpha_PM, 2);
	  //Wavesrc[m][n][i] -= Gamma*Wavesgmean[i]*Wavekw[m][i]/Wavekmean[i];
	  Wavesrc[m][n][i] -= Gamma*Wavesgmean[i];
	}
      }
    }
  }

}

void SinkByWhitecapping_implicit(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, mm, nn, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc, p;
  REAL dthta, trcDM, tmp2, J, dt = prop->dt*Wavewnstep, Cwh, alpha_PM, Nmn, **RHS, DMmn;
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
      for (i= 0; i <  Nc; i++){        
	WaveNtmp[m][n][i] = WaveN[m][n][i] + Wavessrc[m][n][i]*dt*Wavewnstep;
	WaveNtmp[m][n][i] = WaveNtmp[m][n][i]*exp(Wavesrc[m][n][i]*dt*Wavewnstep);\	
      }
    }
  }



  //Second, update the intermediate part of the whitecapping sink.
  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	WaveNold[m][n][i] = WaveNtmp[m][n][i];
	s = Wavekmean[i]*sqrt(WaveEtot[i]);
	if (Wavekmean[i] != 0){
	  alpha = WaveEtot[i]*pow(Wavesgmean[i], 4)/pow(GRAV, 2);
	  Gamma = Cwh*pow(Wavesg[m]/Wavesgmean[i], 2)*pow(alpha/alpha_PM, 2);
	  Wavesrc[m][n][i] = -0.5*Gamma*Wavesgmean[i];
	  WaveNtmp[m][n][i] = WaveNtmp[m][n][i]*exp(Wavesrc[m][n][i]*dt*Wavewnstep);
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
      	if(WaveNtmp[m][n][i] < 0) WaveNtmp[m][n][i] = 0;
	//E(sigma) = INTEGRAL E(sigma, theta)*dtheta = SUM N*sigma*dtheta
      	tmp += WaveNtmp[m][n][i]*dthta;
      }
      //E = INTEGRAL E(sigma)*dsigma
      tmpE+=tmp*Wavesg[m]*Wavedsg[m];      
    }
    WaveEtmp[i] = tmpE;
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
    J = Cwh*pow(Wavesgmean[i], 7)*pow(GRAV, -4)*pow(alpha_PM, -2)*WaveEtmp[i];
    //First, find the trace of A, stored in "tmp", and the RHS
    for (m=0; m < Mw; m++){
      for (n=0; n < Nw; n++){
	Nmn = WaveNold[m][n][i];
	//Nmn = WaveNtmp[m][n][i];
	trcDM += pow(Wavesg[m], 3)*Nmn*Wavedsg[m]*dthta;
	RHS[m][n] = Nmn - 0.25*dt*J*pow(Wavesg[m], 2)*WaveEtot[i]*Nmn;
      }
    }    

    //Second, update N = (I+A)^-1*RHS = RHS + 1/[trace(A)+1]*(-A)*RHS
    for (m=0; m < Mw; m++){
      for (n=0; n < Nw; n++){
	tmp2 = 0.25*dt*J/(1.0 + 0.25*dt*J*trcDM); //Each element in -1/[trace(A0)+1]*A0
	Nmn = RHS[m][n]; //This is equivalent to multiplication of an identity matrix with RHS, stored in a temporary variable Nmn.  
	for (mm = 0; mm < Mw; mm++)
	  for (nn = 0; nn < Nw; nn++){
	    DMmn = WaveNold[m][n][i]*pow(Wavesg[m], 2)*(dthta*Wavesg[mm]*Wavedsg[mm]);
	    Nmn -= tmp2*DMmn*RHS[mm][nn];
	  }

	   
	if (Nmn >= 0.0)	  
	  WaveNtmp[m][n][i] = Nmn;
	else
	  printf("Warning!! Skip negative N at i = %d, m = %d, n = %d, N = %f, RHS = %f\n", i, m, n, Nmn, RHS[m][n]);
      }
    }
  }
  for (m=0; m < Mw; m++)
    for (n=0; n < Nw; n++)
      ISendRecvCellData2D(WaveNtmp[m][n], grid, myproc, comm);      
  
  for (m = 0; m < Mw; m++)
    SunFree(RHS[m], Nw*sizeof(REAL), "SinkByWhitecapping_implicit");
  SunFree(RHS, Mw*sizeof(REAL), "SinkByWhitecapping_implicit");
}


void SinkByMud(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, l, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc, Nl;
  REAL dpth, kd, St_thickness, dm, gamma, cg, Dw, conc, fac, Cb, kwd, fw, dryden;

  dm = Wavebtm_mud_thickness;  
  if (prop->computeSediments){
    dryden=0;
    for(m=0;m<Nsize;m++)
      dryden+=Drydensity[m][0]/Nsize/Gsedi[m];
    conc = dryden*0.000001;
  } else {
    conc = Wavebtm_conc;
  }
  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	dpth = grid->dv[i]+phys->h[i];
	kd = Wavekw[m][i]*dpth;
	if(Wavesg[m] > 0){
	  St_thickness = sqrt(2*Wavebtm_vis/Wavesg[m]);
 
	  if (prop->computeSediments)
	    if (Wavebtm_sedi_erosion)
            {
              dryden=0;
              for(m=0;m<Nsize;m++)
                dryden+=Drydensity[m][0]/Nsize;
	      dm = Thicknesslayer[i][0]/dryden;
	    }
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
	if (Wavebtm_sedi_erosion)
	  gamma = pow(1.0+conc*(sprop->spwght-1.0), -1);
	else
	  gamma = pow(1.0+conc*(2.65-1.0), -1);
	  
	if(kd > SMALL){
	  Dw = St_thickness*gamma*pow(Wavekw[m][i], 2)/(sinh(2*kd)+2*kd)*fac;
          
	  cg = sqrt(pow(Wavecgx[m][n][i], 2) + pow(Wavecgy[m][n][i], 2));
	  Wavesrc[m][n][i] -= 2*Dw*cg;
	}
      }
    }
  }
  
}

void SinkByDrag(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc;
  REAL Cb, fw, dpth, Cf, kwd; 

  
  Cf = 0.15;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){       
	dpth = grid->dv[i]+phys->h[i];	
	if(dpth <= Wavedepth_fw_cutoff)
	  fw = Wavefw_drag;
	else
	  fw = Wavefw[i];	
	Cb = fw*GRAV/sqrt(2)*Waveub[i];	
	kwd = Wavekw[m][i]*(dpth);
	if(kwd > SMALL){
	  Wavesrc[m][n][i]-=Cb*pow(Wavesg[m], 2)*pow(GRAV*sinh(kwd), -2);
	}else{
	  Wavesrc[m][n][i] = 0.0;
	}
      }
    }
  }
  
  
}

void SinkByBreaking(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, Nc=grid->Nc, Nw=WaveNw, Mw=WaveMw;
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
	  Hrms = WaveHs[i]/4.0;
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
	    Dtot = alphaBJ*Qb*Wavesgmean[i]*pow(Hmax, 2)/(8*PI);
	    Wavesrc[m][n][i]-=alphaBJ*Qb*Wavesgmean[i]/pow(beta, 2)/(8*PI);
	  }	  
	}
	  
      }
    }
  }
}


void SourceByTriad(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){
  
  int i, m, n, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc;
  int indx2, indx3, indx1;
  REAL sg1, sg2, Ur, dpth, alpha_EB, J1, J2, dsw, beta;
  REAL sgmax = Wavesgmax, sgmin =Wavesgmin;
  REAL Splus, Sminus, cg, tmp, tmp0;
  REAL E, E05, E2, G, k2;



  dsw = (log(sgmax)-log(sgmin))/(double)Mw;
  alpha_EB = 0.2;

  for (m=0; m < Mw; m++){
    
    indx1 = m + floor(log(0.5)/dsw);
    indx2 = m + floor(log(2)/dsw) + 1;
        
    sg1 = Wavesg[m]/2.0;
    sg2 = Wavesg[m]*2.0;

    for (n=0; n < Nw; n++){
      for(i = 0; i<Nc; i++){
	dpth = grid->dv[i]+phys->h[i];
	if (dpth > 0.1){
	  Ur = GRAV/(8*sqrt(2)*pow(PI, 2))*WaveHs[i]*pow(WaveT01[i], 2)/pow(dpth, 2);
	  if (Ur >=0 && Ur <= 1){
	    if (Ur < 0.02)
	      beta = 0;
	    else
	      beta = -PI/2 + PI/2*tanh(0.2/Ur);
	    
	    E = WaveN[m][n][i]*Wavesg[m];
	    if (indx1 < 0){
	      E05 = 0;
	      J1 = 0;
	    }else{
	      E05 = WaveN[indx1][n][i]*sg1;
	      J1 = pow(Wavekw[indx1][i], 2)*(GRAV*dpth+2*pow(Wavecph[indx1][i], 2))
	      /(Wavekw[m][i]*dpth*(GRAV*dpth + 2/15*GRAV*pow(dpth, 3)*pow(Wavekw[m][i], 2) - 2/5*pow(Wavesg[m]*dpth, 2)));
	    }
	    if (indx2 >= Mw){
	      E2 = WaveN[Mw-1][n][i]*pow(sg2/sgmax, -Wavetail_pow)*sg2;
	      G = pow(sg2, 2)*dpth/GRAV;
	      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
	      tmp = G + 1/tmp;
	      k2 = sg2*sqrt(tmp/(GRAV*dpth));
	      J2 = pow(Wavekw[m][i], 2)*(GRAV*dpth+2*pow(Wavecph[m][i], 2))
	      	/(k2*dpth*(GRAV*dpth + 2/15*GRAV*pow(dpth, 3)*pow(k2, 2) - 2/5*pow(sg2*dpth, 2)));
	    }else{
	      E2 = WaveN[indx2][n][i]*sg2;
	      k2 = Wavekw[indx2][i];
	      J2 = pow(Wavekw[m][i], 2)*(GRAV*dpth+2*pow(Wavecph[m][i], 2))
		/(k2*dpth*(GRAV*dpth + 2/15*GRAV*pow(dpth, 3)*pow(k2, 2) - 2/5*pow(sg2*dpth, 2)));
	    }
	    
	    
	    cg = sqrt(pow(Wavecgx[m][n][i], 2) + pow(Wavecgy[m][n][i], 2));
	    tmp0 = alpha_EB*2*PI*Wavecph[m][i]*cg*pow(J1, 2)*fabs(sin(beta));
	    tmp = tmp0*(pow(E05, 2)-2*E*E05); 
	    if (tmp > 0) Splus = tmp;
	    else Splus = 0.0;
	    tmp0 = alpha_EB*2*PI*Wavecph[m][i]*cg*pow(J2, 2)*fabs(sin(beta));
	    tmp = tmp0*(pow(E, 2)-2*E*E2);
	    if (tmp > 0) Sminus = -2*tmp;
	    else Sminus = 0.0; 
	    Wavessrc[m][n][i] += (Sminus + Splus)/Wavesg[m];
	  }
	    
	}
	
      }
    }
  }
}
void SourceByQuad(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc;
  int nf, nb;
  REAL dsw, lambda = 0.25, sg_pls, Cnl4, sgmax = Wavesgmax, sgmin =Wavesgmin;
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
    sg_pls = Wavesg[m]*(1.0+lambda);
    for (n=0; n < Nw; n++){      
      for(i = 0; i<Nc; i++){
	E1 = WaveN[m][n][i]*Wavesg[m];
	if (m == Mw-1){
	  sgmax = Wavesg[m];
	  E2 = WaveN[m][n][i]*pow(sg_pls/sgmax, -Wavetail_pow)*sg_pls;
	  dEdsg = (WaveN[m][n][i]*Wavesg[m] - WaveN[m-1][n][i]*Wavesg[m-1])
	    /(Wavesg[m]-Wavesg[m-1]);
	  E3 = WaveN[m][n][i]*Wavesg[m] - dEdsg*lambda*Wavesg[m];
	}else if (m == 0){
          dEdsg = (WaveN[m+1][n][i]*Wavesg[m+1] - WaveN[m][n][i]*Wavesg[m])
	    /(Wavesg[m+1]-Wavesg[m]);
	  E2 = WaveN[m][n][i]*Wavesg[m] + dEdsg*lambda*Wavesg[m];
	  E3 = 0;
	}else{
	  dEdsg = (WaveN[m+1][n][i]*Wavesg[m+1] - WaveN[m-1][n][i]*Wavesg[m-1])
	    /(Wavesg[m+1]-Wavesg[m-1]);
	  E2 = WaveN[m][n][i]*Wavesg[m] + dEdsg*lambda*Wavesg[m];
	  E3 = WaveN[m][n][i]*Wavesg[m] - dEdsg*lambda*Wavesg[m];
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
	dEdthta = (WaveN[m][nf][i] - WaveN[m][nb][i])/(2.0*dtw);

	E1pls = E2 + dEdthta*ang1;
	E1mns = E3 + dEdthta*ang2;

	E2pls = E1 + dEdsg*(2.0*lambda + pow(lambda, 2.0))*Wavesg[m]
	  + dEdthta*ang1;
	E2mns = E1 + dEdsg*(-pow(lambda, 2.0))*Wavesg[m] + dEdthta*ang2;

	E3pls = E1 + dEdsg*(-pow(lambda, 2.0))*Wavesg[m] + dEdthta*ang1;
	E3mns = E1 + dEdsg*(-2.0*lambda + pow(lambda, 2.0))*Wavesg[m]
	  + dEdthta*ang2;
	
	coef = Cnl4*pow(2.0*PI, 2)*pow(GRAV, -4)*pow(Wavesg[m]/(2*PI), 11);
	S1 = coef*(pow(E1, 2.0)*(E1pls*pow(1.0+lambda, -4) + E1mns*pow(1.0-lambda, -4))
		   -2*E1*E1pls*E1mns*pow(1.0-lambda*lambda, -4));
	S2 = coef*(pow(E2, 2.0)*(E2pls*pow(1.0+lambda, -4) + E2mns*pow(1.0-lambda, -4))
		   -2*E2*E2pls*E2mns*pow(1.0-lambda*lambda, -4));
	S3 = coef*(pow(E3, 2.0)*(E3pls*pow(1.0+lambda, -4) + E3mns*pow(1.0-lambda, -4))
		   -2*E3*E3pls*E3mns*pow(1.0-lambda*lambda, -4));
	Snl4_1 = 2.0*S1 - S2 - S3;

	E1pls = E2 + dEdthta*ang3;
	E1mns = E3 + dEdthta*ang4;

	E2pls = E1 + dEdsg*(2.0*lambda + pow(lambda, 2.0))*Wavesg[m]
	  + dEdthta*ang3;
	E2mns = E1 + dEdsg*(-pow(lambda, 2))*Wavesg[m] + dEdthta*ang4;


	E3pls = E1 + dEdsg*(-pow(lambda, 2.0))*Wavesg[m] + dEdthta*ang3;
	E3mns = E1 + dEdsg*(-2.0*lambda + pow(lambda, 2.0))*Wavesg[m]
	  + dEdthta*ang4;

	S1 = coef*(pow(E1, 2.0)*(E1pls*pow(1.0+lambda, -4) + E1mns*pow(1.0-lambda, -4))
		   -2*E1*E1pls*E1mns*pow(1.0-lambda*lambda, -4));
	S2 = coef*(pow(E2, 2.0)*(E2pls*pow(1.0+lambda, -4) + E2mns*pow(1.0-lambda, -4))
		   -2*E2*E2pls*E2mns*pow(1.0-lambda*lambda, -4));
	S3 = coef*(pow(E3, 2.0)*(E3pls*pow(1.0+lambda, -4) + E3mns*pow(1.0-lambda, -4))
		   -2*E3*E3pls*E3mns*pow(1.0-lambda*lambda, -4));

	Snl4_2 = 2.0*S1 - S2 - S3;

	Snl4 = Snl4_1 + Snl4_2;

	kp = 0.75*Wavekmean[i];
        if(kp <= 0.5) kp = 0.5;
	dpth = phys->h[i] + grid->dv[i];
	if (dpth <= 0)
	  R = 0;
	else	  
	  R = 1+Csh1/(kp*dpth)*(1-Csh2*kp*dpth)*exp(Csh3*kp*dpth);

	Wavessrc[m][n][i]+=R*Snl4/Wavesg[m];

      }
    }
  }

}


void ObtainTotalEnergy(gridT *grid, physT *phys, propT *prop,MPI_Comm comm, int myproc){
  
  int i, m, n, Nc=grid->Nc, Nw=WaveNw, Mw=WaveMw;
  REAL dsw, dthta, tmp, Nend, mm, dpth, max_wh, bf;
  REAL tmpE, tmpsg, tmpEsg;

  dsw = (Wavesgmax-Wavesgmin)/(double)Mw;
  dthta = 2*PI/(double)Nw;
  mm = Wavetail_pow;

  for(i=0; i< Nc; i++){
    tmpE = 0;
    tmpEsg = 0;
    Nend = 0;
    for(m = 0; m < Mw;m++){
      tmp = 0;
      tmpsg = 0;
      for(n=0; n< Nw; n++){
      	if(WaveN[m][n][i] < 0) WaveN[m][n][i] = 0;
	//E(sigma) = INTEGRAL E(sigma, theta)*dtheta = SUM N*sigma*dtheta
      	tmp += WaveN[m][n][i]*dthta;
	//Esigma(sigma) = INTEGRAL E(sigma, theta)*sigma*dtheta = SUM N*sigma*sigma*dtheta
	tmpsg += WaveN[m][n][i]*dthta;
      }
      //E = INTEGRAL E(sigma)*dsigma
      tmpE+=tmp*Wavesg[m]*Wavedsg[m];
      //E = INTEGRAL Esigma(sigma)*dsigma
      tmpEsg +=tmp*Wavesg[m]*Wavesg[m]*Wavedsg[m];
      if (m == Mw-1)
	Nend = tmp;      
    }
    WaveEtot[i] = tmpE;
    WaveT01[i] = tmpE/tmpEsg*2*PI; //The tail is not considered while obtaining T01    

    if(Wavetail_opt){
      WaveEtail[i] = Nend*pow(Wavesgmax,mm)*1.0/(2.0-mm)
	*(pow(Wavesg99, 2.0-mm)-pow(Wavesgmax, 2.0-mm));
      WaveNtail[i] = WaveEtail[i]/Wavesgtail;
      WaveEtot[i]+=WaveEtail[i];
    }

    WaveHs[i] = 4*sqrt(WaveEtot[i]);

    if (Wavedepth_brk_cutoff){
      dpth = grid->dv[i]+phys->h[i];
      max_wh = Wavedepth_brk_indx*dpth;
      bf = tmpE/pow(max_wh, 2.0);
      if (bf > 1.0){
	WaveEtot[i] = pow(max_wh, 2.0);
	WaveHs[i] = 4.0*max_wh;
	for (m = 0; m < Mw; m++){
	  for (n = 0; n < Nw; n++){
	    WaveN[m][n][i] /= bf;
	  }
	}
      }
    }
  }  
}


void ObtainMeanKSG(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc){

  int i, m, n, Nc=grid->Nc, Nw=WaveNw, Mw=WaveMw;
  REAL dthta, tmp, invE, tmpsinh, dpth, m0, tmpthta, tail, kd;

  dthta = 2*PI/(double)Nw;

  for(i=0; i< Nc; i++){
    Wavekmean[i] = 0;
    Wavesgmean[i] = 0;
    Waveub[i] = 0;
    Waveab[i] = 0;
    Wavethtamean[i] = 0;
    WaveT0[i] = 0;
    dpth = grid->dv[i]+phys->h[i];  
    if (WaveEtot[i] > SMALL){
      invE = 1/WaveEtot[i];
      m0 = 0;
      for(m = 0; m < Mw;m++){
	tmp = 0;
	for(n=0; n< Nw; n++){
	  //E(sigma) = INTEGRAL E(sigma, theta) dthta
	  //= INTERGRAL N(sigma, theta)*sigma dthta
	  tmp += WaveN[m][n][i]*dthta;
	}
	Wavesgmean[i]+=tmp*Wavedsg[m];

	tmp *= Wavesg[m];
	Wavekmean[i] += tmp*Wavedsg[m]*1/sqrt(Wavekw[m][i]);
	m0 += tmp*Wavedsg[m];	

	if (dpth > 0){
	  kd = Wavekw[m][i]*(dpth);
	  tmpsinh = pow(sinh(kd), 2);
	  tmpsinh = 1/tmpsinh;
          Waveab[i] += tmpsinh*tmp*Wavedsg[m];
	  Waveub[i] += tmpsinh*tmp*pow(Wavesg[m], 2)*Wavedsg[m];
	}

      }

      if(Wavetail_opt){       
	Wavesgmean[i]+=WaveNtail[i];
	Wavekmean[i]+=WaveEtail[i]*1/sqrt(Wavektail[i]);
	m0 += WaveEtail[i]*(Wavesg99-Wavesgmax);
	if (dpth > 0){
	  kd = Wavektail[i]*(dpth);
	  tmpsinh = pow(sinh(kd), 2);
	  tmpsinh = 1/tmpsinh;
          Waveab[i] += tmpsinh*WaveEtail[i];
	  Waveub[i] += tmpsinh*WaveEtail[i]*pow(Wavesgtail, 2.0);
	}	
      }
            
      Wavesgmean[i]*= invE;
      Wavesgmean[i] = 1/Wavesgmean[i];      
      
      Wavekmean[i]*= invE;
      Wavekmean[i] = pow(Wavekmean[i], -2);
      
      WaveT0[i] = m0;
      
      Waveab[i] = Waveab[i]*2;
      Waveab[i] = sqrt(Waveab[i]);

      //Waveub[i] = Waveub[i]*2;
      //Waveub[i] = sqrt(Waveub[i]);
      Waveub[i] = Waveab[i]*Wavesgmean[i];
      

      for(n=0; n< Nw; n++){
	tmp = 0;
	for(m = 0; m < Mw;m++){
	  tmp+=WaveN[m][n][i]*Wavesg[m]*Wavedsg[m];
	}
	Wavethtamean[i]+=tmp*Wavethtaw[n]*dthta;
      }
      Wavethtamean[i]*=invE;
    }
  }  

}


//Obtain the wave number field corresponding to each intrinsic frequency using the
//Taylor series expansion, as well as find the phase speed by k/omega. 
void ObtainKField(gridT *grid, physT *phys, propT *prop) 
{
  int i,m;
  REAL dpth, G, tmp, cph;


  for (m=0; m < WaveMw; m++){
    for(i = 0; i < grid->Nc; i++){
      if(phys->h[i] < -grid->dv[i])
	dpth = 0.1;
      else
	dpth = grid->dv[i]+phys->h[i];
      G = pow(Wavesg[m], 2)*dpth/GRAV;
      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
      tmp = G + 1/tmp;
      Wavecph[m][i] = 1.0/sqrt(tmp/(GRAV*dpth));
      Wavekw[m][i] = Wavesg[m]/Wavecph[m][i];
    }
  }
  if (Wavetail_opt)
    for(i = 0; i < grid->Nc; i++){
      if(phys->h[i] < -grid->dv[i])
	dpth = 0.1;
      else
	dpth = grid->dv[i]+phys->h[i];
      G = pow(1.5*Wavesgtail, 2)*dpth/GRAV;
      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
      tmp = G + 1/tmp;
      cph= 1.0/sqrt(tmp/(GRAV*dpth));
      Wavektail[i] = Wavesgtail/cph;    
    }
 
 
}
void ObtainEdgeKField(gridT *grid, physT *phys, propT *prop) 
{
  int j, m, Ne = grid->Ne, Mw = WaveMw;
  REAL dpth, G, tmp, cph;

  for(j = 0; j < Ne; j++){
    dpth = phys->de[j];
    if (dpth < 0.1) dpth = 0.1;    
    for (m=0; m < Mw; m++){      
      G = pow(Wavesg[m], 2)*dpth/GRAV;
      tmp = 1 + 0.6522*G + 0.4622*pow(G, 2) + 0.0864*pow(G, 4) + 0.0675*pow(G, 5);
      tmp = G + 1/tmp;
      cph = 1.0/sqrt(tmp/(GRAV*dpth));
      Wavekw_edge[j][m] = Wavesg[m]/cph;
    }
  } 
}

void ObtainCenterCgField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int i, m, n, Mw=WaveMw,Nw=WaveNw, Nc=grid->Nc;
  REAL tmp, depth, kw, kd;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(i = 0; i < Nc; i++){
	if(phys->h[i] < -grid->dv[i])
	  depth = 0.1;
	else
	  depth = phys->h[i] + grid->dv[i];
	
	kw = Wavekw[m][i];
	kd = kw*depth;
	
	if(kd > 0.1)
	  tmp = 0.5*(1+2*kd/sinh(2*kd))*Wavesg[m]/kw;
	else
	  tmp = sqrt(GRAV*depth);

	Wavecgx[m][n][i] = tmp*cos(Wavethtaw[n]) + phys->uc[i][grid->ctop[i]];
	Wavecgy[m][n][i] = tmp*sin(Wavethtaw[n]) + phys->vc[i][grid->ctop[i]];

      }
      ISendRecvCellData2D(Wavecgx[m][n], grid, myproc, comm);
      ISendRecvCellData2D(Wavecgy[m][n], grid, myproc, comm);
    }
  }
    
}

void ObtainEdgeCgField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int j, m, n, nc1, nc2, Mw=WaveMw,Nw=WaveNw, Ne=grid->Ne;
  REAL tmp, depth, kw;

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      for(j = 0; j<Ne; j++){
	nc1 = grid->grad[2*j];
	nc2 = grid->grad[2*j+1];
	
	if(nc1 == -1){
	  Wavecg[m][n][j] = Wavecgx[m][n][nc2]*grid->n1[j]+Wavecgy[m][n][nc2]*grid->n2[j];	  
	}else if(nc2 == -1){
	  Wavecg[m][n][j] = Wavecgx[m][n][nc1]*grid->n1[j]+Wavecgy[m][n][nc1]*grid->n2[j];
	}else{
	  Wavecg[m][n][j] = InterpCgToFace(m, n, j, grid);
	}
      }
      ISendRecvEdgeData2D(Wavecg[m][n], grid, myproc, comm);
     
    }
  }
    
}

static REAL InterpCgToFace(int m, int n, int j, gridT *grid){
  int nc1, nc2;
  REAL def1, def2, Dj, Cgnc1, Cgnc2;

  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  Dj = grid->dg[j];
  def1 = grid->def[nc1*NFACES+grid->gradf[2*j]];
  def2 = grid->def[nc2*NFACES+grid->gradf[2*j+1]];
  Cgnc1 = Wavecgx[m][n][nc1]*grid->n1[j] + Wavecgy[m][n][nc1]*grid->n2[j];
  Cgnc2 = Wavecgx[m][n][nc2]*grid->n1[j] + Wavecgy[m][n][nc2]*grid->n2[j];

  //  if (def1 == 0 || def2 == 0)
  // return UpWind(Wavecg[m][n][j], Cgnc1, Cgnc2);
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
      Waveab_edge[j] = Waveab[nc2];
      Wavekmean_edge[j] = Wavekmean[nc2];
      Wavesgmean_edge[j] = Wavesgmean[nc2];
      Wavethtamean_edge[j] = Wavethtamean[nc2];      
    }else if(nc2 == -1){
      Waveab_edge[j] = Waveab[nc1];
      Wavekmean_edge[j] = Wavekmean[nc1];
      Wavesgmean_edge[j] = Wavesgmean[nc1];
      Wavethtamean_edge[j] = Wavethtamean[nc1];
    }else{
      Waveab_edge[j] = InterpWavePropToFace(j, Waveab, grid);
      Wavekmean_edge[j] = InterpWavePropToFace(j, Wavekmean, grid);
      Wavesgmean_edge[j] = InterpWavePropToFace(j, Wavesgmean, grid);
      Wavethtamean_edge[j] = InterpWavePropToFace(j, Wavethtamean, grid);
    }
  }
  ISendRecvEdgeData2D(Waveab_edge, grid, myproc, comm);
  ISendRecvEdgeData2D(Wavekmean_edge, grid, myproc, comm);
  ISendRecvEdgeData2D(Wavesgmean_edge, grid, myproc, comm);
  ISendRecvEdgeData2D(Wavethtamean_edge, grid, myproc, comm);

}

static REAL InterpWavePropToFace(int j, REAL *value, gridT *grid){
  int nc1, nc2;
  REAL def1, def2, Cgnc1, Cgnc2;

  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  def1 = grid->def[nc1*NFACES+grid->gradf[2*j]];
  def2 = grid->def[nc2*NFACES+grid->gradf[2*j+1]];
  //  if (def1 == 0 || def2 == 0)
  // return UpWind(Wavecg[m][n][j], Cgnc1, Cgnc2);
  //else
  return (value[nc1]*def2+value[nc2]*def1)/(def1+def2);

}


void ObtainCenterCsCtField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{

  int i, m, n, nf, ne, normal, k, nc1, nc2;
  REAL Ac, df, kw, tmp0, tmp1, tmp00, tmp11, cost, sint, cost2, sint2, tanhkd;
  REAL dudx, dudy, dvdx, dvdy, dddx, dddy, depth;
  int Nc=grid->Nc, Mw=WaveMw, Nw=WaveNw;
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

    for(nf=0; nf<NFACES; nf++){
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
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
	cost = cos(Wavethtaw[n]);
	sint = sin(Wavethtaw[n]);
	for (m = 0; m < Mw; m++){
	  kw = Wavekw[m][i];
	  tanhkd = tanh(kw*depth);
	  tmp1 = 0.5*sqrt(GRAV)*pow(kw, 1.5)*(1-pow(tanhkd, 2))/sqrt(tanhkd);	
	  if (m >= 1){
	    
	    tmp11 = tmp1*(phys->dhdt[i]+phys->uc[i][grid->ctop[i]]*dddx+
			  +phys->vc[i][grid->ctop[i]]*dddy)  
	    //	    tmp11 = tmp1*phys->dhdt[i]
	      - sqrt(pow(Wavecgx[m][n][i], 2)+pow(Wavecgy[m][n][i], 2))
	      * kw*(cost*(dudx*cost+dudy*sint)
		    + sint*(dvdx*cost+dvdy*sint));
	    
	    tmp00 = tmp0*(phys->dhdt[i]+phys->uc[i][grid->ctop[i]]*dddx+
			  +phys->vc[i][grid->ctop[i]]*dddy)
	      //tmp00 = tmp0*phys->dhdt[i]
		       - sqrt(pow(Wavecgx[m-1][n][i], 2)+pow(Wavecgy[m-1][n][i], 2))
	      * kw*(cost*(dudx*cost+dudy*sint)
		    + sint*(dvdx*cost+dvdy*sint));
	    Wavecs[m][n][i] = 0.5*(tmp11 + tmp00);
	  }else{
	    	    Wavecs[m][n][i] = tmp1*(phys->dhdt[i]+phys->uc[i][grid->ctop[i]]*dddx+
	    		      +phys->vc[i][grid->ctop[i]]*dddy)
	    // Wavecs[m][n][i] = tmp1*phys->dhdt[i]
	      - sqrt(pow(Wavecgx[m][n][i], 2)+pow(Wavecgy[m][n][i], 2))
	      * kw*(cost*(dudx*cost+dudy*sint)
		    + sint*(dvdx*cost+dvdy*sint));
	  }
	  tmp0 = tmp1;
	}
      }
      for (n = 0; n < Nw; n++){
	Wavecs[0][n][i] = 2*Wavecs[1][n][i] - Wavecs[2][n][i];
	Wavecs[Mw][n][i] = 2*Wavecs[Mw-1][n][i] - Wavecs[Mw-2][n][i];
      }
    }

    if(depth > 0.1){
      for (m = 0; m < Mw; m++){
	kw = Wavekw[m][i];
	tanhkd = tanh(kw*depth);
	tmp1 = 0.5*sqrt(GRAV)*pow(kw, 1.5)*(1-pow(tanhkd, 2))/sqrt(tanhkd);
	for (n = 0; n < Nw; n++){
	  //Edge value used for theta velocity
	  cost2 = cos(Wavethtaw[n]+0.5*dthta);
	  sint2 = sin(Wavethtaw[n]+0.5*dthta);
	  
	  Wavect[m][n][i] = -1/kw*(tmp1*(-sint2*dddx + cost2*dddy) 
				     + kw*(cost2*(-dudx*sint2+dudy*cost2)
					   + sint2*(-dvdx*sint2+dvdy*cost2)));
	}
      }
    }
  }
    
}

void UpdateWave(gridT *grid, physT *phys, propT *prop, sediT *sedi, spropT *sprop, 
MPI_Comm comm, int myproc, int numprocs){
  
  int i;

  FetchModel=MPI_GetValue(DATAFILE,"Fetchmodel","UpdateWave",myproc);
  if(FetchModel){
    if(prop->n == 1+prop->nstart){
      FetchReadWaveProperties(myproc);
      FetchOpenWaveFiles(prop->mergeArrays,myproc); 
      FetchAllocateWave(grid,myproc);
      FetchInitializeWave(grid,phys,myproc);
    }
    if(ConstantWind==0){
      FetchWindSpeedandDirection(grid,prop,myproc); // add new functions in Boundaries.c
      FetchCalculateFetch(grid,myproc);
      FetchCalculateHwsigTwsig(grid,phys,myproc);
      FetchCalculateWaveexcur(grid,phys,myproc);
      FetchCalculateFw(grid,myproc);
    }

    OutputWave(grid,phys,prop,myproc,numprocs,blowup,comm); 

    //if(prop->n==prop->nstart+prop->nsteps)
     //FetchFreeWave(grid,myproc);
  } else {
    if(prop->n == 1+prop->nstart){
      InitializeWaveProperties(prop, myproc);
      AllocateWaveVariables(grid, prop);
      OpenWaveFiles(myproc);
      InitializeWaveVariables(grid, prop, myproc, comm);

      if (RESTART)
        ReadWaveVariables(grid, prop, myproc, comm);
    
      
      if(Wavewind_forcing){
        ObtainKrigingCoef(grid, myproc, numprocs);
        for(i = 0; i < Wavenstation; i++)
	  InputWind(i, prop, myproc, numprocs);
      }
    }

    if(Wavewind_forcing){
      WindField(prop, grid, phys, myproc, numprocs);
      if(Wavewind_shear)
        WindSurfaceShear(grid, phys, prop, comm, myproc);
    }
    if ((prop->n-1) % Wavewnstep == 0){    
      ObtainKField(grid, phys, prop);
      ObtainEdgeKField(grid, phys, prop); 
      ObtainCenterCgField(grid, phys, prop, comm, myproc);
      ObtainEdgeCgField(grid, phys, prop, comm, myproc);
      ObtainCenterCsCtField(grid, phys, prop, comm, myproc);
      if (Wavewind_forcing)
        SourceByWind(grid, phys, prop, comm, myproc);
      if (Waveform_drag)
        SinkByDrag(grid, phys, prop, comm, myproc);
      if (Wavebtm_mud)
        SinkByMud(grid, phys, prop, sedi, sprop, comm, myproc);
      if (WaveNLquad)
        SourceByQuad(grid, phys, prop, comm, myproc);
      if (WaveNLtriad)
        SourceByTriad(grid, phys, prop, comm, myproc);
      if (WaveBRKdepth)
        SinkByBreaking(grid, phys, prop, comm, myproc);  
      if (Waveimplicit_whitecap)
        SinkByWhitecapping_implicit(grid, phys, prop, comm, myproc);
      else{
        SinkByWhitecapping(grid, phys, prop, comm, myproc);
        UpdateActionDensitySinkSource(grid, phys, prop, comm, myproc, numprocs);
      }
    
      if (Waveimplicit_advection)
        ImplicitUpdateGeographic(grid, phys, prop, comm, myproc, numprocs);
      else
        ExplicitUpdateGeographic(grid, phys, prop, comm, myproc, numprocs); 
    
      UpdateActionDensitySpectral(grid, phys, prop, comm, myproc);
      ObtainTotalEnergy(grid, phys, prop,comm, myproc);
      ObtainMeanKSG(grid, phys, prop, comm, myproc);
      ObtainWaveVelocity(grid, phys, prop, comm, myproc);
      if(Waverad_stress){
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
  
  int i, m, n, Nc=grid->Nc, Mw=WaveMw, Nw=WaveNw, MN=Max(WaveMw, WaveNw);
  REAL dsw, fct, dang;
  REAL x[MN], a[MN], b[MN], c[MN], d[MN];

  dsw = (Wavesgmax-Wavesgmin)/(double)Mw;



  //First update accident density in the frequency space
  for(i=0; i < Nc; i++){
    for(n=0; n< Nw; n++){

      //Compute the source term for the frequency discretization
      //based on the Crank-Nicolson scheme
      for(m=0; m < Mw; m++){
	fct = prop->dt*Wavewnstep/(2*Wavedsg[m]);
	a[m] = fct*Wavesp[m];
	b[m] = 1+(fct*Wavesm[m]-fct*Wavesp[m+1]);
	c[m] = -fct*Wavesm[m+1];
	if(m == 0)
	  d[m] = b[m]*WaveN[m][n][i]+c[m]*WaveN[m+1][n][i];
	else if(m == Mw-1){
	  d[m] = a[m]*WaveN[m-1][n][i]+b[m]*WaveN[m][n][i];
	  d[Mw-1] += c[Mw-1]*WaveN[m][n][i]*pow(1+0.5*Wavedsg[m]/Wavesgmax, -Wavetail_pow);
	}else
	  d[m] = a[m]*WaveN[m-1][n][i]+b[m]*WaveN[m][n][i]
                              + c[m]*WaveN[m+1][n][i];
      }

      //Update upwind velocity operators for the current time step
      for(m=0; m < Mw+1; m++){
	Wavesp[m] = 0.5*(Wavecs[m][n][i]+fabs(Wavecs[m][n][i]));
	Wavesm[m] = 0.5*(Wavecs[m][n][i]-fabs(Wavecs[m][n][i]));	  
      }

      //Generate matrix operator for the matrix solver
      for(m=0; m < Mw; m++){
	fct = prop->dt*Wavewnstep/(2*Wavedsg[m]);
	a[m] = -fct*Wavesp[m];
	b[m] = 1-(fct*Wavesm[m]-fct*Wavesp[m+1]);
	c[m] = fct*Wavesm[m+1];
	x[m] = 0;
      }
      a[0]= 0;
      d[Mw-1] -= c[Mw-1]*WaveN[Mw-1][n][i]*pow(1+0.5*Wavedsg[m]/Wavesgmax, -Wavetail_pow);
      c[Mw-1] = 0;
      
      TriSolve(a, b, c, d, (&x[0]), Mw);

      //Copy to the action density field
      for(m=0; m < Mw; m++){
      	WaveN[m][n][i] = x[m];
      }
    }
  }
  

 
  dang = 360/Nw;
  fct = prop->dt*Wavewnstep/(2*dang);

  for(i=0; i < Nc; i++){
    for(m=0; m< Mw; m++){
      for(n=0; n < Nw; n++){
	if (n == 0){
	  a[n] = fct*Wavetp[Nw-1];
	  b[n] = 1+(fct*Wavetm[Nw-1]-fct*Wavetp[n]);
	}else{
	  a[n] = fct*Wavetp[n-1];
	  b[n] = 1+(fct*Wavetm[n-1]-fct*Wavetp[n]);
	}
	c[n] = -fct*Wavetm[n];
	if(n == 0){
	  d[n] = b[n]*WaveN[m][n][i]+c[m]*WaveN[m][n+1][i];
	  d[n] += a[n]*WaveN[m][Nw-1][i]; //Source due to periodicity
  	}else if(n == Nw-1){
  	  d[n] = a[n]*WaveN[m][n-1][i]+b[m]*WaveN[m][n][i];
	  d[n] += c[n]*WaveN[m][1][i];   //Source due to periodicity
	}else
  	  d[n] = a[n]*WaveN[m][n-1][i]+b[m]*WaveN[m][n][i]
	    + c[n]*WaveN[m][n+1][i];
      }

      for(n=0; n < Nw; n++){
	Wavetp[n] = 0.5*(Wavect[m][n][i]+fabs(Wavect[m][n][i]));
	Wavetm[n] = 0.5*(Wavect[m][n][i]-fabs(Wavect[m][n][i]));	  
      }

      for(n=0; n < Nw; n++){
	if (n == 0){
	  a[n] = -fct*Wavetp[Nw-1];
	  b[n] = 1-(fct*Wavetm[Nw-1]-fct*Wavetp[n]);
	  d[n] -= a[n]*WaveN[m][Nw-1][i]; //Source due to periodicity
	}else{
	  a[n] = -fct*Wavetp[n-1];
	  b[n] = 1-(fct*Wavetm[n-1]-fct*Wavetp[n]);
	  d[n] -= c[n]*WaveN[m][1][i];   //Source due to periodicity
	}
	c[n] = fct*Wavetm[n];
	x[n] = 0;
	
      }
      TriSolve(a, b, c, d, (&x[0]), Nw);
      //Copy to the action density field
      for(n=0; n < Nw; n++){
        WaveN[m][n][i] = x[n];
      }
    }    
  }

  for (m=0; m < Mw; m++){
    for (n=0; n < Nw; n++){
      ISendRecvCellData2D(WaveN[m][n], grid, myproc, comm);
    }
  }
  
  
}

void UpdateActionDensitySinkSource(gridT *grid, physT *phys, propT *prop,MPI_Comm comm, int myproc, int numprocs)
{
  int i, m, n, Nc = grid->Nc, Mw=WaveMw, Nw=WaveNw;
  REAL dt = prop->dt;

  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      //      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //	i = grid->cellp[iptr];
      for (i= 0; i <  Nc; i++){
	WaveNtmp[m][n][i] = WaveN[m][n][i] + Wavessrc[m][n][i]*dt*Wavewnstep;
	WaveNtmp[m][n][i] = WaveNtmp[m][n][i]*exp(Wavesrc[m][n][i]*dt*Wavewnstep);
      }
      ISendRecvCellData2D(WaveNtmp[m][n], grid, myproc, comm);      
    }
  }
}


void ExplicitUpdateGeographic(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs)
{
  int i, iptr, m, nf, ne, n, nc1, nc2, Nc = grid->Nc, Mw=WaveMw, Nw=WaveNw, normal, flag = 0;
  REAL dt = prop->dt, source, Ac, df, dg;
  
  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //i = grid->cellp[iptr];
      for(i = 0; i < Nc; i++){
	Ac = grid->Ac[i];

	source = 0;
	flag = 0;
	for(nf=0;nf<NFACES;nf++) {
	  ne = grid->face[i*NFACES+nf];
	  normal = grid->normal[i*NFACES+nf];
	  df = grid->df[ne];
	  dg = grid->dg[ne];
	  nc1 = grid->grad[2*ne];
	  nc2 = grid->grad[2*ne+1];
	  if(nc1==-1){
	    source -= 0.5*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]+fabs(Wavecg[m][n][ne]))*
	    				   WaveNtmp[m][n][nc2]);
	    flag = 1;
	  }else if(nc2==-1){
	    source -= 0.5*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]-fabs(Wavecg[m][n][ne]))*
	    				   WaveNtmp[m][n][nc1]);
	    flag = 1;
	  }else{
	    source -= 0.5*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]+fabs(Wavecg[m][n][ne]))*
					   WaveNtmp[m][n][nc2]+
					   (Wavecg[m][n][ne]-fabs(Wavecg[m][n][ne]))*
					   WaveNtmp[m][n][nc1]);
	  }  
	}
	
	WaveN[m][n][i] = WaveNtmp[m][n][i] + source; 
       
	
      }
      ISendRecvCellData2D(WaveN[m][n], grid, myproc, comm);      
    }
  }
}

void ImplicitUpdateGeographic(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs)
{  
  int i, iptr, m, nf, ne, n, nc1, nc2, Nc = grid->Nc, Mw=WaveMw, Nw=WaveNw, normal, flag = 0;
  REAL dt = prop->dt, source, Ac, df, dg;
  
  for (m=0; m<Mw; m++){
    for (n=0; n<Nw; n++){
      //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      //i = grid->cellp[iptr];
      for(i = 0; i < Nc; i++){
	Ac = grid->Ac[i];

	source = 0;
	flag = 0;
	for(nf=0;nf<NFACES;nf++) {
	  ne = grid->face[i*NFACES+nf];
	  normal = grid->normal[i*NFACES+nf];
	  df = grid->df[ne];
	  dg = grid->dg[ne];
	  nc1 = grid->grad[2*ne];
	  nc2 = grid->grad[2*ne+1];
	  if(nc1==-1){
	    source -= 0.25*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]+fabs(Wavecg[m][n][ne]))*
	    				   WaveNtmp[m][n][nc2]);
	    flag = 1;
	  }else if(nc2==-1){
	    source -= 0.25*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]-fabs(Wavecg[m][n][ne]))*
	    				   WaveNtmp[m][n][nc1]);
	    flag = 1;
	  }else{
	    source -= 0.25*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]+fabs(Wavecg[m][n][ne]))*
					   WaveNtmp[m][n][nc2]+
					   (Wavecg[m][n][ne]-fabs(Wavecg[m][n][ne]))*
					   WaveNtmp[m][n][nc1]);
	  }  
	}
	
	WaveN[m][n][i] = WaveNtmp[m][n][i] + source; 
	WaveNold[m][n][i] = WaveN[m][n][i];
      }
      ISendRecvCellData2D(WaveN[m][n], grid, myproc, comm);      
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
  int Nc=grid->Nc, Ns = Wavenstation, Ne = grid->Ne;
  int i, j;
  int N;
  REAL r1, r2, sx1, sx2,sy1,sy2, dg1, dg2, rtime_rel, Cd, dg1xy, dg2xy;

  rtime_rel = prop->rtime-(double)prop->nstart*prop->dt;
  N =  floor(rtime_rel/Wavewind_dt);
  r1 = (double)((int)rtime_rel % (int)Wavewind_dt);
  r1/=Wavewind_dt;
  r2 = 1-r1;
  
  if (N+1 < WaveNwind){
    for (i=0; i< Nc; i++){
      Wavewind_spfx[i] = 0;
      Wavewind_spfy[i] = 0;
      Wavewind_dgf[i] = 0;
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
        dg1xy = 0.5*PI - Wavewind_dg[j][N]+PI;
	dg2xy = 0.5*PI - Wavewind_dg[j][N+1]+PI;
  
	sx1 += Waveklambda[i][j]*Wavewind_sp[j][N]*cos(dg1xy);
	sx2 += Waveklambda[i][j]*Wavewind_sp[j][N+1]*cos(dg2xy);
	sy1 += Waveklambda[i][j]*Wavewind_sp[j][N]*sin(dg1xy);
	sy2 += Waveklambda[i][j]*Wavewind_sp[j][N+1]*sin(dg2xy);
	dg1 += Waveklambda[i][j]*dg1xy;
	dg2 += Waveklambda[i][j]*dg2xy;
      }
      Wavewind_spfx[i] = sx1*r2 + sx2*r1;
      Wavewind_spfy[i] = sy1*r2 + sy2*r1;
      Wavewind_spf[i] = sqrt(pow(Wavewind_spfx[i], 2) + pow(Wavewind_spfy[i], 2));

      Wavewind_dgf[i] = asin(Wavewind_spfy[i]/Wavewind_spf[i]);
      if(Wavewind_spfx[i] < 0)
	Wavewind_dgf[i] = PI - Wavewind_dgf[i];

      if(Wavewind_spf[i] < 7.5)
       	Cd = 1.2875*0.001;
      else
       	Cd = (0.8 + 0.065*Wavewind_spf[i])*0.001;	

      Waveuscx[i] = sqrt(Cd)*Wavewind_spf[i]*cos(Wavewind_dgf[i]);
      Waveuscy[i] = sqrt(Cd)*Wavewind_spf[i]*sin(Wavewind_dgf[i]);

     

    }
  }else 
    if(N+1 == WaveNwind){
      for (i=0; i< Nc; i++){
	Wavewind_spfx[i] = 0;
	Wavewind_spfy[i] = 0;
	Wavewind_dgf[i] = 0;
	for (j=0; j< Ns; j++){
	  dg1xy = 0.5*PI - Wavewind_dg[j][N]+PI;
	  Wavewind_spfx[i] += Waveklambda[i][j]*Wavewind_sp[j][N]*cos(dg1xy);
	  Wavewind_spfy[i] += Waveklambda[i][j]*Wavewind_sp[j][N]*sin(dg1xy);
	  Wavewind_dgf[i] += Waveklambda[i][j]*dg1xy;
	}
	Wavewind_spf[i] = sqrt(pow(Wavewind_spfx[i], 2) + pow(Wavewind_spfy[i], 2));

	Wavewind_dgf[i] = asin(Wavewind_spfy[i]/Wavewind_spf[i]);
	if(Wavewind_spfx[i] < 0)
	  Wavewind_dgf[i] = PI - Wavewind_dgf[i];

	if(Wavewind_spf[i] < 7.5)
	  Cd = 1.2875*0.001;
	else
	  Cd = (0.8 + 0.065*Wavewind_spf[i])*0.001;
	
	Waveuscx[i] = sqrt(Cd)*Wavewind_spf[i]*cos(Wavewind_dgf[i]);
	Waveuscy[i] = sqrt(Cd)*Wavewind_spf[i]*sin(Wavewind_dgf[i]);

      }
    }
  if (Wavewind_shear == 1){
    for (j=0; j< Ne; j++)            
      phys->tau_T[j] = InterpWindShearToFace(j, phys, grid);
  }
}

void WaveTurbMixing(gridT *grid, physT *phys, sediT *sedi, propT *prop, MPI_Comm comm, int myproc){

  int Nc = grid->Nc;
  int i,k;
  REAL lz, lz0, z, depth, fzprime;

  for(i = 0; i < Nc; i++){
    Wavefphi[i] = 1.22+0.22*cos(2*WaveCr[i]);
    if (Waveab[i] > 0.02){
      depth = grid->dv[i]+phys->h[i];
      lz0 = log10(sedi->kb[i]/30*Wavesgmean[i]/Waveub[i]);
      z = 0;
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	z -= 0.5*grid->dzz[i][k];
	if (z < 0.0){
	  lz = log(fabs(z)*Wavesgmean[i]/fabs(Waveub[i]));
	  Wavefz[i][k] = -0.0488 + 0.02917*lz + 0.01703*pow(lz, 2)
	    + (1.125*(lz0+5)+0.125*pow((lz0+5), 4))*(-0.0102-0.00253*lz+0.00273*pow(lz, 2));
	  fzprime = 0.02917 + 2.0*0.01703*lz
	    + (1.125*(lz0+5)+0.125*pow((lz0+5), 4))*(-0.00253+2.0*0.00273*lz);	  
	  if(Wavefz[i][k] < 0) Wavefz[i][k] = 0; 
	  if(fzprime <= 0) Wavefz[i][k] = 0;
	}else{
	  Wavefz[i][k] = 0;
	}
	z -= 0.5*grid->dzz[i][k];
      }
    }else{
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	Wavefz[i][k] = 0; 
      }
    }
      
  }

   
}


void ObtainKrigingCoef(gridT *grid, int myproc, int numprocs)
{
  int Nc = grid->Nc, Ns=Wavenstation;
  int i, j, k, jj;

  REAL covMax = 0.78, Dmax = 100000, D, r;
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
      D = sqrt(pow(Wavexw[i]-Wavexw[j], 2)+
	       pow(Waveyw[i]-Waveyw[j], 2));
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
      D = sqrt(pow(grid->xv[i]-Wavexw[j], 2)+
	       pow(grid->yv[i]-Waveyw[j], 2));
      b[j] = semivariogram(covMax, Dmax, D);      
    }
    for (j = 0; j < Ns; j++ ){
      Waveklambda[i][j] = 0;
      for (k = 0; k < Ns; k++ ){
	Waveklambda[i][j]+=H[j][k]*b[k];
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
  x = WaveN[m0][n0];
  r = WaveNtmp[m0][n0];
  p = WaveNold[m0][n0];
  



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
    if(eps0 < epsW)
      printf("Step %d, norm of action density (%d, %d) source at = %e is already small\n", prop->n, m0, n0, eps0);
    else
      if(n==niters) printf("Warning... Step %d, action density (%d, %d) iteration not converging after %d steps! RES=%e > %.2e\n",
			   prop->n, m0, n0, n, eps, SMALL);
      else printf("Step %d, BiCGSolve action density (%d, %d) converged after %d iterations, rsdl = %e < %e\n",
		  prop->n, m0, n0, n, eps, epsW);
    
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
    for(nf=0;nf<NFACES;nf++) {
      ne = grid->face[i*NFACES+nf];
      normal = grid->normal[i*NFACES+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1)
	y[i] += 0.25*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]+fabs(Wavecg[m][n][ne]))*x[nc2]);
      if(nc2==-1)
	y[i] += 0.25*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]-fabs(Wavecg[m][n][ne]))*x[nc1]);
      else
	y[i] += 0.25*dt*Wavewnstep*df*normal/Ac*((Wavecg[m][n][ne]+fabs(Wavecg[m][n][ne]))*x[nc2]+
				     (Wavecg[m][n][ne]-fabs(Wavecg[m][n][ne]))*x[nc1]);         
	       
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
      WavedivScx[i][k] = 0;
      WavedivScy[i][k] = 0;
      volume = 0;
      sumUw = 0;
      for(nf=0; nf<NFACES; nf++){
	ne = grid->face[i*NFACES+nf];
	normal = grid->normal[i*NFACES+nf];
	df = grid->df[ne];
	dg = grid->dg[ne];

	nc1 = grid->grad[2*ne];
	nc2 = grid->grad[2*ne+1];
	if (nc1 == -1){
	  //0.5 is due to phase average, i.e. phase average of cos(phi)^2 

	  WavedivScx[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	    *Waveux[nc2][k]*grid->dzz[nc2][k];
	  WavedivScy[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	    *Waveuy[nc2][k]*grid->dzz[nc2][k];
	  volume += (WaveUw[ne][k]+fabs(WaveUw[ne][k]))*grid->dzz[nc2][k];
	  sumUw += (WaveUw[ne][k]+fabs(WaveUw[ne][k]));

	}else if (nc2 == -1){

	  WavedivScx[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]-fabs(WaveUw[ne][k]))
	    *Waveux[nc1][k]*grid->dzz[nc1][k];
	  WavedivScy[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]-fabs(WaveUw[ne][k]))
	    *Waveuy[nc1][k]*grid->dzz[nc1][k];
	  volume += (WaveUw[ne][k]-fabs(WaveUw[ne][k]))*grid->dzz[nc1][k];
	  sumUw += (WaveUw[ne][k]-fabs(WaveUw[ne][k]));

	}else{
	  Nk_shallow = grid->Nk[nc1];
	  nc_deep = nc2;
	  if (grid->Nk[nc2] < Nk_shallow){
	    Nk_shallow = grid->Nk[nc2];
	    nc_deep = nc1;
	  }
	  //if (k < Nk_shallow){

	  WavedivScx[i][k] -= (0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	      *Waveux[nc2][k]*grid->dzz[nc2][k]
	      +0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]-fabs(WaveUw[ne][k]))
				 *Waveux[nc1][k]*grid->dzz[nc1][k]);

	  WavedivScy[i][k] -= (0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	      *Waveuy[nc2][k]*grid->dzz[nc2][k]
	      +0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]-fabs(WaveUw[ne][k]))
				 *Waveuy[nc1][k]*grid->dzz[nc1][k]);

	    volume += (WaveUw[ne][k]+fabs(WaveUw[ne][k]))*grid->dzz[nc2][k]
	      +(WaveUw[ne][k]-fabs(WaveUw[ne][k]))*grid->dzz[nc1][k];

	    sumUw += (WaveUw[ne][k]+fabs(WaveUw[ne][k]))+(WaveUw[ne][k]-fabs(WaveUw[ne][k]));
	    /*}else if (nc_deep == nc2){
	    WavedivScx[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	      *Waveux[nc2][k]*grid->dzz[nc2][k];
	    WavedivScy[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	      *Waveuy[nc2][k]*grid->dzz[nc2][k];
	    volume += (WaveUw[ne][k]+fabs(WaveUw[ne][k]))*grid->dzz[nc2][k];
	    sumUw += (WaveUw[ne][k]+fabs(WaveUw[ne][k]));
	  }else{
	    WavedivScx[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	      *Waveux[nc1][k]*grid->dzz[nc1][k];
	    WavedivScy[i][k] -= 0.5*0.5*dt*Wavewnstep*df*normal/Ac*(WaveUw[ne][k]+fabs(WaveUw[ne][k]))
	      *Waveuy[nc1][k]*grid->dzz[nc1][k];
	    volume += (WaveUw[ne][k]+fabs(WaveUw[ne][k]))*grid->dzz[nc2][k];
	    sumUw += (WaveUw[ne][k]+fabs(WaveUw[ne][k]));
	    }*/

	}

      }
     
      if (volume > 0.00001){
	WavedivScx[i][k] = WavedivScx[i][k]*sumUw/volume;
	WavedivScy[i][k] = WavedivScy[i][k]*sumUw/volume;
      }else{
	WavedivScx[i][k] = 0.0;
	WavedivScy[i][k] = 0.0;
      }
      if (WavedivScx[i][k] != WavedivScx[i][k]){
	printf("Bad Scx proc=%d Scx[i=%d][k=%d]=%f\n", myproc, i, k, WavedivScx[i][k]);
      }
      if (WavedivScy[i][k] != WavedivScy[i][k]){
	printf("Bad Scx proc=%d Scx[i=%d][k=%d]=%f\n", myproc, i, k, WavedivScy[i][k]);
      }
    }
  }
  // So far, radiation stress due to vertically oscillating motion is not included yet.
  // It will be included when interpolating to the cell edge.
  ISendRecvCellData3D(WavedivScx, grid, myproc, comm);    
  ISendRecvCellData3D(WavedivScy, grid, myproc, comm);        
}

static REAL InterpWindShearToFace(int j, physT *phys, gridT *grid){
  int nc1, nc2;
  REAL def1, def2, Dj, Senc1, Senc2, U, Ustar, usex, usey, cd;

  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  Dj = grid->dg[j];
  cd = sqrt(0.00128);
  
  if (nc1 == -1){
    U = sqrt(pow(Wavewind_spfx[nc2]-phys->uc[nc2][grid->ctop[nc2]], 2.0)
	    +pow(Wavewind_spfy[nc2]-phys->vc[nc2][grid->ctop[nc2]], 2.0));
    if (U <= 7.5)
      cd = 0.00128;
    else
      cd = (0.8+0.065*U)*0.001;

    return 0.00125*cd*U*((Wavewind_spfx[nc2]-phys->uc[nc2][grid->ctop[nc2]])*grid->n1[j] 
			 +(Wavewind_spfy[nc2]-phys->vc[nc2][grid->ctop[nc2]])*grid->n2[j]);
  }else if (nc2 == -1){
    U = sqrt(pow(Wavewind_spfx[nc1]-phys->uc[nc1][grid->ctop[nc1]], 2.0)
	    +pow(Wavewind_spfy[nc1]-phys->uc[nc1][grid->ctop[nc1]], 2.0));
    if (U <= 7.5)
      cd = 0.00128;
    else
      cd = (0.8+0.065*U)*0.001;
    return 0.00125*cd*U*((Wavewind_spfx[nc1]-phys->uc[nc1][grid->ctop[nc1]])*grid->n1[j]
			+(Wavewind_spfy[nc1]-phys->vc[nc1][grid->ctop[nc1]])*grid->n2[j]);

  }else{
    
    def1 = grid->def[nc1*NFACES+grid->gradf[2*j]];
    def2 = grid->def[nc2*NFACES+grid->gradf[2*j+1]];
  
    usex = ((Wavewind_spfx[nc2]-phys->uc[nc2][grid->ctop[nc2]])*def1+
	    (Wavewind_spfx[nc1]-phys->uc[nc1][grid->ctop[nc1]])*def2)/(def1+def2);
    usey = ((Wavewind_spfy[nc2]-phys->vc[nc2][grid->ctop[nc2]])*def1+
	    (Wavewind_spfy[nc1]-phys->vc[nc1][grid->ctop[nc1]])*def2)/(def1+def2);
    
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
      phys->utmp[j][k]+=(1.0-fab)*WavedivSe[j][k];
      
    }
  }


  for(jptr = grid->edgedist[0]; jptr < grid->edgedist[1]; jptr++){
    j = grid->edgep[jptr];
    //  for(j = 0; j<Ne; j++){
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    dz = 0;

    for (k = 0; k < grid->etop[j]; k++) WavedivSe[j][k] = 0.0;
    flag = 0;
    
    for(k = grid->etop[j]; k < grid->Nke[j]; k++){
      if(nc1 == -1){
	WavedivSe[j][k] = WavedivScx[nc2][k]*grid->n1[j]+WavedivScy[nc2][k]*grid->n2[j];
	//Add radiation stress due to vertically oscillating motion
	//0.5 is due to phase average
	WavedivSe[j][k] += prop->dt * Wavewnstep*0.5*(0 - pow(Waveuz[nc2][k], 2))/grid->dg[j];
	if (k == grid->etop[j]){
	  dz = grid->dzz[nc2][k];
	  if (dz > 0.01){
	    WavedivSe[j][k] -= prop->dt * Wavewnstep * 1/16/dz*(0 - WaveEtot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}
	if(k == grid->etop[j] + 1){
	  if (flag == 0){
	    dz += grid->dzz[nc2][k];
	    WavedivSe[j][k] -= prop->dt * Wavewnstep * 1/16/dz*(0 - WaveEtot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}
      
      }else if(nc2 == -1){
	WavedivSe[j][k] = WavedivScx[nc1][k]*grid->n1[j]+WavedivScy[nc1][k]*grid->n2[j];
	WavedivSe[j][k] += prop->dt * Wavewnstep*0.5*(pow(Waveuz[nc1][k], 2) - 0)/grid->dg[j];
	if (k == grid->etop[j]){
	  dz = grid->dzz[nc1][k];
	  if (dz > 0.01){
	    WavedivSe[j][k] -= prop->dt * Wavewnstep * 1/16/dz*(WaveEtot[nc1]*GRAV - 0)/grid->dg[j];
	    flag = 1;
	  }
	}
	if(k == grid->etop[j] + 1){
	  if (flag == 0){
	    dz += grid->dzz[nc1][k];
	    WavedivSe[j][k] -= prop->dt * Wavewnstep * 1/16/dz*(WaveEtot[nc1]*GRAV - 0)/grid->dg[j];
	    flag = 1;
	  }
	}
      }else{
	weight1 = 1 - (grid->def[nc1*NFACES + grid->gradf[2*j]]/grid->dg[j]);
	//weight1 = 0.5;
	WavedivSe[j][k] = 
	  (1.0-weight1)*(WavedivScx[nc1][k] * grid->n1[j] + WavedivScy[nc1][k] * grid->n2[j]) + 
	  weight1*(WavedivScx[nc2][k] * grid->n1[j] + WavedivScy[nc2][k] * grid->n2[j]); 

	WavedivSe[j][k] += prop->dt * Wavewnstep * 0.5 * (pow(Waveuz[nc1][k], 2) - pow(Waveuz[nc2][k], 2))/grid->dg[j];

	if (k == grid->etop[j]){
	  dz = 0.5*(grid->dzz[nc2][k] + grid->dzz[nc1][k]);
	  if (dz > 0.01){
	    WavedivSe[j][k] -= prop->dt * Wavewnstep * 1/16/dz*(WaveEtot[nc1]*GRAV - WaveEtot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}  

	if (k == grid->etop[j] + 1){
	  if (flag == 0){
	    dz = 0.5*(grid->dzz[nc2][k] + grid->dzz[nc1][k]);
	    WavedivSe[j][k] -= prop->dt * Wavewnstep * 1/16/dz*(WaveEtot[nc1]*GRAV - WaveEtot[nc2]*GRAV)/grid->dg[j];
	    flag = 1;
	  }
	}      
      }

      if(WavedivSe[j][k] != WavedivSe[j][k] )
	  printf("----------------proc=%d, divSe[j=%d][k=%d] = %f\n", myproc, j, k, WavedivSe[j][k]);
    }
    for (k = grid->Nke[j]; k < grid->Nkc[j]; k++){
      WavedivSe[j][k] = WavedivSe[j][grid->Nke[j]-1];
    }

    j = grid->edgep[jptr];
    for(k=grid->etop[j];k<grid->Nkc[j];k++){
      phys->utmp[j][k]+=fab*WavedivSe[j][k];
      
    }

  }
  
  ISendRecvEdgeData3D(WavedivSe, grid, myproc, comm);        
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
      phys->utmp[j][k]+=WavedivSe[j][k];
      
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
    kD = Wavekmean[i]*depth;
    kw  = Wavekmean[i];

    if (kD > 0.1){
      invsinh = 1/sinh(kD);
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	z -= 0.5*grid->dzz[i][k];
	kz = Wavekmean[i]*(z+depth);
	if (kD < 100)
	  tmp = WaveHs[i]/8.0*Wavesgmean[i]*cosh(kz)*invsinh;
	else
	  tmp = WaveHs[i]/8.0*Wavesgmean[i]*exp(kz-kD);
	
	Waveuz[i][k] = tmp*tanh(kz);
	Waveux[i][k] = fabs(tmp*cos(Wavethtamean[i]));
	Waveuy[i][k] = fabs(tmp*sin(Wavethtamean[i]));
	z -= 0.5*grid->dzz[i][k];
	if (Waveux[i][k] != Waveux[i][k])
	  printf("Bad ux: p = %d i=%d k=%d ux=%f\n", myproc, i, k, Waveux[i][k]);
	
	if (Waveuz[i][k] != Waveuz[i][k])
	  printf("Bad uy: p = %d i=%d k=%d uy=%f\n", myproc, i, k, Waveuy[i][k]);
	
	if (Waveuz[i][k] != Waveuz[i][k])
	  printf("Bad uz: p = %d i=%d k=%d uz=%f\n", myproc, i, k, Waveuz[i][k]);

      }
    }else{
      for (k=grid->ctop[i]; k < grid->Nk[i]; k++){
	Waveux[i][k] = 0;
	Waveuy[i][k] = 0;
	Waveuz[i][k] = 0;
      }
    }
      
  }
  ISendRecvCellData3D(Waveux, grid, myproc, comm);    
  ISendRecvCellData3D(Waveuy, grid, myproc, comm);    
  ISendRecvCellData3D(Waveuz, grid, myproc, comm);    
}

void ObtainEdgeUwField(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc)
{
  int j, k, Ne=grid->Ne, nc1, nc2, jptr;
  REAL weight1;

  for/*
 * Function: ReadWaveProperties
 * Usage: ReadWaveProperties(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * Based on wave.dat, load in the important parameters for 
 * wave model
 *
 */
void FetchReadWaveProperties(int myproc)
{
  ConstantWind = MPI_GetValue(DATAFILE,"constantwind","ReadWaveProperties",myproc);
}


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

  
  Hwsig = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Twsig = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Fw = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  MPI_GetFile(filename,DATAFILE,"InputCoarseDomain","AllocateWave",myproc); 
  Numdomain = MPI_GetSize(filename,"AllocateWave",myproc);
  Waveexcur = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Fetch = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Coarsedomain = (REAL *)SunMalloc(Numdomain*2*sizeof(REAL), "AllocateWave");

  Uwind = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Winddir = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");

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
    Hwsig[i]=0;
    Twsig[i]=0;
    Fw[i]=0;
    Fetch[i]=0;
    Coarsedomain[i*2]=0;
    Coarsedomain[i*2+1]=0;   
    Uwind[i]=0;
    Winddir[i]=0;
    Waveexcur[i]=0;
  } 

  MPI_GetFile(filename,DATAFILE,"InputCoarseDomain","CalculateFetch",myproc);     
  ifile=MPI_FOpen(filename,"r","CalculateFetch",myproc);
  for(i=0;i<Numdomain;i++){
    Coarsedomain[2*i]=getfield(ifile,filename);
    Coarsedomain[2*i+1]=getfield(ifile,filename);
      // make sure coarsedomain include all the voro points!!!!
  }
  fclose(ifile);

  uwind= MPI_GetValue(DATAFILE,"uwind","InitializeWave",myproc);
  winddir= MPI_GetValue(DATAFILE,"winddir","InitializeWave",myproc);
  
  // give initial value for Uwind and winddir
  if(ConstantWind==0 || ConstantWind==2){  
    for(i=0;i<grid->Nc;i++){
      Uwind[i]=ReturnWindSpeed(grid->xv[i],grid->yv[i]);
      Winddir[i]=ReturnWindDirection(grid->xv[i],grid->yv[i]);
    }
  } else if(ConstantWind==1){
    for(i=0;i<grid->Nc;i++){
      Uwind[i]=uwind;
      Winddir[i]=winddir;
    }
  } else {
    MPI_GetFile(str,DATAFILE,"InitWindFile","InitializeWind",myproc);     
    ifile=MPI_FOpen(str,"r","InitializeWind",myproc);
    for(i=0;i<grid->Nc;i++){ 
      getfield(ifile,str);
      getfield(ifile,str);
      Uwind[i]=(REAL)getfield(ifile,str);
      Winddir[i]=(REAL)getfield(ifile,str);
    }
    fclose(ifile);
  }
  
  CalculateFetch(grid,myproc);

  CalculateHwsigTwsig(grid,phys,myproc);

  CalculateWaveexcur(grid,phys,myproc);

  CalculateFw(grid,myproc);
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
  free(Hwsig);
  free(Twsig);
  free(Uwind);
  free(Winddir);
  free(Fw);
  free(Fetch);
  free(Coarsedomain);
  free(Waveexcur);
}

/*
 * Function: CalculateFetch
 * usage: calculate fetch based on a coarse domain
 * --------------------------------------------------------
 * calculate intersection point and distance
 *
 */
void FetchCalculateFetch(gridT *grid, int myproc)
{
  int i,j;
  REAL tmp,xinter,yinter,distmp,dismin;
  for(i=0;i<grid->Nc;i++){
    Fetch[i]=10000000; // assign a big number
    dismin=10000000; // remember min distance
    for(j=0;j<Numdomain-1;j++){
      tmp=(Coarsedomain[2*j]-Coarsedomain[2*j+2])*sin(Winddir[i])-(Coarsedomain[2*j+1]-Coarsedomain[2*j+3])*cos(Winddir[i]);
      // if tmp!=0 two lines are not parallel and have an intersection point
      if(tmp!=0){
        xinter=Coarsedomain[2*j]*Coarsedomain[2*j+3]-Coarsedomain[2*j+1]*Coarsedomain[2*j+2];
        xinter=xinter*cos(Winddir[i]);
        xinter=xinter+(-grid->yv[i]*cos(Winddir[i])+sin(Winddir[i])*grid->xv[i])*(Coarsedomain[2*j]-Coarsedomain[2*j+2]);
        xinter=xinter/tmp;
        // test whether is in the middle part
        if(xinter>=Min(Coarsedomain[2*j],Coarsedomain[2*j+2]) && xinter<=Max(Coarsedomain[2*j],Coarsedomain[2*j+2])){
          if((Coarsedomain[2*j]-Coarsedomain[2*j+2])!=0){
            yinter=((Coarsedomain[2*j+1]-Coarsedomain[2*j+3])*xinter+(Coarsedomain[2*j]*Coarsedomain[2*j+3]-Coarsedomain[2*j+1]*Coarsedomain[2*j+2]))/(Coarsedomain[2*j]-Coarsedomain[2*j+2]);
          } else {
            yinter=sin(Winddir[i])/cos(Winddir[i])*(xinter-grid->xv[i])+grid->yv[i];
          }
          distmp=sqrt(pow((grid->xv[i]-xinter),2)+pow((grid->yv[i]-yinter),2));
          if(distmp<dismin)
            dismin=distmp;
          if(distmp<Fetch[i] && (xinter-grid->xv[i])*cos(Winddir[i])>=0 && (yinter-grid->yv[i])*sin(Winddir[i])>=0)
            Fetch[i]=distmp;
        }
      }    
    }
    // output error
    if(Fetch[i]==10000000){
      Fetch[i]=dismin;
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
     ua=0.71*pow(Uwind[i],1.23);
     tmp1=grav*Fetch[i]/ua/ua;
     tmp2=grav*(phys->h[i]+grid->dv[i])/ua/ua;
     // for hwsig
     tmp=0;
     tmp=0.00565*pow(tmp1,0.5);
     tmp=tmp/tanh(0.530*pow(tmp2,0.75));
     tmp=tanh(tmp)*tanh(0.530*pow(tmp2,0.75))*0.283;
     Hwsig[i]=tmp*ua*ua/grav;
     // for twsig
     tmp=0;
     tmp=0.0379*pow(tmp1,0.33333333333);
     tmp=tmp/tanh(0.833*pow(tmp2,0.375));
     tmp=tanh(tmp)*tanh(0.833*pow(tmp2,0.75))*7.54;
     Twsig[i]=tmp*ua/grav; 
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
    omega=2*PI/Twsig[i];
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
    Waveexcur[i]=Hwsig[i]/2/sinh(k*(phys->h[i]+grid->dv[i]));
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
    rew=Waveexcur[i]*Waveexcur[i]*2*PI/Twsig[i]/nu;
    if(rew!=0){
      if(rew<3e5){
        Fw[i]=2*pow(rew,-0.5);
      } else if(rew>=3e5 && rew<=6e5) {
        Fw[i]=0.024*pow(rew,-0.123);
      } else {
        Fw[i]=0.005;
      }
    } 
    if(rew==0) // you can set limit
      Fw[i]=0;
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
  
  if(ConstantWind==0){
    if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) { 
       Write2DData(Fetch,prop->mergeArrays,FetchFID,"Error outputting fetch data!\n", grid,numprocs,myproc,comm);
       Write2DData(Hwsig,prop->mergeArrays,HwsigFID,"Error outputting significant wave height data!\n", grid,numprocs,myproc,comm);
       Write2DData(Twsig,prop->mergeArrays,TwsigFID,"Error outputting significant wave period data!\n", grid,numprocs,myproc,comm);
    }
    if(prop->n==prop->nsteps+prop->nstart) {
      fclose(FetchFID);
      fclose(HwsigFID);
      fclose(TwsigFID);
    }
  } else {
    if(prop->n==1+prop->nstart){
      Write2DData(Fetch,prop->mergeArrays,FetchFID,"Error outputting fetch data!\n", grid,numprocs,myproc,comm);
      Write2DData(Hwsig,prop->mergeArrays,HwsigFID,"Error outputting significant wave height data!\n", grid,numprocs,myproc,comm);
      Write2DData(Twsig,prop->mergeArrays,TwsigFID,"Error outputting significant wave period data!\n", grid,numprocs,myproc,comm);
      fclose(FetchFID);
      fclose(HwsigFID);
      fclose(TwsigFID);
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
  HwsigFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"TwsigFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  TwsigFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);
  

  MPI_GetFile(filename,DATAFILE,"FetchFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  FetchFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);
}(jptr = grid->edgedist[0]; jptr < grid->edgedist[1]; jptr++){
    j = grid->edgep[jptr];
    //  for(j = 0; j<Ne; j++){
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j + 1];
    
    weight1 = 0.5;
    for (k = 0; k < grid->etop[j]; k++) WaveUw[j][k] = 0.0;
    for (k = grid->etop[j]; k < grid->Nke[j]; k++){
      if (nc1 == -1)
    	WaveUw[j][k] = Waveux[nc2][k]*grid->n1[j] + Waveuy[nc2][k]*grid->n2[j];
      else if (nc2 == -1)
    	WaveUw[j][k] = Waveux[nc1][k]*grid->n1[j] + Waveuy[nc1][k]*grid->n2[j];
      else
	//weight1 = 1 - (grid->def[nc1*NFACES + grid->gradf[2*j]]/grid->dg[j]);
    	WaveUw[j][k] = weight1 * (Waveux[nc1][k]*grid->n1[j] + Waveuy[nc1][k]*grid->n2[j])
    	  + (1 - weight1) * (Waveux[nc2][k]*grid->n1[j] + Waveuy[nc2][k]*grid->n2[j]);
    }
    for (k = grid->Nke[j]; k < grid->Nkc[j]; k++)
      WaveUw[j][k] = WaveUw[j][grid->Nke[j] - 1];

    
  }
  ISendRecvEdgeData3D(WaveUw, grid, myproc, comm);    
}


static void ReadWaveVariables(gridT *grid, propT *prop, int myproc, MPI_Comm comm) {

  int i, j;
  int m, n, s, l;

  //fread(&(prop->nstart), sizeof(int),1,StartWaveFID);

  for(m=0; m<WaveMw; m++)
    for(n=0; n<WaveNw; n++){
      fread(WaveN[m][n],sizeof(REAL),grid->Nc,StartWaveFID);
      //      ISendRecvCellData2D(WaveN[m][n], grid, myproc, comm);    
    }
  


  fclose(StartWaveFID);

}

static void OpenWaveFiles(int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];


  MPI_GetFile(filename,DATAFILE,"WaveHeightFile","OpenWaveFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  WaveHeightFID = MPI_FOpen(str,"w","OpenFile",myproc);

  MPI_GetFile(filename,DATAFILE,"WindSpeedFile","OpenWaveFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  WindSpeedFID = MPI_FOpen(str,"w","OpenFile",myproc);

  MPI_GetFile(filename,DATAFILE,"WindDirectionFile","OpenWaveFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  WindDirectionFID = MPI_FOpen(str,"w","OpenFile",myproc);

  MPI_GetFile(filename,DATAFILE,"WaveVelocityFile","OpenWaveFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  WaveVelocityFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"StoreWaveFile","OpenWaveFiles",myproc);
  sprintf(str,"%s.%d",filename,myproc);
  StoreWaveFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  if(RESTART) {
    MPI_GetFile(filename,DATAFILE,"StartWaveFile","OpenWaveFiles",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    StartWaveFID = MPI_FOpen(str,"r","OpenWaveFiles",myproc);
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

    if(myproc==0 && VERBOSE>=1) 
      if(!blowup) 
	printf("Outputting wave data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      else
	printf("Outputting blowup wave data at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);

    //    if (prop->wave){
    if(ASCII) 
      for(i=0;i<grid->Nc;i++)
	fprintf(WaveHeightFID,"%f\n",WaveHs[i]);
    else {
      nwritten=fwrite(WaveHs, sizeof(REAL),grid->Nc,WaveHeightFID);
      if(nwritten!=grid->Nc) {
	printf("Error outputting wave height data!\n");
	exit(EXIT_WRITING);
      }
    }
    fflush(WaveHeightFID);
   
     
    if(ASCII) 
      for(i=0;i<grid->Nc;i++)
	fprintf(WindSpeedFID,"%f\n",Wavewind_spf[i]);
    
    else{
      nwritten=fwrite(Wavewind_spf, sizeof(REAL),grid->Nc,WindSpeedFID);
      if(nwritten!=grid->Nc) {
	printf("Error outputting wind speed data!\n");
	exit(EXIT_WRITING);
      }
    }
    fflush(WindSpeedFID);
   
    if(ASCII) 
      for(i=0;i<grid->Nc;i++)
	fprintf(WindDirectionFID,"%f\n",Wavewind_dgf[i]);
      
    else {
      nwritten=fwrite(Wavewind_dgf, sizeof(REAL),grid->Nc,WindDirectionFID);
      if(nwritten!=grid->Nc) {
	printf("Error outputting wind direction data!\n");
	exit(EXIT_WRITING);
      }
    }
    fflush(WindDirectionFID);
    
    if(ASCII) 
      for(i=0;i<grid->Nc;i++)
	fprintf(WaveVelocityFID,"%f\n",Waveub[i]);
    else {
	nwritten=fwrite(Waveub, sizeof(REAL),grid->Nc,WaveVelocityFID);
	if(nwritten!=grid->Nc) {
	  printf("Error outputting wave velocity data!\n");
	  exit(EXIT_WRITING);
	}
	fflush(WaveVelocityFID);
    }
    

  }

  //if(prop->n==prop->nsteps+prop->nstart || blowup) {
  if(!(prop->n%(prop->ntout)) || prop->n==prop->nsteps+prop->nstart || blowup){  
 
    if(VERBOSE>=1 && myproc==0) printf("Writing to wave rstore...\n");
    fseek( StoreWaveFID, 0, SEEK_SET );
    //nwritten=fwrite(&(prop->n),sizeof(int),1,StoreWaveFID);
    

    if (prop->wave)
      for(m=0; m<WaveMw; m++)
	for(n=0; n<WaveNw; n++)
	  fwrite(WaveN[m][n],sizeof(REAL),grid->Nc,StoreWaveFID);


    fflush(StoreWaveFID);

  }
  
  if(prop->n==prop->nsteps+prop->nstart) {
    fclose(WaveHeightFID);
    fclose(WaveVelocityFID);
    fclose(WindSpeedFID);
    fclose(WindDirectionFID);
    fclose(StoreWaveFID);
  }

}

/*---------------------------------------------------------------------------*/
/* fecth model equation, not finish to merge into main code */

/*
 * Function: ReadWaveProperties
 * Usage: ReadWaveProperties(grid,phys,prop,myproc);
 * ----------------------------------------------------
 * Based on wave.dat, load in the important parameters for 
 * wave model
 *
 */
void FetchReadWaveProperties(int myproc)
{
  ConstantWind = MPI_GetValue(DATAFILE,"constantwind","ReadWaveProperties",myproc);
}


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

  
  Hwsig = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Twsig = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Fw = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  MPI_GetFile(filename,DATAFILE,"InputCoarseDomain","AllocateWave",myproc); 
  Numdomain = MPI_GetSize(filename,"AllocateWave",myproc);
  Waveexcur = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Fetch = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Coarsedomain = (REAL *)SunMalloc(Numdomain*2*sizeof(REAL), "AllocateWave");

  Uwind = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");
  Winddir = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "AllocateWave");

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
    Hwsig[i]=0;
    Twsig[i]=0;
    Fw[i]=0;
    Fetch[i]=0;
    Coarsedomain[i*2]=0;
    Coarsedomain[i*2+1]=0;   
    Uwind[i]=0;
    Winddir[i]=0;
    Waveexcur[i]=0;
  } 

  MPI_GetFile(filename,DATAFILE,"InputCoarseDomain","CalculateFetch",myproc);     
  ifile=MPI_FOpen(filename,"r","CalculateFetch",myproc);
  for(i=0;i<Numdomain;i++){
    Coarsedomain[2*i]=getfield(ifile,filename);
    Coarsedomain[2*i+1]=getfield(ifile,filename);
      // make sure coarsedomain include all the voro points!!!!
  }
  fclose(ifile);

  uwind= MPI_GetValue(DATAFILE,"uwind","InitializeWave",myproc);
  winddir= MPI_GetValue(DATAFILE,"winddir","InitializeWave",myproc);
  
  // give initial value for Uwind and winddir
  if(ConstantWind==0 || ConstantWind==2){  
    for(i=0;i<grid->Nc;i++){
      Uwind[i]=ReturnWindSpeed(grid->xv[i],grid->yv[i]);
      Winddir[i]=ReturnWindDirection(grid->xv[i],grid->yv[i]);
    }
  } else if(ConstantWind==1){
    for(i=0;i<grid->Nc;i++){
      Uwind[i]=uwind;
      Winddir[i]=winddir;
    }
  } else {
    MPI_GetFile(str,DATAFILE,"InitWindFile","InitializeWind",myproc);     
    ifile=MPI_FOpen(str,"r","InitializeWind",myproc);
    for(i=0;i<grid->Nc;i++){ 
      getfield(ifile,str);
      getfield(ifile,str);
      Uwind[i]=(REAL)getfield(ifile,str);
      Winddir[i]=(REAL)getfield(ifile,str);
    }
    fclose(ifile);
  }
  
  CalculateFetch(grid,myproc);

  CalculateHwsigTwsig(grid,phys,myproc);

  CalculateWaveexcur(grid,phys,myproc);

  CalculateFw(grid,myproc);
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
  free(Hwsig);
  free(Twsig);
  free(Uwind);
  free(Winddir);
  free(Fw);
  free(Fetch);
  free(Coarsedomain);
  free(Waveexcur);
}

/*
 * Function: CalculateFetch
 * usage: calculate fetch based on a coarse domain
 * --------------------------------------------------------
 * calculate intersection point and distance
 *
 */
void FetchCalculateFetch(gridT *grid, int myproc)
{
  int i,j;
  REAL tmp,xinter,yinter,distmp,dismin;
  for(i=0;i<grid->Nc;i++){
    Fetch[i]=10000000; // assign a big number
    dismin=10000000; // remember min distance
    for(j=0;j<Numdomain-1;j++){
      tmp=(Coarsedomain[2*j]-Coarsedomain[2*j+2])*sin(Winddir[i])-(Coarsedomain[2*j+1]-Coarsedomain[2*j+3])*cos(Winddir[i]);
      // if tmp!=0 two lines are not parallel and have an intersection point
      if(tmp!=0){
        xinter=Coarsedomain[2*j]*Coarsedomain[2*j+3]-Coarsedomain[2*j+1]*Coarsedomain[2*j+2];
        xinter=xinter*cos(Winddir[i]);
        xinter=xinter+(-grid->yv[i]*cos(Winddir[i])+sin(Winddir[i])*grid->xv[i])*(Coarsedomain[2*j]-Coarsedomain[2*j+2]);
        xinter=xinter/tmp;
        // test whether is in the middle part
        if(xinter>=Min(Coarsedomain[2*j],Coarsedomain[2*j+2]) && xinter<=Max(Coarsedomain[2*j],Coarsedomain[2*j+2])){
          if((Coarsedomain[2*j]-Coarsedomain[2*j+2])!=0){
            yinter=((Coarsedomain[2*j+1]-Coarsedomain[2*j+3])*xinter+(Coarsedomain[2*j]*Coarsedomain[2*j+3]-Coarsedomain[2*j+1]*Coarsedomain[2*j+2]))/(Coarsedomain[2*j]-Coarsedomain[2*j+2]);
          } else {
            yinter=sin(Winddir[i])/cos(Winddir[i])*(xinter-grid->xv[i])+grid->yv[i];
          }
          distmp=sqrt(pow((grid->xv[i]-xinter),2)+pow((grid->yv[i]-yinter),2));
          if(distmp<dismin)
            dismin=distmp;
          if(distmp<Fetch[i] && (xinter-grid->xv[i])*cos(Winddir[i])>=0 && (yinter-grid->yv[i])*sin(Winddir[i])>=0)
            Fetch[i]=distmp;
        }
      }    
    }
    // output error
    if(Fetch[i]==10000000){
      Fetch[i]=dismin;
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
     ua=0.71*pow(Uwind[i],1.23);
     tmp1=grav*Fetch[i]/ua/ua;
     tmp2=grav*(phys->h[i]+grid->dv[i])/ua/ua;
     // for hwsig
     tmp=0;
     tmp=0.00565*pow(tmp1,0.5);
     tmp=tmp/tanh(0.530*pow(tmp2,0.75));
     tmp=tanh(tmp)*tanh(0.530*pow(tmp2,0.75))*0.283;
     Hwsig[i]=tmp*ua*ua/grav;
     // for twsig
     tmp=0;
     tmp=0.0379*pow(tmp1,0.33333333333);
     tmp=tmp/tanh(0.833*pow(tmp2,0.375));
     tmp=tanh(tmp)*tanh(0.833*pow(tmp2,0.75))*7.54;
     Twsig[i]=tmp*ua/grav; 
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
    omega=2*PI/Twsig[i];
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
    Waveexcur[i]=Hwsig[i]/2/sinh(k*(phys->h[i]+grid->dv[i]));
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
    rew=Waveexcur[i]*Waveexcur[i]*2*PI/Twsig[i]/nu;
    if(rew!=0){
      if(rew<3e5){
        Fw[i]=2*pow(rew,-0.5);
      } else if(rew>=3e5 && rew<=6e5) {
        Fw[i]=0.024*pow(rew,-0.123);
      } else {
        Fw[i]=0.005;
      }
    } 
    if(rew==0) // you can set limit
      Fw[i]=0;
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
  
  if(ConstantWind==0){
    if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) { 
       Write2DData(Fetch,prop->mergeArrays,wprop->FetchFID,"Error outputting fetch data!\n", grid,numprocs,myproc,comm);
       Write2DData(Hwsig,prop->mergeArrays,wprop->HwsigFID,"Error outputting significant wave height data!\n", grid,numprocs,myproc,comm);
       Write2DData(Twsig,prop->mergeArrays,wprop->TwsigFID,"Error outputting significant wave period data!\n", grid,numprocs,myproc,comm);
    }
    if(prop->n==prop->nsteps+prop->nstart) {
      fclose(FetchFID);
      fclose(HwsigFID);
      fclose(TwsigFID);
    }
  } else {
    if(prop->n==1+prop->nstart){
      Write2DData(Fetch,prop->mergeArrays,FetchFID,"Error outputting fetch data!\n", grid,numprocs,myproc,comm);
      Write2DData(Hwsig,prop->mergeArrays,HwsigFID,"Error outputting significant wave height data!\n", grid,numprocs,myproc,comm);
      Write2DData(Twsig,prop->mergeArrays,TwsigFID,"Error outputting significant wave period data!\n", grid,numprocs,myproc,comm);
      fclose(FetchFID);
      fclose(HwsigFID);
      fclose(TwsigFID);
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
  HwsigFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"TwsigFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  TwsigFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);
  

  MPI_GetFile(filename,DATAFILE,"FetchFile","OpenWaveFiles",myproc);
  if(merge)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  FetchFID = MPI_FOpen(str,"w","OpenWaveFiles",myproc);
}
