/*
 * File: wave.h
 * Author : Yun Zhang & Yi-Ju Chou
 * Institution: Stanford University
 * --------------------------------
 * Header file for wave.c
 * including all the variables for wave model for both fetch model
 * and the wave model from Yi-Ju Chou
 *
 */
#ifndef _wave_h
#define _wave_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

// Yi-Ju's model
REAL *Wavesg;
REAL *Wavedsg;
REAL *Wavethtaw;
REAL **Wavekw;
REAL *Wavektail;
REAL ***Wavecgx;
REAL ***Wavecgy;
REAL ***Wavecg;
REAL **Wavecph;
REAL ***Wavecs;
REAL ***Wavect;
REAL ***WaveN;
REAL ***WaveNold;
REAL ***WaveNtmp;
REAL ***Wavesrc;
REAL *WaveEtot;
REAL *WaveEtmp;
REAL *Wavekmean;
REAL *Wavesgmean;
REAL *WaveT0;
REAL *WaveHs;
REAL *WaveT01;
REAL *Wavesg_PM;
REAL *Wavethtamean;
REAL *WaveCr;
REAL *Wavefw;
REAL *Waveuscx;
REAL *Waveuscy;
REAL *Waveuse;
REAL *WaveUwind;
REAL *WaveEtail;
REAL *WaveNtail;

REAL *Wavea;
REAL *Waveb;
REAL *Wavec;

REAL *Wavesp;
REAL *Wavesm;
REAL ***Wavessrc;

REAL *Wavetp;
REAL *Wavetm;
REAL ***Wavetsrc;

REAL **Wavewind_sp;
REAL **Wavewind_dg;
REAL **Waveklambda;

REAL *Wavewind_spfx;
REAL *Wavewind_spfy;
REAL *Wavewind_spf;
REAL *Wavewind_dgf;

REAL *WaveHw;
REAL *Waveub;
REAL *Waveab;
REAL *Wavetmp;
REAL **Wavefz;
REAL *Wavefphi;

REAL **Waveux;
REAL **Waveuy;
REAL **Waveuz;  
REAL **WaveUw;
REAL **WavedivScx;
REAL **WavedivScy;
REAL **WavedivSe;
  
REAL *Waveab_edge;
REAL *Wavesgmean_edge;
REAL *Wavethtamean_edge;
REAL *Wavekmean_edge;

REAL **Wavekw_edge;


// wave property variables
int WaveMw, WaveNw;
int Wavewnstep;
REAL Wavesgmin;
REAL Wavesgmax;
REAL Wavesg99;
REAL Wavesgtail;
REAL Wavewind_dt;
int Waveimplicit_whitecap;
int Waveimplicit_advection;
int Wavewind_forcing;
int WaveNwind, Wavenwind;
int Wavenstation;
int Wavewind_shear;
int Waverad_stress;
int Waveform_drag;
int Wavebtm_mud;
int Wavebtm_sedi_erosion;
REAL Wavedepth_fw_cutoff;
REAL Wavefw_drag;
REAL *Wavexw, *Waveyw;
REAL **Wavelambda;
REAL Wavebtm_conc;
REAL Wavebtm_vis;
REAL Wavebtm_mud_thickness;
int Wavedepth_brk_cutoff;
REAL Wavedepth_brk_indx;

int Wavetail_opt;
REAL Wavetail_pow;

int WaveNLtriad;
int WaveNLquad;
int WaveBRKdepth;

FILE *WaveHeightFID,  *WaveVelocityFID,  *WindSpeedFID,  *WindDirectionFID,  *StartWaveFID,  *StoreWaveFID;

// Fetch model variables
// FOR FETCH MODEL
REAL *Hwsig, // the significant wave height for each cell, Hwsig[cell] = Hw in Yi-Ju
     *Twsig, // the significant wave period for each cell, Twsig[cell] not exist in Yi-Ju
     *Fw, // the friction factor for wave bottom shear stress, Fw[cell]
     *Fetch, // the Fetch for each cell, Fetch[cell]
     *Winddir, // wind direction for each cell, Winddir[cell]
     *Coarsedomain, // input coarse domain Coarsedomain[2*cell], x y
     *Waveexcur,
     *Uwind; // the wind speed, usually U10[cell]
     
int FetchModel,
    ConstantWind,
    Numdomain;

FILE *HwsigFID,*TwsigFID,*FetchFID;

void UpdateWave(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int numprocs);

#endif
