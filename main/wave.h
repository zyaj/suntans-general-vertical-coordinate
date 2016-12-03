/*
 * File: wave.h
 * Author: Yi-Ju Chou & Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for wave.c
 * including all the variables for wave model
 *
 */
#ifndef _wave_h
#define _wave_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"
 
#define RHOair 1.225

typedef struct _waveT {
  REAL *sg;
  REAL *dsg;
  REAL *thtaw;
  REAL **kw;
  REAL *ktail;
  REAL ***cgx;
  REAL ***cgy;
  REAL ***cg;
  REAL **cph;
  REAL ***cs;
  REAL ***ct;
  REAL ***N;
  REAL ***Nold;
  REAL ***Ntmp;
  REAL ***src;
  REAL *Etot;
  REAL *Etmp;
  REAL *kmean;
  REAL *sgmean;
  REAL *T0;
  REAL *Hs;
  REAL *T01;
  REAL *sg_PM;
  REAL *thtamean;
  REAL *Cr;
  REAL *fw;
  REAL *uscx;
  REAL *uscy;
  REAL *use;
  REAL *Uwind;
  REAL *Etail;
  REAL *Ntail;

  REAL *a;
  REAL *b;
  REAL *c;

  REAL *sp;
  REAL *sm;
  REAL ***ssrc;

  REAL *tp;
  REAL *tm;
  REAL ***tsrc;

  REAL **wind_sp;
  REAL **wind_dg;
  REAL **klambda;

  REAL *wind_spfx;
  REAL *wind_spfy;
  REAL *wind_spf;
  REAL *wind_dgf;

  REAL *Hw;
  REAL *ub;
  REAL *ab;
  REAL *tmp;
  REAL **fz;
  REAL *fphi;

  REAL **ux;
  REAL **uy;
  REAL **uz;  
  REAL **Uw;
  REAL **divScx;
  REAL **divScy;
  REAL **divSe;
  
  REAL *ab_edge;
  REAL *sgmean_edge;
  REAL *thtamean_edge;
  REAL *kmean_edge;

  REAL **kw_edge;

  // fetch model part
  REAL *Hwsig; // the significant wave height for each cell, Hwsig[cell] = Hw in Yi-Ju
  REAL *Twsig; // the significant wave period for each cell, Twsig[cell] not exist in Yi-Ju
  REAL *Fw; // the friction factor for wave bottom shear stress, Fw[cell]
  REAL *Fetch; // the Fetch for each cell, Fetch[cell]
  REAL *Winddir; // wind direction for each cell, Winddir[cell]
  REAL *Coarsedomain; // input coarse domain Coarsedomain[2*cell], x y
  REAL *Waveexcur;
  //REAL *Uwind; // the wind speed, usually U10[cell]

} waveT;


/*
 * Main wave property struct.
 *
 */
typedef struct _wpropT {

  int Mw, Nw;
  int wnstep;
  REAL sgmin;
  REAL sgmax;
  REAL sg99;
  REAL sgtail;
  REAL wind_dt;
  int implicit_whitecap;
  int implicit_advection;
  int wind_forcing;
  int Nwind, nwind;
  int nstation;
  int wind_shear;
  int rad_stress;
  int form_drag;
  int btm_mud;
  int btm_sedi_erosion;
  REAL depth_fw_cutoff;
  REAL fw_drag;
  REAL *xw, *yw;
  REAL **lambda;
  REAL btm_conc;
  REAL btm_vis;
  REAL btm_mud_thickness;
  int depth_brk_cutoff;
  REAL depth_brk_indx;

  int tail_opt;
  REAL tail_pow;

  int NLtriad;
  int NLquad;
  int BRKdepth;
  
  //fetch model
  int FetchModel;
  int ConstantWind;
  int Numdomain;

  FILE *WaveHeightFID,  *WaveVelocityFID,  *WindSpeedFID,  *WindDirectionFID,  *StartWaveFID,  *StoreWaveFID;
  // fetch model
  FILE *HwsigFID,*TwsigFID,*FetchFID;

} wpropT;

waveT *wave;
wpropT *wprop;


void UpdateWave(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int blowup, int myproc, int numprocs);
void UpdateWaveCr(gridT *grid, physT *phys, int myproc);
void UpdateWaveFw(gridT *grid, int myproc); 
#endif
