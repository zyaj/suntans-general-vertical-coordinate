/*
 * Boundaries test file.
 *
 */
#include "boundaries.h"
#include "sediments.h"
#include "wave.h"
static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag);

/*
 * Function: OpenBoundaryFluxes
 * Usage: OpenBoundaryFluxes(q,ubnew,ubn,grid,phys,prop);
 * ----------------------------------------------------
 * This will update the boundary flux at the edgedist[2] to edgedist[3] edges.
 * 
 * Note that phys->uold,vold contain the velocity at time step n-1 and 
 * phys->uc,vc contain it at time step n.
 *
 * The radiative open boundary condition does not work yet!!!  For this reason c[k] is
 * set to 0
 *
 */
void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop) {
  int j, jptr, ib, k, forced;
  REAL *uboundary = phys->a, **u = phys->uc, **v = phys->vc, **uold = phys->uold, **vold = phys->vold;
  REAL z, c0, c1, C0, C1, dt=prop->dt, u0, u0new, uc0, vc0, uc0old, vc0old, ub0;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->n1[j]+phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->n2[j];
  }
}

/*
 * Function: BoundaryScalars
 * Usage: BoundaryScalars(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the scalars at the open boundaries.
 * 
 */
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, iptr, i, ib, k;
  REAL z;

  // At the upstream boundary
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        phys->boundary_T[jptr-grid->edgedist[2]][k]=0;
        phys->boundary_s[jptr-grid->edgedist[2]][k]=0;
      }
  }

  // At the ocean boundary
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=0;k<grid->ctop[i];k++) {
      phys->s[i][k]=0;
      phys->T[i][k]=0;
    } 

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->s[i][k]=29;
      phys->T[i][k]=0;
    }
  }
}

/*
 * Function: BoundaryVelocities
 * Usage: BoundaryVelocities(grid,phys,prop,myproc);
 * -------------------------------------------------
 * This will set the values of u,v,w, and h at the boundaries.
 * 
 */
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, ib, k, iptr, i, winter;
  REAL fac, z, cb, flowflux, flowflux1, flowflux2, area, u0=0.075,time,timeshift,timebc;
  // inflow flux boundary
  fac=1;
  flowflux1=0.0;//0.001;
  flowflux2 = 0.0;
  timebc=prop->n*prop->dt/3600.0/24.0-0.2064814815+23;
  winter=0; 

  // setup inflow boundary condition
  if(winter){
    if(timebc>=28 && timebc<=33)
    {
      flowflux1=fac*0.075*sin((timebc-28)/5*PI);
      flowflux2=0.075*sin((timebc-28)/5*PI);
    }

    // estimated time phase difference between the inflow condition 
    // from the boundary to the location
    timeshift=0.04166666667;
    if(timebc>=38-timeshift && timebc<41-timeshift)
      flowflux1=0.125*fac*sin((timebc-38+timeshift)/3*PI);

    if(timebc>=38-timeshift && timebc<41-timeshift)
      flowflux2=0.25*fac*sin((timebc-38+timeshift)/3*PI);
  } 

  // setup boundary velocity
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) 
  {
    area=0;
    j = grid->edgep[jptr];
    ib = grid->grad[2*j];
    
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      area += grid->dzf[j][k]*grid->df[j];
      
    if((grid->Nke[j]-grid->etop[j])==1)
      if(grid->dzf[j][grid->etop[j]]<0.1)
        area = 0.1*grid->df[j]; 
    
    if(prop->n==0)
      area = 0.1*grid->df[j];
    // At the upstream boundary
    phys->boundary_h[jptr-grid->edgedist[2]]=phys->h[ib];
    
    // castro creek
    if(grid->yv[ib]<1000)
      flowflux=flowflux1;
    // wildcat creek
    else
      flowflux=flowflux2;

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->boundary_u[jptr-grid->edgedist[2]][k]=flowflux*grid->n1[j]/area;
      phys->boundary_v[jptr-grid->edgedist[2]][k]=flowflux*grid->n2[j]/area;
      phys->boundary_w[jptr-grid->edgedist[2]][k]=0; 
    }
  }

  // setup boundary free surface height
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    if(winter){
      timeshift=8640; // winter
      time=timeshift+prop->rtime;
      phys->h[i]=-5+0.3795*sin(0.0001413*time+0.5103)+0.4021*sin(0.000139*time+2.253)+0.3691*sin(0.00007292*time-2.109)+0.1363*sin(0.0000685*time-0.1151)+0.1112*sin(0.00006597*time+1.546)-0.01473*sin(0.000134*time+3.634)-0.00951*sin(0.0001474*time-3726)+0.1519*sin(0.0001468*time-3.681);
    }else{
      timeshift=1740; // summer
      time=timeshift+prop->rtime;
      phys->h[i]=-5+0.6277*sin(0.0001404*time-1.686)+0.3866*sin(7.276e-05*time+0.9212)+0.1549*sin(0.0001445*time+0.4462)+0.2457*sin(0.0001451*time-2.783)+0.1473*sin(0.0001376*time-0.5581)+1.262*sin(1.672e-08*time+3.059)+0.141*sin(6.62e-05*time+5.959)+0.2942*sin(6.705e-5*time+1.676);  
    }
  }
}

static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag) {
  int k;

  if(boundary_flag==open) {
    phys->boundary_h[boundary_index]=phys->h[ib];
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[boundary_index][k]=phys->uc[ib][k];
      phys->boundary_v[boundary_index][k]=phys->vc[ib][k];
      phys->boundary_w[boundary_index][k]=0.5*(phys->w[ib][k]+phys->w[ib][k+1]);
    }
  } else {
    phys->boundary_h[boundary_index]=prop->amp*fabs(cos(prop->omega*prop->rtime));
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[boundary_index][k]=phys->u[j][k]*grid->n1[j];
      phys->boundary_v[boundary_index][k]=phys->u[j][k]*grid->n2[j];
      phys->boundary_w[boundary_index][k]=0.5*(phys->w[ib][k]+phys->w[ib][k+1]);
    }
  }
}
	
/*
 * Function: WindStress
 * Usage: WindStress(grid,phys,prop,myproc);
 * -----------------------------------------
 * Set the wind stress.
 *
 */
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc) {
  int j, jptr;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    phys->tau_T[j]=grid->n2[j]*prop->tau_T;
    phys->tau_B[j]=0;
  }
}

//added part
/*
 * Function: BoundarySediment
 * Usage: BoundarySediment(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the suspended sediment concentration
 * at the open boundaries.
 * 
 */
void BoundarySediment(gridT *grid, physT *phys, propT *prop) {
  int jptr, j, ib, k,nosize,i,iptr;
  REAL z;

  // At the upstream boundary
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j=grid->edgep[jptr];
    ib=grid->grad[2*j];
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        sediments->boundary_sediC[nosize][jptr-grid->edgedist[2]][k]=0;
      }
    }
  }

  // At the ocean boundary
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      for(k=0;k<grid->ctop[i];k++) {
        sediments->SediC[nosize][i][k]=0;
      } 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        sediments->SediC[nosize][i][k]=0;
      }
    }
  }
}
/*
 * Function: WindSpeedandDirection
 * usage: calculate Wind field when Wind is not constant
 * -------------------------------------------------
 * calculate Uwind and Winddir
 *
 */
void FetchWindSpeedandDirection(gridT *grid, propT *prop, int myproc){
   int i;
   for(i=0;i<grid->Nc;i++){
     wave->Uwind[i]=0;
     wave->Winddir[i]=0;
   }
}

/*
 * Function: UserDefineFunction
 * usage: Define user function to output results we want
 * -------------------------------------------------
 *
 */

void UserDefinedFunction(gridT *grid, physT *phys, propT *prop,int myproc){}
void InitBoundaryData(propT *prop, gridT *grid, int myproc,MPI_Comm comm){}
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc,MPI_Comm comm){}
