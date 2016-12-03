/*
 * File: scalars.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * ----------------------------------------
 * This file contains the scalar transport function.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "scalars.h"
#include "util.h"
#include "tvd.h"
#include "subgrid.h"
#include "initialization.h"

#define SMALL_CONSISTENCY 1e-5

REAL smin_value, smax_value;

/*
 * Function: UpdateScalars
 * Usage: UpdateScalars(grid,phys,prop,wnew,scalar,Cn,kappa,kappaH,kappa_tv,theta);
 * ---------------------------------------------------------------------------
 * Update the scalar quantity stored in the array denoted by scal using the
 * theta method for vertical advection and vertical diffusion and Adams-Bashforth
 * for horizontal advection and diffusion.
 *
 * Cn must store the AB terms from time step n-1 for this scalar
 * kappa denotes the vertical scalar diffusivity
 * kappaH denotes the horizontal scalar diffusivity
 * kappa_tv denotes the vertical turbulent scalar diffusivity
 *
 */
void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **scal, REAL **boundary_scal, REAL **Cn, 
    REAL kappa, REAL kappaH, REAL **kappa_tv, REAL theta,
    REAL **src1, REAL **src2, REAL *Ftop, REAL *Fbot, int alpha_top, int alpha_bot,
    MPI_Comm comm, int myproc, int checkflag, int TVDscheme) 
{
  int i, iptr, j, jptr, ib, k, nf, ktop;
  int Nc=grid->Nc, normal, nc1, nc2, ne;
  REAL df, dg, Ac, dt=prop->dt, fab, *a, *b, *c, *d, *ap, *am, *bd, *uflux, dznew, mass, *sp, *temp, fluxtest,fluxtest1,*fluxm,*fluxp;
  REAL smin, smax, div_local, div_da, alpha,fac1,fac2,fac3,sum,sum1,sum2;
  int k1, k2, kmin, imin, kmax, imax, mincount, maxcount, allmincount, allmaxcount, flag;

  prop->TVD = TVDscheme;
  // These are used mostly debugging to turn on/off vertical and horizontal TVD.
  prop->horiTVD = 1;
  prop->vertTVD = 1;

  ap = phys->ap;
  am = phys->am;
  bd = phys->bp;
  temp = phys->bm;
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;

  fluxm=phys->fluxmNk;
  fluxp=phys->fluxpNk;

  // Never use AB2
  if(1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
        Cn[i][k]=0;
  } else
    fab=1.5;

  for(i=0;i<Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=scal[i][k];

  // prop->im
  if(prop->n==1 || prop->wetdry) {
    fac1=theta;
    fac2=1-theta;
    fac3=0;
  } else {
    if(prop->im==0)
    {
      fac1=theta;
      fac2=1-theta;
      fac3=0;
    } else if(prop->im==1){
      fac1=3.0/4.0;
      fac2=0;
      fac3=1.0/4.0;
    } else {
      fac1=5.0/4.0;
      fac2=-1;
      fac3=3.0/4.0;
    }
  }
  
  // Add on boundary fluxes, using stmp2 as the temporary storage
  // variable
  //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      phys->stmp2[i][k]=0;
  }

  if(boundary_scal) {
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];

      ib = grid->grad[2*j];

      // Set the value of stmp2 adjacent to the boundary to the value of the boundary.
      // This will be used to add the boundary flux when stmp2 is used again below.
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++)
        phys->stmp2[ib][k]=boundary_scal[jptr-grid->edgedist[2]][k];
    }
  }

  // Compute the scalar on the vertical faces (for horiz. advection)
  if(prop->TVD && prop->horiTVD)
    HorizontalFaceScalars(grid,phys,prop,scal,boundary_scal,prop->TVD,comm,myproc); 

  //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    Ac = grid->Ac[i];
    fluxtest=0;
    fluxtest1=0;

    if(grid->ctop[i]>=grid->ctopold[i]) {
      ktop=grid->ctop[i];
      dznew=grid->dzz[i][ktop];
    } else {
      ktop=grid->ctopold[i];
      dznew=0;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++) 
        dznew+=grid->dzz[i][k];      
    }
    // These are the advective components of the tridiagonal
    // at the new time step.
    if(!(prop->TVD && prop->vertTVD))
      for(k=0;k<grid->Nk[i]+1;k++) 
      {
        //add subgrid part
        if(prop->subgrid)
          alpha=subgrid->Acveff[i][k]/Ac;
        else
          alpha=1.0;

        if(!prop->im)
        {
          ap[k] = alpha*0.5*(wnew[i][k]+fabs(wnew[i][k]));
          am[k] = alpha*0.5*(wnew[i][k]-fabs(wnew[i][k]));
        } else {
          ap[k] = alpha*0.5*(phys->wtmp4[i][k]+fabs(phys->wtmp4[i][k]));
          am[k] = alpha*0.5*(phys->wtmp4[i][k]-fabs(phys->wtmp4[i][k]));
        }
      }
    else  // Compute the ap/am for TVD schemes
    {
      if(!prop->im)
        GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
          wnew,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
      else
        GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
          phys->wtmp4,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);    
   
      for(k=0;k<grid->Nk[i]+1;k++) 
      {
        //add subgrid part
        if(prop->subgrid)
          alpha=subgrid->Acveff[i][k]/Ac;
        else
          alpha=1.0;
        ap[k]= ap[k]*alpha;
        am[k]= am[k]*alpha;
      }
         
    }
    
    for(k=ktop+1;k<grid->Nk[i];k++) 
    {
      // add subgrid part
      if(prop->subgrid)
        alpha=subgrid->Acceff[i][k]/Ac;
      else
        alpha=1.0;
    
      a[k-ktop]=theta*dt*am[k];
      b[k-ktop]=alpha*grid->dzz[i][k]+theta*dt*(ap[k]-am[k+1]);
      c[k-ktop]=-theta*dt*ap[k+1];
    }

    // Top cell advection
    if(prop->subgrid)
      alpha=subgrid->Acceff[i][ktop]/Ac;
    else
      alpha=1.0;

    a[0]=0;
    b[0]=alpha*dznew-theta*dt*am[ktop+1];
    c[0]=-theta*dt*ap[ktop+1];
    
    // Bottom cell no-flux boundary condition for advection
    b[(grid->Nk[i]-1)-ktop]+=c[(grid->Nk[i]-1)-ktop];

    // Implicit vertical diffusion terms
    for(k=ktop+1;k<grid->Nk[i];k++)
    {
      if(prop->subgrid)
        alpha=subgrid->Acveff[i][k]/Ac;
      else
        alpha=1.0;

      bd[k]=alpha*(2.0*kappa+kappa_tv[i][k-1]+kappa_tv[i][k])/
        (grid->dzz[i][k-1]+grid->dzz[i][k]);

    }

    for(k=ktop+1;k<grid->Nk[i]-1;k++) 
    {
      a[k-ktop]-=theta*dt*bd[k];
      b[k-ktop]+=theta*dt*(bd[k]+bd[k+1]);
      c[k-ktop]-=theta*dt*bd[k+1];
    }

    if(src1)
      for(k=ktop;k<grid->Nk[i];k++)
      {
        if(prop->subgrid)
          alpha=subgrid->Acceff[i][k]/Ac;
        else
          alpha=1.0;

        b[k-ktop]+=alpha*grid->dzz[i][k]*src1[i][k]*theta*dt;
        
      }
    // Diffusive fluxes only when more than 1 layer
    if(ktop<grid->Nk[i]-1) {
      // Top cell diffusion
      b[0]+=theta*dt*(bd[ktop+1]+2*alpha_top*bd[ktop+1]);
      c[0]-=theta*dt*bd[ktop+1];

      // Bottom cell diffusion
      a[(grid->Nk[i]-1)-ktop]-=theta*dt*bd[grid->Nk[i]-1];
      b[(grid->Nk[i]-1)-ktop]+=theta*dt*(bd[grid->Nk[i]-1]+2*alpha_bot*bd[grid->Nk[i]-1]);
    }
  
    // Explicit part into source term d[] 
    for(k=ktop+1;k<grid->Nk[i];k++) 
    {
      if(prop->subgrid)
        alpha=subgrid->Acceffold[i][k]/Ac;
      else 
        alpha=1.0;

      d[k-ktop]=alpha*grid->dzzold[i][k]*phys->stmp[i][k];
    }

    if(src1)
    {
      for(k=ktop+1;k<grid->Nk[i];k++) 
      {
        if(prop->subgrid)
          alpha=subgrid->Acceffold[i][k]/Ac;
        else 
          alpha=1.0;
        d[k-ktop]-=alpha*src1[i][k]*(1-theta)*dt*grid->dzzold[i][k]*phys->stmp[i][k];
      }
    }
    d[0]=0;
    if(grid->ctopold[i]<=grid->ctop[i]) {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
      {
        if(prop->subgrid)
          alpha=subgrid->Acceffold[i][k]/Ac;
        else
          alpha=1.0;
        d[0]+=alpha*grid->dzzold[i][k]*phys->stmp[i][k];
      }

      if(src1)
        for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
        {
          if(prop->subgrid)
            alpha=subgrid->Acceffold[i][k]/Ac;
          else
            alpha=1.0;

          d[0]-=alpha*src1[i][k]*(1-theta)*dt*grid->dzzold[i][k]*phys->stmp[i][k];
        }

    } else {
      if(prop->subgrid)
        alpha=subgrid->Acceffold[i][ktop]/Ac;
      else
        alpha=1.0;      
      
      d[0]=alpha*grid->dzzold[i][ktop]*phys->stmp[i][ktop];
      if(src1)
        d[0]-=alpha*src1[i][ktop]*(1-theta)*dt*grid->dzzold[i][ktop]*phys->stmp[i][k];
    }

    // These are the advective components of the tridiagonal
    // that use the new velocity
    if(!(prop->TVD && prop->vertTVD))
      for(k=0;k<grid->Nk[i]+1;k++) {
        //add subgrid part
        if(prop->subgrid)
          alpha=subgrid->Acveffold[i][k]/Ac;
        else
          alpha=1.0;
        if(!prop->im) // if prop->im!=0, ap and am should be the same as above no need to calculate
        {
          ap[k] = alpha*0.5*(phys->wtmp2[i][k]+fabs(phys->wtmp2[i][k]));
          am[k] = alpha*0.5*(phys->wtmp2[i][k]-fabs(phys->wtmp2[i][k]));
        } 
      }
    else // Compute the ap/am for TVD schemes
    {
      if(!prop->im)
        GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
            phys->wtmp2,grid->dzzold,phys->stmp,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
      else
        GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
          phys->wtmp4,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);    
   
      for(k=0;k<grid->Nk[i]+1;k++) 
      {
        //add subgrid part
        if(prop->subgrid)
          alpha=subgrid->Acveffold[i][k]/Ac;
        else
          alpha=1.0;
        ap[k] = alpha*ap[k];
        am[k] = alpha*am[k];
      }     
    }
      
    // Explicit advection and diffusion
    if(!prop->subgrid){
      for(k=ktop+1;k<grid->Nk[i]-1;k++) 
        d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
            (ap[k]-am[k+1])*phys->stmp[i][k]-
            ap[k+1]*phys->stmp[i][k+1])-
          (1-theta)*dt*(bd[k]*phys->stmp[i][k-1]
              -(bd[k]+bd[k+1])*phys->stmp[i][k]
              +bd[k+1]*phys->stmp[i][k+1]);
    } else {
      // recalculate bd by acveffold
      for(k=ktop+1;k<grid->Nk[i];k++)
      {
        alpha=subgrid->Acveffold[i][k]/Ac;
        bd[k]=alpha*(2.0*kappa+kappa_tv[i][k-1]+kappa_tv[i][k])/
          (grid->dzz[i][k-1]+grid->dzz[i][k]);
      }
      for(k=ktop+1;k<grid->Nk[i]-1;k++) 
        d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
            (ap[k]-am[k+1])*phys->stmp[i][k]-
            ap[k+1]*phys->stmp[i][k+1])-
          (1-theta)*dt*(bd[k]*phys->stmp[i][k-1]
              -(bd[k]+bd[k+1])*phys->stmp[i][k]
              +bd[k+1]*phys->stmp[i][k+1]);
    }


    if(ktop<grid->Nk[i]-1) {
      //Flux through bottom of top cell
      k=ktop;
      d[0]=d[0]-(1-theta)*dt*(-am[k+1]*phys->stmp[i][k]-
          ap[k+1]*phys->stmp[i][k+1])+
        (1-theta)*dt*(-(2*alpha_top*bd[k+1]+bd[k+1])*phys->stmp[i][k]+
            bd[k+1]*phys->stmp[i][k+1]);
 
      // subgrid ??? 
      if(Ftop) d[0]+=dt*(1-alpha_top+2*alpha_top*bd[k+1])*Ftop[i];

      // Through top of bottom cell
      k=grid->Nk[i]-1;
      d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
          ap[k]*phys->stmp[i][k])-
        (1-theta)*dt*(bd[k]*phys->stmp[i][k-1]-
            (bd[k]+2*alpha_bot*bd[k])*phys->stmp[i][k]);
      // subgrid ???
      if(Fbot) d[k-ktop]+=dt*(-1+alpha_bot+2*alpha_bot*bd[k])*Fbot[i];
    }

    // First add on the source term from the previous time step.
    if(grid->ctop[i]<=grid->ctopold[i]) {
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++) 
        d[0]+=(1-fab)*Cn[i][grid->ctopold[i]]/(1+fabs(grid->ctop[i]-grid->ctopold[i]));
      for(k=grid->ctopold[i]+1;k<grid->Nk[i];k++) 
        d[k-grid->ctopold[i]]+=(1-fab)*Cn[i][k];
    } else {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++) 
        d[0]+=(1-fab)*Cn[i][k];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        d[k-grid->ctop[i]]+=(1-fab)*Cn[i][k];
    }

    for(k=0;k<grid->ctop[i];k++)
      Cn[i][k]=0;

    if(src2)
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
      { 
        if(prop->subgrid)
          alpha=subgrid->Acceffold[i][k]/Ac;
        else 
          alpha=1.0; 

        Cn[i][k-ktop]=alpha*dt*src2[i][k]*grid->dzzold[i][k];
      }
    else
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        Cn[i][k]=0; 

    // Now create the source term for the current time step
    for(k=0;k<grid->Nk[i];k++){
      ap[k]=0;
      fluxm[k]=0;
      fluxp[k]=0;
    }
    
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) 
      {
        nc2=nc1;
        if(boundary_scal && (grid->mark[ne]==2 || grid->mark[ne]==3))
          sp=phys->stmp2[nc1];
        else
          sp=phys->stmp[nc1];
      } else 
        sp=phys->stmp[nc2];

      if(!(prop->TVD && prop->horiTVD)) {
        for(k=0;k<grid->Nke[ne];k++) {
          
          temp[k]=UpWind((fac2*phys->utmp2[ne][k]+fac1*phys->u[ne][k]),phys->stmp[nc1][k],sp[k]);
          //if(fabs(temp[k]-1)>1e-8)
            //printf("n %d ne %d i %d k %d temp %e u %e %d\n", prop->n,ne,i,k,temp[k],grid->n1[ne]*(fac2*phys->utmp2[ne][k]+fac1*phys->u[ne][k]),prop->im);
          // new edit
          if(k<grid->ctopold[nc1])
            temp[k]=sp[k];
          if(k<grid->ctopold[nc2])
            temp[k]=phys->stmp[nc1][k];
        }
      } else {
        for(k=0;k<grid->Nke[ne];k++) 
          if(phys->utmp2[ne][k]>0)
            temp[k]=phys->SfHp[ne][k];
          else
            temp[k]=phys->SfHm[ne][k];	    
      }
      // this part don't need to be changed
      for(k=0;k<grid->Nke[ne];k++){
        ap[k] += dt*df*normal/Ac*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])
          *temp[k]*grid->dzf[ne][k];

        if((fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*normal>0)
          fluxp[k]+=dt*df*normal/Ac*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*temp[k]*grid->dzf[ne][k];
        else
          fluxm[k]+=dt*df*normal/Ac*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*temp[k]*grid->dzf[ne][k];
        //if(prop->n==2716 && i==6471){
        //  printf("ne %d fac1 %e u1 %e fac2 %e u2 %e fac3 %e active %d, temp %e, dzf %e sum %e\n", ne,fac1,phys->u[ne][k],fac2,phys->utmp2[ne][k],fac3,phys->active[i],temp[k]
        //    ,grid->dzf[ne][k],dt*df*normal/Ac*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*temp[k]*grid->dzf[ne][k]);
        fluxtest+=dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])
          *grid->dzf[ne][k];
        fluxtest1+=dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])
          *temp[k]*grid->dzf[ne][k];
        /*if(prop->n==6898 && i==6375 && myproc==1 && theta!=1){
          printf("i %d ne %d nke %d etop %d nc1 %d xv1 %e nk1 %d ctopold1 %d ctop1 %d h1 %e nc2 %d xv2 %e nk2 %d ctopold2 %d ctop2 %d h2 %e k %d temp %e u %e dzf %e\n",i,
            ne,grid->Nke[ne],grid->etop[ne],nc1,grid->xv[nc1],grid->Nk[nc1],grid->ctopold[nc1],grid->ctop[nc1],phys->h[nc1],nc2,grid->xv[nc2],grid->Nk[nc2],grid->ctopold[nc2],grid->ctop[nc2],phys->h[nc2],
            k,temp[k],phys->u[ne][k]*grid->n1[ne],grid->dzf[ne][k]);
          printf("fluxtest %e flux %e fluxtest1 %e fluxtmp %e diff %e tempdifff %e\n",fluxtest, dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])
            *grid->dzf[ne][k],fluxtest1,dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*temp[k]*grid->dzf[ne][k],dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])
            *(temp[k]-1)*grid->dzf[ne][k],temp[k]-1);
        }*/

        /*if(i==5175 && prop->n==10024){
          printf("n %d proc %d 222 ne %d mark %d sum0 %e flux %e u %e uold %e dzf %e df %e normal %d utmp3 %e fac1 %e fac2 %e sf %e sc %e b[0] %e\n",prop->n,myproc,ne,grid->mark[ne],fluxtest,dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*grid->dzf[ne][k],
            phys->u[ne][0],phys->utmp2[ne][0],grid->dzf[ne][0],df,normal,phys->utmp3[ne][0],fac1,fac2,temp[k],phys->stmp[i][k],subgrid->Veff[i]/grid->Ac[i]);
          printf("acc*dzz %e accold*dzold%e subgrid->Veff %e subgrid->Veffold %e Vdiff %e active %d\n", subgrid->Acceff[i][k]*grid->dzz[i][k],subgrid->Acceffold[i][k]*grid->dzzold[i][k],subgrid->Veff[i],subgrid->Veffold[i],subgrid->Veffold[i]-subgrid->Veff[i],phys->active[i]); 
        }*/
        // printf("flux sum %e veffdiff %e fluxtest %e\n", dt*df*normal*(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])*grid->dzf[ne][k],subgrid->Veffold[i]-subgrid->Veff[i],fluxtest);
      }
    }
    /*if(prop->n==6898 && i==6375 && myproc==1 && theta!=1)
    for(k=ktop;k<grid->Nk[i];k++)
      printf("k %d V %e Vold %e flux %e\n",k,grid->dzz[i][k]*subgrid->Acceff[i][k]/Ac,subgrid->Acceffold[i][k]*grid->dzzold[i][k]/Ac,ap[k]);


    if(prop->n==6898 && i==6375 && myproc==1 && theta!=1)
      printf("Veff %e Veffold %e flux %e diff %e fluxtemp %e\n", subgrid->Veff[i],subgrid->Veffold[i],fluxtest,subgrid->Veff[i]-subgrid->Veffold[i]+fluxtest,fluxtest1);

    sum1=0;
    sum2=0;
    for(k=ktop;k<grid->Nk[i];k++)
    {
      sum1+=a[k-ktop]+b[k-ktop]+c[k-ktop];
      sum2+=d[k-ktop];
      if(prop->n==6898 && i==6375 && myproc==1 && theta!=1)
        printf("k %d sum1 %e sum2 %e d %e a+b+c %e diff %e\n",k,sum1,sum2,d[k-ktop],a[k-ktop]+b[k-ktop]+c[k-ktop],a[k-ktop]+b[k-ktop]+c[k-ktop]-d[k-ktop]);
    }
    if(prop->n==6898 && i==6375 && myproc==1 && theta!=1)
      printf("Veff %e Veffold %e flux %e diff %e sum1 %e sum2 %e\n", subgrid->Veff[i]/Ac,subgrid->Veffold[i]/Ac,fluxtest/Ac,subgrid->Veff[i]/Ac-subgrid->Veffold[i]/Ac
        +fluxtest/Ac,sum1,sum2);*/


    for(k=ktop+1;k<grid->Nk[i];k++) 
      Cn[i][k-ktop]-=ap[k];

    for(k=0;k<=ktop;k++) 
      Cn[i][0]-=ap[k];

    // Add on the source from the current time step to the rhs.
    for(k=0;k<grid->Nk[i]-ktop;k++) 
      d[k]+=fab*Cn[i][k];

    // Add on the volume correction if h was < -d
    /*
       if(grid->ctop[i]==grid->Nk[i]-1)
       d[grid->Nk[i]-ktop-1]+=phys->hcorr[i]*phys->stmp[i][grid->ctop[i]];
       */

    for(k=ktop;k<grid->Nk[i];k++)
      ap[k]=Cn[i][k-ktop];
    for(k=0;k<=ktop;k++)
      Cn[i][k]=0;
    for(k=ktop+1;k<grid->Nk[i];k++)
      Cn[i][k]=ap[k];
    for(k=grid->ctop[i];k<=ktop;k++)
      Cn[i][k]=ap[ktop]/(1+fabs(grid->ctop[i]-ktop));

    if(grid->Nk[i]-ktop>1) 
      TriSolve(a,b,c,d,&(scal[i][ktop]),grid->Nk[i]-ktop);
    else if(prop->n>1) {
      if(b[0]>0 && phys->active[i])
        scal[i][ktop]=d[0]/b[0];
      else 
        scal[i][ktop]=phys->stmp[i][ktop];
    }

    if(grid->Nk[i]-ktop>1 && !phys->active[i]){
      for(k=ktop;k<grid->Nk[i];k++)
        scal[i][k]=phys->stmp[i][k];
    }

    if(prop->subgrid)
      if(ktop!=grid->Nk[i]-1)
      {
        if((fac1*wnew[i][ktop+1]*subgrid->Acveff[i][ktop+1]+fac2*phys->wtmp2[i][ktop+1]*subgrid->Acveffold[i][ktop+1])>0)
          fluxp[ktop]+=(fac1*wnew[i][ktop+1]*subgrid->Acveff[i][ktop+1]+fac2*phys->wtmp2[i][ktop+1]*subgrid->Acveffold[i][ktop+1]);
        else
          fluxm[ktop]+=(fac1*wnew[i][ktop+1]*subgrid->Acveff[i][ktop+1]+fac2*phys->wtmp2[i][ktop+1]*subgrid->Acveffold[i][ktop+1]);                
        
        for(k=ktop+1;k<grid->Nk[i];k++)
        {
          if((fac1*wnew[i][k+1]*subgrid->Acveff[i][k+1]+fac2*phys->wtmp2[i][k+1]*subgrid->Acveffold[i][k+1])>0)
            fluxp[k]+=(fac1*wnew[i][k+1]*subgrid->Acveff[i][k+1]+fac2*phys->wtmp2[i][k+1]*subgrid->Acveffold[i][k+1]);
          else
            fluxm[k]+=(fac1*wnew[i][k+1]*subgrid->Acveff[i][k+1]+fac2*phys->wtmp2[i][k+1]*subgrid->Acveffold[i][k+1]);                
          
          if((fac1*wnew[i][k]*subgrid->Acveff[i][k]+fac2*phys->wtmp2[i][k]*subgrid->Acveffold[i][k])>0)
            fluxp[k]+=(fac1*wnew[i][k]*subgrid->Acveff[i][k]+fac2*phys->wtmp2[i][k]*subgrid->Acveffold[i][k]);
          else
            fluxm[k]+=(fac1*wnew[i][k]*subgrid->Acveff[i][k]+fac2*phys->wtmp2[i][k]*subgrid->Acveffold[i][k]);               
        }

        for(k=ktop;k<grid->Nk[i];k++)
          if(fabs(fluxm[k])>subgrid->Acceff[i][k]*grid->dzz[i][k] || fabs(fluxp[k])>subgrid->Acceff[i][k]*grid->dzz[i][k])
            scal[i][k]=phys->stmp[i][k];
      }

    for(k=0;k<grid->ctop[i];k++)
      scal[i][k]=0;

    for(k=grid->ctop[i];k<grid->ctopold[i];k++) 
      scal[i][k]=scal[i][ktop];
  }

  // Code to check divergence change CHECKCONSISTENCY to 1 in suntans.h
  if(CHECKCONSISTENCY && checkflag) {

    if(prop->n==1+prop->nstart) {
      smin=INFTY;
      smax=-INFTY;
      for(i=0;i<grid->Nc;i++) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(phys->stmp[i][k]>smax) { 
            smax=phys->stmp[i][k]; 
            imax=i; 
            kmax=k; 
          }
          if(phys->stmp[i][k]<smin) { 
            smin=phys->stmp[i][k]; 
            imin=i; 
            kmin=k; 
          }
        }
      }
      MPI_Reduce(&smin,&smin_value,1,MPI_DOUBLE,MPI_MIN,0,comm);
      MPI_Reduce(&smax,&smax_value,1,MPI_DOUBLE,MPI_MAX,0,comm);
      MPI_Bcast(&smin_value,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&smax_value,1,MPI_DOUBLE,0,comm);

      if(myproc==0)
        printf("Minimum scalar: %.2f, maximum: %.2f\n",smin_value,smax_value);
    }      

    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      flag=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        if(grid->mark[grid->face[i*grid->maxfaces+nf]]==2 || 
            grid->mark[grid->face[i*grid->maxfaces+nf]]==3) {
          flag=1;
          break;
        }
      }

      if(!flag) {
        div_da=0;

        for(k=0;k<grid->Nk[i];k++) {
          if(prop->subgrid) 
          {
            div_da+=(subgrid->Acceff[i][k]*grid->dzz[i][k]-subgrid->Acceffold[i][k]
              *grid->dzzold[i][k])/prop->dt;
          } else
            div_da+=grid->Ac[i]*(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;

          div_local=0;
          for(nf=0;nf<grid->nfaces[i];nf++) {
            ne=grid->face[i*grid->maxfaces+nf];
            div_local+=(fac1*phys->u[ne][k]+fac2*phys->utmp2[ne][k]+fac3*phys->utmp3[ne][k])
              *grid->dzf[ne][k]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne];
          }
          div_da+=div_local;
          if(!prop->subgrid)
          {
            div_local+=grid->Ac[i]*(fac1*(wnew[i][k]-wnew[i][k+1])+
                fac2*(phys->wtmp2[i][k]-phys->wtmp2[i][k+1])+fac3*(phys->wtmp3[i][k]-phys->wtmp3[i][k+1]));
          } else {
            div_local+=(fac1*(wnew[i][k]*subgrid->Acveff[i][k]-subgrid->Acveff[i][k+1]*wnew[i][k+1])+
                fac2*(subgrid->Acveffold[i][k]*phys->wtmp2[i][k]
                -subgrid->Acveffold[i][k]*phys->wtmp2[i][k+1])
                +fac3*(subgrid->Acveffold2[i][k]*phys->wtmp3[i][k]
                -subgrid->Acveffold2[i][k]*phys->wtmp3[i][k+1]));                          
          }
          if(fabs(div_local)>1e-10){
            printf("div_local %e i %d k %d\n",div_local,i,k);
            printf("ctop %d xv %e yv %e h %e d %e\n", grid->ctop[i],grid->xv[i],grid->yv[i],phys->h[i],grid->dv[i]);
          }
          if(k>=grid->ctop[i]) {
            if(fabs(div_local)>SMALL_CONSISTENCY && grid->dzz[imin][0]>DRYCELLHEIGHT) 
              printf("Step: %d, proc: %d, locally-divergent at %d, %d, div=%e\n",
                  prop->n,myproc,i,k,div_local);
          }
        }
        if(fabs(div_da)>SMALL_CONSISTENCY && phys->h[i]+grid->dv[i]>DRYCELLHEIGHT)
          printf("i %d  Step: %d, proc: %d, Depth-Ave divergent at i=%d, div=%e\n",
              i,prop->n,myproc,i,div_da);
      }
    }

    mincount=0;
    maxcount=0;
    smin=INFTY;
    smax=-INFTY;
    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      flag=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        if(grid->mark[grid->face[i*grid->maxfaces+nf]]==2 || grid->mark[grid->face[i*grid->maxfaces+nf]]==3) 
        {
          flag=1;
          break;
        }
      }

      if(!flag) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(scal[i][k]>smax) { 
            smax=scal[i][k]; 
            imax=i; 
            kmax=k; 
          }
          if(scal[i][k]<smin) { 
            smin=scal[i][k]; 
            imin=i; 
            kmin=k; 
          }

          if(scal[i][k]>smax_value+SMALL_CONSISTENCY && grid->dzz[i][k]>DRYCELLHEIGHT)
            maxcount++;
          if(scal[i][k]<smin_value-SMALL_CONSISTENCY && grid->dzz[i][k]>DRYCELLHEIGHT)
            mincount++;
        }
      }
    }
    MPI_Reduce(&mincount,&allmincount,1,MPI_INT,MPI_SUM,0,comm);
    MPI_Reduce(&maxcount,&allmaxcount,1,MPI_INT,MPI_SUM,0,comm);

    if(mincount!=0 || maxcount!=0) 
      printf("Not CWC, step: %d, proc: %d, smin = %e at i=%d,H=%e, smax = %e at i=%d,H=%e\n",
          prop->n,myproc,
          smin,imin,phys->h[imin]+grid->dv[imin],
          smax,imax,phys->h[imax]+grid->dv[imax]);

    if(myproc==0 && (allmincount !=0 || allmaxcount !=0))
      printf("Total number of CWC violations (all procs): s<s_min: %d, s>s_max: %d\n",
          allmincount,allmaxcount);
  }
  }
