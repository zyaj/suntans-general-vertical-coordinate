/*
 * Boundaries test file.
 *
 */
#include "boundaries.h"
#include "sediments.h"
#include "initialization.h"
// #include "wave.h"
#include "util.h"
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
      ub[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->n1[j]+
	phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->n2[j];
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
  int i,iptr,jptr, j, ib, k;
  REAL z;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      z=0;
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	z-=grid->dzz[ib][k]/2;
	phys->boundary_T[jptr-grid->edgedist[2]][k]=ReturnTemperature(0,0,z,0);
	z-=grid->dzz[ib][k]/2;
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
  int iptr, i, jptr, j, ib, k, ne, nc, nf;
  REAL z, u0, v0, w0, A0, dotN[4], nxb, nyb, ustar, urms;

  REAL t_relax = 4*3600;

  // for stage set boundaries, set the free surface to zero
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    phys->h[i]=0;
    for(k=0;k<grid->Nk[i];k++) {
      phys->uc[i][k]=0.88;
      phys->vc[i][k]=0;
      phys->w[i][k]=0;
    }
  }

  // set velocity for open boundaries
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    z=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      z-=grid->dz[k]/2;
      
      phys->boundary_u[jptr-grid->edgedist[2]][k] = 0.88*(1-exp(-4*prop->rtime/t_relax));
      phys->boundary_v[jptr-grid->edgedist[2]][k]=0;
      phys->boundary_w[jptr-grid->edgedist[2]][k]=0;
      
      z-=grid->dz[k]/2;
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

void InitBoundaryData(propT *prop, gridT *grid, int myproc, MPI_Comm comm){}
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc, MPI_Comm comm){}
void BoundarySediment(gridT *grid, physT *phys, propT *prop) {}

/*
 * Function: WindSpeedandDirection
 * usage: calculate Wind field when Wind is not constant
 * -------------------------------------------------
 * calculate Uwind and Winddir
 *
 */
void FetchWindSpeedandDirection(gridT *grid, propT *prop, int myproc){}

/*
 * Function: ReturnCellWindStress
 * usage: calculate Wind field when Wind is not constant
 * -------------------------------------------------
 * calculate Uwind and Winddir
 *
 */
REAL ReturnCellWindStress(REAL x, REAL y, int n, REAL dt, int xyswitch){}



/*
 * Function: UserDefineFunction
 * usage: Define user function to output results we want
 * -------------------------------------------------
 *
 */

// void UserDefinedFunction(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm, int numprocs)
 void UserDefinedFunction(gridT *grid, physT *phys, propT *prop, int myproc)
{
//    REAL **p_prime,p_uw,e_flux, A_f,Q,u_0,**up_prime,**wp_prime,*ubar;
//    int i,j,k,jptr,nc1,nc2,nc_uw;

//   // if(prop->n%40==0)
//   if(!(prop->n%prop->ntout))
//   {
//    // initialize p_prime as 0 in all cell centers
//    p_prime=phys->user_def_nc;
//    for(i=0;i<grid->Nc;i++)
//    {
//      for(k=0;k<grid->Nk[i];k++)
//        p_prime[i][k]=0;
//    }

//    // calculate p_prime;
//    for(i=0;i<grid->Nc;i++)
//    {
//      // p_prime[i][grid->ctop[i]]=0.5*phys->RHO0*phys->rho[i][grid->ctop[i]]*grid->dzz[i][grid->ctop[i]]*prop->grav;
//     p_prime[i][grid->ctop[i]]=0.5*RHO0*phys->rho[i][grid->ctop[i]]*grid->dzz[i][grid->ctop[i]]*prop->grav;

//      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
//      {
//        p_prime[i][k]=p_prime[i][k-1]+0.5*RHO0*prop->grav*(phys->rho[i][k]*grid->dzz[i][k]+phys->rho[i][k-1]*grid->dzz[i][k-1]);
//      }
//    }

//    // calculate up_prime and wp_prime for each cell
//    // first initialize them with zero (along with a depth averaged horizontal velocity, ubar)
//    up_prime=phys->user_def_nc;
//    ubar=phys->user_def_plan;
//    wp_prime=phys->user_def_nc;
//    for(i=0;i<grid->Nc;i++)
//    {
//      ubar[i] = 0;
//      for(k=0;k<grid->Nk[i];k++)
//        up_prime[i][k]=0;
//        wp_prime[i][k]=0;
//    }

//    // find a depth averaged horizontal velocity to subtract from u to get u_prime
//    for(i=0;i<grid->Nc;i++)
//    {
//       A_f=0;
//       Q=0;
//       for(k=grid->etop[i];k<grid->Nke[i];k++)
//       {
//          A_f+=grid->dzf[i][k]*grid->df[i];
//          Q+=phys->u[i][k]*grid->n1[i]*grid->dzf[i][k]*grid->df[i];
//       }
//       ubar[i]=Q/A_f;
//    }

//    // calculate up_prime and wp_prime 
//    for(j=0;j<grid->Nc;j++)
//    {
//       for(k=grid->etop[j];k<grid->Nke[j];k++)
//       {
//          up_prime[j][k] = (phys->uc[j][k]-ubar[j])*p_prime[j][k];
//          wp_prime[j][k] = phys->wc[j][k]*p_prime[j][k];
//       }
//    }

//    // write these to a file
//     FILE *ufp = fopen( "up.dat" , "w" );
//     FILE *wfp = fopen( "wp.dat" , "w" );

//     Write3DData(up_prime,phys->htmp,prop->mergeArrays,ufp,
//     "Error outputting up-data!\n",grid,numprocs,myproc,comm);
//     Write3DData(wp_prime,phys->htmp,prop->mergeArrays,wfp,
//     "Error outputting wp-data!\n",grid,numprocs,myproc,comm);




//    // // calculate leewave energy flux
//    // for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) 
//    // {
//    //    j = grid->edgep[jptr];
//    //    nc1 = grid->grad[2*j];
//    //    nc2 = grid->grad[2*j+1];
//    //    if(nc1==-1) nc1=nc2;
//    //    if(nc2==-1) nc2=nc1;
//    //    e_flux=0;
//    //    A_f=0;
//    //    Q=0;
//    //    for(k=grid->etop[j];k<grid->Nke[j];k++)
//    //    {
//    //       A_f+=grid->dzf[j][k]*grid->df[j];
//    //       Q+=phys->u[j][k]*grid->n1[j]*grid->dzf[j][k]*grid->df[j];
//    //       p_uw=UpWind(phys->u[j][k],p_prime[nc1][k],p_prime[nc2][k]);
//    //       e_flux+=phys->u[j][k]*grid->n1[j]*grid->dzf[j][k]*grid->df[j]*p_uw;
//    //    }
//    //    u_0=Q/A_f;
//    //    printf("proc %d xe %e Af %e Q %e u0 %e eflux %e\n",myproc,grid->xe[j],A_f,Q,u_0,e_flux);
//    // }

//    // // Now go for the vertical energy flux
//    // // pick a depth about halfway up
//    // int k_half = grid->Nkmax/2;
//    // REAL depth_of_e_flux = DepthFromDZ(grid, phys, i/2, k_half);;

//    // // set a counter and initialize the sum of w_prime*p_prime and the vertical upwind p_prime
//    // // also offer a variable to use comparing sum(abs(w * p)) with sum(w*p).
//    // int count=0;
//    // REAL sum_w_p=0, sum_abs_w_p=0, p_w_uw, e_flux_vertical=0, sum_vertical_area=0, e_flux_per_area;
   

//    // // loop through the grid at this depth, computing w_prime*p_prime at each vertical face
//    // for(i=0; i<grid->Nc; i++){


//    //  // do not need to filter out an area... so the following commented out if clause is void
//    //  // but if you wanted to keep from including cells with no IGW activity in the average, 
//    //  // this if clause checks that this cell 
//    //  // has a non-negligable perturbation velocity (in all directions!)
//    //  //if(fabs(phys->u[i][k_half] - u_0) > 0.0001 && fabs(phys->w[i][k_half]) > 0.0001){

//    //    //find p_prime upwind of the vertical face
//    //    p_w_uw = UpWind(phys->w[i][k_half],p_prime[i][k_half],p_prime[i][k_half-1]);

//    //    //compute sum_w_p (times the area of the vertical face!) and increment count
//    //    sum_w_p += phys->w[i][k_half] * p_w_uw;
//    //    sum_abs_w_p += fabs(phys->w[i][k_half] * p_w_uw);
//    //    sum_vertical_area += grid->Ac[i];
//    //    e_flux_vertical += phys->w[i][k_half] * p_w_uw * grid->Ac[i];
//    //    // sum_w_p += phys->w[i][k_half] * p_prime[i][k_half];
//    //    count+=1;
//    //  // this bracket for the commented out if clause
//    //  //}
//    // }

//    // // compute the average of w_prime*p_prime
//    // e_flux_per_area = sum_w_p/count;

//    // printf("proc %d vertical eflux at depth %e  over area of %e meters is %e watts \n",myproc,depth_of_e_flux, sum_vertical_area, e_flux_vertical);
//    // printf("average vertical eflux per area is %e watts/m^2 \n",e_flux_per_area);

   

//    // // make an array of e_flux information to output
//    // // REAL eflux_output[4] = {(REAL) myproc,depth_of_e_flux, sum_vertical_area, e_flux_vertical};
//    // // make a pointer to this array
//    // // REAL *eflux_out_ptr= &eflux_output;

//    // // // make a file pointer
//    // // FILE *fp;
//    // // // open the file eflux.dat to append
//    // // fp = fopen( "eflux.dat" , "a" );   
//    // // // write the array
//    // // fprintf(fp, "%d %e %e %e \n", myproc, depth_of_e_flux, sum_vertical_area, e_flux_vertical);
//    // // // fwrite(eflux_out_ptr,sizeof(REAL),4,fp);    
//    // // //close the file
//    // // fclose(fp);
 // }
}
