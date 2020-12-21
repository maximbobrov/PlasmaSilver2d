#include "mcm.h"
#define POINTS_IN_BETWEEN  5
namespace monte{
void find_intercept(double xp, double yp, double vxp, double vyp,
                    int k, double *xint, double *yint, int *side);
int  emitE(int isp, double xp, double yp);
/***************************************************************/

void move(int isp)
{
  register int n, i, j;
  register double delx, dely, w1, w2, w3, w4;
  
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) {
      ax[i][j] = sp_ex[isp][i][j]*ax_scale[isp];
      ay[i][j] = sp_ey[isp][i][j]*ay_scale[isp];
    }
  for (n=0; n<np[isp]; n++) {
    i = x[isp][n];   delx = x[isp][n] -i;
    j = y[isp][n];   dely = y[isp][n] -j;
    
    vx[isp][n]+= ax[i][j]*(w1=(1-delx)*(1-dely)) +ax[i+1][j]*(w2=delx*(1-dely))
      +ax[i][j+1]*(w3=(1-delx)*dely) +ax[i+1][j+1]*(w4=delx*dely);
    
    vy[isp][n]+= ay[i][j]*w1 +ay[i+1][j]*w2 +ay[i][j+1]*w3 +ay[i+1][j+1]*w4;
    
    x[isp][n] += vx[isp][n];
    y[isp][n] += vy[isp][n];
  }
}

/***************************************************************/

void move_internal(int isp)
{
  register int n, i, j, k;
  register double delx, dely, w1, w2, w3, w4;
  register int ix, iy, xold, yold;
  double xint, yint, vx_int, vy_int;
  int side;
  
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) {
      ax[i][j] = sp_ex[isp][i][j]*ax_scale[isp];
      ay[i][j] = sp_ey[isp][i][j]*ay_scale[isp];
    }
  for (n=0; n<np[isp]; n++) {
    i = x[isp][n];   delx = x[isp][n] -i;
    j = y[isp][n];   dely = y[isp][n] -j;
    
    vx[isp][n]+= ax[i][j]*(w1=(1-delx)*(1-dely)) +ax[i+1][j]*(w2=delx*(1-dely))
      +ax[i][j+1]*(w3=(1-delx)*dely) +ax[i+1][j+1]*(w4=delx*dely);
    
    vy[isp][n]+= ay[i][j]*w1 +ay[i+1][j]*w2 +ay[i][j+1]*w3 +ay[i+1][j+1]*w4;
    
    if(vx[isp][n]>1.0 || vy[isp][n]>1.0) {
      vx_int=vx[isp][n]/POINTS_IN_BETWEEN;
      vy_int=vy[isp][n]/POINTS_IN_BETWEEN;
      xold = x[isp][n];
      yold = y[isp][n];
      for(k=0; k<POINTS_IN_BETWEEN; k++) {
	x[isp][n] += vx_int;
	y[isp][n] += vy_int;
	ix= x[isp][n];
	iy= y[isp][n];
	if(ix<0 || ix>ncx || iy<0 || iy>ncy){
	  x[isp][n] = xold +vx[isp][n];
	  y[isp][n] = yold +vy[isp][n];
	  break;
	}
	else if(cell_mask[ix][iy]) {
	  stc_np[strc[cell_k[ix][iy]].type][isp]++;
	  find_intercept(x[isp][n], y[isp][n], vx[isp][n], vy[isp][n],
			 cell_k[ix][iy], &xint, &yint, &side);
	  ix= xint;
	  iy= yint;
	  if(side == DOWN || side == UP) {
	    delx= xint -ix;
	    sp_sigma[isp][ix][iy]  += 1-delx;
	    sp_sigma[isp][ix+1][iy]+= delx;
	  }
	  else { 
	    dely= yint -iy;
	    sp_sigma[isp][ix][iy]  += 1-dely;
	    sp_sigma[isp][ix][iy+1]+= dely;
	  }
	  
	  /* Inject a secondary */
	  if(sec_flag[isp] && frand() < grid_sec[isp])
        emitE(sec_sp[isp], xint, yint); //, side);
	  
	  /* Remove the particle */
	  x[isp][n]  =  x[isp][np[isp]-1];
	  y[isp][n]  =  y[isp][np[isp]-1];
	  vx[isp][n] = vx[isp][np[isp]-1];
	  vy[isp][n] = vy[isp][np[isp]-1];
	  vz[isp][n] = vz[isp][np[isp]-1];
	  n--;
	  np[isp]--;
	  break;
	}
      }
    }
    else {
      x[isp][n] += vx[isp][n];
      y[isp][n] += vy[isp][n];
    }
  }
}

/***************************************************************/

void mag_move(int isp)
{
  register int n, i, j;
  double vxtemp, vytemp, vztemp;
  double delx, dely, axtemp, aytemp;
  double w1, w2, w3, w4;
  
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) {
      ax[i][j] = sp_ex[isp][i][j]*ax_scale[isp];
      ay[i][j] = sp_ey[isp][i][j]*ay_scale[isp];
    }
  
  for (n=np[isp]-1; n>=0; n--) {
    i = x[isp][n];   delx = x[isp][n] -i;
    j = y[isp][n];   dely = y[isp][n] -j;
    
    axtemp = ax[i][j]*(w1=(1-delx)*(1-dely)) +ax[i+1][j]*(w2=delx*(1-dely))
      +ax[i][j+1]*(w3=(1-delx)*dely) +ax[i+1][j+1]*(w4=delx*dely);
    aytemp = ay[i][j]*w1 +ay[i+1][j]*w2 +ay[i][j+1]*w3 +ay[i+1][j+1]*w4;
    
    vxtemp = vx[isp][n] +axtemp;
    vytemp = vy[isp][n] +aytemp;
    vztemp = vz[isp][n];
    
    vx[isp][n]= vxtemp + vytemp*ttz[isp] - vztemp*tty[isp];
    vy[isp][n]= vytemp + vztemp*ttx[isp] - vxtemp*ttz[isp];
    vz[isp][n]= vztemp + vxtemp*tty[isp] - vytemp*ttx[isp];
    
    vxtemp += vy[isp][n]*ssz[isp] - vz[isp][n]*ssy[isp];
    vytemp += vz[isp][n]*ssx[isp] - vx[isp][n]*ssz[isp];
    vztemp += vx[isp][n]*ssy[isp] - vy[isp][n]*ssx[isp];
    
    vx[isp][n] = vxtemp +axtemp;
    vy[isp][n] = vytemp +aytemp;
    vz[isp][n] = vztemp;
    
    x[isp][n] += dtdx[isp]*vx[isp][n];
    y[isp][n] += dtdy[isp]*vy[isp][n];
  }
}
}
/***************************************************************/
