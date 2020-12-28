#include "mcm.h"

/***************************************************************/
namespace monte{
void gather(int isp)
{
  register int i, j, n;
  register double tempvx, tempvy, tempvz, tempv, delx, dely, w1, w2, w3, w4;

  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) {
      sp_n_0[isp][i][j]  = sp_n_k[isp][i][j] +sp_n_mcc[isp][i][j];
      sp_n_mcc[isp][i][j]= sp_n_k[isp][i][j] = 0.0;
      sp_ex[isp][i][j]   = sp_ey[isp][i][j]  = 0.0;
      sp_vx0[isp][i][j]  = sp_vy0[isp][i][j]= 0.0;
      sp_vz0[isp][i][j]  = sp_vt[isp][i][j] = 0.0;
    }
  
  /* Compute n, vx0, vy0, vz0, vt at grid point */
  for (n=np[isp]-1; n>=0; n--) {
    i = x[isp][n];   delx = x[isp][n] -i;
    j = y[isp][n];   dely = y[isp][n] -j;      
    
    sp_n_k[isp][i][j]    += (w1=(1-delx)*(1-dely));
    sp_n_k[isp][i+1][j]  += (w2=delx*(1-dely));
    sp_n_k[isp][i][j+1]  += (w3=(1-delx)*dely);
    sp_n_k[isp][i+1][j+1]+= (w4=delx*dely);
    
    tempvx= vxscale[isp]*vx[isp][n];
    sp_vx0[isp][i][j]    += w1*tempvx;
    sp_vx0[isp][i+1][j]  += w2*tempvx;
    sp_vx0[isp][i][j+1]  += w3*tempvx;
    sp_vx0[isp][i+1][j+1]+= w4*tempvx;
    
    tempvy= vyscale[isp]*vy[isp][n];
    sp_vy0[isp][i][j]    += w1*tempvy;
    sp_vy0[isp][i+1][j]  += w2*tempvy;
    sp_vy0[isp][i][j+1]  += w3*tempvy;
    sp_vy0[isp][i+1][j+1]+= w4*tempvy;

    tempvz= vzscale[isp]*vz[isp][n];
    sp_vz0[isp][i][j]    += w1*tempvz;
    sp_vz0[isp][i+1][j]  += w2*tempvz;
    sp_vz0[isp][i][j+1]  += w3*tempvz;
    sp_vz0[isp][i+1][j+1]+= w4*tempvz;
   
    tempv = tempvx*tempvx +tempvy*tempvy +tempvz*tempvz;
    sp_vt[isp][i][j]    += w1*tempv;
    sp_vt[isp][i+1][j]  += w2*tempv;
    sp_vt[isp][i][j+1]  += w3*tempv;
    sp_vt[isp][i+1][j+1]+= w4*tempv;
  }
  
  /**** compute v0's and v^2's ****/
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) {
      if(sp_n_k[isp][i][j]>1.0e-10){
	sp_vx0[isp][i][j] /= sp_n_k[isp][i][j];
	sp_vy0[isp][i][j] /= sp_n_k[isp][i][j];
	sp_vz0[isp][i][j] /= sp_n_k[isp][i][j];
	
	sp_vt[isp][i][j]  /= sp_n_k[isp][i][j];
	sp_vt[isp][i][j]  -= sp_vx0[isp][i][j]*sp_vx0[isp][i][j];
	sp_vt[isp][i][j]  -= sp_vy0[isp][i][j]*sp_vy0[isp][i][j];
	sp_vt[isp][i][j]  -= sp_vz0[isp][i][j]*sp_vz0[isp][i][j];
	sp_vt[isp][i][j]   = (sp_vt[isp][i][j]>0.0) ? sqrt(sp_vt[isp][i][j]/3.): 0.0;      
      }
      else {
	sp_vx0[isp][i][j]  = 0.0;
	sp_vy0[isp][i][j]  = 0.0;
	sp_vz0[isp][i][j]  = 0.0;
	sp_vt[isp][i][j]   = 0.0;
      }
    }
  /* Adjust for the boundaries */
  for (j=0; j<= ncy; j++) {
    sp_n_k[isp][0][j] *= 2; sp_n_k[isp][ncx][j] *= 2;
  }
  if(ybc_flag == PERIODIC) {
    for (i=0; i<= ncx; i++) {
      sp_n_k[isp][i][0] = sp_n_k[isp][i][ncy] = sp_n_k[isp][i][0] +sp_n_k[isp][i][ncy];
    }
  }
  else {
    for (i=0; i<= ncx; i++) {
      sp_n_k[isp][i][0] *= 2; sp_n_k[isp][i][ncy] *= 2;
    }
  }
  
  /* To account for the relative weight of each species */
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++)
      sp_n_k[isp][i][j] *= weight[isp];
}

/***************************************************************/

void setrho()
{
  register int i, j, isp, n;
  
  for(isp=0; isp<nsp; isp++) {
    k_count[isp]=0;
    it[isp]=0;
    gather(isp);
    for(i=0; i<=ncx; i++) {
      for(j=0; j<=ncy; j++) {
	sp_n[isp][i][j]= sp_n_k[isp][i][j];
	sp_n_mcc[isp][i][j] =0.0;
      }
    }
    /**** Digitally smooth the species densities in x and y. */
    for (n=0; n<sflag; n++) {
      if(ybc_flag == PERIODIC){
	Periodic_Smooth(sp_n[isp],   ncx, ncy);
	Periodic_Smooth(sp_vx0[isp], ncx, ncy);
	Periodic_Smooth(sp_vy0[isp], ncx, ncy);
	Periodic_Smooth(sp_vz0[isp], ncx, ncy);
	Periodic_Smooth(sp_vt[isp],  ncx, ncy);
      }
      else { 
    TWOD_One_2_One(sp_n[isp],   ncx, ncy);
    TWOD_One_2_One(sp_n[isp],   ncx, ncy);
    TWOD_One_2_One(sp_n[isp],   ncx, ncy);
	TWOD_One_2_One(sp_vx0[isp], ncx, ncy);
	TWOD_One_2_One(sp_vy0[isp], ncx, ncy);
	TWOD_One_2_One(sp_vz0[isp], ncx, ncy);
	TWOD_One_2_One(sp_vt[isp],  ncx, ncy);
      }
    }
  }
}

/***************************************************************/
/* Smoothing theArray using the 1-2-1 method with the proper   */
/* boundary conditions to conserve charge.                     */

void TWOD_One_2_One(double **matrix, int nx, int ny)
{
  register int i, j;
  static int init_smth_flag=1;
  static double **temp;
  
  if(init_smth_flag) {
    temp= (double **)malloc((nx+1)*sizeof(double *));
    for(i=0; i<=nx; i++) temp[i]= (double *)malloc((ny+1)*sizeof(double));
    init_smth_flag=0;
  }
  /********* first do the y pass **********/
  for (i=0; i<=nx; i++) {
    temp[i][0]= (matrix[i][0] +matrix[i][1])/2;
    for(j=1; j< ny; j++) {
      if(grid_mask[i][j] == 1) {
	if(face[i][j] == UP || face[i][j] == UL_CORN || face[i][j] == UR_CORN)
	  temp[i][j]= (matrix[i][j] +matrix[i][j+1])/2;
	else if(face[i][j] == DOWN || face[i][j] == LL_CORN || face[i][j] == LR_CORN)
	  temp[i][j]= (matrix[i][j] +matrix[i][j-1])/2;
	else if(face[i][j] == LEFT || face[i][j] == RIGHT)
	  temp[i][j]= (matrix[i][j-1] +2*matrix[i][j] +matrix[i][j+1])/4;
	else
	  temp[i][j]= matrix[i][j];
      }
      else
	temp[i][j]= (matrix[i][j-1] +2*matrix[i][j] +matrix[i][j+1])/4;
    }
    temp[i][ny]= (matrix[i][ny] +matrix[i][ny-1])/2;
  }
  
  /********* then do the x pass **********/
  for (j=0; j<=ny; j++) {
    matrix[0][j]= (temp[0][j] +temp[1][j])/2;
    for(i=1; i< nx; i++) {
      if(grid_mask[i][j] == 1) {
	if(face[i][j] == LEFT || face[i][j] == UL_CORN || face[i][j] == LL_CORN)
	  matrix[i][j]= (temp[i][j] +temp[i-1][j])/2;
	else if(face[i][j] == RIGHT || face[i][j] == LR_CORN || face[i][j] == UR_CORN)
	  matrix[i][j]= (temp[i][j] +temp[i+1][j])/2;
	else if(face[i][j] == DOWN || face[i][j] == UP)
	  matrix[i][j]= (temp[i-1][j] +2*temp[i][j] +temp[i+1][j])/4;
	else
	  matrix[i][j]= temp[i][j];
      }
      else
	matrix[i][j]= (temp[i-1][j] +2*temp[i][j] +temp[i+1][j])/4;
    } 
    matrix[nx][j]= (temp[nx][j] +temp[nx-1][j])/2;
  }
}

/***************************************************************/
/* Smoothing the array using the 1-2-1 method with the proper   */
/* boundary conditions to conserve charge.                     */

int One_2_One(double *ary, int nmax)
{
  register int i;
  static int nlocal=0;
  static double *temp;
  
  if(nlocal < nmax) {
    free(temp);
    temp = (double *)malloc((nmax+1)*sizeof(double));
    nlocal = nmax;
  }
  
  temp[0]= (ary[0] +ary[1])/2;
  for(i=1; i< nmax; i++)  temp[i]= (ary[i-1] +2*ary[i] +ary[i+1])/4;
  temp[nmax]= (ary[nmax] +ary[nmax-1])/2;
  
  for(i=0; i<= nmax; i++) ary[i]= temp[i];
}
}
/***************************************************************/
