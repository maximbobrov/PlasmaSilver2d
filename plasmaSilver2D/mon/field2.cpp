#include "mcm.h"
namespace monte{
double     lhs_source2();
double     *gam, b0;
double     se;
void      circuit1(), circuit2(), circuit3(), circuit4(), (*circuitptr)();
double     (*bcptr)(), bc1(), bc2(), bc3(), bc4();
int       nstrt;
/*************************************************/

void fields2()
{
  static int init_flag=1;
  int isp, i, j, k;
  double s, bet;
  
  if(init_flag)
  {
    gam  = (double *)malloc(ngx*sizeof(double));
    se = dx*dx/epsilon;

    /* Deciding which circuit solver to use.. */
    if(!lhs_flag) {  /* SHORT CIRCUIT */
      bcptr = bc3;
      circuitptr = circuit3;
      nstrt = 1;
      b0 = -2.0;
    }
    else if(lhs_extc < 1e-30) {   /* When C tends to 0.0 */
      /*  OPEN CIRCUIT */
      if(fabs(lhs_dc)+fabs(lhs_ac) > 0)
	puts("START: Active external source is ignored since C < 1E-30\n");
      
      bcptr = bc1;
      circuitptr = circuit1;
      
      nstrt = 0;
      b0 = -1.0;
    }
    else if(lhs_flag==2) {   /* When current source is applied */
      /*  CURRENT SOURCE */
      bcptr = bc2;
      circuitptr = circuit2;
      nstrt = 0;
      b0 = -1.0;
    }
    else if(lhs_extl < 1e-30 && lhs_extr < 1e-30 && lhs_extc >= 1e5*epsilon*zlength*ylength/xlength) {
      /* When L=R=0.0 and C tends to infinity ----  SHORT CIRCUIT */
      bcptr = bc3;
      circuitptr = circuit3;
      
      nstrt = 1;
      b0 = -2.0;
    }
    else {                   /* The general case with external voltage source */
      /*  GENERAL CASE */
      bcptr = bc4;
      circuitptr = circuit4;
      
      nstrt = 0;
      b0 = 2.25*lhs_extl/dt/dt + 1.5*lhs_extr/dt + 1/lhs_extc;
      b0 = -(1 + dx/b0/ylength/zlength/epsilon);
    }
    init_flag = 0;
    for (j=0; j<=ncy; j++) ey[0][j] = ey[ncx][j] = 0.0;
  }
  
  /**************************************************/
  /*** Get the charge densities..                ****/
  
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) rho[i][j] = rhoback;
  
  for (isp=0; isp< nsp; isp++)
    for (i=0; i<= ncx; i++)
      for (j=0; j<= ncy; j++) rho[i][j] += sp_n[isp][i][j]*q_per_cell[isp];
  
  /****************************************************/

  /****************************************************/
  for (i=0; i<= ncx; i++) {
    for (j=0; j<=ncy; j++) rhok[i][j] = rho[i][j];
    realft(rhok[i]-1, khi, 1);
  }
  /*****************************************************/
  /* Solve a tridiagnal matrix to get Phi[i] at time t */
  
  /* first Ky=0 */
  for (k=0; k<ncy; k++)
    phik[0][k] = phik[ncx][k] = 0.0;

  lhs_conv_chrg = 0.0;
  for(isp=0; isp<nsp; isp++) lhs_conv_chrg += lhs_chrg[isp];
  
  if (nstrt)  phik[0][0] = lhs_source2();
  bet = b0;
  phik[nstrt][0] = (*bcptr)()/bet;
  
  for (i=nstrt+1; i<ncx; i++) {
    gam[i] = 1./bet;
    bet = -2. - gam[i];
    phik[i][0] = -(se*rhok[i][0] +phik[i-1][0])/bet;
  }
  
  for(i=ncx-2; i>=nstrt; i--) phik[i][0] -= gam[i+1]*phik[i+1][0];

  /* now Ky!=0 */
  for (k=1; k < ncy; k++) {
    bet = -2.0*dk[k];
    phik[1][k] = -se*rhok[1][k]/bet;
    for (i=2; i< ncx; i++) {
      gam[i] = 1.0/bet;
      bet = -2.0*dk[k] - gam[i];
      phik[i][k] = -(se*rhok[i][k] +phik[i-1][k])/bet;
    }
    for(i=ncx-2; i>=1; i--) 
      phik[i][k] -= gam[i+1]*phik[i+1][k];
  }
  
  /****************************************************/
  
  for (i=0; i<= ncx; i++) {
    for (j=0; j< ncy; j++) phi[i][j] = phik[i][j];
    realft(phi[i]-1, khi, -1); 
    phi[i][ncy]= phi[i][0];
  }
  
  /********************************************************************/


  /* ..and calculate E(x) from phi(x)... */
  
  s = 0.5/dx;
  for(j=0; j <= ncy; j++) {
    ex[0][j] = (phi[0][j] - phi[1][j])/dx - rho[0][j]*0.5*dx/epsilon;
    ex[ncx][j] = (phi[ncx-1][j] - phi[ncx][j])/dx + rho[ncx][j]*0.5*dx/epsilon;
    sigma[0][j] = ex[0][j]*epsilon;
    for (i=1; i< ncx; i++) ex[i][j] = (phi[i-1][j] - phi[i+1][j])*s;
  }
  
  s = 0.5/dy;
  for(i=0; i <= ncx; i++) {
    ey[i][ncy] = ey[i][0] = (phi[i][ncy-1] -phi[i][1])*s;
    for (j=1; j< ncy; j++) ey[i][j] = (phi[i][j-1] - phi[i][j+1])*s;
  }

  (*circuitptr)();
  lhs_phi0=phik[0][0];

  /**************************************************/
  /***  and finally, calculate the species E field  */
  
  /*ek: comment out since it will be done in pdp2.c */
  /* for(isp=0; isp<nsp; isp++)
    for(i=0; i <= ncx; i++)
      for(j=0; j <= ncy; j++) {
	sp_ex[isp][i][j] += ex[i][j];
	sp_ey[isp][i][j] += ey[i][j];
      }
  */     

}   /* end of FIELDS */

/***************************************************************/

double lhs_source2()
{
  if(lhs_flag) return (lhs_dc + lhs_ac*sin(t*lhs_w0 + lhs_theta0));
  return(0.0);
}

/***************************************************************/

/***************************************************************/
/* When C = 0: OPEN CIRCUIT */

double bc1()
{
  return (-(lhs_sigma_total + lhs_conv_chrg/del_area + 0.5*rhok[0][0]*dx)*dx/epsilon);
}

void circuit1()
{
  lhs_sigma_total_1 = lhs_sigma_total;
  lhs_sigma_total += lhs_conv_chrg/del_area;
}

/***************************************************************/
/* When CURRENT SOURCE is applied */

double bc2()
{
  return (-(lhs_sigma_total + (dt*lhs_source2() + lhs_conv_chrg)/del_area  + 0.5*rhok[0][0]*dx)*dx/epsilon);
}

void circuit2()
{
  lhs_exti  = lhs_source2();
  lhs_sigma_total_1 = lhs_sigma_total;
  lhs_sigma_total += (lhs_exti*dt +lhs_conv_chrg)/del_area;
}

/***************************************************************/
/* When R=L=0, C -> infinity: SHORT CIRCUIT */

double bc3()
{
  return (-rhok[1][0]*dx*dx/epsilon - lhs_source2());
}

void circuit3()
{
  lhs_sigma_total_1 = lhs_sigma_total;
  lhs_sigma_total = (phik[0][0] - phik[1][0])*epsilon/dx - rhok[0][0]*dx*0.5;
  lhs_exti = del_area*(lhs_sigma_total - lhs_sigma_total_1)/dt - lhs_conv_chrg/dt;
}

/***************************************************************/
/* The General case */

double bc4()
{
  static double a0, a1, a2, a3, a4;
  double k;
  
  if(!a0) {
    a0 = 2.25*lhs_extl/dt/dt + 1.5*lhs_extr/dt + 1/lhs_extc;
    a1 = -6*lhs_extl/dt/dt - 2*lhs_extr/dt;
    a2 = 5.5*lhs_extl/dt/dt + .5*lhs_extr/dt;
    a3 = -2*lhs_extl/dt/dt;
    a4 = .25*lhs_extl/dt/dt;
  }
  k = (a1*lhs_extq + a2*lhs_extq_1 + a3*lhs_extq_2 + a4*lhs_extq_3)/a0;
  return (-(0.5*rhok[0][0]*dx + lhs_sigma_total + lhs_conv_chrg/del_area
	    + (lhs_source2()/a0 - k - lhs_extq)/del_area)*dx/epsilon);
}

void circuit4()
{
  static double a0, a1, a2, a3, a4;
  double k;
  
  if(!a0) {
    a0 = 2.25*lhs_extl/dt/dt + 1.5*lhs_extr/dt + 1/lhs_extc;
    a1 = -6*lhs_extl/dt/dt - 2*lhs_extr/dt;
    a2 = 5.5*lhs_extl/dt/dt + .5*lhs_extr/dt;
    a3 = -2*lhs_extl/dt/dt;
    a4 = .25*lhs_extl/dt/dt;
  }
  k = (a1*lhs_extq + a2*lhs_extq_1 + a3*lhs_extq_2 + a4*lhs_extq_3)/a0;
  lhs_extq_3 = lhs_extq_2;
  lhs_extq_2 = lhs_extq_1;
  lhs_extq_1 = lhs_extq;
  lhs_extq = (lhs_source2() - phik[0][0])/a0 - k;
  
  lhs_exti = (lhs_extq - lhs_extq_1)/dt;
  lhs_sigma_total_1 = lhs_sigma_total;
  lhs_sigma_total  += (lhs_conv_chrg + lhs_exti*dt)/del_area;
}

/***************************************************************/
/***************************************************************/
/*                               1 2 1                          */
/* Smoothing the array using the 1 4 1 method with the proper  */
/*                               1 2 1                          */
/* boundary conditions to conserve charge.                     */

void Periodic_Smooth(double **a, int nx, int ny)
{
  int i,j;
  static int init_flag = 1;
  static double **temp_smooth;
  
  if (init_flag) {
    temp_smooth= (double **)malloc((nx+1)*sizeof(double *));
    for(i=0; i<=nx; i++) {
      if(!(temp_smooth[i]= (double *)malloc((ny+1)*sizeof(double)))) {
	puts("null ptr in 1-4-1");
	exit(1);
      }
    }
    init_flag=0;
  }
  
  for (j=1; j<ny; j++) {
    temp_smooth[0][j] = (2.0*a[0][j] + a[0][j+1] + a[0][j-1] + 2.0*a[1][j] + a[1][j+1] + a[1][j-1])/8.0;
    temp_smooth[nx][j] = (2.0*a[nx][j] + a[nx][j+1] + a[nx][j-1] + 2.0*a[nx-1][j] + a[nx-1][j+1] + a[nx-1][j-1])/8.0;
    for (i=1; i<nx; i++)
      temp_smooth[i][j] = (4.0*a[i][j] + 2.0*a[i][j+1] +2.0*a[i][j-1] + 2.0*a[i+1][j] + 2.0*a[i-1][j] + a[i+1][j+1] + a[i-1][j+1] + a[i+1][j-1] + a[i-1][j-1])/16.0;
  }

  temp_smooth[0][ny] = temp_smooth[0][0] = (2.0*a[0][0] + a[0][1] + a[0][ny-1] + 2.0*a[1][0] + a[1][1] + a[1][ny-1])/8.0;
  temp_smooth[nx][ny] = temp_smooth[nx][0] = (2.0*a[nx][0] + a[nx][1] + a[nx][ny-1] + 2.0*a[nx-1][0] + a[nx-1][1] + a[nx-1][ny-1])/8.0;
  for (i=1; i<nx; i++)
    temp_smooth[i][ny] = temp_smooth[i][0] = (4.0*a[i][0] + 2.0*a[i][1] + 2.0*a[i][ny-1] + 2.0*a[i+1][0] + 2.0*a[i-1][0] + a[i+1][1] + a[i+1][ny-1] + a[i-1][1] + a[i-1][ny-1])/16.0;
  
  for(i=0; i<=nx; i++) 
    for(j=0; j<=ny; j++) a[i][j]=temp_smooth[i][j];
}
}
/***************************************************************/
