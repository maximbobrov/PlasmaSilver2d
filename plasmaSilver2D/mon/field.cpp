#include "mcm.h"

#define my_div(x, y)     ((fabs(y) <= 1e-10) ? 0 : x/y)
namespace monte{
double     lhs_source();
void      lhs_circuit(), get_lhs_phi0();
double     lhs_a, lhs_b, lhs_c, lhs_bc_den;

double     rhs_source();
void      rhs_circuit(), get_rhs_phi0();
double     rhs_a, rhs_b, rhs_c, rhs_bc_den;

double     hdx, hdy;
double dadi_dt;
 
/***************************************************************/

void fields()
{
  register int isp, i, j, n, ip1, im1, jp1, jm1;
  double ex_mag, ey_mag, e_mag;
  static int init_field_flag=1, stat_counter=0;

  if(init_field_flag) {
    dadi_dt=0.1*(dx*dx+dy*dy);
    init_fields();
    init_field_flag=0;
    
    for(i=0; i<= ncx; i++)
      for(j=0; j<=ncy; j++) eps_array[i][j] *= EPS0;
  }
  
  /**************************************************/
  /*** Get the charge densities..                ****/
  
  for (i=0; i<= ncx; i++)
    for (j=0; j<= ncy; j++) {
      rho[i][j] = rhoback;
      sigma[i][j] = 0.0;
      source[i][j]= rhoback;
    }
  
  for (isp=0; isp< nsp; isp++)
    for (i=0; i<= ncx; i++)
      for (j=0; j<= ncy; j++) {
    source[i][j]-= q_per_cell[isp]*(sp_n[isp][i][j] +sp_sigma[isp][i][j])/EPS0;
	rho[i][j]   += q_per_cell[isp]*sp_n[isp][i][j];
	sigma[i][j] += q[isp]*sp_sigma[isp][i][j];
      }
  
  /*********************************************************/
  /*** Solve Poisson's eqn with zero boundary condition  ***/
  dadi_dt= dadi(dadi_dt, phi_pois, source, 20, tol_pois, 0,
        (double *)NULL, (double *)NULL, (double *)NULL, (double *)NULL);
  
  for (i=0; i< ngx; i++)
    for (j=0; j< ngy; j++) phi[i][j] = phi_pois[i][j];
  
  /***  Add the solution to adjusted Laplace's eqn with lhs bc ***/
  if(lhs_flag) {
    get_lhs_phi0();
    for (i=0; i< ngx; i++)
      for (j=0; j< ngy; j++) phi[i][j] += lhs_phi0*phi_lap_lhs[i][j];
  }
  
  /***  Add the solution to adjusted Laplace's eqn with rhs bc ***/
  if(rhs_flag) {
    get_rhs_phi0();
    for (i=0; i< ngx; i++)
      for (j=0; j< ngy; j++) phi[i][j] += rhs_phi0*phi_lap_rhs[i][j];
  }

  //for (j=0; j< ngy; j++) phi[ngx-1][j] = -1000;//phi[ngx-2][j];

    /*if(nstrcmax>0)
    for (i=0; i< ngx; i++)
      for (j=0; j< ngy; j++) phi[i][j] += phi_intl[i][j];*/

  /**************************************************/
  /*** Calculate the sigma's..                    ***/
  
  for (i=0; i<=ncx; i++) {
    for (j=0; j<=ncy; j++) {
      if(conductor[i][j]==0 && face[i][j] != NO_FACE) {
	ip1= (i==ncx)?ncx:i+1;
	jp1= (j==ncy)?ncy:j+1;
	im1= (i)?i-1:0;
	jm1= (j)?j-1:0;
	sigma[i][j] = zlength*dy*0.5*(eps_array[i][j]+eps_array[i][jm1])*(phi[i][j]-phi[ip1][j])/dx;
	sigma[i][j]+= zlength*dy*0.5*(eps_array[im1][j]+eps_array[im1][jm1])*(phi[i][j]-phi[im1][j])/dx;
	sigma[i][j]+= zlength*dx*0.5*(eps_array[i][j]+eps_array[im1][j])*(phi[i][j]-phi[i][jp1])/dy;
	sigma[i][j]+= zlength*dx*0.5*(eps_array[i][jm1]+eps_array[im1][jm1])*(phi[i][j] -phi[i][jm1])/dy;
	sigma[i][j]-= zlength*dx*dy*rho[i][j];
      }
      sigma[i][j] /= area[i][j];
    }
  }
  
  /**************************************************/
  /*** External circuit stuff                     ***/
  
  if(lhs_flag) lhs_circuit();
  if(rhs_flag) rhs_circuit();
  
  /**************************************************/
  /***  calculate the E field...                  ***/
  
  for(i=0; i<=ncx; i++) {
    for(j=0; j<=ncy; j++) {
      ip1= (i==ncx)?ncx:i+1;
      jp1= (j==ncy)?ncy:j+1;
      im1= (i)?i-1:0;
      jm1= (j)?j-1:0;
      if(face[i][j] == LEFT) {
	if(conductor[i][j]==0) {
	  ex[i][j] = -2.0*sigma[i][j]/(eps_array[im1][j]+eps_array[im1][jm1]+SNG_MIN);
	  ey[i][j] = 0.0;
	}
	else {
	  ex[i][j] = sigma[i][j];
	  ex[i][j]+= 0.5*(eps_array[i][j]+eps_array[i][jm1])*(phi[ip1][j] -phi[i][j])/dx;
	  ex[i][j]-= (0.5*dx/dy)*eps_array[i][j]*(phi[i][j] -phi[i][jp1])/dy;
	  ex[i][j]-= (0.5*dx/dy)*eps_array[i][jm1]*(phi[i][j] -phi[i][jm1])/dy;
	  ex[i][j]/= -0.5*(eps_array[im1][j]+eps_array[im1][jm1]);
	  
	  ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
	}
      }
      else if(face[i][j] == RIGHT) {
	if(conductor[i][j]==0) {
	  ex[i][j] = 2.0*sigma[i][j]/(eps_array[i][j]+eps_array[i][jm1]+SNG_MIN);
	  ey[i][j] = 0.0;
	}
	else {
	  ex[i][j] = sigma[i][j];
	  ex[i][j]+= 0.5*(eps_array[im1][j]+eps_array[im1][jm1])*(phi[im1][j] -phi[i][j])/dx;
	  ex[i][j]-= (0.5*dx/dy)*eps_array[im1][j]*(phi[i][j] -phi[i][jp1])/dy;
	  ex[i][j]-= (0.5*dx/dy)*eps_array[im1][jm1]*(phi[i][j] -phi[i][jm1])/dy;
	  ex[i][j]/= 0.5*(eps_array[i][j]+eps_array[i][jm1]);
	  
	  ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
	}
      }
      else if(face[i][j] == UP) {
	if(conductor[i][j]==0) {
	  ex[i][j] = 0.0;
	  ey[i][j] = 2.0*sigma[i][j]/(eps_array[i][j]+eps_array[im1][j]+SNG_MIN);
	}
	else {
	  ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	  
	  ey[i][j] = sigma[i][j];
	  ey[i][j]+= 0.5*(eps_array[i][jm1]+eps_array[im1][jm1])*(phi[i][jm1] -phi[i][j])/dy;
	  ey[i][j]-= (0.5*dy/dx)*eps_array[i][jm1]*(phi[i][j] -phi[ip1][j])/dx;
	  ey[i][j]-= (0.5*dy/dx)*eps_array[im1][jm1]*(phi[i][j] -phi[im1][j])/dx;
	  ey[i][j]/= 0.5*(eps_array[i][j]+eps_array[im1][j]);
	}
      }
      else if(face[i][j] == DOWN) {
	if(conductor[i][j]==0) {
	  ex[i][j] = 0.0;
	  ey[i][j] = -2.0*sigma[i][j]/(eps_array[i][jm1]+eps_array[im1][jm1]+SNG_MIN);
	}
	else {
	  ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	  
	  ey[i][j] = sigma[i][j];
	  ey[i][j]+= 0.5*(eps_array[i][j]+eps_array[im1][j])*(phi[i][jp1] -phi[i][j])/dy;
	  ey[i][j]-= (0.5*dy/dx)*eps_array[i][j]*(phi[i][j] -phi[ip1][j])/dx;
	  ey[i][j]-= (0.5*dy/dx)*eps_array[im1][j]*(phi[i][j] -phi[im1][j])/dx;
	  ey[i][j]/= -0.5*(eps_array[i][jm1]+eps_array[im1][jm1]);
	}
      }
      else if(face[i][j]==UL_CORN) {
	if(conductor[i][j]==0) {
	  ex_mag = -(phi[i][j] -phi[im1][j])/dx;
	  ey_mag = -(phi[i][jp1] -phi[i][j])/dy;
	  e_mag  = sqrt(ex_mag*ex_mag + ey_mag*ey_mag);
	  ex[i][j] = fabs(sigma[i][j])*my_div(ex_mag, e_mag)*2.0/(eps_array[im1][j]+eps_array[im1][jm1]);
	  ey[i][j] = fabs(sigma[i][j])*my_div(ey_mag, e_mag)*2.0/(eps_array[i][j]+eps_array[im1][j]);
	}
	else {
	  ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	  ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
	}
      }
      else if(face[i][j]==UR_CORN) {
	if(conductor[i][j]==0) {
	  ex_mag = -(phi[ip1][j] -phi[i][j])/dx;
	  ey_mag = -(phi[i][jp1] -phi[i][j])/dy;
	  e_mag  = sqrt(ex_mag*ex_mag + ey_mag*ey_mag);
	  ex[i][j] = fabs(sigma[i][j])*my_div(ex_mag, e_mag)*2.0/(eps_array[i][j]+eps_array[i][jm1]);
	  ey[i][j] = fabs(sigma[i][j])*my_div(ey_mag, e_mag)*2.0/(eps_array[i][j]+eps_array[im1][j]);
	}
	else {
	  ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	  ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
	}
      }
      else if(face[i][j]==LL_CORN) {
	if(conductor[i][j]==0) {
	  ex_mag = -(phi[i][j] -phi[im1][j])/dx;
	  ey_mag = -(phi[i][j] -phi[i][jm1])/dy;
	  e_mag  = sqrt(ex_mag*ex_mag + ey_mag*ey_mag);
	  ex[i][j] = fabs(sigma[i][j])*my_div(ex_mag, e_mag)*2.0/(eps_array[im1][j]+eps_array[im1][jm1]);
	  ey[i][j] = fabs(sigma[i][j])*my_div(ey_mag, e_mag)*2.0/(eps_array[i][jm1]+eps_array[im1][jm1]);
	}
	else {
	  ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	  ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
	}
      }
      else if(face[i][j]==LR_CORN) {
	if(conductor[i][j]==0) {
	  ex_mag = -(phi[ip1][j] -phi[i][j])/dx;
	  ey_mag = -(phi[i][j] -phi[i][jm1])/dy;
	  e_mag  = sqrt(ex_mag*ex_mag + ey_mag*ey_mag);
	  ex[i][j] = fabs(sigma[i][j])*my_div(ex_mag, e_mag)*2.0/(eps_array[i][j]+eps_array[i][jm1]);
	  ey[i][j] = fabs(sigma[i][j])*my_div(ey_mag, e_mag)*2.0/(eps_array[i][jm1]+eps_array[im1][jm1]);
	}
	else {
	  ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	  ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
	}
      }
      else {
	ex[i][j] = hdx*(phi[im1][j] -phi[ip1][j]);
	ey[i][j] = hdy*(phi[i][jm1] -phi[i][jp1]);
      }
    }
  }
}

/***************************************************************/

double lhs_source()
{
  return (lhs_dc + lhs_ac*sin(t*lhs_w0 + lhs_theta0));
}

/***************************************************************/

double rhs_source()
{
  return (rhs_dc + rhs_ac*sin(t*rhs_w0 + rhs_theta0));
}

/***************************************************************/

void get_lhs_phi0()
{
  register int isp, i, j;
  register double numerator;
  
  lhs_conv_chrg = 0.0;
  for(isp=0; isp<nsp; isp++) lhs_conv_chrg += lhs_chrg[isp];
  
  numerator = lhs_sigma_total;
  if(lhs_flag==CURRENT_D)
    numerator += (lhs_conv_chrg +dt*lhs_source())/del_area;
  else if(lhs_flag==VOLTAGE_D)
    numerator += (lhs_conv_chrg -lhs_extq +lhs_extc*lhs_source())/del_area;
  
  for(j=lhs_nys; j<=lhs_nyf; j++)
    numerator += phi_pois[1][j]*EPS0/dx +0.5*dx*rho[0][j];
  
  if(ylength!=lhs_length) {
    numerator -= (rho[0][lhs_nys] +rho[0][lhs_nyf])*0.5*dx*(dx -2*dy)/(dx+dy);
    numerator += lhs_b*phi_pois[1][lhs_nys] +lhs_c*phi_pois[1][lhs_nyf];
  }
  
  lhs_phi0 = numerator/lhs_bc_den;
}

/***************************************************************/

void get_rhs_phi0()
{
  register int isp, i, j;
  register double numerator;
  
  rhs_conv_chrg = 0.0;
  for(isp=0; isp<nsp; isp++) rhs_conv_chrg += rhs_chrg[isp];
  
  numerator = rhs_sigma_total;
  if(rhs_flag==CURRENT_D)
    numerator += (rhs_conv_chrg +dt*rhs_source())/del_area;
  else if(rhs_flag==VOLTAGE_D)
    numerator += (rhs_conv_chrg -rhs_extq +rhs_extc*rhs_source())/del_area;
  
  for(j=rhs_nys; j<=rhs_nyf; j++)
    numerator += phi_pois[1][j]*EPS0/dx +0.5*dx*rho[0][j];
  
  if(ylength!=rhs_length) {
    numerator -= (rho[ncx][rhs_nys] +rho[ncx][rhs_nyf])*0.5*dx*(dx -2*dy)/(dx+dy);
    numerator += rhs_b*phi_pois[ncx-1][rhs_nys] +rhs_c*phi_pois[ncx-1][rhs_nyf];
  }
  
  rhs_phi0 = numerator/rhs_bc_den;
}

/***************************************************************/

void lhs_circuit()
{
  register int j;
  
  lhs_extq_1 = lhs_extq;
  lhs_sigma_total_1 = lhs_sigma_total;
  if(lhs_flag==VOLTAGE_D) {
    lhs_sigma_total = 0.0;
    for(j=lhs_nys; j<=lhs_nyf; j++) lhs_sigma_total += sigma[0][j];
    lhs_extq = lhs_extq_1 +del_area*(lhs_sigma_total -lhs_sigma_total_1) -lhs_conv_chrg;
    lhs_exti = (lhs_extq -lhs_extq_1)/dt;
  }
  else if(lhs_flag==CURRENT_D) {
    lhs_exti = lhs_source();
    lhs_sigma_total += (lhs_exti*dt +lhs_conv_chrg)/del_area;
  }
}

/***************************************************************/

void rhs_circuit()
{
  register int j;
  
  rhs_extq_1 = rhs_extq;
  rhs_sigma_total_1 = rhs_sigma_total;
  
  if(rhs_flag==VOLTAGE_D) {
    rhs_sigma_total = 0.0;
    for(j=rhs_nys; j<=rhs_nyf; j++) rhs_sigma_total += sigma[ncx][j];
    rhs_extq = rhs_extq_1 +del_area*(rhs_sigma_total -rhs_sigma_total_1) -rhs_conv_chrg;
    rhs_exti = (rhs_extq -rhs_extq_1)/dt;
  }
  else if(rhs_flag==CURRENT_D) {
    rhs_exti = rhs_source();
    rhs_sigma_total += (rhs_exti*dt +rhs_conv_chrg)/del_area;
  }
}

/***************************************************************/
/* Initialize the the boundaries profile                       */

void init_fields()
{
  register int i, j;
  double *phi_lhs, *phi_rhs, *phi_dwn, *phi_up;

  hdx = 0.5/dx;
  hdy = 0.5/dy;
  
  //if(ybc_flag == NEUMANN)
    init_dadi_arrays(NEUMANN, NEUMANN, NEUMANN, NEUMANN);
    //init_dadi_arrays(DIRICHLET, DIRICHLET, NEUMANN, NEUMANN);
  //else if(ybc_flag == DIRICHLET)
    //init_dadi_arrays(DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET);
  
  phi_lhs = (double *)malloc(ngy*sizeof(double));
  phi_rhs = (double *)malloc(ngy*sizeof(double));
  phi_dwn = (double *)malloc(ngx*sizeof(double));
  phi_up  = (double *)malloc(ngx*sizeof(double));
  
  for (i=0; i< ngx; i++)
    for (j=0; j< ngy; j++) source[i][j] = 0.0;
  
  /********************************************************/
  /*** Set boundary conditions at y=0 and y= ly to zero ***/
  for (i=0; i<= ncx; i++) phi_dwn[i] = phi_up[i]= 0.0;
  
  /***********  Fixing the face flags for top    *******/
  if(ybc_flag == DIRICHLET)
    for (i=1; i< ncx; i++) face[i][0]= UP;
  
  /***********  Fixing the face flags for bottom *******/
  if(ybc_flag == DIRICHLET)
    for (i=1; i< ncx; i++) face[i][ncy]= DOWN;
  
  /***************************************************************/
  /*** Set boundary condition at x=lhs; solve Laplace's eqn ******/
  if(lhs_flag==VOLTAGE_D || lhs_flag==CURRENT_D) {    
    for (j=0; j<= ncy; j++) {
      if(j>=lhs_nys-lhs_gap && j<lhs_nys)      phi_lhs[j]= (j-lhs_nys+lhs_gap)/(double)(lhs_gap);
      else if(j>=lhs_nys && j<=lhs_nyf)        phi_lhs[j]= 1.0;
      else if(j>lhs_nyf && j<lhs_nyf+lhs_gap)  phi_lhs[j]= (lhs_nyf+lhs_gap-j)/(double)(lhs_gap);
      else                                     phi_lhs[j]= 0.0;
    }
    for (j=0; j<= ncy; j++) phi_rhs[j]= 0.0; 
    
    dadi_dt= dadi(dadi_dt, phi_lap_lhs, source, 60, 1e-3, 1, phi_lhs, phi_rhs, phi_dwn, phi_up);

    lhs_b = lhs_c = EPS0*(dy-dx)/(dy+dx)/dx;

    lhs_bc_den  = lhs_ny+1;
    if(ylength!=lhs_length) lhs_bc_den += 2*(dy-dx)/(dy+dx) +4*dx*dx/lhs_gap/dy/(dx+dy);
    lhs_bc_den *= EPS0/dx;

    if(lhs_flag==VOLTAGE_D) lhs_bc_den += lhs_extc/del_area;
    
    if(ylength!=lhs_length) {
      lhs_bc_den -= lhs_b*phi_lap_lhs[1][lhs_nys];
      lhs_bc_den -= lhs_c*phi_lap_lhs[1][lhs_nyf];
    }    
    for(j=lhs_nys; j<=lhs_nyf; j++)
      lhs_bc_den -= phi_lap_lhs[1][j]*EPS0/dx;
    
    /***********  Fixing the face flags  *******/
    for(j=1; j< ncy; j++) {
      if (j<lhs_nys-lhs_gap)                   face[0][j]= RIGHT;
      else if (j==lhs_nys-lhs_gap)             face[0][j]= UR_CORN;
      else if (j<lhs_nys && lhs_nys-lhs_gap<j) face[0][j]= NO_FACE;
      else if (j==lhs_nys)                     face[0][j]= LR_CORN;
      else if (lhs_nys<j && j<lhs_nyf)         face[0][j]= RIGHT;
      else if (j==lhs_nyf)                     face[0][j]= UR_CORN;
      else if (lhs_nyf<j && j<lhs_nyf+lhs_gap) face[0][j]= NO_FACE;
      else if (j==lhs_nyf+lhs_gap)             face[0][j]= LR_CORN;
      else                                     face[0][j]= RIGHT;
    }
  }
  
  /***************************************************************/
  /*** Set boundary condition at x=rhs; solve Laplace's eqn ******/
  if(rhs_flag==VOLTAGE_D || rhs_flag==CURRENT_D) {    
    for (j=0; j<= ncy; j++) {
      if(j>=rhs_nys-rhs_gap && j<rhs_nys)      phi_rhs[j]= (j-rhs_nys+rhs_gap)/(double)(rhs_gap);
      else if(j>=rhs_nys && j<=rhs_nyf)        phi_rhs[j]= 1.0;
      else if(j>rhs_nyf && j<rhs_nyf+rhs_gap)  phi_rhs[j]= (rhs_nyf+rhs_gap-j)/(double)(rhs_gap);
      else                                     phi_rhs[j]= 0.0;
    }
    for (j=0; j<= ncy; j++) phi_lhs[j]= 0.0;
    
    dadi_dt= dadi(dadi_dt, phi_lap_rhs, source, 60, 1e-3, 1, phi_lhs, phi_rhs, phi_dwn, phi_up);
    
    rhs_b = rhs_c = EPS0*(dy-dx)/(dy+dx)/dx;

    rhs_bc_den  = rhs_ny+1;
    if(ylength!=rhs_length) rhs_bc_den += 2*(dy-dx)/(dy+dx) +4*dx*dx/rhs_gap/dy/(dx+dy);
    rhs_bc_den *= EPS0/dx;

    if(rhs_flag==VOLTAGE_D) rhs_bc_den += rhs_extc/del_area;
    
    if(ylength!=rhs_length) {
      rhs_bc_den -= rhs_b*phi_lap_rhs[ncx-1][rhs_nys];
      rhs_bc_den -= rhs_c*phi_lap_rhs[ncx-1][rhs_nyf];
    }    
    for(j=rhs_nys; j<=rhs_nyf; j++)
      rhs_bc_den -= phi_lap_rhs[ncx-1][j]*EPS0/dx;
    
    /***********  Fixing the face flags  *******/
    for(j=1; j< ncy; j++) {
      if (j<rhs_nys-rhs_gap)                   face[ncx][j]= LEFT;
      else if (j==rhs_nys-rhs_gap)             face[ncx][j]= UL_CORN;
      else if (j<rhs_nys && rhs_nys-rhs_gap<j) face[ncx][j]= NO_FACE;
      else if (j==rhs_nys)                     face[ncx][j]= LL_CORN;
      else if (rhs_nys<j && j<rhs_nyf)         face[ncx][j]= LEFT;
      else if (j==rhs_nyf)                     face[ncx][j]= UL_CORN;
      else if(rhs_nyf<j && j<rhs_nyf+rhs_gap)  face[ncx][j]= NO_FACE;
      else if(j==rhs_nyf+rhs_gap)              face[ncx][j]= LL_CORN;
      else                                     face[ncx][j]= LEFT;
    }
  }
  
  /***************************************************************/
  free(phi_lhs);  free(phi_rhs);  free(phi_dwn);  free(phi_up);
}
}
/***************************************************************/

