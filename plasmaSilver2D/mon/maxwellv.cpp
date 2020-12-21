#include "mcm.h"

#define  NVTS    3.3
#define  NVVEL   50000
namespace monte{
double F(double v), Fv(double v, double v0, double vth, double r);
double *vinject[NSMAX][SIDES];

/**************************************************************/

void maxwellv(double *vx, double *vy, double *vz, double vth,
          double Rv, double Rphi, double Rthe)
{
  static int nvel, init_flag= 1;
  static double *vsave;
  
  int i, n;
  double aphi, dv, rr, sintheta, costheta;
  
  if (init_flag) {
    nvel= 1./(1-F((double)NVTS));
    dv = sqrt(M_PI)/(4.0*nvel);
    if(nvel > 100000)
      puts("Warning: Your choice of NVTS has made nvel > 1e5"); 
    
    vsave= (double *) malloc(nvel*sizeof(double));
    
    init_flag= 0;
    i=n=0;
    for (n=0; n<nvel; n++) {
      rr=(1.0*n)/nvel;
      while (F(i*dv)< rr) i++;
      vsave[n]=sqrt(2.0)*i*dv;
    }
  }
  n = (nvel-1)*Rv;
  aphi=2*M_PI*Rphi;
  costheta = 1-2*Rthe;
  sintheta = sqrt(1-costheta*costheta);
  *vx = vth*vsave[n]*sintheta*cos(aphi);
  *vy = vth*vsave[n]*sintheta*sin(aphi);
  *vz = vth*vsave[n]*costheta;
}

/**************************************************************/

double F(double v)
{
  return(-2*v*exp(-v*v)/sqrt(M_PI) +erf(v));
}

/**************************************************************/

void init_vmaxwellv(int isp, int side)
{
  int n;
  double df, rr, vdrift, vth, dv, vn;
  
  vdrift= v0x[isp][side];
  vth   = vt[isp][side];
  
  if(vth >1e-5*vdrift) {
    vinject[isp][side] = (double *)malloc(NVVEL*sizeof(double));
    
    if(vdrift>4*vth) vn = vdrift -3.0*vth;
    else             vn = (vth+vdrift)/NVVEL;
    
    df = 1./(2*NVVEL);
    for (n=0; n<NVVEL; n++) {
      rr = df*(1+2*n);
      dv = Fv(vn, vdrift, vth, rr);
      while (fabs(dv) > 1e-5*(vth+vdrift)) {
	if(fabs(dv)>0.5*fabs(vn)) dv *= fabs(vn/dv);
	vn -= dv;
	dv = Fv(vn, vdrift, vth, rr);
      }
      vinject[isp][side][n]= vn;
    }
  }
}

/**************************************************************/

void vmaxwellv(int isp, int side, double *vx, double *vy, double *vz)
{
  double del, vmag, aphi;
  int i;

  if(vt[isp][side] >1e-5*v0x[isp][side]) {
    i = del = frand()*(NVVEL-1);  del -= i;
    *vx = del*vinject[isp][side][i+1] +(1-del)*vinject[isp][side][i];
  }
  else
    *vx = v0x[isp][side];
  
  aphi= 2.0*M_PI*frand();
  vmag= sqrt(2.0*fabs(log(frand())));
  *vy = v0y[isp][side] +vt[isp][side]*vmag*cos(aphi);
  *vz = v0z[isp][side] +vt[isp][side]*vmag*sin(aphi);
}

/**************************************************************/

double Fv(double v, double v0, double vth, double r)
{
  double temp1, temp2, temp3;

  temp1  = (v-v0)/vth/sqrt(2.0);
  temp2  = v0/vth/sqrt(2.0);
  temp3  = (erf(temp1) +(1-r)*erf(temp2) -r)*exp(temp1*temp1)*sqrt(0.5*M_PI)*vth*v0/v
    -vth*vth/v +(1-r)*exp(temp1*temp1-temp2*temp2)*vth*vth/v;
  
  return(temp3);
}
}
/**************************************************************/
