#include "mcm.h"
/************************************************************************/
/* frand() returns values 0 through 1.                                  */
/* From "Random number generators: good ones are hard to find", S. Park */
/* and K. Miller, Communications of ACM, October 1988, pp 1192-1201.    */
/* This is from page 1195, and is to work on any system for which       */
/* maxint is 2**31-1 or larger. Due earlier to Schrage, as cited by P&M.*/
/*                                                                      */
/* Note: OK values for iseed are 1 through 2**31-2. Give it 0 or 2*31-1 */
/* and it will return the same values thereafter!                       */
/*                                                                      */
/* C version 6/91, Bruce Langdon.                                       */
/*                                                                      */
/* Algorithm replaces seed by mod(a*seed,m). First represent            */
/* seed = q*hi + lo.  Then                                              */
/* a*seed = a*q*hi + lo = (m - r)*hi + a*lo = (a*lo - r*hi) + m*hi,     */
/* and new seed = a*lo - r*hi unless negative; if so, then add m.       */
namespace monte{
double frand()
{
  long a = 16807, m = 2147483647, q = 127773, r = 2836;
  long hi, lo;
  double fnumb;
  /* static long seed=31207321; */

  hi = seed/q;
  lo = seed - q*hi;
  seed = a*lo - r*hi;
  /* "seed" will always be a legal integer of 32 bits (including sign). */
  if(seed <= 0) seed = seed + m;
  fnumb = seed/2147483646.0;

  return(fnumb);
}
/***************************************************************/
double revers(int num, int n)
{
  double power, rev;
  int inum, iquot, irem;

  rev = 0.;
  inum = num;
  power = 1.;

  do {
    iquot = inum/n;
    irem = inum - n*iquot;
    power /= n;
    rev += irem*power;
    inum = iquot;
  } while (inum > 0);

  return (rev);
}

void load(int isp, double initn, int loader, int fill_region, double xleft, double xright,
      double ylow, double yhigh, double v0xi, double v0yi, double v0zi, double vti)
{
  static int i=0;
  register int n, nx, ny;
  double fnp, x0, xlen, y0, ylen;
  
  if(fill_region) {
    x0  = 0.0;
    xlen= fncx;
    y0  = 0.0;
    ylen= fncy;
  }
  else {
    x0  = fncx*xleft/xlength;
    xlen= fncx*(xright -xleft)/xlength;
    y0  = fncy*ylow/ylength;
    ylen= fncy*(yhigh -ylow)/ylength;
  }
  
  np[isp]= 0.001 * initn*xlen*ylen*dx*dy*zlength/nc2p/weight[isp];
  if (np[isp] > maxnp[isp]) {
    printf("LOAD: too many particles, species %d\n", isp);
    exit(1);
  }
  
  /******************************************************/
  /* Loading the velocities and positions depending on
     the given flag.  This Loader can load particles
     using a quiet start method, uniform XY or randomly */
  
  switch (loader) {
  case QUIET_START:
    for(n=0; n<np[isp]; n++) {
      vx[isp][n] = vy[isp][n] = vz[isp][n] = 0.0;
      if (vti>0.0)
	maxwellv(&vx[isp][n], &vy[isp][n], &vz[isp][n], vti,
		 frand(), revers(n+i,2), revers(n+i,7));
      
      vx[isp][n] =(vx[isp][n] +v0xi)/vxscale[isp];
      vy[isp][n] =(vy[isp][n] +v0yi)/vyscale[isp];
      vz[isp][n] =(vz[isp][n] +v0zi)/vzscale[isp];
      
      x[isp][n]= x0 +xlen*revers(n+i, 5);
      y[isp][n]= y0 +ylen*revers(n+i, 3);
    }
    break;
  case UNIFORM_XY:
    fnp= np[isp];
    nx = sqrt(fnp)/fncy+.5;
    ny = sqrt(fnp)/fncx+.5;
    np[isp]= nx*ny*fncx*fncy;
    
    for(n=0; n<np[isp]; n++) {
      vx[isp][n] = vy[isp][n] = vz[isp][n] = 0.0;
      if (vti>0.0)
	maxwellv(&vx[isp][n], &vy[isp][n], &vz[isp][n], vti,
		 revers(n+i,2), revers(n+i,7), revers(n+i,3));
      
      vx[isp][n] =(vx[isp][n] +v0xi)/vxscale[isp];
      vy[isp][n] =(vy[isp][n] +v0yi)/vyscale[isp];
      vz[isp][n] =(vz[isp][n] +v0zi)/vzscale[isp];
      
      x[isp][n]= (n%(ncx*nx) +.5)/nx;
      y[isp][n]= (n/(ncx*nx) +.5)/ny;
    }
    break;
  case RANDOM:
    for(n=0; n<np[isp]; n++) {
      vx[isp][n] = vy[isp][n] = vz[isp][n] = 0.0;
      if (vti>0.0)
	maxwellv(&vx[isp][n], &vy[isp][n], &vz[isp][n], vti,
		 frand(), frand(), frand());
      
      vx[isp][n] =(vx[isp][n] +v0xi)/vxscale[isp];
      vy[isp][n] =(vy[isp][n] +v0yi)/vyscale[isp];
      vz[isp][n] =(vz[isp][n] +v0zi)/vzscale[isp];
      
      x[isp][n]= x0 +xlen*frand();
      y[isp][n]= y0 +ylen*frand();
    }
    break;
  default:
    puts("LOAD: Bad value for loader flag");
    exit(-1);
    break;
  }
  
  /***************************************************/
  /* if one particle, stick it in center and give it */
  /* drift + random thermal part                     */
  if (np[isp]==1) {
    x[isp][0]=x0 +xlen/2.0;
    y[isp][0]=y0 +ylen/2.0;
    vx[isp][0]=(v0xi+2.0*(frand()-0.5)*vti)/vxscale[isp];
    vy[isp][0]=(v0yi+2.0*(frand()-0.5)*vti)/vyscale[isp];
    vz[isp][0]=(v0zi+2.0*(frand()-0.5)*vti)/vzscale[isp];
  }
  i+= np[isp];
}
}
/***************************************************************/
