#include "mcm.h"
#include <string.h>

namespace monte{
char  name[NSMAX][50], strc_name[NSTRCTYPEMAX][50];

double nc2p, fncx, fncy, xlength, ylength, zlength, rhoback, df,
epsilon, dx, dy, dt, gtemp, bmag, btheta, bphi, pressure,
lhs_conv_chrg, lhs_exti, lhs_phi0, lhs_fnys, lhs_fnyf,
rhs_conv_chrg, rhs_exti, rhs_phi0, rhs_fnys, rhs_fnyf,
lhs_sigma_total, lhs_sigma_total_1, lhs_extq, lhs_extq_1,
rhs_sigma_total, rhs_sigma_total_1, rhs_extq, rhs_extq_1,
lhs_length, lhs_dc, lhs_ac, lhs_w0, lhs_theta0, lhs_extc,
rhs_length, rhs_dc, rhs_ac, rhs_w0, rhs_theta0, rhs_extc,
tol_pois, del_area;

double vt[NSMAX][SIDES], v0x[NSMAX][SIDES], v0y[NSMAX][SIDES], v0z[NSMAX][SIDES],
enter[NSMAX][SIDES], q[NSMAX], m[NSMAX], qm[NSMAX], q_per_cell[NSMAX],
Exscale[NSMAX], Eyscale[NSMAX], Ezscale[NSMAX], dtdx[NSMAX], dtdy[NSMAX],
ttx[NSMAX], tty[NSMAX], ttz[NSMAX], ssx[NSMAX], ssy[NSMAX], ssz[NSMAX],
vxscale[NSMAX], vyscale[NSMAX], vzscale[NSMAX], ax_scale[NSMAX], ay_scale[NSMAX],
grid_sec[NSMAX], rhs_sec[NSMAX], lhs_sec[NSMAX], wall_sec[NSMAX],
lhs_chrg[NSMAX], rhs_chrg[NSMAX], weight[NSMAX];

int   nsp, ncx, ngx, ncy, ngy, sflag, interval, hist_hi, nfft, freq_hi,
thist_hi, ecollisional, icollisional, ionsp,
wall_flag, wall_read_flag, lhs_flag, lhs_nys, lhs_ny, lhs_nyf, lhs_gap,
rhs_flag, rhs_nys, rhs_ny, rhs_nyf, rhs_gap,
sp_k[NSMAX], np[NSMAX], inject[NSMAX], maxnp[NSMAX],
k_count[NSMAX], it[NSMAX], sec_flag[NSMAX], sec_sp[NSMAX];

double **x, **y, **vx, **vy, **vz, *x_array, *y_array, **source, **eps_array, **area,
**phi, **phi_lap_lhs, **phi_lap_rhs, **phi_pois, **phi_ave, **phi_ave_show,
***sp_n, ***sp_n_0, ***sp_n_k, ***sp_n_mcc, ***sp_n_ave, ***sp_n_ave_show,
**rho, **ex, **ey, **ax, **ay, ***sp_ex, ***sp_ey, **sigma, ***sp_sigma,
***sp_vx0, ***sp_vy0, ***sp_vz0, ***sp_vt, **phi_intl,**Py_,**Py0_,**RHS_p,*extra_e;

double **sp_n_sm;
int saveTime=0;
int fileNumber=0;
pzSolver* pz_solver;

int   nbin[SIDES][NSMAX];
double emin[SIDES][NSMAX], de[SIDES][NSMAX];
double ****fe, ***e_array, ****fth, ***th_array, ***wall_q,
***th_ave1, ***th_ave2, ***th_ave_show, ***e_ave_show,
***fetot_center, ***fetot_edge, ***fthtot_center, ***fthtot_edge;

double **elasrate, **extrate, **ionrate, **chrgxrate;

double *t_array, *Local_t_array, *f_array, **np_hist, **kes_hist,
*phi_lhs_hist, *phi_lhs_fft, *com_phi_lhs_hist, *com_phs_lhs_hist, *com_phc_lhs_hist,
*phi_rhs_hist, *phi_rhs_fft, *com_phi_rhs_hist, *com_phs_rhs_hist, *com_phc_rhs_hist,
*cur_lhs_hist, *cur_lhs_fft, *com_cur_lhs_hist,
*cur_rhs_hist, *cur_rhs_fft, *com_cur_rhs_hist,
*pow_lhs_hist, *pow_lhs_fft, *com_pow_lhs_hist,
*pow_rhs_hist, *pow_rhs_fft, *com_pow_rhs_hist,
*phi_mid_hist, *phi_mid_fft, *com_phi_mid_hist;

int   stc_np[NSTRCTYPEMAX][NSMAX], nstrctype, nstrcmax, **cell_mask, **grid_mask,
**cell_k, **face, **conductor;

double ***stc_np_hist;

IntStructType *strc;
double t;

long   int seed;

int    xbc_flag, ybc_flag, khi, khi_plus1;
double  **rhok, **phik, **show_rhok, **show_phik, **phi_mid_k_t, **phi_mid_k_f,
*k_array, *dk, lhs_extl, lhs_extr, lhs_extq_2, lhs_extq_3;
void   (*fields_ptr)();
void   (*move_ptr)(int isp);

}

namespace monte{
void solve_PY();
void solve_mcm();
FILE *InputDeck;

/***********************************************************/
/***************************************************************/

void species(int isp)
{
    char aachar[80];
    int loader, fill_region;
    double xleft, xright, ylow, yhigh, initn, v0xi, v0yi, v0zi, vti, nvl, nvr;

    weight[isp]=1.0;
    /* Read in SPECIES parameters... */
    fscanf(InputDeck, "%s", name[isp]);

    while (fscanf(InputDeck, "%lf %lf %d %d %lf",
                  &q[isp], &m[isp], &sp_k[isp], &maxnp[isp], &weight[isp]) <4)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %d %lf %lf %lf %lf", &loader,
                  &fill_region, &xleft, &xright, &ylow, &yhigh) <6)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%lf %lf %lf %lf %lf", &initn, &v0xi, &v0yi, &v0zi, &vti) <5)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%lf %lf %lf %lf %lf", &nvl, &v0x[isp][LEFT],
                  &v0y[isp][LEFT], &v0z[isp][LEFT], &vt[isp][LEFT]) <5)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%lf %lf %lf %lf %lf", &nvr, &v0x[isp][RIGHT],
                  &v0y[isp][RIGHT], &v0z[isp][RIGHT], &vt[isp][RIGHT]) <5)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %lf %lf %d %lf %lf",
                  &nbin[LEFT][isp], &emin[LEFT][isp], &de[LEFT][isp],
                  &nbin[RIGHT][isp],&emin[RIGHT][isp],&de[RIGHT][isp]) <6)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %lf %lf %lf %lf", &sec_flag[isp],
                  &wall_sec[isp], &lhs_sec[isp], &rhs_sec[isp], &grid_sec[isp]) <5)
        fscanf(InputDeck, "%s", aachar);

    /*****************************************************/
    /* The drift and temperatures in the input file (given
     in eVs) must be converted to velocities (m/s) */

    v0xi = sqrt(2.0*1.602e-19*v0xi/m[isp]);
    v0yi = sqrt(2.0*1.602e-19*v0yi/m[isp]);
    v0zi = sqrt(2.0*1.602e-19*v0zi/m[isp]);
    vti  = sqrt(1.602e-19*vti/m[isp]);

    v0x[isp][LEFT] = sqrt(2.0*1.602e-19*v0x[isp][LEFT]/m[isp]);
    v0y[isp][LEFT] = sqrt(2.0*1.602e-19*v0y[isp][LEFT]/m[isp]);
    v0z[isp][LEFT] = sqrt(2.0*1.602e-19*v0z[isp][LEFT]/m[isp]);
    vt[isp][LEFT]  = sqrt(1.602e-19*vt[isp][LEFT]/m[isp]);

    v0x[isp][RIGHT] = sqrt(2.0*1.602e-19*v0x[isp][RIGHT]/m[isp]);
    v0y[isp][RIGHT] = sqrt(2.0*1.602e-19*v0y[isp][RIGHT]/m[isp]);
    v0z[isp][RIGHT] = sqrt(2.0*1.602e-19*v0z[isp][RIGHT]/m[isp]);
    vt[isp][RIGHT]  = sqrt(1.602e-19*vt[isp][RIGHT]/m[isp]);

    /* Normalization and other species related variables */
    if(sec_flag[isp]) sec_sp[isp] = sec_flag[isp]-1;

    qm[isp] = q[isp]/m[isp];
    q[isp] *= nc2p;
    q_per_cell[isp] = q[isp]/(zlength*dx*dy);
    dtdx[isp] = sp_k[isp]*dt/dx;
    dtdy[isp] = sp_k[isp]*dt/dy;

    if(bmag > 0.0) {
        vxscale[isp]= 1.0;
        vyscale[isp]= 1.0;
        vzscale[isp]= 1.0;
        ax_scale[isp] = 0.5*qm[isp]*dt;
        ay_scale[isp] = 0.5*qm[isp]*dt;
    }
    else {
        vxscale[isp]= dx/(sp_k[isp]*dt);
        vyscale[isp]= dy/(sp_k[isp]*dt);
        vzscale[isp]= zlength/(sp_k[isp]*dt);
        ax_scale[isp] = qm[isp]*sp_k[isp]*dt*dt/dx;
        ay_scale[isp] = qm[isp]*sp_k[isp]*dt*dt/dy;
    }

    /* Allocate particle arrays */
    if (isp==0)
    {extra_e = (double *)malloc(maxnp[isp]*sizeof(double));
        for (int i=0;i<maxnp[isp];i++)
        {extra_e[i]=0.0;}
    }

    x[isp] = (double *)malloc(maxnp[isp]*sizeof(double));
    y[isp] = (double *)malloc(maxnp[isp]*sizeof(double));
    vx[isp]= (double *)malloc(maxnp[isp]*sizeof(double));
    vy[isp]= (double *)malloc(maxnp[isp]*sizeof(double));
    vz[isp]= (double *)malloc(maxnp[isp]*sizeof(double));

    /* Load particles if the density is non-zero */
    xleft = max(0.0, xleft);
    ylow  = max(0.0, ylow);
    xright= min(xright, xlength);
    yhigh = min(yhigh, ylength);
    if(initn>0.0)
        load(isp, initn, loader, fill_region, xleft, xright, ylow, yhigh, v0xi, v0yi, v0zi, vti);

    /* Calculate the number of particles to be injected from each side */
    enter[isp][LEFT]  = nvl*dt*sp_k[isp]*zlength*lhs_length/nc2p/weight[isp];
    if(enter[isp][LEFT]> 0.0) {
        init_vmaxwellv(isp, LEFT);  inject[isp] = 1;
    }

    enter[isp][RIGHT] = nvr*dt*sp_k[isp]*zlength*rhs_length/nc2p/weight[isp];
    if(enter[isp][RIGHT]> 0.0) {
        init_vmaxwellv(isp, RIGHT); inject[isp] = 1;
    }
    printf("isp=%d, enterl=%lf, enterr=%lf, inject=%d\n",
           isp, enter[isp][LEFT], enter[isp][RIGHT], inject[isp]);
    printf("v0xi = %lf, v0yi=%lf, v0zi=%lf vti=%lf\n",
           v0xi, v0yi, v0zi, vti);
    printf("v0xl = %lf, v0yl=%lf, v0zl=%lf vtl=%lf\n",
           v0x[isp][LEFT], v0y[isp][LEFT], v0z[isp][LEFT], vt[isp][LEFT]);
    printf("v0xr = %lf, v0yr=%lf, v0zr=%lf vtr=%lf\n\n",
           v0x[isp][RIGHT], v0y[isp][RIGHT], v0z[isp][RIGHT], vt[isp][RIGHT]);
}


void start()
{
    register int i, j, k, n, isp, ix, iy;
    double temp, ttmag;
    char aachar[80], c;
    pz_solver= new pzSolver();

    /***********************************************/
    /* Open files and read in "general" parameters */

    makethefile();
    InputDeck = fopen("argon.inp", "r");

    /*********************/
    /*   DEFAULTS        */

    wall_flag      = 0;
    wall_read_flag = 0;
    ecollisional   = 0;
    icollisional   = 0;
    ionsp          = 0;
    epsilon        = 1.0;
    ybc_flag       = NEUMANN;

    /* Read lines until we get to numbers */
    while (fscanf(InputDeck,"%d %d %d %lf %lf %lf %lf %lf %lf", &nsp, &ncx, &ncy,
                  &nc2p, &dt, &xlength, &ylength, &zlength, &epsilon) <9)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %lf %lf %lf %lf %lf %lf %d", &lhs_flag, &lhs_length,
                  &lhs_dc, &lhs_ac, &lhs_w0, &lhs_theta0, &lhs_extc, &lhs_gap) <8)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %lf %lf %lf %lf %lf %lf %d", &rhs_flag, &rhs_length,
                  &rhs_dc, &rhs_ac, &rhs_w0, &rhs_theta0, &rhs_extc, &rhs_gap) <8)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%lf %lf %lf %lf %d %d %lf %d", &rhoback, &bmag,
                  &btheta, &bphi, &nfft, &sflag, &tol_pois, &nstrcmax) <8)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %d %d %lf %lf", &ecollisional, &icollisional,
                  &ionsp, &pressure, &gtemp) <5)
        fscanf(InputDeck, "%s", aachar);

    while (fscanf(InputDeck, "%d %lf %d %d", &ybc_flag,&temp,&wall_flag,&wall_read_flag) <4)
        fscanf(InputDeck, "%s", aachar);
    ybc_flag       = NEUMANN;
    xbc_flag       = NEUMANN;
    /*****************************************/

    btheta     *= 2*M_PI/360.0;
    bphi       *= 2*M_PI/360.0;

    lhs_theta0 *= 2*M_PI/360.0;
    rhs_theta0 *= 2*M_PI/360.0;
    lhs_w0     *= 2*M_PI;
    rhs_w0     *= 2*M_PI;
    epsilon    *= EPS0;

    lhs_extq    = 0.0;
    rhs_extq    = 0.0;
    seed        = (long) getpid(); /* This is the seed used in the random number generator */
    /* seed=31207321; */

    ngx = ncx+1;
    ngy = ncy+1;
    fncx= ncx;
    fncy= ncy;

    dx = xlength/fncx;
    dy = ylength/fncy;

    /*****************************************/
    /*****   Fix field parameters     ********/

    fields_ptr= fields;
    del_area  = zlength*dy;

    lhs_ny = lhs_length*ncy/ylength;   /* Left hand side  */
    if(!(ncy%2) && lhs_ny%2) lhs_ny++;
    if(ncy%2 && !(lhs_ny%2)) lhs_ny++;
    lhs_fnys = lhs_nys = (ncy-lhs_ny)/2;
    lhs_fnyf = lhs_nyf = ncy-(ncy-lhs_ny)/2;
    if(lhs_gap==0 && lhs_length != ylength) lhs_gap=1;

    rhs_ny = rhs_length*ncy/ylength;   /* Right hand side */
    if(!(ncy%2) && rhs_ny%2) rhs_ny++;
    if(ncy%2 && !(rhs_ny%2)) rhs_ny++;
    rhs_fnys = rhs_nys = (ncy-rhs_ny)/2;
    rhs_fnyf = rhs_nyf = ncy-(ncy-rhs_ny)/2;
    if(rhs_gap==0 && rhs_length != ylength) rhs_gap=1;

    /******* If y is periodic  ********/
    if(ybc_flag == PERIODIC) {
        puts("Periodic in y.\n");
        fields_ptr= fields2;
        del_area  = ylength*zlength;

        rhs_flag=rhs_gap=0;       /* set rhs to ground */
        rhs_length=ylength;
        rhs_dc=rhs_ac=rhs_w0=rhs_theta0=rhs_extc=0.0;
        rhs_fnyf = rhs_nyf= ncy;
        rhs_fnys = rhs_nys= 0;

        nstrcmax=0;
        lhs_extq_2=lhs_extq_3=0.0;

        lhs_length=ylength;       /* adjust lhs */
        lhs_gap=0;
        lhs_fnyf = lhs_nyf= ncy;
        lhs_fnys = lhs_nys= 0;
        lhs_extl=lhs_extr = 0.0;

        /**  In the periodic case ncy= 2^n **/
        for(i=0; i<= 15; i++)
            if(ncy == (1<<i)) {
                i=0;  break;
            }
        if(i) {
            puts("START: ncy is not an integer power of 2");
            exit(1);
        }
    }

    /******************************************************/
    /* Check for errors in the "general" input parameters */

    move_ptr = move;
    if(bmag > 0.0) {
        move_ptr = mag_move;
    }
    /*
  else if(nstrcmax>0)
    move_ptr = move_internal;
  */
    if (ncx<2 || ncy<2) {
        puts("START: dont be cute with ncx,ncy<2\n");
        exit(1);
    }
    if (lhs_length>ylength || rhs_length>ylength) {
        puts("START: dont be cute with electrode length > y length\n");
        exit(1);
    }
    /**end djc**/
    if (dt <= 0.0) {
        puts("START: dt <= 0!!\n");
        exit(1);
    }

    /* Check to see if nfft is an integer power of 2 */
    if(nfft) {
        for(i=0; i<= 15; i++)
            if(nfft == (1<<i)) {
                i=0;  break;
            }
        if(i) {
            puts("START: nfft is not an integer power of 2");
            exit(1);
        }
    }

    interval=1;

    /***************************************/
    /* Allocate space for field parameters */

    x_array= (double *)malloc(ngx*sizeof(double));
    for(i=0; i<=ncx; i++) x_array[i]= i*dx;

    y_array= (double *)malloc(ngy*sizeof(double));
    for(i=0; i<=ncy; i++) y_array[i]= i*dy;

    /***********  Species Stuff  ********/

    sp_n    = (double ***)malloc(nsp*sizeof(double **));
    sp_n_0  = (double ***)malloc(nsp*sizeof(double **));
    sp_n_k  = (double ***)malloc(nsp*sizeof(double **));
    sp_n_mcc= (double ***)malloc(nsp*sizeof(double **));
    sp_vx0  = (double ***)malloc(nsp*sizeof(double **));
    sp_vy0  = (double ***)malloc(nsp*sizeof(double **));
    sp_vz0  = (double ***)malloc(nsp*sizeof(double **));
    sp_vt   = (double ***)malloc(nsp*sizeof(double **));
    sp_ex   = (double ***)malloc(nsp*sizeof(double **));
    sp_ey   = (double ***)malloc(nsp*sizeof(double **));
    sp_sigma= (double ***)malloc(nsp*sizeof(double **));

    for(isp=0; isp<nsp; isp++) {
        sp_n[isp]    = (double **)malloc(ngx*sizeof(double *));
        sp_n_0[isp]  = (double **)malloc(ngx*sizeof(double *));
        sp_n_k[isp]  = (double **)malloc(ngx*sizeof(double *));
        sp_n_mcc[isp]= (double **)malloc(ngx*sizeof(double *));
        sp_vx0[isp]  = (double **)malloc(ngx*sizeof(double *));
        sp_vy0[isp]  = (double **)malloc(ngx*sizeof(double *));
        sp_vz0[isp]  = (double **)malloc(ngx*sizeof(double *));
        sp_vt[isp]   = (double **)malloc(ngx*sizeof(double *));
        sp_ex[isp]   = (double **)malloc(ngx*sizeof(double *));
        sp_ey[isp]   = (double **)malloc(ngx*sizeof(double *));
        sp_sigma[isp]= (double **)malloc(ngx*sizeof(double *));
        for(i=0; i<=ncx; i++) {
            sp_n[isp][i]    = (double *)malloc(ngy*sizeof(double));
            sp_n_0[isp][i]  = (double *)malloc(ngy*sizeof(double));
            sp_n_k[isp][i]  = (double *)malloc(ngy*sizeof(double));
            sp_n_mcc[isp][i]= (double *)malloc(ngy*sizeof(double));
            sp_vx0[isp][i]  = (double *)malloc(ngy*sizeof(double));
            sp_vy0[isp][i]  = (double *)malloc(ngy*sizeof(double));
            sp_vz0[isp][i]  = (double *)malloc(ngy*sizeof(double));
            sp_vt[isp][i]   = (double *)malloc(ngy*sizeof(double));
            sp_ex[isp][i]   = (double *)malloc(ngy*sizeof(double));
            sp_ey[isp][i]   = (double *)malloc(ngy*sizeof(double));
            sp_sigma[isp][i]= (double *)malloc(ngy*sizeof(double));
            for(j=0; j<=ncy; j++) {
                sp_n[isp][i][j]    = sp_n_0[isp][i][j]= sp_n_k[isp][i][j]= 0.0;
                sp_n_mcc[isp][i][j]= sp_vx0[isp][i][j]= sp_vy0[isp][i][j]= 0.0;
                sp_vz0[isp][i][j]= sp_vt[isp][i][j]= 0.0;
                sp_ex[isp][i][j] = sp_ey[isp][i][j] = sp_sigma[isp][i][j]= 0.0;
            }
        }
    }


    sp_n_sm = (double **)malloc(ngx*sizeof(double *));
    for(i=0; i<=ncx; i++) {
        sp_n_sm[i] = (double *)malloc(ngy*sizeof(double));
        for(j=0; j<=ncy; j++) {
            sp_n_sm[i][j] = sp_n[1][i][j];
        }
    }

    /***********  More field parameters  ********/

    rho     = (double **)malloc(ngx*sizeof(double *));
    source  = (double **)malloc(ngx*sizeof(double *));
    phi     = (double **)malloc(ngx*sizeof(double *));
    phi_pois= (double **)malloc(ngx*sizeof(double *));
    phi_intl= (double **)malloc(ngx*sizeof(double *));
    ex      = (double **)malloc(ngx*sizeof(double *));
    ey      = (double **)malloc(ngx*sizeof(double *));
    ax      = (double **)malloc(ngx*sizeof(double *));
    ay      = (double **)malloc(ngx*sizeof(double *));
    sigma   = (double **)malloc(ngx*sizeof(double *));
    eps_array=(double **)malloc(ngx*sizeof(double *));
    Py_=(double **)malloc(ngx*sizeof(double *));
    Py0_=(double **)malloc(ngx*sizeof(double *));
    RHS_p=(double **)malloc(ngx*sizeof(double *));


    for(i=0; i<=ncx; i++) {
        rho[i]     = (double *)malloc(ngy*sizeof(double));
        source[i]  = (double *)malloc(ngy*sizeof(double));
        phi[i]     = (double *)malloc(ngy*sizeof(double));
        phi_pois[i]= (double *)malloc(ngy*sizeof(double));
        phi_intl[i]= (double *)malloc(ngy*sizeof(double));
        ex[i]      = (double *)malloc(ngy*sizeof(double));
        ey[i]      = (double *)malloc(ngy*sizeof(double));
        ax[i]      = (double *)malloc(ngy*sizeof(double));
        ay[i]      = (double *)malloc(ngy*sizeof(double));
        sigma[i]   = (double *)malloc(ngy*sizeof(double));
        eps_array[i]=(double *)malloc(ngy*sizeof(double));
        Py_[i]=(double *)malloc(ngy*sizeof(double));
        Py0_[i]=(double *)malloc(ngy*sizeof(double));
        RHS_p[i]=(double *)malloc(ngy*sizeof(double));

        for(j=0; j<=ncy; j++) {
            phi[i][j] = phi_pois[i][j]= phi_intl[i][j]= sigma[i][j] = 0.0;
            eps_array[i][j]= 1.0;
            Py_[i][j]=Py0_[i][j]=0.0;
            RHS_p[i][j]=0.0;
        }
    }

    /**djc**new**/
    if(ybc_flag == PERIODIC) {
        khi=ncy/2;
        khi_plus1=ncy/2 +1;
        rhok      = (double **)malloc(ngx*sizeof(double *));
        phik      = (double **)malloc(ngx*sizeof(double *));
        show_phik = (double **)malloc(ngx*sizeof(double *));
        show_rhok = (double **)malloc(ngx*sizeof(double *));
        for(i=0; i<=ncx; i++) {
            rhok[i]      = (double *)malloc(ngy*sizeof(double));
            phik[i]      = (double *)malloc(ngy*sizeof(double));
            show_phik[i] = (double *)malloc(khi_plus1*sizeof(double));
            show_rhok[i] = (double *)malloc(khi_plus1*sizeof(double));
            for(j=0;j<khi_plus1;j++) show_rhok[i][j] = show_phik[i][j] = 0.0;
        }

        dk = (double  *)malloc(ngy*sizeof(double));
        for (k=0; k< khi; k++)
            dk[2*k] = dk[2*k+1] = 1.0 +2.0*(dx/dy*sin(M_PI*k/ncy))*(dx/dy*sin(M_PI*k/ncy));
        dk[1]=1.0 +2.0*(dx/dy)*(dx/dy);

        k_array = (double  *)malloc(khi_plus1*sizeof(double));
        for (k=0; k<=khi; k++) k_array[k]= (double) k;

        phi_mid_k_t = (double **)malloc(khi_plus1*sizeof(double *));
        phi_mid_k_f = (double **)malloc(khi_plus1*sizeof(double *));
        for(i=0;i<khi_plus1;i++) {
            phi_mid_k_t[i] = (double *)malloc(2*nfft*sizeof(double));
            phi_mid_k_f[i] = (double *)malloc(2*nfft*sizeof(double));
            for(j=0;j<2*nfft;j++) phi_mid_k_f[i][j] = 0.0;
        }
    }
    /**end djc**/

    if(lhs_flag) {
        phi_lap_lhs= (double **)malloc(ngx*sizeof(double *));
        for(i=0; i<=ncx; i++) {
            phi_lap_lhs[i]= (double *)malloc(ngy*sizeof(double));
            for(j=0; j<=ncy; j++) phi_lap_lhs[i][j]= 0.0;
        }
    }
    if(rhs_flag) {
        phi_lap_rhs= (double **)malloc(ngx*sizeof(double *));
        for(i=0; i<=ncx; i++) {
            phi_lap_rhs[i]= (double *)malloc(ngy*sizeof(double));
            for(j=0; j<=ncy; j++) phi_lap_rhs[i][j]= 0.0;
        }
    }

    /****************************/
    /* Allocate particle arrays */

    x = (double **)malloc(nsp*sizeof(double *));
    y = (double **)malloc(nsp*sizeof(double *));
    vx= (double **)malloc(nsp*sizeof(double *));
    vy= (double **)malloc(nsp*sizeof(double *));
    vz= (double **)malloc(nsp*sizeof(double *));

    /********************************************/
    /* Set up the parameters for each "species" */

    for (isp=0; isp<nsp; ++isp) {
        species(isp);

        Exscale[isp] = 0.5*vxscale[isp]*vxscale[isp]*m[isp]/(1.602e-19);
        Eyscale[isp] = 0.5*vyscale[isp]*vyscale[isp]*m[isp]/(1.602e-19);
        Ezscale[isp] = 0.5*vzscale[isp]*vzscale[isp]*m[isp]/(1.602e-19);

        if(bmag > 0.0) {
            ttmag=tan(0.5*sp_k[isp]*dt*qm[isp]*bmag);
            ttx[isp]=ttmag*sin(btheta)*cos(bphi);
            tty[isp]=ttmag*sin(btheta)*sin(bphi);
            ttz[isp]=ttmag*cos(btheta);
            ssx[isp]=2.0*ttx[isp]/(1.0+ttmag*ttmag);
            ssy[isp]=2.0*tty[isp]/(1.0+ttmag*ttmag);
            ssz[isp]=2.0*ttz[isp]/(1.0+ttmag*ttmag);
        }
    }

    /*********************************************/
    /*************  Reading structures  **********/

    nstrctype=0;
    strc = (IntStructType *)malloc(nstrcmax*sizeof(IntStructType));
    for(k=0; k<nstrcmax; k++) {
        while (fscanf(InputDeck,"%d %d %d %d %lf %lf %s",
                      &strc[k].xl, &strc[k].yl, &strc[k].xr, &strc[k].yu,
                      &strc[k].matrl, &strc[k].volt, strc[k].name) <7)
            fscanf(InputDeck, "%s", aachar);

        if(k==0){
            strcpy(strc_name[nstrctype], strc[k].name);
            strc[k].type=nstrctype;
            nstrctype=1;
        }
        else {
            for(j=0; j<nstrctype; j++) {
                if(!strcmp(strc_name[j], strc[k].name)) break;
            }
            if(j<nstrctype) strc[k].type = j;
            else {
                strcpy(strc_name[nstrctype], strc[k].name);
                strc[k].type=nstrctype;
                nstrctype++;
            }
        }
    }
    fclose(InputDeck);

    stc_np_hist = (double ***)malloc(nstrctype*sizeof(double **));
    for(j=0; j<nstrctype; j++){
        stc_np_hist[j] = (double **)malloc(nsp*sizeof(double *));
        for (i=0; i<nsp; i++) {
            stc_np_hist[j][i] = (double *)malloc(HISTMAX*sizeof(double));
        }
    }

    grid_mask= (int **)malloc(ngx*sizeof(int *));
    cell_mask= (int **)malloc(ngx*sizeof(int *));
    cell_k   = (int **)malloc(ngx*sizeof(int *));
    conductor= (int **)malloc(ngx*sizeof(int *));
    face     = (int **)malloc(ngx*sizeof(int *));
    area     = (double **)malloc(ngx*sizeof(double *));
    for(i=0; i<=ncx; i++) {
        grid_mask[i]= (int *)malloc(ngy*sizeof(int));
        cell_mask[i]= (int *)malloc(ngy*sizeof(int));
        cell_k[i]   = (int *)malloc(ngy*sizeof(int));
        conductor[i]= (int *)malloc(ngy*sizeof(int));
        face[i]     = (int *)malloc(ngy*sizeof(int));
        area[i]     = (double *)malloc(ngy*sizeof(double));
        for(j=0; j<=ncy; j++) {
            grid_mask[i][j]= 0;          /* No structure */
            cell_mask[i][j]= 0;          /* No structure */
            cell_k[i][j]   = 0;          /* No structure */
            conductor[i][j]= 1;          /* No conductor */
            face[i][j]     = NO_FACE;    /* No face      */
            area[i][j]     = dy*zlength; /* surface area */
        }
    }

    /****************************************/
    /* Setting cell flags due to structures */

    for(i=0; i<= ncx; i++)
        for(j=0; j<= ncy; j++)
            for(k=0; k<nstrcmax; k++) {
                if(j>= strc[k].yl && j < strc[k].yu && i>= strc[k].xl && i < strc[k].xr) {
                    cell_mask[i][j]= 1;
                    cell_k[i][j]   = k;
                    eps_array[i][j]= strc[k].matrl;
                }
            }

    for(i=0; i< ncx; i++)
        for(j=0; j< ncy; j++) {
            if(cell_mask[i][j]) {
                grid_mask[i][j] = grid_mask[i+1][j] = grid_mask[i][j+1] = grid_mask[i+1][j+1] = 1;
                if(strc[cell_k[i][j]].matrl == 0.0) {
                    conductor[i][j]= conductor[i+1][j]= conductor[i][j+1]= conductor[i+1][j+1]= 0;
                    phi_intl[i][j] = phi_intl[i+1][j] = phi_intl[i][j+1] = phi_intl[i+1][j+1] = strc[cell_k[i][j]].volt;
                }
            }
        }

    for(j=0; j<=ncy; j++) eps_array[ncx][j] = eps_array[ncx-1][j];
    for(i=0; i<=ncx; i++) eps_array[i][ncy] = eps_array[i][ncy-1];

    for(i=1; i< ncx; i++)
        for(j=1; j< ncy; j++) {
            if (grid_mask[i][j]) {
                if((eps_array[i-1][j]==1.0 && eps_array[i][j]==1.0 && eps_array[i-1][j-1]==1.0) ||
                        (eps_array[i-1][j]!=0.0 && eps_array[i][j]!=0.0 && eps_array[i-1][j-1]!=0.0 && eps_array[i][j-1]==0.0)) {
                    face[i][j]=UL_CORN;
                    area[i][j]= 0.5*(dx+dy)*zlength;
                }
                else if((eps_array[i-1][j]==1.0 && eps_array[i][j]==1.0 && eps_array[i][j-1]==1.0) ||
                        (eps_array[i-1][j]!=0.0 && eps_array[i][j]!=0.0 && eps_array[i][j-1]!=0.0 && eps_array[i-1][j-1]==0.0)) {
                    face[i][j]=UR_CORN;
                    area[i][j]= 0.5*(dx+dy)*zlength;
                }
                else if((eps_array[i-1][j]==1.0 && eps_array[i][j]==1.0 && eps_array[i-1][j-1]!=1.0 && eps_array[i][j-1]!=1.0) ||
                        (eps_array[i-1][j]!=0.0 && eps_array[i][j]!=0.0 && eps_array[i-1][j-1]==0.0 && eps_array[i][j-1]==0.0)) {
                    face[i][j]=UP;
                    area[i][j]= dx*zlength;
                }
                else if((eps_array[i-1][j-1]==1.0 && eps_array[i][j-1]==1.0 && eps_array[i-1][j]==1.0) ||
                        (eps_array[i-1][j]!=0.0 && eps_array[i-1][j-1]!=0.0 && eps_array[i][j-1]!=0.0 && eps_array[i][j]==0.0)) {
                    face[i][j]=LL_CORN;
                    area[i][j]= 0.5*(dx+dy)*zlength;
                }
                else if((eps_array[i-1][j-1]==1.0 && eps_array[i][j-1]==1.0 && eps_array[i][j]==1.0) ||
                        (eps_array[i][j]!=0.0 && eps_array[i-1][j-1]!=0.0 && eps_array[i][j-1]!=0.0 && eps_array[i-1][j]==0.0)) {
                    face[i][j]=LR_CORN;
                    area[i][j]= 0.5*(dx+dy)*zlength;
                }
                else if((eps_array[i-1][j-1]==1.0 && eps_array[i][j-1]==1.0 && eps_array[i][j]!=1.0 && eps_array[i-1][j]!=1.0) ||
                        (eps_array[i-1][j-1]!=0.0 && eps_array[i][j-1]!=0.0 && eps_array[i][j]==0.0 && eps_array[i-1][j]==0.0)) {
                    face[i][j]=DOWN;
                    area[i][j]= dx*zlength;
                }
                else if((eps_array[i-1][j-1]==1.0 && eps_array[i-1][j]==1.0 && eps_array[i][j]!=1.0 && eps_array[i][j-1]!=1.0) ||
                        (eps_array[i-1][j-1]!=0.0 && eps_array[i-1][j]!=0.0 && eps_array[i][j]==0.0 && eps_array[i][j-1]==0.0)) {
                    face[i][j]=LEFT;
                    area[i][j]= dy*zlength;
                }
                else if((eps_array[i-1][j-1]!=1.0 && eps_array[i-1][j]!=1.0 && eps_array[i][j]==1.0 && eps_array[i][j-1]==1.0) ||
                        (eps_array[i-1][j-1]==0.0 && eps_array[i-1][j]==0.0 && eps_array[i][j]!=0.0 && eps_array[i][j-1]!=0.0)) {
                    face[i][j]=RIGHT;
                    area[i][j]= dy*zlength;
                }
            }
        }
    /*********************************************/
    /* Fixing the face flags for left and right  */

    for(j=0; j<=ncy; j++) {
        /*conductor[0][j] = 0;    face[0][j]  = RIGHT; area[0][j]  = dy*zlength;
        conductor[ncx][j]= 0;   face[ncx][j]= LEFT;  area[ncx][j]= dy*zlength;*/
    }

    /* Fixing the face flags for top and bottom  */
    for(i=1; i<ncx; i++) {
        if(ybc_flag == DIRICHLET) {
            conductor[i][ncy]= 0; face[i][ncy]= DOWN; area[i][ncy]= dx*zlength;
            conductor[i][0]  = 0; face[i][0]  = UP;   area[i][0]  = dx*zlength;
        }
        else if(ybc_flag == NEUMANN) {
            if(grid_mask[i][ncy]){
                if((eps_array[i-1][ncy-1]!=0.0 && eps_array[i][ncy-1]==0.0) ||
                        (eps_array[i-1][ncy-1]==1.0 && eps_array[i][ncy-1]!=1.0)) {
                    face[i][ncy]=LEFT; area[i][ncy]= dy*zlength;
                }
                else if((eps_array[i-1][ncy-1]==0.0 && eps_array[i][ncy-1]!=0.0) ||
                        (eps_array[i-1][ncy-1]!=1.0 && eps_array[i][ncy-1]==1.0)) {
                    face[i][ncy]=RIGHT; area[i][ncy]= dy*zlength;
                }
            }
            if(grid_mask[i][0]){
                if((eps_array[i-1][0]!=0.0 && eps_array[i][0]==0.0) ||
                        (eps_array[i-1][0]==1.0 && eps_array[i][0]!=1.0)) {
                    face[i][0]=LEFT; area[i][0]= dy*zlength;
                }
                else if((eps_array[i-1][0]==0.0 && eps_array[i][0]!=0.0) ||
                        (eps_array[i-1][0]!=1.0 && eps_array[i][0]==1.0)) {
                    face[i][0]=RIGHT; area[i][ncy]= dy*zlength;
                }
            }
        }
    }

    /***************************************************/
    /* If particles are inside structures remove them  */

    if (nstrcmax>0) {
        for (isp=0; isp<nsp; ++isp) {
            for(n=0; n < np[isp]; n++) {
                ix= x[isp][n];
                iy= y[isp][n];
                if(cell_mask[ix][iy]) {
                    x[isp][n]  =  x[isp][np[isp]-1];
                    y[isp][n]  =  y[isp][np[isp]-1];
                    vx[isp][n] = vx[isp][np[isp]-1];
                    vy[isp][n] = vy[isp][np[isp]-1];
                    vz[isp][n] = vz[isp][np[isp]-1];
                    n--;
                    np[isp]--;
                }
            }
        }
    }

    /*************************************************/
    /* Setting up the energy and angular dist. stuff */

    if(nsp && wall_flag) {
        e_array = (double ***) malloc(SIDES*sizeof(double **));
        th_array= (double ***) malloc(SIDES*sizeof(double **));
        fe      = (double ****)malloc(SIDES*sizeof(double ***));
        fth     = (double ****)malloc(SIDES*sizeof(double ***));
        fetot_center = (double ***) malloc(SIDES*sizeof(double **));
        fetot_edge   = (double ***) malloc(SIDES*sizeof(double **));
        fthtot_center= (double ***) malloc(SIDES*sizeof(double **));
        fthtot_edge  = (double ***) malloc(SIDES*sizeof(double **));
        wall_q  = (double ***) malloc(SIDES*sizeof(double **));
        th_ave1 = (double ***) malloc(SIDES*sizeof(double **));
        th_ave2 = (double ***) malloc(SIDES*sizeof(double **));
        e_ave_show = (double ***)malloc(SIDES*sizeof(double **));
        th_ave_show= (double ***)malloc(SIDES*sizeof(double **));

        for(k=0; k<SIDES; k++) {
            e_array[k] = (double **) malloc(nsp*sizeof(double *));
            th_array[k]= (double **) malloc(nsp*sizeof(double *));
            fe[k]      = (double ***)malloc(nsp*sizeof(double **));
            fth[k]     = (double ***)malloc(nsp*sizeof(double **));
            fetot_center[k] = (double **) malloc(nsp*sizeof(double *));
            fetot_edge[k]   = (double **) malloc(nsp*sizeof(double *));
            fthtot_center[k]= (double **) malloc(nsp*sizeof(double *));
            fthtot_edge[k]  = (double **) malloc(nsp*sizeof(double *));
            wall_q[k]  = (double **) malloc(nsp*sizeof(double *));
            th_ave1[k] = (double **) malloc(nsp*sizeof(double *));
            th_ave2[k] = (double **) malloc(nsp*sizeof(double *));
            e_ave_show[k] = (double **)malloc(nsp*sizeof(double *));
            th_ave_show[k]= (double **)malloc(nsp*sizeof(double *));

            /* For each species allocate the diagnostic arrays */
            for(isp=0; isp<nsp; isp++) {

                /* Energy and theta arrays */
                nbin[k][isp]++;
                e_array[k][isp] = (double *)malloc(nbin[k][isp]*sizeof(double));
                th_array[k][isp]= (double *)malloc(nbin[k][isp]*sizeof(double));
                fetot_center[k][isp] = (double *)malloc(nbin[k][isp]*sizeof(double));
                fetot_edge[k][isp]   = (double *)malloc(nbin[k][isp]*sizeof(double));
                fthtot_center[k][isp]= (double *)malloc(nbin[k][isp]*sizeof(double));
                fthtot_edge[k][isp]  = (double *)malloc(nbin[k][isp]*sizeof(double));
                for(i=0; i <nbin[k][isp]; i++) {
                    e_array[k][isp][i] = i*de[k][isp];
                    th_array[k][isp][i]= i*(0.5*M_PI/(nbin[k][isp]-1));
                }

                /* Average angles */
                wall_q[k][isp] = (double *)malloc(ngy*sizeof(double));
                th_ave1[k][isp] = (double *)malloc(ngy*sizeof(double));
                th_ave2[k][isp] = (double *)malloc(ngy*sizeof(double));
                e_ave_show[k][isp] = (double *)malloc(ngy*sizeof(double));
                th_ave_show[k][isp]= (double *)malloc(ngy*sizeof(double));

                /* Energy and theta diagnostic arrays */
                fe[k][isp] = (double **)malloc(ngy*sizeof(double *));
                fth[k][isp]= (double **)malloc(ngy*sizeof(double *));

                for(j=0; j<ngy; j++) {
                    wall_q[k][isp][j]  = 0.0;
                    th_ave1[k][isp][j] = 0.0;
                    th_ave2[k][isp][j] = 0.0;

                    fe[k][isp][j] = (double *)malloc(nbin[k][isp]*sizeof(double));
                    fth[k][isp][j]= (double *)malloc(nbin[k][isp]*sizeof(double));
                    for(i=0; i <nbin[k][isp]; i++)
                        fe[k][isp][j][i] = fth[k][isp][j][i] = 0.0;
                }
            }
        }
    }

}

void saveInTecplot()
{
    fileNumber++;
    printf("saving file, t=%e\n",t);
    char str[80];
    sprintf (str, "output/out_%d_%e.dat",fileNumber, t);
    FILE *file = fopen(str, "w");

    fprintf(file, "TITLE =\"  \" \n VARIABLES=\"x\" \n \"y\"\n \"z\" \n \"phi\" \n \"ex\" \n \"ey\" \n \"ne\" \n \"ni\" \n  ZONE T=\"  \" \n");
    fprintf(file," I=%d J=%d K=%d F=POINT \n", ncx, ncy, 1);

    for( int j = 0; j < ncy; j++ ){
        for( int i = 0; i < ncx; i++ ){
            double x = i * dx;
            double y = j * dy;
            double z = 0;

            fprintf(file,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",
                    x, y, z, phi[i][j], ex[i][j], ey[i][j], sp_n[0][i][j], sp_n[1][i][j]);
        }
    }

    fclose(file);
}

void updateEforPz()
{
    for (int i=0;i<pz_solver->m_p_num;i++)
    {
        double x,ym, y0,yp;
        x=pz_solver->m_p[i].r.x;
        ym=pz_solver->m_p[i].r.y+pz_solver->m_p[i].dl*0.22;
        y0=pz_solver->m_p[i].r.y;
        yp=pz_solver->m_p[i].r.y+pz_solver->m_p[i].dl*0.5;
        int xidx = x/monte::dx;
        int yidx = yp/monte::dy;
        vec2 Ep0 = pz_solver->getEdepol(x,y0);
        pz_solver->m_p[i].E =(1000 * ey[ xidx][ yidx] + Ep0.y)*100;


        //printf("pppp=%e \n",pz_solver->m_p[i].E);//<<endl;
        /*vec2 Em = m_Esolver->getE(x,ym);
                vec2 Epm = m_pzSolver->getEdepol(x,ym);
                vec2 Eem = m_elecSolver->getEe(x,ym);

                vec2 E0 = m_Esolver->getE(x,y0);
                vec2 Ep0 = m_pzSolver->getEdepol(x,y0);
                vec2 Ee0 = m_elecSolver->getEe(x,y0);

                vec2 Ep = m_Esolver->getE(x,yp);
                vec2 Epp = m_pzSolver->getEdepol(x,yp);
                vec2 Eep = m_elecSolver->getEe(x,yp);
                m_pzSolver->m_p[i].E_elec = (Eem.y +Ee0.y+Eep.y)/3;//0.95 * m_pzSolver->m_p[i].E_elec + 0.05 * Ee.y;
                m_pzSolver->m_p[i].E = ((Em.y + Epm.y+E0.y + Ep0.y+ Ep.y + Epp.y)/3.0 + m_pzSolver->m_p[i].E_elec);*/
        //m_pzSolver->m_p[i].E*=3;
        /*if(i%10 == 0)
                      //printf("IIIIII=%d EEEEEE=%e\n", i, m_pzSolver->m_p[i].E);*/
    }
}
/***************************************************************/

void solve_mcm()
{
    register int i, j, n, isp;
    double frac;

    t += dt;
    saveTime++;
    if(frand()>0.95)
        printf("t=%e\n",t);
    for(isp=0; isp<nsp; isp++) {
        if(!(k_count[isp]%sp_k[isp])) {
            it[isp]++;
            (*move_ptr)(isp);
            boundary(isp);
            if(wall_flag) wall_diagnos(isp);
            mcc(isp);
            diagnos();
            gather(isp);
            k_count[isp]=0;
        }
        //solve_PY();
        k_count[isp]++;
        frac = ((double)k_count[isp])/sp_k[isp];
        for(i=0; i<=ncx; i++)
            for(j=0; j<=ncy; j++)
                sp_n[isp][i][j]= (1- frac)*sp_n_0[isp][i][j]
                        +frac*sp_n_k[isp][i][j] +sp_n_mcc[isp][i][j];

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
                /*TWOD_One_2_One(sp_n[isp],   ncx, ncy);
                TWOD_One_2_One(sp_n[isp],   ncx, ncy);
                TWOD_One_2_One(sp_n[isp],   ncx, ncy);
                TWOD_One_2_One(sp_n[isp],   ncx, ncy);*/
                TWOD_One_2_One(sp_vx0[isp], ncx, ncy);
                TWOD_One_2_One(sp_vy0[isp], ncx, ncy);
                TWOD_One_2_One(sp_vz0[isp], ncx, ncy);
                TWOD_One_2_One(sp_vt[isp],  ncx, ncy);
            }
        }
    }

    for(i=0; i<=ncx; i++)
        for(j=0; j<=ncy; j++)
            sp_n_sm[i][j]= sp_n[1][i][j];



    /**************************************************/
    /** Now do the field solve and finally, calculate */
    /** the species E field  */

    (*fields_ptr)();
    for(isp=0; isp<nsp; isp++)
        for(i=0; i <= ncx; i++)
            for(j=0; j <= ncy; j++) {
                sp_ex[isp][i][j] += ex[i][j];
                sp_ey[isp][i][j] += ey[i][j];
            }

    /*if(saveTime==100){
        saveTime=0;
        saveInTecplot();
    }*/
}



double calc_poly(poly p, double r, double x) //возвращает значение функции f(x) = x^2-2
{
    double sum=0.0;
    double xn=x;
    for (int i=0;i<p.order;++i)
    {
        sum+=xn*p.C[i];
        xn*=x;
    }

    return sum-r;
}

float calc_d_poly(poly p, double x) //возвращает значение производной
{
    double sum=p.C[0];
    double xn=x;
    for (int i=0;i<p.order-1;++i)
    {
        sum+=xn*p.C[i+1]*(i+2);
        xn*=x;
    }
    return sum;
}

float calc_d2_poly(poly p, double x) // значение второй производной
{
    double sum=2.0*p.C[1];
    double xn=x;
    for (int i=0;i<p.order-2;++i)
    {
        sum+=xn*p.C[i+1]*(i+2)*(i+3);
        xn*=x;
    }
    return sum;
}

double solve_poly(poly p,double _x0, double rhs,int itn)
{
    double x0,x,xn,eps;// вычисляемые приближения для корня
    double a, b,deltax;// границы отрезка и необходимая точность

    // deltax=deltax0;
    x=_x0;
    /* a=x0-deltax;
    b=x0+deltax;

    int nn=0;
    double sucs=f(a)*f(b);
    while ((sucs>0)&&(nn<10)) // если знаки функции на краях отрезка одинаковы
    {
        deltax=deltax*2;
        a=x0-deltax;
        b=x0+deltax;
        nn++;
        sucs=f(a)*f(b);
    }

    if (sucs>0)
    {
        printf("DIVERGENCEE in POLY \n",)
        return -1e30
    }
*/
    for (int i=0; i<itn; ++i)
    {
        x = x-calc_poly(p,rhs,x)/calc_d_poly(p,x); // считаем первое приближение
    }

    // printf("x = %lf f(x)=%e \n",x*1000.0,calc_poly(p,rhs,x));

    return x;

}


int jacobi_polynomial(INPUT_PARAM par, poly pol,double** field,double **rhs, int itn)
{
    int i,j,n;
    double res=0;

    //static poly pol_new[N_X][N_Y];

    poly poly_new;

    double a,b_p,b_m,c_p,c_m;

    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    b_p=par.bp;//-1.0/(dx*dx);
    b_m=par.bm;//-1.0/(dx*dx);
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);

    static double rhs_[600][600]; //careful here :)

    //field_x*a+....+pol*fieldx^..=rhs

    poly_new=pol;
    poly_new.C[0]+=a;


    int j0 =2; //dielectric between this
    int j1= 5; //

    for (i=1; i<ncx-1; i+=30)
    {
        printf("i=%d P=%e ey=%e \n",i,Py_[i][2],ey[i][2]);
    }

    for(n=0;n<itn;n++)
    {
        for (i=1; i<ncx-1; i++)
        {
            //field[i][0]=1;
            for (j=j0+1; j<j1-1; j++)
            {
                //p[i][j]=0;

                rhs_[i][j]=rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]);

                field[i][j]=field[i][j]*0.7+0.3*solve_poly(poly_new,field[i][j],rhs_[i][j],4);//(rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))/a;
            }
        }

        //0--fixed value, 1--fixed gradient,2 --cyclic, 3 --init
        for (int j=j0; j< j1; j++)
        {
            if (par.w_bc_type==0)//fixed_value
                field[0][j]=par.w_bc_val;
            if (par.w_bc_type==1)//fixed_gradient
                field[0][j]=field[1][j];//-dx*par.w_bc_val;
            if (par.w_bc_type==2)//cyclic
                field[0][j]=field[ncx-2][j];
            // if (par.w_bc_type==3)// init

            if (par.e_bc_type==0)//fixed_value
                field[ncx-1][j]=par.e_bc_val;
            if (par.e_bc_type==1)//fixed_gradient
                field[ncx-1][j]=field[ncx-2][j]+dx*par.e_bc_val;
            if (par.e_bc_type==2)//cyclic
                field[ncx-1][j]=field[1][j];
            // if (par.e_bc_type==3)// init
        }

        for (int i=0; i<ncx; i++ )
        {
            if (par.n_bc_type==0)//fixed_value
                field[i][j1-1]=par.n_bc_val;
            if (par.n_bc_type==1)//fixed_gradient
                field[i][j1-1]=field[i][j1-2]+dy*par.n_bc_val;
            if (par.n_bc_type==2)//cyclic
                field[i][j1-1]=field[i][j0+1];
            if (par.n_bc_type==5)//zero div
                field[i][j1-1]=field[i+1][j1-1]-dx*rhs[i][j1-1];
            //  if (par.n_bc_type==3)// init

            if (par.s_bc_type==0)//fixed_value
                field[i][j0]=par.s_bc_val;
            if (par.s_bc_type==1)//fixed_gradient
                field[i][j0]=field[i][j0+1]-dy*par.s_bc_val;
            if (par.s_bc_type==2)//cyclic
                field[i][j0]=field[i][j1-1];
            if (par.s_bc_type==5)//zero div
                field[i][j0]=field[i+1][j0]-dx*rhs[i][j0];
            // j=N_Y_DIEL-1;
            // Py_[i][j]=Py_[i][j-1] -(Px_[i+1][j]-Px_[i][j])/(1.0*dx)*(1.0*dy);
            //Py_[i][j]=0.5*(Py_[i][j+1]+Py_[i][j-1]);

            /* j=1;
            Py_[i][j-1]=Py_[i][j+1] +(Px_[i+1][j]-Px_[i-1][j])/(2.0*dx)*(2.0*dy);
            Py_[i][j]=0.5*(Py_[i][j+1]+Py_[i][j-1]);*/
        }
    }

    return 0;
}


void solve_PY()
{
    INPUT_PARAM par_ferr;
    double kap=1.38e-10;
    par_ferr.a=(1.0/dt)+(kap*2.0/(dx*dx) + kap*2.0/(dy*dy));
    par_ferr.bp=-kap/(dx*dx);
    par_ferr.bm=-kap/(dx*dx);
    par_ferr.cp=-kap/(dy*dy);
    par_ferr.cm=-kap/(dy*dy);

    par_ferr.w_bc_type=1;
    par_ferr.e_bc_type=1;
    par_ferr.n_bc_type=1;
    par_ferr.s_bc_type=1;

    par_ferr.w_bc_val=0.0;
    par_ferr.e_bc_val=0.0;
    par_ferr.n_bc_val=0.0;
    par_ferr.s_bc_val=0.0;

    poly p;

    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;

    p.order=5;
    p.C[0]=2*alp*(T-T0);//+1/EPS0;//(T-T0); //x
    p.C[1]=0.0;         //xx
    p.C[2]=-4.0*bet;   //xxx
    p.C[3]=0.0;        //x^4
    p.C[4]=6.0*gam;

    for (int i=1; i<ncx-1; i++)
    {
        for (int j=1; j<ncy-1; j++)
        {

            RHS_p[i][j]= 5000.0*ey[i][j]+Py0_[i][j]/dt;

        }
    }

    //  printf("emin=%e emax=%e \n",emin,emax);

    jacobi_polynomial( par_ferr, p,Py_,RHS_p, 5);

    for (int i=1; i<ncx-1; i++)
    {
        for (int j=1; j<ncy-1; j++)
        {
            Py0_[i][j]=Py_[i][j];
        }
    }

}




void solve_landau(double d_t) //landau dipole theory in fixed external field
{
    poly p;
    double alp,bet,gam,T,T0,rh;
    alp=3.324e5;
    bet=6.381e8;
    gam=7.89e9;
    T=300;
    T0=381;
    rh=0.0;

    double x1=0;
    double x0=0;
    //double d_t=1e-11;
    double E0=-1e5;
    double E;
    for (int i=0;i<10000;i++)
    {
        rh=E0+x0/d_t;

        p.order=5;
        p.C[0]=2*alp*(T-T0)+1.0/d_t +1/EPS0; //x
        p.C[1]=0.0;         //xx
        p.C[2]=-4.0*bet;   //xxx
        p.C[3]=0.0;        //x^4
        p.C[4]=6.0*gam;

        x1=solve_poly(p,x1, rh,10);
        x0=x1;
        if (i%10==0)
            printf("t=%e P=%e E=%e \n",i*d_t,x1,(E0-x1/EPS0));
    }

    printf("now exact:--- \n");
    rh=E0;

    p.order=5;
    p.C[0]=2*alp*(T-T0) +1/EPS0; //x
    p.C[1]=0.0;         //xx
    p.C[2]=-4.0*bet;   //xxx
    p.C[3]=0.0;        //x^4
    p.C[4]=6.0*gam;

    x1=solve_poly(p,x1, rh,100);
    printf("P_exact=%e E=%e poly=%e \n",x1,(E0-x1/EPS0),calc_poly(p,rh,x1));


}


}
