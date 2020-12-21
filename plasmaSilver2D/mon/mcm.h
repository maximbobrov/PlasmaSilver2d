#ifndef MCM_H
#define MCM_H

#include <math.h>
#include <stdio.h>
#include <unistd.h>


#define double_CHAR "double"
#define XG_double_DOUBLE



#define EPS0            8.8542e-12      /* (F/m)  */
#define NperTORR        8.3221e20
#define HISTMAX         512

#define NSMAX           12

#define QUIET_START     0
#define UNIFORM_XY      1
#define RANDOM          2

#define GROUNDED        0
#define VOLTAGE_D       1
#define CURRENT_D       2

#define DIRICHLET       0
#define PERIODIC        1
#define NEUMANN         2

#define LEFT            0
#define RIGHT           1
#define UP              2
#define DOWN            3
#define UL_CORN         4
#define UR_CORN         5
#define LL_CORN         6
#define LR_CORN         7
#define NO_FACE         8

#define SIDES           2

#ifndef max
#define max(x, y)       (((x) > (y)) ? (x) : (y))
#endif

#ifndef min
#define min(x, y)       (((x) < (y)) ? (x) : (y))
#endif

#define SNG_MIN         1E-5

#ifndef True
#define True            1
#endif

#ifndef False
#define False           0
#endif

/**************************************/
/* defining internal structures       */

#define NSTRCTYPEMAX    20

namespace monte{

typedef struct structure IntStructType;
struct structure {
  char  name[20];
  int   yu;
  int   yl;
  int   xl;
  int   xr;
  double matrl;
  double volt;
  int   type;
};

/**************************************/


extern char  name[NSMAX][50], strc_name[NSTRCTYPEMAX][50];

extern double nc2p, fncx, fncy, xlength, ylength, zlength, rhoback, df,
      epsilon, dx, dy, dt, gtemp, bmag, btheta, bphi, pressure,
      lhs_conv_chrg, lhs_exti, lhs_phi0, lhs_fnys, lhs_fnyf,
      rhs_conv_chrg, rhs_exti, rhs_phi0, rhs_fnys, rhs_fnyf,
      lhs_sigma_total, lhs_sigma_total_1, lhs_extq, lhs_extq_1,
      rhs_sigma_total, rhs_sigma_total_1, rhs_extq, rhs_extq_1,
      lhs_length, lhs_dc, lhs_ac, lhs_w0, lhs_theta0, lhs_extc,
      rhs_length, rhs_dc, rhs_ac, rhs_w0, rhs_theta0, rhs_extc,
      tol_pois, del_area;

extern double vt[NSMAX][SIDES], v0x[NSMAX][SIDES], v0y[NSMAX][SIDES], v0z[NSMAX][SIDES],
      enter[NSMAX][SIDES], q[NSMAX], m[NSMAX], qm[NSMAX], q_per_cell[NSMAX],
      Exscale[NSMAX], Eyscale[NSMAX], Ezscale[NSMAX], dtdx[NSMAX], dtdy[NSMAX],
      ttx[NSMAX], tty[NSMAX], ttz[NSMAX], ssx[NSMAX], ssy[NSMAX], ssz[NSMAX],
      vxscale[NSMAX], vyscale[NSMAX], vzscale[NSMAX], ax_scale[NSMAX], ay_scale[NSMAX],
      grid_sec[NSMAX], rhs_sec[NSMAX], lhs_sec[NSMAX], wall_sec[NSMAX],
      lhs_chrg[NSMAX], rhs_chrg[NSMAX], weight[NSMAX];

extern int   nsp, ncx, ngx, ncy, ngy, sflag, interval, hist_hi, nfft, freq_hi,
      thist_hi, ecollisional, icollisional, ionsp,
      wall_flag, wall_read_flag, lhs_flag, lhs_nys, lhs_ny, lhs_nyf, lhs_gap,
      rhs_flag, rhs_nys, rhs_ny, rhs_nyf, rhs_gap,
      sp_k[NSMAX], np[NSMAX], inject[NSMAX], maxnp[NSMAX],
      k_count[NSMAX], it[NSMAX], sec_flag[NSMAX], sec_sp[NSMAX];

extern double **x, **y, **vx, **vy, **vz, *x_array, *y_array, **source, **eps_array, **area,
      **phi, **phi_lap_lhs, **phi_lap_rhs, **phi_pois, **phi_ave, **phi_ave_show,
      ***sp_n, ***sp_n_0, ***sp_n_k, ***sp_n_mcc, ***sp_n_ave, ***sp_n_ave_show,
      **rho, **ex, **ey, **ax, **ay, ***sp_ex, ***sp_ey, **sigma, ***sp_sigma,
      ***sp_vx0, ***sp_vy0, ***sp_vz0, ***sp_vt, **phi_intl;

extern int   nbin[SIDES][NSMAX];
extern double emin[SIDES][NSMAX], de[SIDES][NSMAX];
extern double ****fe, ***e_array, ****fth, ***th_array, ***wall_q,
      ***th_ave1, ***th_ave2, ***th_ave_show, ***e_ave_show,
      ***fetot_center, ***fetot_edge, ***fthtot_center, ***fthtot_edge;

extern double **elasrate, **extrate, **ionrate, **chrgxrate;

extern double *t_array, *Local_t_array, *f_array, **np_hist, **kes_hist,
      *phi_lhs_hist, *phi_lhs_fft, *com_phi_lhs_hist, *com_phs_lhs_hist, *com_phc_lhs_hist,
      *phi_rhs_hist, *phi_rhs_fft, *com_phi_rhs_hist, *com_phs_rhs_hist, *com_phc_rhs_hist,
      *cur_lhs_hist, *cur_lhs_fft, *com_cur_lhs_hist,
      *cur_rhs_hist, *cur_rhs_fft, *com_cur_rhs_hist,
      *pow_lhs_hist, *pow_lhs_fft, *com_pow_lhs_hist,
      *pow_rhs_hist, *pow_rhs_fft, *com_pow_rhs_hist,
      *phi_mid_hist, *phi_mid_fft, *com_phi_mid_hist;

extern int   stc_np[NSTRCTYPEMAX][NSMAX], nstrctype, nstrcmax, **cell_mask, **grid_mask,
      **cell_k, **face, **conductor;

extern double ***stc_np_hist;

extern IntStructType *strc;
extern double t;

extern long   int seed;

extern int    xbc_flag, ybc_flag, khi, khi_plus1;
extern double  **rhok, **phik, **show_rhok, **show_phik, **phi_mid_k_t, **phi_mid_k_f,
       *k_array, *dk, lhs_extl, lhs_extr, lhs_extq_2, lhs_extq_3;
extern void   (*fields_ptr)();
extern void   (*move_ptr)(int isp);



/*******************************************/
/*     Global functions                    */

void      mcc(int isp);
void      boundary(int isp), wall_diagnos(int isp);
void      init_dadi_arrays(int x0_flag, int xl_flag, int y0_flag, int ly_flag);
double    dadi(double dadi_dt, double **u_in, double **s, int itermax, double tol_test,
           int bound_flag, double *u_x0, double *u_xlx, double *u_y0, double *u_yly);
void      four1(double *data, int nn, int isign), realft(double *data, int n, int isign);
void      init_fields(), fields();
void      gather(int isp), setrho(), TWOD_One_2_One(double **matrix, int nx, int ny);
double     revers(int num, int n), frand();

void      load(int isp, double initn, int loader, int fill_region, double xleft, double xright,
           double ylow, double yhigh, double v0xi, double v0yi, double v0zi, double vti);
void      maxwellv(double *vx, double *vy, double *vz, double vth, double Rv, double Rphi, double Rthe);
void      vmaxwellv(int isp, int side, double *vx, double *vy, double *vz);
void      init_vmaxwellv(int isp, int side);

void      move(int isp), mag_move(int isp), move_internal(int isp);
void solve_mcm();
void start();

void fields2(), Periodic_Smooth(double **a, int nx, int ny);


}

#endif // MCM_H
