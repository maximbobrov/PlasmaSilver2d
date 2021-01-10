#include "mcm.h"

/**********************************************************************
  Dynamic ADI solution with Direchlet 0
  boundary conditions for the equation:
  
  dxdxu + dydyu =  s
  
  The function is based on a Peaceman Rachford Douglas
  advance of the parabolic equation:
  
  dtu = dxdxu + dydyu -s
  
  But with the time step dynamically chosen so as to speed convergence
  to the dtu=0 steady state solution which is what we want.  It is the
  user's responsiblity to put the initial guess of the solution stored
  in u when the function is called.  u=0 or u=u(previous time step)
  are possible choices.
  
  The function sends back the finishing iteration number, the finishing
  normalized maximum error and location, and the failure code.
  
  The failure codes are:
    ifail=0, success actually
    ifail=1, had to discard too many times
    ifail=2, used up the iterations
    
  Send in a negative tol_test if you want to freeze the time step to
  the initial choice.  Send in a 0 adel_t to let program pick an
  initial del_t.
  **********************************************************************/
/*
#define   DADI_DEBUG
*/

namespace monte{
int         x_start, x_end, y_start, y_end;

/*  Various copies of the 'answer' we're working on */
double **u, **uwork, **ustor, **ustar;

/*  The arrays used internal to DADI which contain the coefficients
    for each tridiagonal matrix solution */
double  *a_x, *b_x, *c_x, *r_x, *v_x, *gam_x;
double  *a_y, *b_y, *c_y, *r_y, *v_y, *gam_y;

/* The actual coefficients of the differencing on the mesh,
   and the influence of the spatial variation of epsilon:  */
double **a_xgeom, **b_xgeom, **c_xgeom, **a_ygeom, **b_ygeom, **c_ygeom;

void        adi(double **uadi, double **s, double del_t, int bound_flag,
                double *u_x0, double *u_xlx, double *u_y0, double *u_yly);
void        tridag(double *a, double *b, double *c, double *r,
                   double *utri, double *gam, int n);
void        ptridag(double *a, double *b, double *c, double *r,
                    double *utri, double *gam, int n);

/**********************************************************************/
void init_dadi_arrays(int x0_flag, int xl_flag, int y0_flag, int yl_flag)
{
    register int i, j;

    a_xgeom= (double **) malloc(ngx*sizeof(double *));
    b_xgeom= (double **) malloc(ngx*sizeof(double *));
    c_xgeom= (double **) malloc(ngx*sizeof(double *));
    a_ygeom= (double **) malloc(ngx*sizeof(double *));
    b_ygeom= (double **) malloc(ngx*sizeof(double *));
    c_ygeom= (double **) malloc(ngx*sizeof(double *));
    u      = (double **) malloc(ngx*sizeof(double *));
    uwork  = (double **) malloc(ngx*sizeof(double *));
    ustor  = (double **) malloc(ngx*sizeof(double *));
    ustar  = (double **) malloc(ngx*sizeof(double *));

    for(i=0; i<ngx; i++) {
        a_xgeom[i]= (double *) malloc(ngy*sizeof(double));
        b_xgeom[i]= (double *) malloc(ngy*sizeof(double));
        c_xgeom[i]= (double *) malloc(ngy*sizeof(double));
        a_ygeom[i]= (double *) malloc(ngy*sizeof(double));
        b_ygeom[i]= (double *) malloc(ngy*sizeof(double));
        c_ygeom[i]= (double *) malloc(ngy*sizeof(double));
        u[i]      = (double *) malloc(ngy*sizeof(double));
        uwork[i]  = (double *) malloc(ngy*sizeof(double));
        ustor[i]  = (double *) malloc(ngy*sizeof(double));
        ustar[i]  = (double *) malloc(ngy*sizeof(double));
    }

    for(j=0; j<ngy; j++) {
        uwork[0][j]  = ustor[0][j]  = u[0][j]  = 0.0;
        uwork[ncx][j]= ustor[ncx][j]= u[ncx][j]= 0.0;
    }
    for(i=0; i<ngx; i++) {
        uwork[i][0]  = ustor[i][0]  = u[i][0]  = 0.0;
        uwork[i][ncy]= ustor[i][ncy]= u[i][ncy]= 0.0;
    }

    for(i=1; i<ncx; i++)
        for(j=1; j<ncy; j++) {
            a_xgeom[i][j] = conductor[i][j]*0.5*(eps_array[i-1][j]+eps_array[i-1][j-1])/dx/dx;
            c_xgeom[i][j] = conductor[i][j]*0.5*(eps_array[i][j]+eps_array[i][j-1])/dx/dx;
            b_xgeom[i][j] = a_xgeom[i][j] +c_xgeom[i][j];

            a_ygeom[i][j] = conductor[i][j]*0.5*(eps_array[i][j-1]+eps_array[i-1][j-1])/dy/dy;
            c_ygeom[i][j] = conductor[i][j]*0.5*(eps_array[i-1][j]+eps_array[i][j])/dy/dy;
            b_ygeom[i][j] = a_ygeom[i][j] +c_ygeom[i][j];
        }

    /***** Fixing x0 boundary conditions  **********/
    if(x0_flag== DIRICHLET) {
        x_start=1;
        for(j=0; j<ngy; j++) {
            a_xgeom[0][j] = 0.0;
            b_xgeom[0][j] = 0.0;
            c_xgeom[0][j] = 0.0;

            a_ygeom[0][j] = 0.0;
            b_ygeom[0][j] = 0.0;
            c_ygeom[0][j] = 0.0;
        }
    }
    else if(x0_flag== NEUMANN) {
        x_start=0;
        for(j=1; j<ncy; j++) {
            a_xgeom[0][j] = 0.0;
            c_xgeom[0][j] = conductor[0][j]*(eps_array[0][j-1]+eps_array[0][j])/dx/dx;
            b_xgeom[0][j] = a_xgeom[0][j] +c_xgeom[0][j];

            a_ygeom[0][j] = conductor[0][j]*eps_array[0][j-1]/dy/dy;
            c_ygeom[0][j] = conductor[0][j]*eps_array[0][j]/dy/dy;
            b_ygeom[0][j] = a_ygeom[0][j] +c_ygeom[0][j];
        }
    }

    /***** Fixing xl boundary conditions  **********/
    if(xl_flag== DIRICHLET) {
        x_end=ncx-1;
        for(j=0; j<ngy; j++) {
            a_xgeom[ncx][j] = 0.0;
            b_xgeom[ncx][j] = 0.0;
            c_xgeom[ncx][j] = 0.0;

            a_ygeom[ncx][j] = 0.0;
            b_ygeom[ncx][j] = 0.0;
            c_ygeom[ncx][j] = 0.0;
        }
    }
    else if(xl_flag== NEUMANN) {
        x_end=ncx;
        for(j=1; j<ncy; j++) {
            //printf("eps [%d] = %f  eps = %f  cond = %d\n", j, eps_array[ncx-1][j], eps_array[ncx-1][j-1],conductor[ncx][j]);
            a_xgeom[ncx][j] = conductor[ncx][j]*(eps_array[ncx-1][j]+eps_array[ncx-1][j-1])/dx/dx;
            c_xgeom[ncx][j] = 0.0;
            b_xgeom[ncx][j] = a_xgeom[ncx][j] +c_xgeom[ncx][j];

            a_ygeom[ncx][j] = conductor[ncx][j]*eps_array[ncx-1][j-1]/dy/dy;
            c_ygeom[ncx][j] = conductor[ncx][j]*eps_array[ncx-1][j]/dy/dy;
            b_ygeom[ncx][j] = a_xgeom[ncx][j] +c_xgeom[ncx][j];
        }
    }

    /***** Fixing y0 boundary conditions  **********/
    if(y0_flag== DIRICHLET) {
        y_start=1;
        for(i=1; i<ncx; i++) {
            a_xgeom[i][0] = 0.0;
            b_xgeom[i][0] = 0.0;
            c_xgeom[i][0] = 0.0;

            a_ygeom[i][0] = 0.0;
            b_ygeom[i][0] = 0.0;
            c_ygeom[i][0] = 0.0;
        }
    }
    else if(y0_flag== NEUMANN) {
        y_start=0;
        for(i=1; i<ncx; i++) {
            a_xgeom[i][0] = conductor[i][0]*eps_array[i-1][0]/dx/dx;
            c_xgeom[i][0] = conductor[i][0]*eps_array[i][0]/dx/dx;
            b_xgeom[i][0] = a_xgeom[i][0] +c_xgeom[i][0];

            a_ygeom[i][0] = 0.0;
            c_ygeom[i][0] = conductor[i][0]*(eps_array[i-1][0]+eps_array[i][0])/dy/dy;
            b_ygeom[i][0] = a_ygeom[i][0] +c_ygeom[i][0];
        }
    }

    /***** Fixing yl boundary conditions  **********/
    if(yl_flag== DIRICHLET) {
        y_end=ncy-1;
        for(i=1; i<ncx; i++) {
            a_xgeom[i][ncy] = 0.0;
            b_xgeom[i][ncy] = 0.0;
            c_xgeom[i][ncy] = 0.0;

            a_ygeom[i][ncy] = 0.0;
            b_ygeom[i][ncy] = 0.0;
            c_ygeom[i][ncy] = 0.0;
        }
    }
    else if(yl_flag== NEUMANN) {
        y_end=ncy;
        for(i=1; i<ncx; i++) {
            a_xgeom[i][ncy] = conductor[i][ncy]*eps_array[i-1][ncy-1]/dx/dx;
            c_xgeom[i][ncy] = conductor[i][ncy]*eps_array[i][ncy-1]/dx/dx;
            b_xgeom[i][ncy] = a_xgeom[i][ncy] +c_xgeom[i][ncy];

            a_ygeom[i][ncy] = conductor[i][ncy]*(eps_array[i][ncy-1]+eps_array[i-1][ncy-1])/dy/dy;
            c_ygeom[i][ncy] = 0.0;
            b_ygeom[i][ncy] = a_ygeom[i][ncy] +c_ygeom[i][ncy];
        }
    }


    if(x0_flag==NEUMANN && y0_flag== NEUMANN) {
        a_xgeom[0][0] = 0.0;
        c_xgeom[0][0] = conductor[0][0]*eps_array[0][0]/dx/dx;
        b_xgeom[0][0] = a_xgeom[0][0] +c_xgeom[0][0];

        a_ygeom[0][0] = 0.0;
        c_ygeom[0][0] = conductor[0][0]*eps_array[0][0]/dy/dy;
        b_ygeom[0][0] = a_ygeom[0][0] +c_ygeom[0][0];
    }

    if(x0_flag==NEUMANN && yl_flag== NEUMANN) {
        a_xgeom[0][ncy] = 0.0;
        c_xgeom[0][ncy] = conductor[0][ncy]*eps_array[0][ncy-1]/dx/dx;
        b_xgeom[0][ncy] = a_xgeom[0][ncy] +c_xgeom[0][ncy];

        a_ygeom[0][ncy] = conductor[0][ncy]*eps_array[0][ncy-1]/dy/dy;
        c_ygeom[0][ncy] = 0.0;
        b_ygeom[0][ncy] = a_ygeom[0][ncy] +c_ygeom[0][ncy];
    }

    if(xl_flag==NEUMANN && yl_flag== NEUMANN) {
        a_xgeom[ncx][ncy] = conductor[ncx][ncy]*eps_array[ncx-1][ncy-1]/dx/dx;
        c_xgeom[ncx][ncy] = 0.0;
        b_xgeom[ncx][ncy] = a_xgeom[ncx][ncy] +c_xgeom[ncx][ncy];

        a_ygeom[ncx][ncy] = conductor[ncx][ncy]*eps_array[ncx-1][ncy-1]/dy/dy;
        c_ygeom[ncx][ncy] = 0.0;
        b_ygeom[ncx][ncy] = a_ygeom[ncx][ncy] +c_ygeom[ncx][ncy];
    }

    if(xl_flag==NEUMANN && y0_flag== NEUMANN) {
        a_xgeom[ncx][0] = conductor[ncx][0]*eps_array[ncx-1][0]/dx/dx;
        c_xgeom[ncx][0] = 0.0;
        b_xgeom[ncx][0] = a_xgeom[ncx][0] +c_xgeom[ncx][0];

        a_ygeom[ncx][0] = 0.0;
        c_ygeom[ncx][0] = conductor[ncx][0]*eps_array[ncx-1][0]/dy/dy;
        b_ygeom[ncx][0] = a_ygeom[ncx][0] +c_ygeom[ncx][0];
    }

    /************************************************************/
    /* Set up variables for tridiagonal inversion along x and y */

    a_x = (double *) malloc(ngx*sizeof(double));
    b_x = (double *) malloc(ngx*sizeof(double));
    c_x = (double *) malloc(ngx*sizeof(double));
    r_x = (double *) malloc(ngx*sizeof(double));
    v_x = (double *) malloc(ngx*sizeof(double));
    gam_x=(double *) malloc(ngx*sizeof(double));

    a_y = (double *) malloc(ngy*sizeof(double));
    b_y = (double *) malloc(ngy*sizeof(double));
    c_y = (double *) malloc(ngy*sizeof(double));
    r_y = (double *) malloc(ngy*sizeof(double));
    v_y = (double *) malloc(ngy*sizeof(double));
    gam_y=(double *) malloc(ngy*sizeof(double));

    for(i=0; i<ngx; i++) {
        a_x[i] = 0.0;
        b_x[i] = 1.0;
        c_x[i] = 0.0;
        r_x[i] = 0.0;
    }
    for(j=0; j<ngy; j++) {
        a_y[j] = 0.0;
        b_y[j] = 1.0;
        c_y[j] = 0.0;
        r_y[j] = 0.0;
    }
}

/**********************************************************************/

double dadi(double dadi_dt, double **u_in, double **s, int itermax,
            double tol_test, int bound_flag, double *u_x0,
            double *u_xlx, double *u_y0, double *u_yly)
{
    register int i, j, ip1, im1, jp1, jm1;
    int iter, ndiscard;
    double del_t, del_td, tptop, tpbot, ratio;
    double rnorm, rsum, res, errchk, dxdxutrm, dydyutrm;
#ifdef DADI_DEBUG
    static double iter_ave=0.0, res_ave=0.0;
#endif

    /***************************************************************/
    /* Start initial step-size as previous step size divided by n. */

    del_t = dadi_dt/16.0;
    del_td= 2.0*del_t;
    ndiscard=0;

    /* Residual normalization.  */
    rnorm=0.0;
    for(i=x_start; i<=x_end; i++)
        for(j=y_start; j<=y_end; j++) rnorm += s[i][j]*s[i][j]*conductor[i][j];

    if (fabs(rnorm) < 1e-30) {
        if(!bound_flag) return(dadi_dt);
        else            rnorm = 1.0;
    }
    rnorm=sqrt(rnorm);

    /*******************************************/
    /* copy u_in to u for working purposes.    */

    for(i=0; i<ngx; i++)
        for(j=0; j<ngy; j++) u[i][j] = u_in[i][j];

    /**********************************************************/
    /* Fixing the boundary conditions for u and other arrays. */

    if(bound_flag) {
        for(i=0; i<ngx; i++) {
            uwork[i][0]  = ustor[i][0]  = u[i][0]  = u_y0[i];
            uwork[i][ncy]= ustor[i][ncy]= u[i][ncy]= u_yly[i];
        }
        for(j=0; j<ngy; j++) {
            uwork[0][j]  = ustor[0][j]  = u[0][j]  = u_x0[j];
            uwork[ncx][j]= ustor[ncx][j]= u[ncx][j]= u_xlx[j];
        }
    }

    /********************/
    /* Begin iteration. */

    for(iter=0; iter<itermax; iter++) {
        /*************************************************/
        /* Copy u into the work array and storage array. */

        for(i=0; i<ngx; i++)
            for(j=0; j<ngy; j++) uwork[i][j] = ustor[i][j] = u[i][j];

        /************************************/
        /* Two advances of u via ADI at del_t. */

        adi(u, s, del_t, bound_flag, u_x0, u_xlx, u_y0, u_yly);
        adi(u, s, del_t, bound_flag, u_x0, u_xlx, u_y0, u_yly);

        /*****************************************/
        /* One advance of uwork via ADI at 2*del_t. */

        adi(uwork, s, del_td, bound_flag, u_x0, u_xlx, u_y0, u_yly);

        /*******************************************************/
        /* Calculate test parameter and normalized error.
       For Dirichlet BCs, no need to worry about boundary
       points since u,uwork, and ustor should be the same. */

        tptop= tpbot= rsum = 0.0;
        for(i=x_start; i<=x_end; i++) {
            im1= (i)?i-1:0;
            ip1= (i==ncx)?ncx:i+1;
            for(j=y_start; j<=y_end; j++) {
                jm1= (j)?j-1:0;
                jp1= (j==ncy)?ncy:j+1;

                /* Test paramter sums. */
                tptop += (u[i][j]-uwork[i][j])*(u[i][j]-uwork[i][j]);
                tpbot += (u[i][j]-ustor[i][j])*(u[i][j]-ustor[i][j]);

                /* Residual terms. */
                dxdxutrm= a_xgeom[i][j]*u[im1][j] -b_xgeom[i][j]*u[i][j]
                        +c_xgeom[i][j]*u[ip1][j];
                dydyutrm= a_ygeom[i][j]*u[i][jm1] -b_ygeom[i][j]*u[i][j]
                        +c_ygeom[i][j]*u[i][jp1];

                /* Residual sums.  Only include points outside of structures.  */
                errchk = conductor[i][j]*(dxdxutrm +dydyutrm -s[i][j]);
                rsum += errchk*errchk;
            }
        }
        /* Calculate normalized residual. */
        res = sqrt(rsum)/rnorm;

        /* If the residual is less than the tolerance, SUCCESS! */
        if ((res < tol_test) && (iter)) {
#ifdef DADI_DEBUG
            res_ave += res;
            iter_ave+= iter;
            printf("dadi: SUCCESS it= %d iter=%d, res=%lg, iter_ave=%lg, res_ave=%lg, del_t=%lg\n\n",
                   it[0], iter, res, (it[0])?iter_ave/it[0]:1, (it[0])?res_ave/it[0]:1, del_t);
#endif
            for(i=0;i<ngx;i++)
                for(j=0;j<ngy;j++) u_in[i][j] = u[i][j];
            return(del_t);
        }

        /* Determine ratio used to find the time step change.  If tpbot
       is zero but tptop is finite, consider this a case of a large
       ratio and act accordingly.  DWH does about the same thing
       except he does NOT discard the solution if tpbot=0. */

        if (tpbot >  0.0) ratio=tptop/tpbot;
        if (tpbot <1e-30) ratio=1.0;
#ifdef DADI_DEBUG
        printf("dadi: iter=%d, res=%lg, tol=%lg, del_t=%lg\n", iter, res, tol_test, del_t);
        printf("      rat=%lg, tptop=%lg, tpbot=%lg\n", ratio, tptop, tpbot);
#endif
        /* Get next time step. */
        if (ratio < 0.02)      del_t *= 8.000;
        else if (ratio < 0.05) del_t *= 4.000;
        else if (ratio < 0.10) del_t *= 2.000;
        else if (ratio < 0.30) del_t *= 0.800;
        else if (ratio < 0.40) del_t *= 0.500;
        else if (ratio < 0.60) del_t *= 0.250;

        /* Ratio is too large. */
        else {
            ndiscard++;
#ifdef DADI_DEBUG
            printf("    iter=%d, res=%lg, ndiscard= %d\n", iter, res, ndiscard);
#endif     
            /* Check if too many discards. */
            if (ndiscard > 20) {
#ifdef DADI_DEBUG
                printf("  it= %d ndiscard>20\n\n", it[0]);
#endif
                for(i=0; i<ngx; i++)
                    for(j=0; j<ngy; j++) u_in[i][j] = u[i][j];
                return(del_t);
            }
            /* Discard by replacing u with what we started with. */
            for(i=0; i<ngx; i++)
                for(j=0; j<ngy; j++) u[i][j] = ustor[i][j];

            /* Reduce del_t. */
            del_t /= 16.0;
        }
        del_td=2*del_t;
    }
    /* Fail if used up the maximum iterations. */
#ifdef DADI_DEBUG
    printf("        FAILED iter>= %d, res=%lg, tol=%lg, del_t=%lg\n\n",
           itermax, res, tol_test, del_t);
#endif
    for(i=0; i<ngx; i++)
        for(j=0; j<ngy; j++) u_in[i][j] = u[i][j];
    return(del_t);
}

/**********************************************************************
  Single Peaceman Rachford Douglas pass with Direchlet 0 c boundary
  conditions for the equation:
  
  dtu = dxdxu + dydyu -s, where s is constant in time.
  
  The Crank-Nicolson finite difference approximation to the
  above equation leads to the fractional step or
  ADI equations:
  
  u*(i,j)-(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)]
  = un(i,j)+(del_t/2dydy)[un(i,j+1)-2un(i,j)+un(i,j-1)] - (del_t/2)s(i,j)
  
  un+1(i,j)-(del_t/2dydy)[un+1(i,j+1)-2un+1(i,j)+un+1(i,j-1)]
  = u*(i,j)+(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)] - (del_t/2)s(i,j)
  
  **********************************************************************/

void adi(double **uadi, double **s, double del_t, int bound_flag,
         double *u_x0, double *u_xlx, double *u_y0, double *u_yly)
{
    register int i, j, ip1, im1, jp1, jm1;
    double dthi;

    dthi = -2.0/del_t;

    /***************************************/
    /* Do x pass.  Set up variables for    */
    /* tridiagonal inversion along x.      */

    for(j=y_start; j<=y_end; j++) {
        if(bound_flag) {
            r_x[0]  = u_x0[j];
            r_x[ncx]= u_xlx[j];
        }
        else {
            r_x[0]  = 0.0;
            r_x[ncx]= 0.0;
        }
        jm1= (j)?j-1:0;
        jp1= (j==ncy)?ncy:j+1;
        for(i=x_start; i<=x_end; i++) {
            a_x[i] = a_xgeom[i][j];
            b_x[i] = dthi -b_xgeom[i][j];
            c_x[i] = c_xgeom[i][j];
            r_x[i] = (dthi*uadi[i][j] +s[i][j])*conductor[i][j]
                    -( a_ygeom[i][j]*uadi[i][jm1]
                       -b_ygeom[i][j]*uadi[i][j]
                       +c_ygeom[i][j]*uadi[i][jp1]) +  /*sin(t * 6.28e7) **/ dthi*phi_intl[i][j]*(1-bound_flag) /* * (1 + 8*((sin(t*6.28/1e-6))))*/;
        }

        /* Solve tridiagonal system. */
        tridag(a_x, b_x, c_x, r_x, v_x, gam_x, ngx);

        /* Copy solution into ustar. */
        for(i=0; i<ngx; i++) ustar[i][j] =v_x[i];
    }

    /***************************************/
    /* Do y pass.  Set up variables for    */
    /* tridiagonal inversion along y.      */

    for(i=x_start; i<=x_end; i++) {

        if(bound_flag) {
            r_y[0]  = u_y0[i];
            r_y[ncy]= u_yly[i];
        }
        else {
            r_y[0]  = 0.0;
            r_y[ncy]= 0.0;
        }
        im1= (i)?i-1:0;
        ip1= (i==ncx)?ncx:i+1;
        for(j=y_start; j<=y_end; j++) {
            a_y[j] = a_ygeom[i][j];
            b_y[j] = dthi -b_ygeom[i][j];
            c_y[j] = c_ygeom[i][j];
            r_y[j] = (dthi*ustar[i][j] +s[i][j])*conductor[i][j]
                    -( a_xgeom[i][j]*ustar[im1][j]
                       -b_xgeom[i][j]*ustar[i][j]
                       +c_xgeom[i][j]*ustar[ip1][j]) + /*sin(t * 6.28e7) **/ dthi*phi_intl[i][j]*(1-bound_flag) /* * (1 + 8*((sin(t*6.28/1e-6))))*/;
        }

        /* Solve tridiagonal system. */
        tridag(a_y, b_y, c_y, r_y, v_y, gam_y, ngy);

        /* Copy solution into uadi. */
        for(j=0; j<ngy; j++) uadi[i][j] =v_y[j];
    }
    /***************************************/
}

/***********************************************************************
  Tridiagonal field solver:
  
  | b0 c0 0                      | | u0 |     | r0 |
  | a1 b1 c1                     | | u1 |     | r1 |
  |      ............            | | .  |  =  | .  |
  |               an-2 bn-2 cn-2 | |un-2|     |rn-2|
  |               0    an-1 bn-1 | |un-1|     |rn-1|
  
  **********************************************************************/

void tridag(double *a, double *b, double *c, double *r,
            double *utri, double *gam, int n)
{
    register int i;
    double bet;

    /*******************************************/
    /* Decomposition and forward substitution. */

    bet = b[0];
    utri[0]= r[0]/bet;

    for (i=1; i<n; i++) {
        gam[i] = c[i-1]/bet;
        bet = b[i] - a[i]*gam[i];
        utri[i] = (r[i] -a[i]*utri[i-1])/bet;
    }

    /**********************/
    /* Back substitution. */
    for(i=n-2; i>=0; i--) utri[i] -= gam[i+1]*utri[i+1];
}

/******************************************************
 Periodic Tridiagnol field solver from Dennis Hewett
  via Greg DiPeso.
  
  | b0 c0 0                  a0  | | u0 |     | r0 |
  | a1 b1 c1                     | | u1 |     | r1 |
  |      ............            | | .  |  =  | .  |
  |               an-2 bn-2 cn-2 | |un-2|     |rn-2|
  |cn-1           0    an-1 bn-1 | |un  |     |rn-1|
  
  **********************************************************************/

void ptridag(double *a, double *b, double *c, double *r,
             double *utri, double *gam, int n)
{
    register int i, n_1;
    double bet;

    /*******************************************/
    /* Decomposition and forward substitution. */

    n_1 = n-1;
    for (i=2; i<n_1; i++) {
        c[i-1] /= b[i-1];
        r[i-1] /= b[i-1];
        a[i-1] /= b[i-1];

        b[i] -= a[i]*c[i-1];
        r[i] -= a[i]*q[i-1];
        a[i] *= -a[i-1];

        b[n_1] -= c[n_1]*a[i-1];
        r[n_1] -= c[n_1]*r[i-1];
        c[n_1] *= -c[i-1];
    }
    r[n-2]   /= b[n-2];
    c[n-2]    = (a[n-2]+c[n-2])/b[n-2];
    b[n-1]   -= (a[n-1]+c[n-1])*c[n-2];
    r[n-1]   -= (a[n-1]+c[n-1])*q[n-2];
    utri[n-1] = r[n-1]/b[n-1];
    utri[n-2] = r[n-2] -c[n-2]*utri[n];

    /**********************/
    /* Back substitution. */
    for(i=n-2; i>=1; i--) utri[i] = r[i] -c[i]*utri[i+1] -a[i]*utri[n-1];
    utri[0] = utri[n-1];
}
}
/******************************************************/

