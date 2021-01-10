#include "mcm.h"

#define REFLECTION_RATIO    1.0
namespace monte{
void boundary_diag(int k, int isp, int i);
void find_intercept(double xp, double yp, double vxp, double vyp,
                    int k, double *xint, double *yint, int *side);
int  emitE(int isp, double xp, double yp);

/****************************************************************/
/* Routine to adjust (initialize, re-pack and inject) particles */
/* to the desired boundary conditions                           */

void boundary(int isp)
{
    register int i, nnp, ninjected;
    static int init_flag=1;
    static double extra[NSMAX][SIDES];
    double del_t;

    if(init_flag) {
        for(i=0; i<nsp; i++)
            extra[i][LEFT] = extra[i][RIGHT] = 0.5;
        init_flag=0;
    }

    rhs_chrg[isp] = lhs_chrg[isp]= 0.0;

    /* Injection from walls */
    if(inject[isp]) {
        extra[isp][LEFT]+= enter[isp][LEFT];
        ninjected = extra[isp][LEFT];
        while(extra[isp][LEFT] >= 1.0) {
            del_t= ((int)extra[isp][LEFT])/(ninjected +1.);
            nnp = np[isp];
            vmaxwellv(isp, LEFT, &vx[isp][nnp], &vy[isp][nnp], &vz[isp][nnp]);
            vx[isp][nnp] /= vxscale[isp];
            vy[isp][nnp] /= vyscale[isp];
            vz[isp][nnp] /= vzscale[isp];
            x[isp][nnp] = del_t*vx[isp][nnp];
            y[isp][nnp] = (lhs_length*frand()+0.5*(ylength -lhs_length))/dy;
            np[isp]++;
            lhs_chrg[isp] -= weight[isp]*q[isp];
            extra[isp][LEFT] -=1.0;
        }

        extra[isp][RIGHT]+= enter[isp][RIGHT];
        ninjected = extra[isp][RIGHT];
        while(extra[isp][RIGHT] >= 1.0) {
            del_t= ((int)extra[isp][RIGHT])/(ninjected +1.);
            nnp = np[isp];
            vmaxwellv(isp, RIGHT, &vx[isp][nnp], &vy[isp][nnp], &vz[isp][nnp]);
            vx[isp][nnp] /=-vxscale[isp];
            vy[isp][nnp] /= vyscale[isp];
            vz[isp][nnp] /= vzscale[isp];
            x[isp][nnp] = fncx +del_t*vx[isp][nnp];
            y[isp][nnp] = (rhs_length*frand()+0.5*(ylength -rhs_length))/dy;
            np[isp]++;
            rhs_chrg[isp] -= weight[isp]*q[isp];
            extra[isp][RIGHT] -=1.0;
        }
    }

    if (np[isp]>0) {

        /****************************************/
        /***  Then do the other boundaries   ****/
        for(i=0; i < np[isp]; i++) {

            /*****  Particle has crossed the RHS boundary  ****/
            if (x[isp][i] >= fncx) {
                if(rhs_fnys <= y[isp][i] && y[isp][i] <= rhs_fnyf) {
                    rhs_chrg[isp] += weight[isp]*q[isp];

                    /* Inject a secondary */
                    if(sec_flag[isp] && frand() < rhs_sec[isp])
                        if(emitE(sec_sp[isp], fncx-SNG_MIN, y[isp][i]))
                            rhs_chrg[sec_sp[isp]] -= weight[isp]*q[sec_sp[isp]];
                }
                else if(rhs_fnyf+rhs_gap <= y[isp][i] || y[isp][i] <= rhs_fnys-rhs_gap) {

                    /**** Inject a secondary *****/
                    if(sec_flag[isp] && frand() <wall_sec[isp])
                        emitE(sec_sp[isp], fncx-SNG_MIN, y[isp][i]);
                }
                if(wall_flag) boundary_diag(RIGHT, isp, i);
                /* Remove the particle */
                x[isp][i] =  x[isp][np[isp]-1];
                y[isp][i] =  y[isp][np[isp]-1];
                vx[isp][i]= vx[isp][np[isp]-1];
                vy[isp][i]= vy[isp][np[isp]-1];
                vz[isp][i]= vz[isp][np[isp]-1];
                i--;
                np[isp]--;
            }
            /*****  Particle has crossed the LHS boundary  ****/
            else if (x[isp][i] <= 0.0) {
                if(lhs_fnys <= y[isp][i] && y[isp][i] <= lhs_fnyf) {
                    lhs_chrg[isp] += weight[isp]*q[isp];

                    /**** Inject a secondary *****/
                    if(sec_flag[isp] && frand() < lhs_sec[isp])
                        if(emitE(sec_sp[isp], SNG_MIN, y[isp][i]))
                            lhs_chrg[sec_sp[isp]] -= weight[isp]*q[sec_sp[isp]];
                }
                else if(lhs_fnyf+lhs_gap <= y[isp][i] || y[isp][i] <= lhs_fnys-lhs_gap) {

                    /**** Inject a secondary *****/
                    if(sec_flag[isp] && frand() <wall_sec[isp])
                        emitE(sec_sp[isp], SNG_MIN, y[isp][i]);
                }
                if(wall_flag) boundary_diag(LEFT, isp, i);
                /* Remove the particle */
                x[isp][i]  =  x[isp][np[isp]-1];
                y[isp][i]  =  y[isp][np[isp]-1];
                vx[isp][i] = vx[isp][np[isp]-1];
                vy[isp][i] = vy[isp][np[isp]-1];
                vz[isp][i] = vz[isp][np[isp]-1];
                i--;
                np[isp]--;
            }
            /*****  Particle has crossed the TOP boundary  ****/
            else if (y[isp][i] >= fncy) {

                if(ybc_flag == DIRICHLET)
                {
                    /**** Inject a secondary *****/
                    /*if(sec_flag[isp] && frand() < wall_sec[isp])
        emitE(sec_sp[isp], x[isp][i], fncy-SNG_MIN);*/

                    /* Remove the particle */
                    x[isp][i]  =  x[isp][np[isp]-1];
                    y[isp][i]  =  y[isp][np[isp]-1];
                    vx[isp][i] = vx[isp][np[isp]-1];
                    vy[isp][i] = vy[isp][np[isp]-1];
                    vz[isp][i] = vz[isp][np[isp]-1];
                    i--;
                    np[isp]--;
                }
                else if(ybc_flag == NEUMANN) {
      //if(frand() <= REFLECTION_RATIO)
      {
        vy[isp][i] *= -1.0;
        if(bmag > 0.0) y[isp][i] += dtdy[isp]*vy[isp][i];
        else           y[isp][i] += vy[isp][i];
      }
    }
    else
      while ((y[isp][i]-=fncy) > fncy);
            }
            /*****  Particle has crossed the BOTTOM boundary  ****/
            else if (y[isp][i] <= 0.0) {

                if(ybc_flag == DIRICHLET) {
                    /**** Inject a secondary *****/
                    if(sec_flag[isp] && frand() < wall_sec[isp])
                        emitE(sec_sp[isp], x[isp][i], SNG_MIN);

                    /* Remove the particle */
                    x[isp][i]  =  x[isp][np[isp]-1];
                    y[isp][i]  =  y[isp][np[isp]-1];
                    vx[isp][i] = vx[isp][np[isp]-1];
                    vy[isp][i] = vy[isp][np[isp]-1];
                    vz[isp][i] = vz[isp][np[isp]-1];
                    i--;
                    np[isp]--;
                }
                else if(ybc_flag == NEUMANN) {
                    if(frand() <= REFLECTION_RATIO) {
                        vy[isp][i] *= -1.0;
                        if(bmag > 0.0) y[isp][i] += dtdy[isp]*vy[isp][i];
                        else           y[isp][i] += vy[isp][i];
                    }
                    else
                        while ( (y[isp][i]+=fncy) < 0.0);
                }
                else
                    while ( (y[isp][i]+=fncy) < 0.0);
            }
        }

        /****************************************/
        /* For Stuctures, remove the particles  */
        if (nstrcmax>0) {
            register int ix, iy;
            double xint, yint, delx, dely;
            int side;

            for(i=0; i < np[isp]; i++) {
                ix= x[isp][i];
                iy= y[isp][i];
                if(cell_mask[ix][iy]) {
                    stc_np[strc[cell_k[ix][iy]].type][isp]+= weight[isp];
                    find_intercept(x[isp][i], y[isp][i], vx[isp][i], vy[isp][i],
                                   cell_k[ix][iy], &xint, &yint, &side);
                    ix= xint;
                    iy= yint;
                    if(side == DOWN || side == UP) {
                        delx= xint -ix;
                        sp_sigma[isp][ix][iy]  += weight[isp]*(1-delx);
                        sp_sigma[isp][ix+1][iy]+= weight[isp]*delx;
                    }
                    else {
                        dely= yint -iy;
                        sp_sigma[isp][ix][iy]  += weight[isp]*(1-dely);
                        sp_sigma[isp][ix][iy+1]+= weight[isp]*dely;
                    }

                    /* Inject a secondary */
                    if(sec_flag[isp] && frand() < grid_sec[isp]) {
                        if(side == DOWN)       emitE(sec_sp[isp], xint, yint-SNG_MIN);
                        else if(side == UP)    emitE(sec_sp[isp], xint, yint+SNG_MIN);
                        else if(side == LEFT)  emitE(sec_sp[isp], xint-SNG_MIN, yint);
                        else if(side == RIGHT) emitE(sec_sp[isp], xint+SNG_MIN, yint);
                    }
                    //if(eps_array[ix][iy]/EPS0 < 1)
                    {
                        /* Remove the particle */
                        x[isp][i]  =  x[isp][np[isp]-1];
                        y[isp][i]  =  y[isp][np[isp]-1];
                        vx[isp][i] = vx[isp][np[isp]-1];
                        vy[isp][i] = vy[isp][np[isp]-1];
                        vz[isp][i] = vz[isp][np[isp]-1];
                        i--;
                        np[isp]--;
                    }
                }
            }
        }
    }
    rhs_chrg[isp] /= sp_k[isp];
    lhs_chrg[isp] /= sp_k[isp];
}

/**************************************************************/

void find_intercept(double xp, double yp, double vxp, double vyp,
                    int k, double *xint, double *yint, int *side)
{
    /***** particle arrived in the box from below **/
    if(vyp > 0.0) {
        *xint = (strc[k].yl - yp)*vxp/vyp +xp;

        if(*xint < strc[k].xl) {
            *xint = (double)strc[k].xl;
            *yint = (strc[k].xl - xp)*vyp/vxp +yp;
            *side= LEFT;
        }
        else if(*xint > strc[k].xr) {
            *xint = (double)strc[k].xr;
            *yint = (strc[k].xr - xp)*vyp/vxp +yp;
            *side= RIGHT;
        }
        else {
            *yint = (double)strc[k].yl;
            *side= DOWN;
        }
    }
    /***** picle arrived in the box from above **/
    else if(vyp < 0.0) {
        *xint = (strc[k].yu - yp)*vxp/vyp +xp;

        if(*xint < strc[k].xl) {
            *xint = (double)strc[k].xl;
            *yint = (strc[k].xl - xp)*vyp/vxp +yp;
            *side= LEFT;
        }
        else if(*xint > strc[k].xr) {
            *xint = (double)strc[k].xr;
            *yint = (strc[k].xr - xp)*vyp/vxp +yp;
            *side= RIGHT;
        }
        else {
            *yint = (double)strc[k].yu;
            *side= UP;
        }
    }
    /***** particle arrived in the box from above **/
    else {
        if (vxp > 0.0) {
            *xint = (double)strc[k].xl;
            *yint = yp;
            *side= LEFT;
        }
        else {
            *xint = (double)strc[k].xr;
            *yint = yp;
            *side= RIGHT;
        }
    }
}

/**************************************************************/

int emitE(int isp, double xp, double yp)
{
    int nnp;

    if (np[isp] >= maxnp[isp]) {
        puts("emitE: too many particles. MUST EXIT!");
        exit(-1);
    }
    nnp = np[isp];

    x[isp][nnp] = xp;
    y[isp][nnp] = yp;
    vx[isp][nnp] = 0.0;
    vy[isp][nnp] = 0.0;
    vz[isp][nnp] = 0.0;

    if(ybc_flag != PERIODIC) {
        if(x[isp][nnp]>=0.0 && x[isp][nnp]<fncx && y[isp][nnp]>=0.0 && y[isp][nnp]<fncy) {
            np[isp]++;
            return(1);
        }
        else
            return(0);
    }
    else {
        if(x[isp][nnp]< 0.0 || x[isp][nnp]>=fncx)
            return(0);
        else {
            if (y[isp][nnp] >= fncy)
                while ((y[isp][nnp]-=fncy) > fncy);
            if (y[isp][nnp] <= 0.0)
                while ((y[isp][nnp]+=fncy) < 0.0);
            np[isp]++;
            return(1);
        }
    }
}
/**************************************************************/

void boundary_diag(int k, int isp, int i)
{
    double dum, dum1, vxtemp, vytemp, vztemp;
    int s, s1;

    vxtemp= vxscale[isp]*vx[isp][i];
    vytemp= vyscale[isp]*vy[isp][i];
    vztemp= vzscale[isp]*vz[isp][i];

    dum= y[isp][i];
    s  = dum;

    if(s<ncy && s>=0) {
        dum -=s;
        wall_q[k][isp][s]  += weight[isp]*(1-dum);
        wall_q[k][isp][s+1]+= weight[isp]*dum;

        /****  Calculating the f(y, E) diagnostic  ****/
        dum1 = 0.5*m[isp]*(vxtemp*vxtemp +vytemp*vytemp +vztemp*vztemp)/(1.602e-19);
        dum1 = (dum1 -emin[k][isp])/de[k][isp];
        s1   = dum1;

        if (s1<nbin[k][isp]-1 && dum1>=0) {
            dum1 -= s1;
            fe[k][isp][s][s1]    += (1-dum)*(1-dum1);
            fe[k][isp][s+1][s1]  += dum*(1-dum1);
            fe[k][isp][s][s1+1]  += (1-dum)*dum1;
            fe[k][isp][s+1][s1+1]+= dum*dum1;
        }

        /****  Calculating the f(y, Theta) diagnostic  ****/
        dum1 = atan(sqrt(vytemp*vytemp+vztemp*vztemp)/(fabs(vxtemp) + SNG_MIN));
        dum1*= (nbin[k][isp]-1)/(0.5*M_PI);
        s1   = dum1;

        if (s1<nbin[k][isp]-1 && dum1>=0) {
            dum1 -= s1;
            fth[k][isp][s][s1]    += (1-dum)*(1-dum1);
            fth[k][isp][s+1][s1]  += dum*(1-dum1);
            fth[k][isp][s][s1+1]  += (1-dum)*dum1;
            fth[k][isp][s+1][s1+1]+= dum*dum1;
        }
    }
}

/*********************************************************************/
/* Wall diagnostics such as ion angular and energy       *************/
/* distributions are collected here.                     *************/

void wall_diagnos(int isp)
{
    register int i, j, k;
    int locl_nys[SIDES], locl_nyf[SIDES];
    double e_ave_tot,  fetot_center_tot,  fetot_edge_tot;
    double th_ave_tot, fthtot_center_tot, fthtot_edge_tot, temp;

    locl_nys[LEFT] = lhs_nys;
    locl_nyf[LEFT] = lhs_nyf;
    locl_nys[RIGHT]= rhs_nys;
    locl_nyf[RIGHT]= rhs_nyf;

    for(k=0; k<SIDES; k++) {
        /**************************************************/
        /** Local averages of ion energy dist. function ***/

        fetot_center_tot= fetot_edge_tot =0.0;
        fthtot_center_tot=fthtot_edge_tot=0.0;
        for(i=0; i<nbin[k][isp]; i++) {

            /** Average over the center  **/
            fetot_center[k][isp][i] = fe[k][isp][ncy/2][i];
            fthtot_center[k][isp][i]= fth[k][isp][ncy/2][i];
            for(j=1; j< ncy/8; j++){
                fetot_center[k][isp][i] += fe[k][isp][ncy/2+j][i];
                fetot_center[k][isp][i] += fe[k][isp][ncy/2-j][i];
                fthtot_center[k][isp][i]+= fth[k][isp][ncy/2+j][i];
                fthtot_center[k][isp][i]+= fth[k][isp][ncy/2-j][i];
            }

            /** Average over the edge    **/
            fetot_edge[k][isp][i] = 0.0;
            fthtot_edge[k][isp][i]= 0.0;
            for(j=0; j< ncy/8; j++){
                fetot_edge[k][isp][i] += fe[k][isp][locl_nys[k]+j][i];
                fetot_edge[k][isp][i] += fe[k][isp][locl_nyf[k]-j][i];
                fthtot_edge[k][isp][i]+= fth[k][isp][locl_nys[k]+j][i];
                fthtot_edge[k][isp][i]+= fth[k][isp][locl_nyf[k]-j][i];
            }
            fetot_center_tot += fetot_center[k][isp][i];
            fetot_edge_tot   += fetot_edge[k][isp][i];
            fthtot_center_tot+= fthtot_center[k][isp][i];
            fthtot_edge_tot  += fthtot_edge[k][isp][i];
        }

        if(fetot_center_tot > 1e-20)
            for(i=0; i<nbin[k][isp]; i++) fetot_center[k][isp][i] /= fetot_center_tot;
        if(fetot_edge_tot   > 1e-20)
            for(i=0; i<nbin[k][isp]; i++) fetot_edge[k][isp][i]   /= fetot_edge_tot;
        if(fthtot_center_tot> 1e-20)
            for(i=0; i<nbin[k][isp]; i++) fthtot_center[k][isp][i]/= fthtot_center_tot;
        if(fthtot_edge_tot  > 1e-20)
            for(i=0; i<nbin[k][isp]; i++) fthtot_edge[k][isp][i]  /= fthtot_edge_tot;

        /** Average over all angles at each position **/
        for(j=0; j<=ncy; j++) {
            e_ave_tot = e_ave_show[k][isp][j] = 0.0;
            th_ave_tot = th_ave_show[k][isp][j] = 0.0;

            for(i=0; i<nbin[k][isp]; i++) {
                e_ave_show[k][isp][j] += e_array[k][isp][i]*fe[k][isp][j][i];
                e_ave_tot += fe[k][isp][j][i];

                th_ave_show[k][isp][j]+= th_array[k][isp][i]*fth[k][isp][j][i];
                th_ave_tot+= fth[k][isp][j][i];
            }
            if(e_ave_tot > 1e-20) e_ave_show[k][isp][j] /= e_ave_tot;
            if(th_ave_tot > 1e-20) th_ave_show[k][isp][j] /= th_ave_tot;
        }
    }
}
}
/**************************************************************/
