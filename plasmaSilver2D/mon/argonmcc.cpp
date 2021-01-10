#include "mcm.h"

#define   NEMAX   500
#define   NIMAX   500
namespace monte{
void newvel(double energy, double vel, double *vx, double *vy, double *vz, int e_flag);
double sigma1(double energy), sigma2(double energy), sigma3(double energy);
double sigma4(double energy), sigma5(double energy);
double gden, vgth, extengy0, ionengy0;
int ecolsp, icolsp;

/*********************************************************/
/* Monte Carlo Collisions for electron- and ion-neutrals */

void mcc(int isp)
{

    static double _max_sigmav_e, _max_sigmav_i;
    static double _col_prob_e, _col_prob_i;

    register int i, j, k, index;
    int N, nnp, s1, s2;
    double random, vel, sigma_total;
    double engy, rengy, phi1, delx, dely;
    double cosphi, sinphi, coschi, sinchi;
    double sum_sigma, temp, vneutx, vneuty, vneutz;


    static int init_mcc_flag=1;
    static double ecol_extra, icol_extra;
    static double max_sigmav_e, max_sigmav_i;
    static double col_prob_e, col_prob_i;

    if (init_mcc_flag) {
        gden = NperTORR*pressure/(gtemp+SNG_MIN);  /* Calculating the gas density */

        ecolsp= ecollisional-1;
        icolsp= icollisional-1;
        ionsp--;
        if(ionsp>0) vgth= sqrt(2*1.602e-19*gtemp/m[ionsp]);

        /*****************************************/
        /* Calculating the null collision prob.  */

        if(ecollisional) {
            extengy0 = 11.55;
            ionengy0 = 15.76;

            max_sigmav_e = 0.0;
            for (engy=0; engy<NEMAX; engy += 0.1)
                max_sigmav_e= max(max_sigmav_e,
                                  sqrt(2*1.602e-19*engy/m[ecolsp])
                                  *(sigma1(engy)+sigma2(engy)+sigma3(engy)));

            col_prob_e = 1 -exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
            ecol_extra = 0.5;
        }
        if(icollisional) {
            max_sigmav_i = 0.0;
            for (engy=0; engy<NIMAX; engy += 0.1)
                max_sigmav_i= max(max_sigmav_i,
                                  sqrt(2*1.602e-19*engy/m[icolsp])
                                  *(sigma4(engy)+sigma5(engy)));

            col_prob_i = 1 -exp(-max_sigmav_i*gden*dt*sp_k[icolsp]);
            icol_extra = 0.5;
        }
        init_mcc_flag = 0;
    }

    /*****KOSTIL_START*******/
    int sum_rand0=0;
    int sum_rand1=0;
    int sum_rand2=0;
    int sum_rand3=0;
    int sum_high=0;

    if(ecollisional && isp==ecolsp) {
        for(j=0; j< np[ecolsp]; j++) {
            engy = Exscale[ecolsp]*vx[ecolsp][j]*vx[ecolsp][j];
            engy+= Eyscale[ecolsp]*vy[ecolsp][j]*vy[ecolsp][j];
            engy+= Ezscale[ecolsp]*vz[ecolsp][j]*vz[ecolsp][j];
            if(engy>ionengy0) {
                sum_high++;
                max_sigmav_e=sqrt(2*1.602e-19*engy/m[ecolsp])
                        *(sigma1(engy)+sigma2(engy)+sigma3(engy));

                col_prob_e = 1 -exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
                if(frand() < col_prob_e) {
                    sum_rand0++;
                    vel  = sqrt(2.0*1.602e-19*engy/m[ecolsp]);
                    sigma_total = max_sigmav_e/(vel + SNG_MIN);   /**djc**/
                    random= frand();

                    /*********************************************************************/
                    /* determine the type of collision and calculate new velocities, etc */

                    /*******************************/
                    /* if the collision is elastic */

                    if (random <= (sum_sigma = sigma1(engy))/sigma_total) {
                        sum_rand1++;
                        /* first normalize vel */
                        vx[ecolsp][j] *= vxscale[ecolsp]/vel;
                        vy[ecolsp][j] *= vyscale[ecolsp]/vel;
                        vz[ecolsp][j] *= vzscale[ecolsp]/vel;

                        /* scatter the electron */
                        newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 1);

                        /* normalize back to deminsionless vel */
                        vx[ecolsp][j] /= vxscale[ecolsp];
                        vy[ecolsp][j] /= vyscale[ecolsp];
                        vz[ecolsp][j] /= vzscale[ecolsp];
                    }
                    /**********************************/
                    /* if the collision is excitation */

                    else if (engy >= extengy0 && random <= (sum_sigma +=sigma2(engy))/sigma_total) {
                        sum_rand2++;
                        /* first normalize vel */
                        vx[ecolsp][j] *= vxscale[ecolsp]/vel;
                        vy[ecolsp][j] *= vyscale[ecolsp]/vel;
                        vz[ecolsp][j] *= vzscale[ecolsp]/vel;

                        engy -= extengy0;
                        vel   = sqrt(2.0*1.602e-19*engy/m[ecolsp]);

                        /* scatter the electron */
                        newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);

                        /* normalize back to deminsionless vel */
                        vx[ecolsp][j] /= vxscale[ecolsp];
                        vy[ecolsp][j] /= vyscale[ecolsp];
                        vz[ecolsp][j] /= vzscale[ecolsp];

                    }
                    /***************************************************************************/
                    /* if the collision is ionization, add the released electron and ion first */

                    else if(engy >= ionengy0 && random <= (sum_sigma +=sigma3(engy))/sigma_total) {
                        sum_rand3++;
                        /* first normalize vel */
                        vx[ecolsp][j] *= vxscale[ecolsp]/vel;
                        vy[ecolsp][j] *= vyscale[ecolsp]/vel;
                        vz[ecolsp][j] *= vzscale[ecolsp]/vel;

                        /********************************/
                        /* subtract the ion. energy and */
                        /* partition the remaining energy */

                        engy -= ionengy0;
                        rengy = 10.0*tan(frand()*atan(engy/20.0));
                        engy -= rengy;

                        /********************************/
                        /* scatter the created electron */

                        vel= sqrt(2.0*1.602e-19*rengy/m[ecolsp]);
                        k  = np[ecolsp];
                        vx[ecolsp][k] = vx[ecolsp][j];
                        vy[ecolsp][k] = vy[ecolsp][j];
                        vz[ecolsp][k] = vz[ecolsp][j];
                        newvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], 0);

                        /* normalize back to deminsionless vel */
                        vx[ecolsp][k] /= vxscale[ecolsp];
                        vy[ecolsp][k] /= vyscale[ecolsp];
                        vz[ecolsp][k] /= vzscale[ecolsp];

                        x[ecolsp][k] = x[ecolsp][j];
                        y[ecolsp][k] = y[ecolsp][j];

                        /****************************************/
                        /* assign velocities to the created ion */

                        k = np[ionsp];
                        maxwellv(&vx[ionsp][k], &vy[ionsp][k], &vz[ionsp][k], vgth, frand(), frand(), frand());
                        vx[ionsp][k] /= vxscale[ionsp];
                        vy[ionsp][k] /= vyscale[ionsp];
                        vz[ionsp][k] /= vzscale[ionsp];

                        x[ionsp][k] = x[ecolsp][j];
                        y[ionsp][k] = y[ecolsp][j];

                        /*****************************************/
                        /* finally scatter the incident electron */

                        vel= sqrt(2.0*1.602e-19*engy/m[ecolsp]);
                        newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);

                        /* normalize back to dimensionless vel */
                        vx[ecolsp][j] /= vxscale[ecolsp];
                        vy[ecolsp][j] /= vyscale[ecolsp];
                        vz[ecolsp][j] /= vzscale[ecolsp];

                        if(++np[ionsp] >maxnp[ionsp] || ++np[ecolsp] >maxnp[ecolsp]) {
                            puts("ADJUST(Ionization): too many particles. MUST EXIT!");
                            exit(1);
                        }
                        s1 = x[ecolsp][j];
                        s2 = y[ecolsp][j];
                        delx = x[ecolsp][j] -s1;
                        dely = y[ecolsp][j] -s2;
                        sp_n_mcc[ionsp][s1][s2]    += weight[ecolsp]*(1-delx)*(1-dely);
                        sp_n_mcc[ionsp][s1+1][s2]  += weight[ecolsp]*delx*(1-dely);
                        sp_n_mcc[ionsp][s1][s2+1]  += weight[ecolsp]*(1-delx)*dely;
                        sp_n_mcc[ionsp][s1+1][s2+1]+= weight[ecolsp]*delx*dely;

                    }
                }
            }
        }
    }
    if(sum_high>0)
      printf("sh=%d s0=%d s1=%d  s2=%d  s3=%d\n", sum_high, sum_rand0, sum_rand1, sum_rand2, sum_rand3);
    /*****KOSTIL_END*******/


    /********************************/
    /* Electron collisions with Ar  */

    if(ecollisional && isp==ecolsp) {
        ecol_extra += np[ecolsp]*col_prob_e;
        N = ecol_extra;
        ecol_extra -= N;

        nnp = np[ecolsp];
        for(j=0; j< N; j++) {
            index= nnp*frand();
            nnp--;
            temp = x[ecolsp][nnp];
            x[ecolsp][nnp] = x[ecolsp][index];
            x[ecolsp][index] = temp;

            temp = y[ecolsp][nnp];
            y[ecolsp][nnp] = y[ecolsp][index];
            y[ecolsp][index] = temp;

            temp = vx[ecolsp][nnp];
            vx[ecolsp][nnp] = vx[ecolsp][index];
            vx[ecolsp][index] = temp;

            temp = vy[ecolsp][nnp];
            vy[ecolsp][nnp] = vy[ecolsp][index];
            vy[ecolsp][index] = temp;

            temp = vz[ecolsp][nnp];
            vz[ecolsp][nnp] = vz[ecolsp][index];
            vz[ecolsp][index] = temp;
        }

        for(j=nnp; j<nnp+N; j++) {
            engy = Exscale[ecolsp]*vx[ecolsp][j]*vx[ecolsp][j];
            engy+= Eyscale[ecolsp]*vy[ecolsp][j]*vy[ecolsp][j];
            engy+= Ezscale[ecolsp]*vz[ecolsp][j]*vz[ecolsp][j];
            vel  = sqrt(2.0*1.602e-19*engy/m[ecolsp]);
            sigma_total = max_sigmav_e/(vel + SNG_MIN);   /**djc**/
            random= frand();

            /*********************************************************************/
            /* determine the type of collision and calculate new velocities, etc */

            /*******************************/
            /* if the collision is elastic */

            if (random <= (sum_sigma = sigma1(engy))/sigma_total) {
                /* first normalize vel */
                vx[ecolsp][j] *= vxscale[ecolsp]/vel;
                vy[ecolsp][j] *= vyscale[ecolsp]/vel;
                vz[ecolsp][j] *= vzscale[ecolsp]/vel;

                /* scatter the electron */
                newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 1);

                /* normalize back to deminsionless vel */
                vx[ecolsp][j] /= vxscale[ecolsp];
                vy[ecolsp][j] /= vyscale[ecolsp];
                vz[ecolsp][j] /= vzscale[ecolsp];
            }
            /**********************************/
            /* if the collision is excitation */

            else if (engy >= extengy0 && random <= (sum_sigma +=sigma2(engy))/sigma_total) {
                /* first normalize vel */
                vx[ecolsp][j] *= vxscale[ecolsp]/vel;
                vy[ecolsp][j] *= vyscale[ecolsp]/vel;
                vz[ecolsp][j] *= vzscale[ecolsp]/vel;

                engy -= extengy0;
                vel   = sqrt(2.0*1.602e-19*engy/m[ecolsp]);

                /* scatter the electron */
                newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);

                /* normalize back to deminsionless vel */
                vx[ecolsp][j] /= vxscale[ecolsp];
                vy[ecolsp][j] /= vyscale[ecolsp];
                vz[ecolsp][j] /= vzscale[ecolsp];

            }
            /***************************************************************************/
            /* if the collision is ionization, add the released electron and ion first */

            else if(engy >= ionengy0 && random <= (sum_sigma +=sigma3(engy))/sigma_total) {
                /* first normalize vel */
                vx[ecolsp][j] *= vxscale[ecolsp]/vel;
                vy[ecolsp][j] *= vyscale[ecolsp]/vel;
                vz[ecolsp][j] *= vzscale[ecolsp]/vel;

                /********************************/
                /* subtract the ion. energy and */
                /* partition the remaining energy */

                engy -= ionengy0;
                rengy = 10.0*tan(frand()*atan(engy/20.0));
                engy -= rengy;

                /********************************/
                /* scatter the created electron */

                vel= sqrt(2.0*1.602e-19*rengy/m[ecolsp]);
                k  = np[ecolsp];
                vx[ecolsp][k] = vx[ecolsp][j];
                vy[ecolsp][k] = vy[ecolsp][j];
                vz[ecolsp][k] = vz[ecolsp][j];
                newvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], 0);

                /* normalize back to deminsionless vel */
                vx[ecolsp][k] /= vxscale[ecolsp];
                vy[ecolsp][k] /= vyscale[ecolsp];
                vz[ecolsp][k] /= vzscale[ecolsp];

                x[ecolsp][k] = x[ecolsp][j];
                y[ecolsp][k] = y[ecolsp][j];

                /****************************************/
                /* assign velocities to the created ion */

                k = np[ionsp];
                maxwellv(&vx[ionsp][k], &vy[ionsp][k], &vz[ionsp][k], vgth, frand(), frand(), frand());
                vx[ionsp][k] /= vxscale[ionsp];
                vy[ionsp][k] /= vyscale[ionsp];
                vz[ionsp][k] /= vzscale[ionsp];

                x[ionsp][k] = x[ecolsp][j];
                y[ionsp][k] = y[ecolsp][j];

                /*****************************************/
                /* finally scatter the incident electron */

                vel= sqrt(2.0*1.602e-19*engy/m[ecolsp]);
                newvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], 0);

                /* normalize back to dimensionless vel */
                vx[ecolsp][j] /= vxscale[ecolsp];
                vy[ecolsp][j] /= vyscale[ecolsp];
                vz[ecolsp][j] /= vzscale[ecolsp];

                if(++np[ionsp] >maxnp[ionsp] || ++np[ecolsp] >maxnp[ecolsp]) {
                    puts("ADJUST(Ionization): too many particles. MUST EXIT!");
                    exit(1);
                }
                s1 = x[ecolsp][j];
                s2 = y[ecolsp][j];
                delx = x[ecolsp][j] -s1;
                dely = y[ecolsp][j] -s2;
                sp_n_mcc[ionsp][s1][s2]    += weight[ecolsp]*(1-delx)*(1-dely);
                sp_n_mcc[ionsp][s1+1][s2]  += weight[ecolsp]*delx*(1-dely);
                sp_n_mcc[ionsp][s1][s2+1]  += weight[ecolsp]*(1-delx)*dely;
                sp_n_mcc[ionsp][s1+1][s2+1]+= weight[ecolsp]*delx*dely;

            }
        }
    }

    /**************************************/
    /* Ar+ + Ar -> .....                  */

    if(icollisional && isp==icolsp) {
        icol_extra += np[icolsp]*col_prob_i;
        N = icol_extra;
        icol_extra -= N;

        nnp = np[icolsp];
        for(j=0; j< N; j++) {
            index= nnp*frand();
            nnp--;
            temp = x[icolsp][nnp];
            x[icolsp][nnp] = x[icolsp][index];
            x[icolsp][index] = temp;

            temp = y[icolsp][nnp];
            y[icolsp][nnp] = y[icolsp][index];
            y[icolsp][index] = temp;

            temp = vx[icolsp][nnp];
            vx[icolsp][nnp] = vx[icolsp][index];
            vx[icolsp][index] = temp;

            temp = vy[icolsp][nnp];
            vy[icolsp][nnp] = vy[icolsp][index];
            vy[icolsp][index] = temp;

            temp = vz[icolsp][nnp];
            vz[icolsp][nnp] = vz[icolsp][index];
            vz[icolsp][index] = temp;
        }

        for(j=nnp; j<nnp+N; j++) {
            maxwellv(&vneutx, &vneuty, &vneutz, vgth, frand(), frand(), frand());
            vneutx /= vxscale[ionsp];
            vneuty /= vyscale[ionsp];
            vneutz /= vzscale[ionsp];

            vx[icolsp][j] -= vneutx;
            vy[icolsp][j] -= vneuty;
            vz[icolsp][j] -= vneutz;
            engy = Exscale[icolsp]*vx[icolsp][j]*vx[icolsp][j];
            engy+= Eyscale[icolsp]*vy[icolsp][j]*vy[icolsp][j];
            engy+= Ezscale[icolsp]*vz[icolsp][j]*vz[icolsp][j];
            vel  = sqrt(2.0*1.602e-19*engy/m[icolsp]);
            sigma_total = max_sigmav_i/(vel + SNG_MIN);   /**djc**/
            random= frand();

            /**********************************/
            /* if the collision is scattering */

            if (random <= (sum_sigma =sigma5(engy))/sigma_total) {
                double up1, up2, up3, mag;
                double r11, r12, r13, r21, r22, r23, r31, r32, r33;

                coschi= sqrt(frand());
                sinchi= sqrt(fabs(1. -coschi*coschi));

                phi1  = 2*M_PI*frand();
                cosphi= cos(phi1);
                sinphi= sin(phi1);

                r13 = vx[icolsp][j]*vxscale[icolsp]/vel;
                r23 = vy[icolsp][j]*vyscale[icolsp]/vel;
                r33 = vz[icolsp][j]*vzscale[icolsp]/vel;

                if(r33 == 1.0) { up1= 0;  up2= 1;  up3= 0; }
                else           { up1= 0;  up2= 0;  up3= 1; }

                r12 = r23*up3 -r33*up2;
                r22 = r33*up1 -r13*up3;
                r32 = r13*up2 -r23*up1;
                mag = sqrt(r12*r12 + r22*r22 + r32*r32);

                r12/= mag;
                r22/= mag;
                r32/= mag;

                r11 = r22*r33 -r32*r23;
                r21 = r32*r13 -r12*r33;
                r31 = r12*r23 -r22*r13;

                vx[icolsp][j]= vel*coschi*(r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi)/vxscale[icolsp];
                vy[icolsp][j]= vel*coschi*(r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi)/vyscale[icolsp];
                vz[icolsp][j]= vel*coschi*(r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi)/vzscale[icolsp];
            }
            /***************************************/
            /* if the collision is charge exchange */

            else if (random <= (sum_sigma +=sigma4(engy))/sigma_total) {
                vx[icolsp][j] = vy[icolsp][j] = vz[icolsp][j] = 0.0;

            }
            vx[icolsp][j] += vneutx;
            vy[icolsp][j] += vneuty;
            vz[icolsp][j] += vneutz;
        }
    }
}

/**************************************************************/

void newvel(double energy, double vel, double *vx, double *vy, double *vz, int e_flag)
{
    double phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
    double mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;

    if(energy < 1e-30)  coschi = 1;
    else  coschi = (energy +2 -2*pow(energy +1,frand()))/energy;
    sinchi= sqrt(fabs(1. - coschi*coschi));

    phi1  = 2*M_PI*frand();
    cosphi= cos(phi1);
    sinphi= sin(phi1);

    if(e_flag) vel *= sqrt(1.0 - 2.0*m[ecolsp]*(1-coschi)/m[ionsp]);

    r13 = *vx;
    r23 = *vy;
    r33 = *vz;

    if(r33 == 1.0) { up1= 0;  up2= 1;  up3= 0; }
    else           { up1= 0;  up2= 0;  up3= 1; }

    r12 = r23*up3 -r33*up2;
    r22 = r33*up1 -r13*up3;
    r32 = r13*up2 -r23*up1;
    mag = sqrt(r12*r12 + r22*r22 + r32*r32);

    r12/= mag;
    r22/= mag;
    r32/= mag;

    r11 = r22*r33 -r32*r23;
    r21 = r32*r13 -r12*r33;
    r31 = r12*r23 -r22*r13;

    *vx= vel*(r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi);
    *vy= vel*(r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi);
    *vz= vel*(r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi);
}

/***********************************/
/*  e + Ar -> e + Ar  Elastic      */

double sigma1(double energy)
{
    int i;
    double engy, els_sigma, alpha;
    static int init_flag=1;
    static double *elastic1, *elastic2;

    /******  Initialization   *********/
    if(init_flag) {
        double engy;
        elastic1 = (double *) malloc(101*sizeof(double));
        elastic2 = (double *) malloc(NEMAX*sizeof(double));

        /*****  With no Ramseur min.  ******/
        /*
       for(i=0; i<100; i++) elastic1[i]= 1.4e-19;
       for(i=0; i<NEMAX; i++)
       {
       engy = i;
       if(engy <= 20.) elastic2[i]= 1.4e-19;
       else
       elastic2[i]= 9.07e-19*pow(engy, 1.55)*pow(engy+70.0, 1.10)/pow(14.+engy, 3.25);
       }
       */

        /*****  With the Ramseur min.  *****/
        for(i=0; i<101; i++) {
            engy = .01*i;
            if(engy < 0.2) elastic1[i]= 1./pow(10.0, 19.0 +engy/.11);
            else elastic1[i]= 9.07e-19*pow(engy, 1.55)*pow(engy+70.0, 1.10)/pow(14.+engy, 3.25);
        }
        for(i=0; i<NEMAX; i++) {
            engy = i;
            elastic2[i]= 9.07e-19*pow(engy, 1.55)*pow(engy+70.0, 1.10)/pow(14.+engy, 3.25);
        }
        init_flag =0;
    }
    /****************************/

    if(energy < 1.0) {
        i= 100*energy;
        alpha= 100*energy -i;
        els_sigma = elastic1[i] + alpha*(elastic1[i+1] - elastic1[i]);
    }
    else {
        i= energy;
        if(i < NEMAX-1) {
            alpha= energy -i;
            els_sigma = elastic2[i] + alpha*(elastic2[i+1] - elastic2[i]);
        }
        else
            els_sigma = elastic2[NEMAX-1];
    }
    return(els_sigma);
}

/************************************/
/*  e + Ar -> e + Ar  Excitation    */

double sigma2(double energy)
{
    int i;
    double engy, exc_sigma, alpha;
    static int init_flag=1;
    static double *excit;

    /******  Initialization   *********/
    if(init_flag)
    {
        double engy;
        excit = (double *) malloc(NEMAX*sizeof(double));

        for(i=0; i<NEMAX; i++) {
            engy = i;
            if(engy < 12.0) excit[i] = 0;
            else
                excit[i] = (3.85116e-19*log(engy/3.4015) -4.85227e-19)/engy;
        }
        init_flag =0;
    }
    /****************************/

    i= energy;
    if(i < NEMAX-1) {
        alpha= energy -i;
        exc_sigma = excit[i] + alpha*(excit[i+1] - excit[i]);
    }
    else
        exc_sigma = excit[NEMAX-1];

    return(exc_sigma);
}

/************************************/
/*  e + Ar -> e + e + Ar+  Ion.     */

double sigma3(double energy)
{
    int i;
    double engy, ion_sigma, alpha;
    static int init_flag=1;
    static double *ioniz;

    /******  Initialization   *********/
    if(init_flag) {
        double engy;
        ioniz = (double *) malloc(NEMAX*sizeof(double));

        for(i=0; i<NEMAX; i++) {
            engy = i;
            if(engy < 15.76) ioniz[i] = 0;
            else
                ioniz[i] = (1.3596e-18/engy)*log((engy +120.0/engy)/15.76)*
                        (atan((engy*engy -9.76*engy +2.4)/(20.6*engy +206)) +
                         atan((2*engy -80.0)/(10.3*engy +103.0)));
        }
        init_flag =0;
    }
    /****************************/

    i= energy;
    if(i < NEMAX-1) {
        alpha= energy -i;
        ion_sigma = ioniz[i] + alpha*(ioniz[i+1] - ioniz[i]);
    }
    else
        ion_sigma = ioniz[NEMAX-1];

    return(ion_sigma);
}

/************************************/
/*  Ar + Ar+ -> Ar+ + Ar  Charge X  */

double sigma4(double energy)
{
    if(energy > 4.0) return(2.0e-19 +5.5e-19/sqrt(energy));
    else             return(-2.95e-19*sqrt(energy) +10.65e-19);
}

/***********************************/
/*  Ar + Ar+ -> Ar + Ar+   Scat.   */

double sigma5(double energy)
{
    if(energy > 4.0) return(1.8e-19 +4.0e-19/sqrt(energy));
    else             return(-2.0e-19*sqrt(energy) +7.8e-19);
}

/**************************************************************/

int makethefile()
{
    double e;
    FILE *DMPFile;

    DMPFile = fopen("arxsects.dat", "w");

    for(e=0; e<10; e+=0.01)
        fprintf(DMPFile, "%e %e %e %e %e %e\n", e+1e-30, sigma1(e)+1e-30,
                sigma2(e)+1e-30, sigma3(e)+1e-30, sigma4(e)+1e-30, sigma5(e)+1e-30);

    for(e=10; e<100; e+=0.1)
        fprintf(DMPFile, "%e %e %e %e %e %e\n", e+1e-30, sigma1(e)+1e-30,
                sigma2(e)+1e-30, sigma3(e)+1e-30, sigma4(e)+1e-30, sigma5(e)+1e-30);

    fclose(DMPFile);
}
}
/**************************************************************/
