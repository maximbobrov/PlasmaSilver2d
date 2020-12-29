#ifndef PZSOLVER_H
#define PZSOLVER_H


#define w_x0 0
#define w_x1 0.08


#define w_y0 0.0005
#define w_y1 0.006

#define w_z0 -1.0e-7
#define w_z1 1.0e-7

#define N_Y_DIEL 64

#define dl_pz  0.0015
#define eps_pz 100


#define qe -1.6e-19
#define Me 9.11e-31
#define eps0 8.85e-12

#define pi4 12.5663706144
#define pi2 6.28318530718

struct vec2 {
    vec2(): x(0.0), y(0.0), charge(0.0){}
    vec2(double ix, double iy, double icharge) : x(ix), y(iy), charge(icharge){}
    double x;
    double y;
    double charge;
};

typedef struct
{
    double a,bm,bp,cm,cp;

    int w_bc_type,e_bc_type,n_bc_type,s_bc_type; //0--fixed value, 1--fixed gradient,2 --cyclic; 3 --init
    int w_bc_val,e_bc_val,n_bc_val,s_bc_val;

}INPUT_PARAM;

typedef struct
{
    int order;

    double C[8]; //C0*x + C1*x^2 + C2*x^3..

}poly;

class pzSolver
{
public:
    struct pElem   //linear electrode elem
    {
        double p,p_prev, ds, dl,E,E_elec,E_prev,RHS;//C/m2, m, V
        vec2 r;
        double q, q_ext; //q is dipolar charge; q_ext is the attached charge
    };

public:
    poly m_poly;
    vec2* m_rCentre;
    INPUT_PARAM m_par;

    double m_dx;
    double m_dt;
    double kappa;
    pElem *m_p;
    int m_p_num;
    pzSolver();
    void init();
    void solvePz(int itn);
    void getRHS();
    void get_q();
    void step();
    vec2 getEdepol(double x, double y);
    double getPhidepol(double x, double y);
    void updateCharge();
    void updateGridProp();
    void setWallPos(double a);

   // static vec2 getEField(const vec2& iCenterPos, const vec2& iFarPos);
   // static vec2 getPhiField(const vec2& iCenterPos, const vec2& iFarPos);
};

#endif // PZSOLVER_H
