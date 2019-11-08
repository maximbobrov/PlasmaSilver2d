#ifndef REACTION_H
#define REACTION_H

//#include "simulationdata.h"
#include "crosssection.h"
class simulationData;

class reaction
{
public:
    reaction(simulationData* data);// = nullptr);
    virtual void calc();
    virtual double getDe();
    double** getR();
protected:
    double** m_R;
    simulationData* m_pData;
};

//без нижних подчеркиваний пиздец непонятно, где левая часть реакции, а где правая
class reactionEAr_EAr:public reaction //e+Ar=>e+Ar  ELASTIC
{
public:
    reactionEAr_EAr(simulationData* data);// = nullptr);
    virtual void calc();
private:
    double m_massRatio = 0.136E-04; //m_e/m_Ar
    crossSection * m_cs;
};


class reactionEAr_EArs:public reaction //e+Ar=>e+Ars  EXCITATION
{
public:
    reactionEAr_EArs(simulationData* data);// = nullptr);
    virtual void calc();
private:
    double m_energy =  -11.50 ; //excitation energy eV
    double m_wRat = 12;  // statistical weight ratio of initial state to the excited state
    int m_Detailed = 1; //use detailed balance (if 0-- otherwise)
    crossSection * m_cs;
};


class reactionEAr_2EArp:public reaction //e+Ar=>2e+Ar+  Ionization
{
public:
    reactionEAr_2EArp(simulationData* data);// = nullptr);
    virtual void calc();
     virtual double getDe();
private:
    double m_energy =  -15.80;// eV threshold energy
    crossSection * m_cs;
};



class reactionEAr_2EArp_comsol:public reaction //e+Ar=>2e+Ar+  Ionization
{
public:
    reactionEAr_2EArp_comsol(simulationData* data);// = nullptr);
    virtual void calc();
     virtual double getDe();
private:
    double m_energy =  -15.80;// eV threshold energy
    splineInterp * m_spline;
};

class reactionEArs_2EArp:public reaction //e+Ars=>2e+Ar+  Ionization
{
public:
    reactionEArs_2EArp(simulationData* data);// = nullptr);
    virtual void calc();
private:
    double m_energy =  -4.427;// eV threshold energy
    crossSection * m_cs;
};


#endif // REACTION_H
