#include "reaction.h"
#include "simulationdata.h"
#include <QDebug>

static const int g_EAr_EAr_num =67;

static  double g_EAr_EAr[g_EAr_EAr_num][2] = {
    #include "reactions/e+Ar_e+Ar.txt"
};

static const int g_EAr_EArs_num =34;

static double g_EAr_EArs[g_EAr_EArs_num][2] = {
    #include "reactions/e+Ar_e+Ars.txt"
};

static const int g_EAr_2EArp_num =31;

static double g_EAr_2EArp[g_EAr_2EArp_num][2] = {
    #include "reactions/e+Ar_2e+Ar+.txt"
};

static const int g_EArs_2EArp_num =20;

static double g_EArs_2EArp[g_EArs_2EArp_num][2] = {
    #include "reactions/e+Ars_2e+Ar+.txt"
};


static const int comsol_EAr_2EArp_num =101;

static double comsol_EAr_2EArp[comsol_EAr_2EArp_num][2] = {
    #include "comsol/e+Ar_2e+Ar+.txt"
};




reaction::reaction(simulationData* data)
{
    m_pData=data;
    if (data != nullptr)
    {
        int cx=data->getCellsXNumber();
        int cy=data->getCellsYNumber();
        m_R=new double*[cx];
        for (int i = 0; i < cx; ++i)
        {
            m_R[i] = new double[cy];
        }
        for (int i = 0 ; i < cx ; i++)
        {
            for (int j = 0 ; j < cy ; j++)
            {
                m_R[i][j] = 0;
            }
        }
    }
}

void reaction::calc()
{
    //just a placeholder here
}

double reaction::getDe()
{
    return 0.0;
}

double** reaction::getR()
{
    return m_R;
}

reactionEAr_EAr::reactionEAr_EAr(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EAr_EAr_num,12);
        m_cs->fillSigmas2(g_EAr_EAr,g_EAr_EAr_num); //fill from static arrays (safer then messing up with external txt files)
    }
}

void reactionEAr_EAr::calc()
{
    double** En=m_pData->getFieldEnergy()->arr;
    double** Ne=m_pData->getFieldNe()->arr;
    double** Ar=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar)->arr;

    double N=m_pData->getN();

    for (int i=0;i<m_pData->getCellsXNumber();i++)
    {
        for (int j=0;j<m_pData->getCellsYNumber();j++)
        {
            int a=m_pData->getCellsXNumber();
            m_R[i][j]= 0.0;//m_cs->getSpline(i*100.0/m_pData->getCellsNumber());
            //m_cs->getSpline(En[i])*N*Ar[i]*Ne[i]*0.0;
        }
    }
}

reactionEAr_EArs::reactionEAr_EArs(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EAr_EArs_num,12);
        m_cs->fillSigmas2(g_EAr_EArs,g_EAr_EArs_num);
    }
}

void reactionEAr_EArs::calc()
{
    double** En=m_pData->getFieldEnergy()->arr;
    double** Ne=m_pData->getFieldNe()->arr;
    double** Ar=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar)->arr;

    double N=m_pData->getN();

    for (int i=0;i<m_pData->getCellsXNumber();i++)
    {
        for (int j=0;j<m_pData->getCellsYNumber();j++)
        {
            m_R[i][j]= 0.0;//m_cs->getSpline(i*100.0/m_pData->getCellsNumber());
            //m_cs->getSpline(En[i])*N*Ar[i]*Ne[i]*0.0;
        }
    }
}

reactionEAr_2EArp::reactionEAr_2EArp(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EAr_2EArp_num,12);
        m_cs->fillSigmas2(g_EAr_2EArp,g_EAr_2EArp_num);
    }
}

void reactionEAr_2EArp::calc()
{
    //<<m_pData->getFieldEnergy()->name;
    double** En=m_pData->getFieldEnergy()->arr;
    double** Ne=m_pData->getFieldNe()->arr;
    double** Ar=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar)->arr;

    double N=m_pData->getN();
    for (int i=0;i<m_pData->getCellsXNumber();i++)
    {
        for (int j=0;j<m_pData->getCellsYNumber();j++)
        {
            m_R[i][j] = m_cs->getSpline(i*100.0/m_pData->getCellsXNumber());
            //m_cs->getSpline(En[i])*N*Ar[i]*Ne[i]*0.0;
        }
    }
}

double reactionEAr_2EArp::getDe()
{
    return m_energy;
}

reactionEArs_2EArp::reactionEArs_2EArp(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(36,g_EArs_2EArp_num,12);
        m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}

void reactionEArs_2EArp::calc()
{


    double** En=m_pData->getFieldEnergy()->arr;
    double** Ne=m_pData->getFieldNe()->arr;
    double** Ars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    double N=m_pData->getN();

    for (int i=0;i<m_pData->getCellsXNumber();i++)
    {
        for (int j=0;j<m_pData->getCellsYNumber();j++)
        {
          m_R[i][j] = 0.0;//m_cs->getSpline(i*100.0/m_pData->getCellsNumber());
          //m_cs->getSpline(En[i])*N*Ars[i]*Ne[i]*0.0;
        }
    }
}

reactionEAr_2EArp_comsol::reactionEAr_2EArp_comsol(simulationData *data):reaction(data)
{
    if (data != nullptr)
    {
        m_spline=new splineInterp(20);

        m_spline->fillData(comsol_EAr_2EArp,comsol_EAr_2EArp_num);
        //m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}

void reactionEAr_2EArp_comsol::calc()
{
    //qDebug()<<" calccc!";
    double ** arrTe = m_pData->getArrTe();
    for (int i=0;i<m_pData->getCellsXNumber();i++)
    {
        for (int j=0;j<m_pData->getCellsYNumber();j++)
        {
         m_R[i][j] =m_pData->getFieldNe()->arr[i][j]*0.0000001*m_spline->getSpline(arrTe[i][j]);//i*100.0/m_pData->getCellsNumber());
        //m_cs->getSpline(En[i])*N*Ar[i]*Ne[i]*0.0;

        // qDebug()<<"i="<<i<<" R="<<m_R[i];
        }
    }
}

double reactionEAr_2EArp_comsol::getDe()
{
    return m_energy;
}
