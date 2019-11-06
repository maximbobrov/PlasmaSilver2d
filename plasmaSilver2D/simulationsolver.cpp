#include "simulationsolver.h"
#include "simulationtools.h"
#include "phi_mult.h"
#include "math.h"
#include <QDebug>


simulationSolver::simulationSolver(simulationData* ipData/*= nullptr*/)
{
    m_pData = ipData;
    m_aRHS = nullptr;
    m_field = nullptr;
    m_aNu = nullptr;
}

simulationSolver::~simulationSolver()
{

}

double simulationSolver::solve(int iNumberIteration)
{
    double res = 0.0;
    /*if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dt = m_pData->getDt();
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();
            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = (dt * m_aNu[j] * (m_field->arr[j+1] + m_field->arr[j-1]) / (dz * dz) +  m_aRHS[j] * dt + m_field->arrPrev[j]) / (1.0 + 2.0 * m_aNu[j] * dt / (dz * dz));
            }
            setBc();

        }
        for (int j=1; j < m_field->cellsNumber - 1; j++)
        {
            res += (m_field->arr[j] - m_field->arrPrev[j]) / dt - (m_aNu[j] * (m_field->arr[j+1] - 2.0 * m_field->arr[j] + m_field->arr[j-1]) / (dz * dz) + m_aRHS[j]);
        }
        return res;
    }*/
    return -1;
}


double simulationSolver::getRhs()
{
 return -1;
}

void simulationSolver::getStepEuler()
{
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        for (int j = 0; j < m_field->cellsY; ++j)
        {
         m_field->arrPrev[i][j] = m_field->arr[i][j];
        }
    }
}

void simulationSolver::setBc()
{
    //periodic bc initally
    /*double tmp=m_field->arr[1];
    int last=m_field->cellsNumber - 1;
    m_field->arr[0]=m_field->arr[m_field->cellsNumber - 2];
    m_field->arr[last]=tmp;*/

}

double simulationSolver::init(double value)
{
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        for (int j = 0; j < m_field->cellsY; ++j)
        {
           m_field->arr[i][j] = value;
           m_field->arrPrev[i][j] = value;
        }
    }
}

solverNe::~solverNe()
{

}

solverNe::solverNe(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldNe();
    m_aNu = m_pData->getParameters()->arrDe;
    m_aRHS = new double*[m_field ->cellsX];
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }

}

double solverNe::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();


    m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    double* R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double mult=pParams->p/(pParams->T*8.314);

    //double dz = m_pData->getDz();
   /* for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i]= 0.25*mult*R1_Ar_e[i]*m_field->arr[i]
                   + pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);
                   ;//+ simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDe, m_field ->cellsNumber, dz, i);
    }*/
    return 1;
}

solverEnergy::~solverEnergy()
{

}

solverEnergy::solverEnergy(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldEnergy();
    m_aNu = m_pData->getParameters()->arrDeps;
    m_aRHS = new double*[m_field ->cellsX];
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

double solverEnergy::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField2d* pNe = m_pData->getFieldNe();
    //double dz = m_pData->getDz();


    m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);

    double* R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double de=m_pData->getReactionDe(simulationData::ReactionName::comsol_eAr_2eArp);

    double mult=m_pData->q*de*pParams->p/(pParams->T*8.314);

    /*for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*simulationTools::ddzCentral(pNe->arr, m_field ->cellsNumber, dz, i);
        //(pNe->arr[i+1]-pNe->arr[i]);
        m_aRHS[i]= mult*R1_Ar_e[i]*pNe->arr[i]
                   + pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                   + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i)
                   - electronFlux * pParams->arrE[i];
    }*/
    return 1;
}

solverPhi::~solverPhi()
{

}

solverPhi::solverPhi(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldPhi();
    m_aRHS = new double*[m_field ->cellsX];
    for (int i = 0; i < m_field ->cellsX; ++i)
    {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

double solverPhi::getRhs()
{
    /*double q_e= 1.6*10e-19;
    double eps_0=8.85*10e-12;
    double q_over_eps0=q_e/eps_0;

    double** pNe = m_pData->getFieldNe()->arr;
    double** pArp = m_pData->getFieldHeavySpicies(0)->arr;

    simulationData::simulationParameters* pParams=m_pData->getParameters();
    double mult= (pParams->p*6.022e23)/(pParams->T*8.314);*/
            //m_fHeavy[j]->arr[i] =(m_fNe->arr[i]*pParams->T*8.314)/(pParams->p*6.022e23);

   /* for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i] = q_over_eps0*0.5*((pNe[i]+pNe[i+1])+mult*(pArp[i]+pArp[i+1]));
    }*/


    return 1;
}

void solverPhi::setBc()
{

}

double solverPhi::solve(int iNumberIteration)
{
    double res = 0.0;
   /* if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dzdz=dz*dz;
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();

            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = 0.5* ((m_field->arr[j+1] + m_field->arr[j-1]) +  m_aRHS[j]*dzdz);
            }
            setBc();
        }
        return res;
    }*/

    double dx = m_pData->getDx();
    double dy = m_pData->getDy();

    for (int i=1; i<m_field ->cellsX-1; i++)
    {
        for (int j=1; j<m_field ->cellsY-1; j++)
        {
            m_aRHS[i][j]= 0.0;//q_e*(-n_1[i][j]+n_2[i][j])/eps_0;
        }
    }

    for (int i=0;i<m_field ->cellsX;i++)
    {
        m_field->arr[i][m_field ->cellsY-1]=0.0;
        m_field->arr[i][0]=0.0;
        if ((i<20)&&(i>0)) m_field->arr[i][m_field ->cellsY-1]=-1500;
        if ((i>=60)&&(i<m_field ->cellsX-1)) m_field->arr[i][0]=1500;
    }

    for (int nnn=0;nnn<13;nnn++)
    {
        for (int i=1;i<m_field ->cellsX-1;i++)
        {
            m_field->arr[i][m_field ->cellsY-1]=(m_field->arr[i][m_field ->cellsY-1]+m_field->arr[i-1][m_field ->cellsY-1]+m_field->arr[i+1][m_field ->cellsY-1])/3.0;//-A;//.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
            m_field->arr[i][0]=(m_field->arr[i-1][0]+m_field->arr[i][0]+m_field->arr[i+1][0])/3.0;//A*sin(i*6*M_PI/(N_X-1)+1800.0*alpha*t*dt);
        }
    }

    INPUT_PARAM par;
    par.a=((2.0)/(dx*dx)+2.0/(dy*dy));
    par.bp=-1.0/(dx*dx);
    par.bm=-1.0/(dx*dx);
    par.cp=-1.0/(dy*dy);
    par.cm=-1.0/(dy*dy);

    par.w_bc_type=2;
    par.e_bc_type=2;
    par.n_bc_type=3;
    par.s_bc_type=3;

    par.w_bc_val=0.0;
    par.e_bc_val=0.0;
    par.n_bc_val=0.0;
    par.s_bc_val=0.0;
    par.dx = dx;
    par.dy = dy;

    //jacobi(par,phi_,div_,20);//phi
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);
    multigrid_N(par, m_field->arr, m_aRHS, 8, 3);

    return -1;
}

solverHeavySpicies::solverHeavySpicies(simulationData *pData, int num)
{
    m_pData = pData;
    m_field = m_pData->getFieldHeavySpicies(num);
    m_charge = m_pData->getHeavySpiciesCharge(num);
    m_aNu = m_pData->getParameters()->arrDomega;
    m_aRHS = new double*[m_field ->cellsX];
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

solverHeavySpicies::~solverHeavySpicies()
{

}

double solverHeavySpicies::getRhs()
{

    simulationData::simulationField2d* pNe = m_pData->getFieldNe();
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    //double dz = m_pData->getDz();

    m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    double* R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double mult=pParams->p / ( pParams->T * 8.314 * 6.022e23 * pParams->rho);


    /*for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i]= 0.02* mult*R1_Ar_e[i]*pNe->arr[i]
                -  m_charge  * pParams->arrMuomega[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   -  m_charge  * pParams->arrMuomega[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);
    }*/
    return 1;
}


