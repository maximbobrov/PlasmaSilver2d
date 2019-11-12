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
    if(m_field != nullptr)
    {
        double dx = m_pData->getDx();
        double dy = m_pData->getDy();
        double dt = m_pData->getDt();
        simulationData::simulationParameters* pParams = m_pData->getParameters();
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();
            setBc();

            for (int i = 1; i < m_field->cellsX-1; ++i)
            {
                for (int j = 1; j < m_field->cellsY-1; ++j)
                {
                    double D = m_aNu[i][j];
                    double a = 1.0/dt+D*((2.0)/(dx*dx)+(2.0)/(dy*dy));
                    double bp = -D/(dx*dx);
                    double bm = bp;
                    double cp = -D/(dy*dy);
                    double cm = cp;
                    m_field->arr[i][j]=((m_aRHS[i][j]-(bp*m_field->arr[i+1][j]+bm*m_field->arr[i-1][j]+cp*m_field->arr[i][j+1]+cm*m_field->arr[i][j-1]))
                            /a)*m_mask[i][j]+m_maskValue[i][j];
                }
            }

        }
    }
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
    for (int j = 0; j < m_field->cellsY; ++j)
    {
        m_field->arr[0][j] =0.0;// m_field->arr[m_field->cellsX-2][j];
        m_field->arr[m_field->cellsX-1][j] =0.0;// m_field->arr[1][j];
        m_field->arrPrev[0][j] =0.0;// m_field->arrPrev[m_field->cellsX-2][j];
        m_field->arrPrev[m_field->cellsX-1][j] =0.0;// m_field->arrPrev[1][j];
    }
    for (int i = 0; i < m_field->cellsX; ++i)
    {
        m_field->arr[i][0] =0.0;// m_field->arr[i][m_field->cellsY-2];
        m_field->arr[i][m_field->cellsY-1] = 0.0;//m_field->arr[i][1];
        m_field->arrPrev[i][0] = 0.0;//= m_field->arrPrev[i][m_field->cellsY-2];
        m_field->arrPrev[i][m_field->cellsY-1]= 0.0;// = m_field->arrPrev[i][1];
    }

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
    m_mask = m_pData->getParameters()->arrMaskNe;
    m_maskValue= m_pData->getParameters()->arrMaskNeValue;
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }

}

double solverNe::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();


   // m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    double** R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double mult=pParams->p/(pParams->T*8.314);
    double ** phi  = m_pData->getFieldPhi()->arr;
    double dx=m_pData->getDx();
    double dy=m_pData->getDy();
    double dt=m_pData->getDt();
    double mu=1.0;
    for (int i = 1; i < m_field->cellsX-1; ++i)
    {
        for (int j = 1; j < m_field->cellsY-1; ++j)
        {
            int mm=(m_field->arr[i][j]>1e-5);
            m_aRHS[i][j]= mult*R1_Ar_e[i][j] //*m_field->arr[i][j]+
                    +m_field->arrPrev[i][j]/dt
                    + mm*pParams->arrMue[i][j]*(pParams->arrEx[i][j] * simulationTools::ddxUpWind(m_field->arr, m_field->cellsX, dx, i, j, pParams->arrEx[i][j])
                                                + pParams->arrEy[i][j] * simulationTools::ddyUpWind(m_field->arr, m_field->cellsY, dy, i, j, pParams->arrEy[i][j])
                                                + m_field->arr[i][j] * pParams->arrLaplPhi[i][j])
                    + simulationTools::ddxCentral(m_field->arr,m_field->cellsX, dx, i, j) * simulationTools::ddxCentral(pParams->arrDe,m_field->cellsX, dx, i, j)
                    + simulationTools::ddyCentral(m_field->arr,m_field->cellsY, dy, i, j) * simulationTools::ddyCentral(pParams->arrDe,m_field->cellsY, dy, i, j);
        }
    }
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
    m_mask = m_pData->getParameters()->arrMaskEnergy;
    m_maskValue= m_pData->getParameters()->arrMaskEnergyValue;
    for (int i = 0; i < m_field ->cellsX; ++i) {
        m_aRHS[i] = new double [m_field ->cellsY];
    }
}

double solverEnergy::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField2d* pNe = m_pData->getFieldNe();
    //double dz = m_pData->getDz();


   // m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);

    double** R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double de=m_pData->getReactionDe(simulationData::ReactionName::comsol_eAr_2eArp);

    double mult=fabs(m_pData->q*de*pParams->p/(pParams->T*8.314));
    qDebug()<<"mult="<<mult;
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
    double dx=m_pData->getDx();
    double dy=m_pData->getDy();
    double dt=m_pData->getDt();
    double mu=1.0;
    for (int i = 1; i < m_field->cellsX-1; ++i)
    {
        for (int j = 1; j < m_field->cellsY-1; ++j)
        {
            int mm=(m_field->arr[i][j]>1e-5);
            int mm2=(pNe->arr[i][j]>1e-5);
            double electronFluxX = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEx[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddxCentral(pNe->arr, m_field->cellsX, dx, i, j));
            double electronFluxY = 0.5*(- mm2 * pParams->arrMue[i][j] * pParams->arrEy[i][j] * pNe->arr[i][j] - pParams->arrDe[i][j]*simulationTools::ddyCentral(pNe->arr, m_field->cellsY, dy, i, j));
            m_aRHS[i][j]= mult*R1_Ar_e[i][j]//*pNe->arr[i][j]+
                    +m_field->arrPrev[i][j]/dt
                    + mm*pParams->arrMueps[i][j]*(pParams->arrEx[i][j] * simulationTools::ddxUpWind(m_field->arr, m_field->cellsX, dx, i, j, pParams->arrEx[i][j])
                                                  + pParams->arrEy[i][j] * simulationTools::ddyUpWind(m_field->arr, m_field->cellsY, dy, i, j, pParams->arrEy[i][j])
                                                  + m_field->arr[i][j] * pParams->arrLaplPhi[i][j])
                    + simulationTools::ddxCentral(m_field->arr,m_field->cellsX, dx, i, j) * simulationTools::ddxCentral(pParams->arrDeps,m_field->cellsX, dx, i, j)
                    + simulationTools::ddyCentral(m_field->arr,m_field->cellsY, dy, i, j) * simulationTools::ddyCentral(pParams->arrDeps,m_field->cellsY, dy, i, j)
                    - electronFluxX * pParams->arrEx[i][j] - electronFluxY * pParams->arrEy[i][j];
        }
    }
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
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double q_e= 1.6*10e-19;
    double eps_0=8.85*10e-12;
    double q_over_eps0=q_e/eps_0;

    double** pNe = m_pData->getFieldNe()->arr;
    double** pArp = m_pData->getFieldHeavySpicies(0)->arr;

    double mult= (pParams->p*6.022e23)/(pParams->T*8.314);
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();

    for (int i=1; i<m_field ->cellsX-1; i++)
    {
        for (int j=1; j<m_field ->cellsY-1; j++)
        {
            m_aRHS[i][j]= /*q_over_eps0*0.5*((pNe[i][j]+pNe[i+1][j])+mult*(pArp[i][j]+pArp[i+1][j]))*/
                    - ((m_field->arrPrev[i+1][j]-m_field->arrPrev[i-1][j])*(pParams->arrEps[i+1][j]-pParams->arrEps[i-1][j])/(4*dx*dx)
                    +(m_field->arr[i][j+1]-m_field->arr[i][j])*(pParams->arrEps[i][j+1]-pParams->arrEps[i][j])/(dy*dy))/pParams->arrEps[i][j];
            pParams->arrLaplPhi[i][j] = m_aRHS[i][j];
        }
    }
    return 1;
}

void solverPhi::setBc()
{

}

double solverPhi::solve(int iNumberIteration)
{

    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();

    for (int i = 0; i < iNumberIteration; ++i)
    {

        for (int i=0;i<m_field ->cellsX;i++)
        {
            m_field->arr[i][m_field ->cellsY-1]=0.0;
            m_field->arr[i][0]=0.0;
            //if ((i<20)&&(i>0)) m_field->arr[i][m_field ->cellsY-1]=-1500.0;
            if ((i>=60)&&(i<m_field ->cellsX-1)) m_field->arr[i][0]=+1500;
        }

        getRhs();

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

        par.w_bc_type=0;
        par.e_bc_type=0;
        par.n_bc_type=3;
        par.s_bc_type=3;

        par.w_bc_val=0.0;
        par.e_bc_val=0.0;
        par.n_bc_val=0.0;
        par.s_bc_val=0.0;
        par.dx = dx;
        par.dy = dy;

        for (int i=1;i<10;i++)
        {
            multigrid_N(par, m_field->arr, m_aRHS, pParams->arrMaskPhi, pParams->arrMaskPhiValue, 8, 3);
        }

    }

    return -1;
}

solverHeavySpicies::solverHeavySpicies(simulationData *pData, int num)
{
    m_pData = pData;
    m_field = m_pData->getFieldHeavySpicies(num);
    m_charge = m_pData->getHeavySpiciesCharge(num);
    m_aNu = m_pData->getParameters()->arrDomega;
    m_aRHS = new double*[m_field ->cellsX];
    m_mask = m_pData->getParameters()->arrMaskHeavy;
    m_maskValue= m_pData->getParameters()->arrMaskHeavyValue;
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
    double dx = m_pData->getDx();
    double dy = m_pData->getDy();
    double dt = m_pData->getDt();

  //  m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    double** R1_Ar_e=m_pData->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    double mult=pParams->p / ( pParams->T * 8.314 * 6.022e23 * pParams->rho);


    for (int i=1; i<m_field ->cellsX-1; i++)
    {
        for (int j=1; j<m_field ->cellsY-1; j++)
        {
            m_aRHS[i][j]= mult*R1_Ar_e[i][j]//*pNe->arr[i][j]+
                  + m_field->arrPrev[i][j]/dt +
                         -  m_charge  * pParams->arrMuomega[i][j] * (pParams->arrEx[i][j] * simulationTools::ddxCentral(m_field->arr, m_field ->cellsX, dx, i, j)
                                                                   + pParams->arrEy[i][j] * simulationTools::ddyCentral(m_field->arr, m_field ->cellsY, dy, i, j)
                                                                   + m_field->arr[i][j] * simulationTools::ddxCentral(pParams->arrEx, m_field ->cellsX, dx, i, j)
                                                                   + m_field->arr[i][j] * simulationTools::ddyCentral(pParams->arrEy, m_field ->cellsY, dy, i, j));
        }
    }
    return 1;
}


