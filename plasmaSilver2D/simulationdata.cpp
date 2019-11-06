#include "simulationdata.h"
#include "simulationtools.h"
#include "string.h"
#include "reaction.h"
#include <QDebug>
#include <math.h>

simulationData::simulationData(int iCellsX,int iCellsY)
{
    m_defaultCellsNumber = iCellsX;
    m_fieldNe = new simulationField2d(iCellsX, iCellsY,"electrons");
    m_fieldEnergy = new simulationField2d(iCellsX, iCellsY,"energy");
    m_fieldPhi = new simulationField2d(iCellsX+1, iCellsY+1,"potential");
    m_fieldsHeavySpecies.push_back(new simulationField2d(iCellsX, iCellsY,"Ar+"));
    m_chargeHeavySpecies.push_back(0);
 //   m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ars"));
  //  m_chargeHeavySpecies.push_back(1);
  //  m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar+"));
   // m_chargeHeavySpecies.push_back(1);
    m_numberHeavySpicies = m_fieldsHeavySpecies.size();
    m_params = new simulationData::simulationParameters(iCellsX,iCellsY);

    m_reactions.push_back(new reactionEAr_EAr(this));
    m_reactions.push_back(new reactionEAr_EArs(this));
    m_reactions.push_back(new reactionEAr_2EArp(this));
    m_reactions.push_back(new reactionEArs_2EArp(this));
    m_reactions.push_back(new reactionEAr_2EArp_comsol(this));
}

simulationData::~simulationData()
{

}

void simulationData::setDt(double idt)
{
    m_dt = idt;
}

void simulationData::setDx(double dx)
{
    m_dx = dx;
}

void simulationData::setDy(double dy)
{
    m_dy = dy;
}

void simulationData::setCellsNumber(int iCellsNumber)
{
    m_defaultCellsNumber = iCellsNumber;
}

double simulationData::getDt()
{
    return m_dt;
}

double simulationData::getDx()
{
    return m_dx;
}

double simulationData::getDy()
{
    return m_dy;
}

simulationData::simulationField2d* simulationData::getFieldNe()
{
    return m_fieldNe;
}

double **simulationData::getArrTe()
{
   return m_params->arrTe;
}


simulationData::simulationField2d* simulationData::getFieldHeavySpicies(int num)
{
    return m_fieldsHeavySpecies[num];
}

int simulationData::getHeavySpiciesCharge(int num)
{
    return m_chargeHeavySpecies[num];
}

int simulationData::getNumberHeavySpicies()
{
    return m_numberHeavySpicies;
}

double simulationData::getN()
{
   // pv=nkT
   // n/v=p/kT;
    double pres=101505;
    double T=300;
    double k=1.38e-23;


    return pres/(k*T);
}

double *simulationData::getReactionRate(simulationData::ReactionName reactName)
{
    return m_reactions[reactName]->getR();
}

double simulationData::getReactionDe(simulationData::ReactionName reactName)
{
     return m_reactions[reactName]->getDe();
}

void simulationData::calcReaction(simulationData::ReactionName reactName)
{
    m_reactions[reactName]->calc();
}

simulationData::simulationField2d* simulationData::getFieldEnergy()
{
    return m_fieldEnergy;
}

simulationData::simulationField2d* simulationData::getFieldPhi()
{
    return m_fieldPhi;
}

simulationData::simulationParameters *simulationData::getParameters()
{
    return m_params;
}

int simulationData::getCellsNumber()
{
    return m_defaultCellsNumber;
}

simulationData::simulationField::simulationField(int iCellsNumber, char* iName)
{
    init(iCellsNumber, iName);
}

void simulationData::simulationField::init(int iCellsNumber, char* iName)
{
  cellsNumber =  iCellsNumber;
  arr = new double[iCellsNumber];
  arrPrev = new double[iCellsNumber];
  int len=strlen(iName);
  name = new char[len+1];
  name[len]=0;
  strcpy(name, iName);
}

simulationData::simulationField2d::simulationField2d(int iCellsX, int iCellsY, char* iName)
{

    init(iCellsX, iCellsY, iName);
}

void simulationData::simulationField2d::init(int iCellsX, int iCellsY, char* iName)
{
  cellsX =  iCellsX;
  cellsY =  iCellsY;
  arr = new double*[cellsX];
  arrPrev = new double*[cellsX];
  for (int i = 0; i < cellsX; ++i) {
      arr[i] = new double [cellsY];
      arrPrev[i] = new double [cellsY];
  }
  int len=strlen(iName);
  name = new char[len+1];
  name[len]=0;
  strcpy(name, iName);
}

void simulationData::simulationParameters::init(int iCellsX, int iCellsY)
{
    cellsX =  iCellsX;
    cellsY =  iCellsY;
    arrDe = new double*[cellsX];
    arrDeps  = new double*[cellsX];
    arrDomega  = new double*[cellsX];
    arrMue  = new double*[cellsX];
    arrMueps  = new double*[cellsX];
    arrMuomega  = new double*[cellsX];
    arrE  = new double*[cellsX];
    arrTe  = new double*[cellsX];
    for (int i = 0; i < cellsX; ++i) {
        arrDe[i] = new double [cellsY];
        arrDeps[i] = new double [cellsY];
        arrDomega[i] = new double [cellsY];
        arrMue[i] = new double [cellsY];
        arrMueps[i] = new double [cellsY];
        arrMuomega[i] = new double [cellsY];
        arrE[i] = new double [cellsY];
        arrTe[i] = new double [cellsY];
    }
    T=400; //K
    p=101325;//pa
    mAr=39.948/1000.0;//kg/mol
    rho=p*mAr/(8.314*T);//kg/m^3
    double Na=6.022e23; //1/mol
    N=p*Na/(8.314*T);

}

void simulationData::updateParams()
{
    simulationData::simulationParameters* pParams = m_params;
    simulationData::simulationField2d* pEn = m_fieldEnergy;
    simulationData::simulationField2d* pNe = m_fieldNe;
    simulationData::simulationField2d* pPhi = m_fieldPhi;
    //double dz = m_dz;


   /* for (int j=0; j<pParams->cellsNumber; j++)
    {
        pParams->arrMue[j] =4e24/m_params->N; ; //4e4; m^2/(V*s)
        pParams->arrMueps[j] =5.0 * pParams->arrMue[j] / 3.0;

        pParams->arrTe[j]=(2.0/3.0)*fabs(pEn->arr[j])/(fabs(pNe->arr[j])+1);
        if (pParams->arrTe[j]>60.0) pParams->arrTe[j]=60.0; //some limiter here


        pParams->arrDe[j] = pParams->arrMue[j] * pParams->arrTe[j];
        pParams->arrDeps[j] = pParams->arrMueps[j] * pParams->arrTe[j];
        //pParams->arrE[j] = - simulationTools::ddzCentral(pPhi->arr, pPhi->cellsNumber, dz, j);
        pParams->arrDomega[j] = 0.01;// m^2/s
        pParams->arrMuomega[j]=pParams->arrDomega[j]*q/(k_B_const*pParams->T);

    }*/
}



simulationData::simulationParameters::simulationParameters(int iCellsX, int iCellsY)
{
    init(iCellsX,iCellsY);
}
