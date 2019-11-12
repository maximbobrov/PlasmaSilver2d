#ifndef SIMULATIONSOLVER_H
#define SIMULATIONSOLVER_H
#include "simulationdata.h"


class simulationSolver
{
public:

    simulationSolver(simulationData* pData = nullptr);


    virtual double solve(int numberIteration);
    virtual double getRhs();
    virtual void getStepEuler();
    virtual void setBc();
    virtual double init(double value);

    virtual ~simulationSolver();

protected:
    simulationData* m_pData;
    simulationData::simulationField2d* m_field;
    double** m_aRHS;
    double** m_aNu;
    double** m_mask;
    double** m_maskValue;
};

class solverNe : public simulationSolver
{
public:
    solverNe(simulationData* pData = nullptr);
    virtual double getRhs();

    virtual ~solverNe();
};


class solverEnergy : public simulationSolver
{
public:
    solverEnergy(simulationData* pData = nullptr);
    virtual double getRhs();

    virtual ~solverEnergy();
};

class solverPhi : public simulationSolver
{
public:
    solverPhi(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double solve(int numberIteration);
    virtual void setBc();

    virtual ~solverPhi();
private:
    double q_e = 1.6*10e-19;
    double eps_0 = 8.85*10e-12;
};

class solverHeavySpicies : public simulationSolver
{
public:
    solverHeavySpicies(simulationData* pData = nullptr, int num = 0);
    virtual double getRhs();


    virtual ~solverHeavySpicies();
private:
    int m_charge;
};

#endif // SIMULATIONSOLVER_H
