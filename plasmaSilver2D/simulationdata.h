#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H
#include <iostream>
#include <vector>
#include "reaction.h"

class simulationData
{
public:
    simulationData(int iCellsX,int iCellsY);
    ~simulationData();

    enum SpecieName
    {
       Ar,
       Ar_star,
       Ar_plus
    };

    enum ReactionName
    {
       eAr_eAr,
       eAr_eArs,
       eAr_2eArp,
       eArs_2eArp,
        comsol_eAr_2eArp
    };

    struct simulationField
    {
        char* name;
        double* arr;
        double* arrPrev;
        int cellsNumber;
        simulationField(int cellsNumber, char* name);
        private: void init(int cellsNumber, char* name);
    };

    struct simulationField2d
    {
        char* name;
        double** arr;
        double** arrPrev;
        int cellsX;
        int cellsY;
        simulationField2d(int cellsX, int cellsY, char* name);
        private: void init(int cellsX, int cellsY, char* name);
    };

    struct simulationParameters
    {
        double** arrDe;
        double** arrDeps;
        double** arrDomega;
        double** arrMue;
        double** arrMueps;
        double** arrMuomega;
        double** arrTe;
        double** arrE;
        int cellsX;
        int cellsY;
        double rho; //mixture density
        double p; //pressure
        double T; //temperature
        double mAr; //argon molar mass
        double N; //neutral number denisty
        simulationParameters(int cellsX, int cellsY);
        private: void init(int cellsX, int cellsY);
    };

    const double q=1.6022e-19; //coulumbs elementary charge
    const double k_B_const=1.3806e-23;// boltzmann constant

    void setDt(double dt = 0.000003);
    void setDx(double dx = 0.01);
    void setDy(double dy = 0.01);

    void setCellsNumber(int cellsNumber);

    int getCellsNumber();
    double getDt();
    double getDx();
    double getDy();
    void updateParams();
    double** getArrTe();
    simulationField2d* getFieldNe();
    simulationField2d* getFieldEnergy();
    simulationField2d* getFieldPhi();
    simulationField2d* getFieldHeavySpicies(int num);
    int getHeavySpiciesCharge(int num);
    int getNumberHeavySpicies();
    double getN(); //total number of particles in one m^3
    double* getReactionRate(simulationData::ReactionName reactName);
    double getReactionDe(simulationData::ReactionName reactName);

    void calcReaction(ReactionName reactName);
    simulationParameters *getParameters();

private:
    double m_dt;
    double m_dx;
    double m_dy;
    int m_defaultCellsNumber;
    int m_numberHeavySpicies;
    simulationField2d* m_fieldNe;
    simulationField2d* m_fieldEnergy;
    simulationField2d* m_fieldPhi;
    std::vector<simulationField2d* > m_fieldsHeavySpecies;
    std::vector<int> m_chargeHeavySpecies;
    simulationParameters* m_params;
    std::vector<reaction *> m_reactions;
};

#endif // SIMULATIONDATA_H
