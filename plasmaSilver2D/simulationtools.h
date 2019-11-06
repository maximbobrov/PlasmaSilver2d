#ifndef SIMULATIONTOOLS_H
#define SIMULATIONTOOLS_H

class simulationTools
{
public:
    static double ddz(double *arr, double dz, int i);
    static double ddzCentral(double *arr, int cellsNumber, double dz, int i);
    double d2dz(double *arr, int cellsNumber, double dz, int j);
    static double gauss(double x, double d);
    simulationTools();
};
#endif // SIMULATIONTOOLS_H
