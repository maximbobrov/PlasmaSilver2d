#include "simulationtools.h"
#include "math.h"

simulationTools::simulationTools()
{

}

double simulationTools::ddz(double* arr, double dz, int j)
{
    return j == 0 ? (arr[1] - arr[0]) / (dz) : (arr[j]-arr[j-1])/(dz);
}


double simulationTools::ddzCentral(double *arr, int cellsNumber, double dz, int j)
{
    return j == 0 ? (-3.0*arr[j] + 4.0*arr[j+1] - arr[j+2])/(2.0*dz) : (j == cellsNumber-1 ? (3.0*arr[j] - 4.0*arr[j-1] + arr[j-2])/(2.0*dz) : (arr[j+1]-arr[j-1])/(2.0*dz));
}

double simulationTools::d2dz(double *arr, int cellsNumber, double dz, int j)
{
    return j == 0 ? (arr[j] - 2.0 * arr[j+1] + arr[j+2])/(dz*dz) : (j == cellsNumber-1 ? (arr[j] - 2.0*arr[j-1] + arr[j-2])/(dz*dz) : (arr[j] - 2.0*arr[j-1] + arr[j-2])/(dz*dz));
}

double simulationTools::gauss(double x, double d)
{
    return exp(-x*x/(d*d));
}
