#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "functions.h"

using namespace std;

double computeR2(double *obs, double* sim, int nrPoints)
{
    double R2=0;
    double meanObs=0;
    double RSS=0;
    double TSS=0;
    for (int i=0;i<nrPoints;i++)
    {
        meanObs += obs[i];
    }
    meanObs /= nrPoints;
    //compute RSS and TSS
    for (int i=0;i<nrPoints;i++)
    {
        RSS += (obs[i]-sim[i])*(obs[i]-sim[i]);
        TSS += (obs[i]-meanObs)*(obs[i]-meanObs);
    }
    R2 = 1 - RSS/TSS;
    return R2;
}

double functionTemperatureVsHeight(double* x, double* par)
{
    double y;
    y = par[0] + par[1]*x[0] + par[2]*(1/(1+exp(-par[3]*(x[0] - par[4]))));
    return y;
}

double multilinear(double* x, int xDim, double* par)
{
    double y = 0;
    for (int i=0;i<xDim;i++)
        y += par[i] * x[i];
    y += par[xDim];
    return y;
}

double parabolicFunction(double* x, double* par)
{
    return par[1]*(x[0] - par[0])*(x[0] - par[0]) + par[2] ;
}

int bestFittingMarquardt_nDimension(double (*func)(double*, double*), int nrTrials, int nrMinima,
                                                    double* parametersMin, double* parametersMax, double* parameters, int nrParameters,
                                                    double* parametersDelta, int maxIterationsNr, double myEpsilon,
                                                    double** x, double* y, int nrData, int xDim, bool isWeighted, double* weights)
{
    double bestR2 = -9999;
    double R2;
    double* R2Previous = (double *) calloc(nrMinima, sizeof(double));
    double* ySim = (double *) calloc(nrData, sizeof(double));
    double* bestParameters = (double *) calloc(nrParameters, sizeof(double));
    int i;
    int iRandom = 0;
    int counter = 0;
    for (i=0; i<nrMinima; i++)
    {
        R2Previous[i] = NODATA;
    }
    srand (time(nullptr));
    do
    {
        for (i=0;i<nrParameters;i++)
        {
            parameters[i] = parametersMin[i] + ((double) rand() / (RAND_MAX))*(parametersMax[i]-parametersMin[i]);
        }
        fittingMarquardt_nDimension(func,parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,x,y,nrData,xDim,isWeighted,weights);
        for (i=0;i<nrData;i++)
        {
            double xSim;
            xSim = x[i][0];
            ySim[i]= func(&xSim,parameters); // TODO to be generalized!!
        }
        R2 = computeR2(y,ySim,nrData);
        //printf("%d R2 = %f\n",iRandom,R2);
        if (R2 > bestR2-EPSILON)
        {
            for (int j=0;j<nrMinima-1;j++)
            {
                R2Previous[j] = R2Previous[j+1];
            }
            R2Previous[nrMinima-1] = R2;
            bestR2 = R2;
            for (i=0;i<nrParameters;i++)
            {
                bestParameters[i] = parameters[i];
            }
        }
        iRandom++;
        counter++;
    } while(iRandom<nrTrials && R2<(1 - EPSILON) && fabs(R2Previous[0]-R2Previous[nrMinima-1])>0.0001);

    for (i=0;i<nrParameters;i++)
    {
        parameters[i] = bestParameters[i];
    }
    free(bestParameters);
    free(ySim);
    free(R2Previous);
    return counter;
}

bool fittingMarquardt_nDimension(double (*func)(double*, double*),double* parametersMin, double* parametersMax, double* parameters, int nrParameters,
                      double* parametersDelta, int maxIterationsNr, double myEpsilon,
                      double** x, double* y, int nrData, int xDim,bool isWeighted, double* weights)
{
    // Sum of Squared Erros
    double mySSE, diffSSE, newSSE;
    static double VFACTOR = 10;

    double* paramChange = (double *) calloc(nrParameters, sizeof(double));
    double* newParameters = (double *) calloc(nrParameters, sizeof(double));
    double* lambda = (double *) calloc(nrParameters, sizeof(double));

    for(int i = 0; i < nrParameters; i++)
    {
        lambda[i] = 0.01;       // damping parameter
        paramChange[i] = 0;
    }

    mySSE = normGeneric_nDimension(func, parameters, nrParameters, x, y, nrData,xDim);

    int iterationNr = 0;
    do
    {
        leastSquares_nDimension(func, parameters, nrParameters, x, y, nrData,xDim, lambda, parametersDelta, paramChange,isWeighted,weights);
            // change parameters
        for (int i = 0; i < nrParameters; i++)
        {
            newParameters[i] = parameters[i] + paramChange[i];
            if ((newParameters[i] > parametersMax[i]) && (lambda[i] < 1000))
            {
                newParameters[i] = parametersMax[i];
                if (lambda[i] < 1000)
                    lambda[i] *= VFACTOR;
            }
            if (newParameters[i] < parametersMin[i])
            {
                newParameters[i] = parametersMin[i];
                if (lambda[i] < 1000)
                    lambda[i] *= VFACTOR;
            }
        }

        newSSE = normGeneric_nDimension(func, newParameters, nrParameters, x, y, nrData,xDim);

        if (newSSE == NODATA)
        {
            // free memory
            free(lambda);
            free(paramChange);
            free(newParameters);

            return false;
        }

        diffSSE = mySSE - newSSE ;

        if (diffSSE > 0)
        {
            mySSE = newSSE;
            for (int i = 0; i < nrParameters ; i++)
            {
                parameters[i] = newParameters[i];
                lambda[i] /= VFACTOR;
            }
        }
        else
        {
            for(int i = 0; i < nrParameters; i++)
            {
                lambda[i] *= VFACTOR;
            }
        }

        iterationNr++;
    }
    while (iterationNr <= maxIterationsNr && fabs(diffSSE) > myEpsilon);

    // free memory
    free(lambda);
    free(paramChange);
    free(newParameters);

    return (fabs(diffSSE) <= myEpsilon);
}

void leastSquares_nDimension(double (*func)(double*, double*), double* parameters, int nrParameters,
                  double** x, double* y, int nrData,int xDim, double* lambda,
                  double* parametersDelta, double* parametersChange, bool isWeighted, double* weights)
{
    int i, j, k;
    double pivot, mult, top;

    double* g = (double *) calloc(nrParameters+1, sizeof(double));
    double* z = (double *) calloc(nrParameters+1, sizeof(double));
    double* firstEst = (double *) calloc(nrData+1, sizeof(double));

    double** a = (double **) calloc(nrParameters+1, sizeof(double*));
    double** P = (double **) calloc(nrParameters+1, sizeof(double*));
    double* xPoint = (double *) calloc(xDim, sizeof(double));
    for (i = 0; i < nrParameters+1; i++)
    {
            a[i] = (double *) calloc(nrParameters+1, sizeof(double));
            P[i] = (double *) calloc(nrData+1, sizeof(double));
    }

    // first set of estimates
    for (i = 0; i < nrData; i++)
    {
        for (k=0; k<xDim; k++)
        {
            xPoint[k] = x[i][k];
        }
        firstEst[i] = func(xPoint,parameters);
        //firstEst[i] = estimateFunction_nDimensionExternalFunction(idFunction, parameters, nrParameters, xPoint,xDim);
    }

    // change parameters and compute derivatives
    for (i = 0; i < nrParameters; i++)
    {
        parameters[i] += parametersDelta[i];
        for (j = 0; j < nrData; j++)
        {
            for (k=0; k<xDim; k++)
            {
                xPoint[k] = x[j][k];
            }
            double newEst = func(xPoint,parameters);
            //double newEst = estimateFunction_nDimensionExternalFunction(idFunction, parameters, nrParameters, xPoint,xDim);
            P[i][j] = (newEst - firstEst[j]) / MAXVALUE(parametersDelta[i], EPSILON) ;
        }
        parameters[i] -= parametersDelta[i];
    }

    for (i = 0; i < nrParameters; i++)
    {
        for (j = i; j < nrParameters; j++)
        {
            a[i][j] = 0;
            for (k = 0; k < nrData; k++)
            {
                if (isWeighted)
                {
                    a[i][j] += (weights[k]*(P[i][k] * P[j][k]));
                }
                else
                {
                    a[i][j] += (P[i][k] * P[j][k]);
                }
            }
        }
        z[i] = sqrt(a[i][i]) + EPSILON; //?
    }

    for (i = 0; i < nrParameters; i++)
    {
        g[i] = 0.;
        for (k = 0 ; k<nrData ; k++)
        {
            g[i] += P[i][k] * (y[k] - firstEst[k]);
        }
        g[i] /= z[i];
        for (j = i; j < nrParameters; j++)
        {
            a[i][j] /= (z[i] * z[j]);
        }
    }

    for (i = 0; i < (nrParameters+1); i++)
    {
        a[i][i] += lambda[i];
        for (j = i+1; j < nrParameters; j++)
        {
            a[j][i] = a[i][j];
        }
    }

    for (j = 0; j < (nrParameters - 1); j++)
    {
        pivot = a[j][j];
        for (i = j + 1 ; i < nrParameters; i++)
        {
            mult = a[i][j] / pivot;
            for (k = j + 1; k < nrParameters; k++)
            {
                a[i][k] -= mult * a[j][k];
            }
            g[i] -= mult * g[j];
        }
    }

    parametersChange[nrParameters - 1] = g[nrParameters - 1] / a[nrParameters - 1][nrParameters - 1];

    for (i = nrParameters - 2; i >= 0; i--)
    {
        top = g[i];
        for (k = i + 1; k < nrParameters; k++)
        {
            top -= a[i][k] * parametersChange[k];
        }
        parametersChange[i] = top / a[i][i];
    }

    for (i = 0; i < nrParameters; i++)
    {
        parametersChange[i] /= z[i];
    }

    // free memory
    for (i = 0; i < nrParameters+1; i++)
    {
        free(a[i]);
        free(P[i]);
    }
    free(a);
    free(P);
    free(g);
    free(z);
    free(firstEst);
    free(xPoint);
}

double normGeneric_nDimension(double (*func)(double*, double*), double *parameters,int nrParameters, double** x, double *y, int nrData, int xDim)
{
    double estimate, error;
    double norm = 0;
    double* xPoint = (double *) calloc(xDim, sizeof(double));
    for (int i = 0; i < nrData; i++)
    {
        for (int j=0; j<xDim; j++)
        {
            xPoint[j] = x[i][j];
        }
        //estimate = estimateFunction_nDimensionExternalFunction(idFunction, parameters, nrParameters, xPoint,xDim);
        estimate = func(xPoint,parameters);
        if (estimate == NODATA)
        {
            free(xPoint);
            return NODATA;
        }
        error = y[i] - estimate;
        norm += error * error;
    }
    free(xPoint);
    return norm;
}


/* ////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////// */

int bestFittingMarquardt_nDimension(int nrTrials,int nrMinima, double* parametersMin, double* parametersMax, double* parameters, int nrParameters,
                                      double* parametersDelta, int maxIterationsNr, double myEpsilon, int idFunction,
                                      double** x, double* y, int nrData, int xDim,bool isWeighted, double* weights)
{
    double bestR2 = -9999;
    double R2;
    double* R2Previous = (double *) calloc(nrMinima, sizeof(double));
    double* ySim = (double *) calloc(nrData, sizeof(double));
    double* bestParameters = (double *) calloc(nrParameters, sizeof(double));
    double* xPoint = (double *) calloc(xDim, sizeof(double));
    int i;
    int iRandom = 0;
    int counter = 0;
    for (i=0; i<nrMinima; i++)
    {
        R2Previous[i] = NODATA;
    }
    srand (time(nullptr));
    do
    {
        for (i=0;i<nrParameters;i++)
        {
            parameters[i] = parametersMin[i] + ((double) rand() / (RAND_MAX))*(parametersMax[i]-parametersMin[i]);
        }
        fittingMarquardt_nDimension(parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,idFunction,x,y,nrData,xDim,isWeighted,weights);
        for (i=0;i<nrData;i++)
        {
            for (int k=0; k<xDim; k++)
            {
                xPoint[k] = x[i][k];
            }
            ySim[i]= estimateFunction_nDimension(idFunction,parameters,nrParameters,xPoint,xDim);
        }
        R2 = computeR2(y,ySim,nrData);
        if (R2 > bestR2-EPSILON)
        {
            for (int j=0;j<nrMinima-1;j++)
            {
                R2Previous[j] = R2Previous[j+1];
            }
            R2Previous[nrMinima-1] = R2;
            bestR2 = R2;
            for (i=0;i<nrParameters;i++)
            {
                bestParameters[i] = parameters[i];
            }
        }
        iRandom++;
        counter++;
    } while(iRandom<nrTrials && R2<(1 - EPSILON) && fabs(R2Previous[0]-R2Previous[nrMinima-1])>0.0001);

    for (i=0;i<nrParameters;i++)
    {
        parameters[i] = bestParameters[i];
    }
    free(bestParameters);
    free(ySim);
    free(R2Previous);
    free(xPoint);
    return counter;
}

bool fittingMarquardt_nDimension(double* parametersMin, double* parametersMax, double* parameters, int nrParameters,
                      double* parametersDelta, int maxIterationsNr, double myEpsilon, int idFunction,
                      double** x, double* y, int nrData, int xDim,bool isWeighted, double* weights)
{
    // Sum of Squared Erros
    double mySSE, diffSSE, newSSE;
    static double VFACTOR = 10;

    double* paramChange = (double *) calloc(nrParameters, sizeof(double));
    double* newParameters = (double *) calloc(nrParameters, sizeof(double));
    double* lambda = (double *) calloc(nrParameters, sizeof(double));

    for(int i = 0; i < nrParameters; i++)
    {
        lambda[i] = 0.01;       // damping parameter
        paramChange[i] = 0;
    }

    mySSE = normGeneric_nDimension(idFunction, parameters, nrParameters, x, y, nrData,xDim);

    int iterationNr = 0;
    do
    {
        leastSquares_nDimension(idFunction, parameters, nrParameters, x, y, nrData,xDim, lambda, parametersDelta, paramChange,isWeighted,weights);
            // change parameters
        for (int i = 0; i < nrParameters; i++)
        {
            newParameters[i] = parameters[i] + paramChange[i];
            if ((newParameters[i] > parametersMax[i]) && (lambda[i] < 1000))
            {
                newParameters[i] = parametersMax[i];
                if (lambda[i] < 1000)
                    lambda[i] *= VFACTOR;
            }
            if (newParameters[i] < parametersMin[i])
            {
                newParameters[i] = parametersMin[i];
                if (lambda[i] < 1000)
                    lambda[i] *= VFACTOR;
            }
        }

        newSSE = normGeneric_nDimension(idFunction, newParameters, nrParameters, x, y, nrData,xDim);

        if (newSSE == NODATA)
        {
            // free memory
            free(lambda);
            free(paramChange);
            free(newParameters);

            return false;
        }

        diffSSE = mySSE - newSSE ;

        if (diffSSE > 0)
        {
            mySSE = newSSE;
            for (int i = 0; i < nrParameters ; i++)
            {
                parameters[i] = newParameters[i];
                lambda[i] /= VFACTOR;
            }
        }
        else
        {
            for(int i = 0; i < nrParameters; i++)
            {
                lambda[i] *= VFACTOR;
            }
        }

        iterationNr++;
    }
    while (iterationNr <= maxIterationsNr && fabs(diffSSE) > myEpsilon);

    // free memory
    free(lambda);
    free(paramChange);
    free(newParameters);

    return (fabs(diffSSE) <= myEpsilon);
}

void leastSquares_nDimension(int idFunction, double* parameters, int nrParameters,
                  double** x, double* y, int nrData,int xDim, double* lambda,
                  double* parametersDelta, double* parametersChange, bool isWeighted, double* weights)
{
    int i, j, k;
    double pivot, mult, top;

    double* g = (double *) calloc(nrParameters+1, sizeof(double));
    double* z = (double *) calloc(nrParameters+1, sizeof(double));
    double* firstEst = (double *) calloc(nrData+1, sizeof(double));

    double** a = (double **) calloc(nrParameters+1, sizeof(double*));
    double** P = (double **) calloc(nrParameters+1, sizeof(double*));
    double* xPoint = (double *) calloc(xDim, sizeof(double));
    for (i = 0; i < nrParameters+1; i++)
    {
            a[i] = (double *) calloc(nrParameters+1, sizeof(double));
            P[i] = (double *) calloc(nrData+1, sizeof(double));
    }

    // first set of estimates
    for (i = 0; i < nrData; i++)
    {
        for (k=0; k<xDim; k++)
        {
            xPoint[k] = x[i][k];
        }
        firstEst[i] = estimateFunction_nDimension(idFunction, parameters, nrParameters, xPoint,xDim);
    }

    // change parameters and compute derivatives
    for (i = 0; i < nrParameters; i++)
    {
        parameters[i] += parametersDelta[i];
        for (j = 0; j < nrData; j++)
        {
            for (k=0; k<xDim; k++)
            {
                xPoint[k] = x[j][k];
            }
            double newEst = estimateFunction_nDimension(idFunction, parameters, nrParameters, xPoint,xDim);
            P[i][j] = (newEst - firstEst[j]) / MAXVALUE(parametersDelta[i], EPSILON) ;
        }
        parameters[i] -= parametersDelta[i];
    }

    for (i = 0; i < nrParameters; i++)
    {
        for (j = i; j < nrParameters; j++)
        {
            a[i][j] = 0;
            for (k = 0; k < nrData; k++)
            {
                if (isWeighted)
                {
                    a[i][j] += (weights[k]*(P[i][k] * P[j][k]));
                }
                else
                {
                    a[i][j] += (P[i][k] * P[j][k]);
                }
            }
        }
        z[i] = sqrt(a[i][i]) + EPSILON; //?
    }

    for (i = 0; i < nrParameters; i++)
    {
        g[i] = 0.;
        for (k = 0 ; k<nrData ; k++)
        {
            g[i] += P[i][k] * (y[k] - firstEst[k]);
        }
        g[i] /= z[i];
        for (j = i; j < nrParameters; j++)
        {
            a[i][j] /= (z[i] * z[j]);
        }
    }

    for (i = 0; i < (nrParameters+1); i++)
    {
        a[i][i] += lambda[i];
        for (j = i+1; j < nrParameters; j++)
        {
            a[j][i] = a[i][j];
        }
    }

    for (j = 0; j < (nrParameters - 1); j++)
    {
        pivot = a[j][j];
        for (i = j + 1 ; i < nrParameters; i++)
        {
            mult = a[i][j] / pivot;
            for (k = j + 1; k < nrParameters; k++)
            {
                a[i][k] -= mult * a[j][k];
            }
            g[i] -= mult * g[j];
        }
    }

    parametersChange[nrParameters - 1] = g[nrParameters - 1] / a[nrParameters - 1][nrParameters - 1];

    for (i = nrParameters - 2; i >= 0; i--)
    {
        top = g[i];
        for (k = i + 1; k < nrParameters; k++)
        {
            top -= a[i][k] * parametersChange[k];
        }
        parametersChange[i] = top / a[i][i];
    }

    for (i = 0; i < nrParameters; i++)
    {
        parametersChange[i] /= z[i];
    }

    // free memory
    for (i = 0; i < nrParameters+1; i++)
    {
        free(a[i]);
        free(P[i]);
    }
    free(a);
    free(P);
    free(g);
    free(z);
    free(firstEst);
    free(xPoint);
}


double estimateFunction_nDimension(int idFunction, double *parameters, int nrParameters, double* xPoint, int xDim)
{
    switch (idFunction)
    {
        case FUNCTION_CODE_PARABOLIC :
            return parabolicFunction(xPoint, parameters);

        case FUNCTION_CODE_MULTILINEAR :
            return multilinear(xPoint,xDim,parameters);

        case FUNCTION_CODE_TEMPVSHEIGHT :
            return functionTemperatureVsHeight(xPoint,parameters);

        default:
            return NODATA ;
    }
}


double normGeneric_nDimension(int idFunction, double *parameters,int nrParameters, double** x, double *y, int nrData, int xDim)
{
    double estimate, error;
    double norm = 0;
    double* xPoint = (double *) calloc(xDim, sizeof(double));
    for (int i = 0; i < nrData; i++)
    {
        for (int j=0; j<xDim; j++)
        {
            xPoint[j] = x[i][j];
        }
        estimate = estimateFunction_nDimension(idFunction, parameters, nrParameters, xPoint,xDim);
        if (estimate == NODATA)
        {
            free(xPoint);
            return NODATA;
        }
        error = y[i] - estimate;
        norm += error * error;
    }
    free(xPoint);
    return norm;
}

int main()
{
    int nrParameters =5;
    int nrData =10;
    int xDim = 1;
    int maxIterationsNr = 10000;
    int nrMinima = 10;
    double myEpsilon = EPSILON;
    int idFunction = FUNCTION_CODE_TEMPVSHEIGHT;
    double* parametersMin = (double *) calloc(nrParameters, sizeof(double));
    double* parametersMax = (double *) calloc(nrParameters, sizeof(double));
    double* parameters = (double *) calloc(nrParameters, sizeof(double));
    double* parametersDelta = (double *) calloc(nrParameters, sizeof(double));
    double** x = (double **) calloc(nrData, sizeof(double*));
    double* y = (double *) calloc(nrData, sizeof(double));
    double* weights = (double *) calloc(nrData, sizeof(double));
    bool isWeighted = false;



    for (int i=0;i<nrData;i++)
    {
        x[i] = (double *) calloc(xDim, sizeof(double));
    }

    for (int i=0;i<nrParameters;i++)
    {
        parametersMin[i]= -1000;
        parametersMax[i]= 1000;
        parametersDelta[i] = 0.00001;
        parameters[i]= 0;
    }

    parametersMin[0]= -20;
    parametersMax[0]= 45;
    parametersMin[1]= -0.05;
    parametersMax[1]= 0.001;
    parametersMin[2]= -0.01;
    parametersMax[2]= 100;
    parametersMin[3]= -0.01;
    parametersMax[3]= 1;
    parametersMin[4]= -10;
    parametersMax[4]= 1000;
/*
    parametersDelta[0]= 0.01;
    parametersDelta[1]= 0.0001;
    parametersDelta[2]= 0.01;
    parametersDelta[3]= 0.0001;
    parametersDelta[4]= 0.01;
*/
    parameters[0]= 20;
    parameters[1]= -0.01;
    parameters[2]= 5;
    parameters[3]=0.2;
    parameters[4]= 50;



    /* test multilinear
    x[0][0] = 2;
    x[1][0] = 3;
    x[2][0] = 5;
    x[3][0] = 7;
    x[4][0] = 8;
    x[0][1] = 1;
    x[1][1] = 5;
    x[2][1] = 3;
    x[3][1] = 6;
    x[4][1] = 7;

    y[0] = 3;
    y[1] = 2;
    y[2] = 4;
    y[3] = 5;
    y[4] = 8;
    */

    /* test parabola
    x[0][0] = 0;
    x[1][0] = 1;
    x[2][0] = 2;
    x[3][0] = 3;
    x[4][0] = 4;
    y[0] = -0.;
    y[1] = 1;
    y[2] = 4.;
    y[3] = 9;
    y[4] = 16.;
    */

    x[0][0] = 0;
    x[1][0] = 150;
    x[2][0] = 300;
    x[3][0] = 450;
    x[4][0] = 600;
    x[5][0] = 750;
    x[6][0] = 900;
    x[7][0] = 1050;
    x[8][0] = 1200;
    x[9][0] = 1350;

    y[0] = 18.1;
    y[1] = 23.5;
    y[2] = 22.;
    y[3] = 20.5;
    y[4] = 19;
    y[5] = 14.5;
    y[6] = 16.0;
    y[7] = 14.5;
    y[8] = 13.;
    y[9] = 11.5;
    int nrSteps;
    //nrSteps = bestFittingMarquardt_nDimension(maxIterationsNr,nrMinima,parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,idFunction,x,y,nrData,xDim,isWeighted,weights);
    for (int jj=0; jj<1000; jj++)
    {
        nrSteps = bestFittingMarquardt_nDimension(&functionTemperatureVsHeight, maxIterationsNr,nrMinima,parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,x,y,nrData,xDim,isWeighted,weights);
        printf("%d nrSteps %d\n",jj,nrSteps);
    }
    for (int i=0;i<nrParameters;i++)
    {
        //printf("%f\t",parameters[i]);
    }
    printf("\n");
    /*
    parameters[0]= 20;
    parameters[1]= -0.01;
    parameters[2]= 5;
    parameters[3]=0.2;
    parameters[4]= 50;
    */
    for (int i=0;i<1500;i=i+150)
    {
        double xtry;
        xtry = i*1.0;
        printf(" %f  %f\n",xtry,functionTemperatureVsHeight(&xtry,parameters));
    }

    free(x);
    free(y);
    free(weights);
    free(parameters);
    free(parametersDelta);
    free(parametersMax);
    free(parametersMin);

    //printf("Test parabola\n");
    //double* xParabola = (double *) calloc(1, sizeof(double));
    //double* parParabola = (double *) calloc(3, sizeof(double));;
    //xParabola[0] = 2;
    //parParabola[1] = 1 ;
    //provaFunzione(&parabolicFunction,xParabola,parParabola);


    printf("The end\n");
    return 0;
}

double provaFunzione(double (*func)(double*, double*), double* x, double* param)
{
    double y=0;
    //double x,param;
    y = ((*func)(x,param));
    return y;
}
