#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "functions.h"

using namespace std;

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

/*
void leastSquares_nDimension(int idFunction, double* parameters, int nrParameters,
                  double** x, double* y, int nrData,int xDim, double* lambda,
                  double* parametersDelta, double* parametersChange)
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
        for (int jj=0; jj<xDim; jj++)
        {
            xPoint[jj] = x[i][jj];
        }
        firstEst[i] = estimateFunction_nDimension(idFunction, parameters, nrParameters, xPoint,xDim);
    }

    // change parameters and compute derivatives
    for (i = 0; i < nrParameters; i++)
    {
        parameters[i] += parametersDelta[i];
        for (j = 0; j < nrData; j++)
        {
            for (int jj=0; jj<xDim; jj++)
            {
                xPoint[jj] = x[j][jj];
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
                a[i][j] += P[i][k] * P[j][k];
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
*/
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
    int nrParameters =3;
    int nrData =5;
    int xDim = 2;
    int maxIterationsNr = 10000;
    double myEpsilon = EPSILON;
    int idFunction = FUNCTION_CODE_MULTILINEAR;
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
        parametersMin[i]= -10;
        parametersMax[i]= 10;
        parametersDelta[i] = 0.00001;
        parameters[i]= 0;
    }
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

    fittingMarquardt_nDimension(parametersMin,parametersMax,parameters,nrParameters,parametersDelta,maxIterationsNr,myEpsilon,idFunction,x,y,nrData,xDim,isWeighted,weights);
    for (int i=0;i<nrParameters;i++)
    {
        printf("parameter %d  =  %f;\n",i,parameters[i]);
    }
    printf("\n");
    for (int i=0;i<5;i++)
    {
        //printf(" %d%f\n",i,multilinear(x[i],xDim,parameters));
    }
    free(x);
    free(y);
    free(weights);
    free(parameters);
    free(parametersDelta);
    free(parametersMax);
    free(parametersMin);
    printf("The end\n");
    return 0;
}
