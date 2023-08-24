#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define NODATA -9999
#define EPSILON 0.00001
#define MINVALUE(a, b) (((a) < (b))? (a) : (b))
#define MAXVALUE(a, b) (((a) > (b))? (a) : (b))

enum estimatedFunction {FUNCTION_CODE_SPHERICAL, FUNCTION_CODE_LINEAR, FUNCTION_CODE_PARABOLIC,
                   FUNCTION_CODE_EXPONENTIAL, FUNCTION_CODE_LOGARITMIC,
                   FUNCTION_CODE_TWOPARAMETERSPOLYNOMIAL, FUNCTION_CODE_FOURIER_2_HARMONICS,
                   FUNCTION_CODE_FOURIER_GENERAL_HARMONICS,
                   FUNCTION_CODE_MODIFIED_VAN_GENUCHTEN, FUNCTION_CODE_MODIFIED_VAN_GENUCHTEN_RESTRICTED,
                   FUNCTION_CODE_MULTILINEAR};


double parabolicFunction(double* x, double* par);
double multilinear(double* x, int xDim, double* par);
bool fittingMarquardt_nDimension(double* parametersMin, double* parametersMax, double* parameters, int nrParameters,
                      double* parametersDelta, int maxIterationsNr, double myEpsilon, int idFunction,
                      double** x, double* y, int nrData, int xDim, bool isWeighted, double *weights);

double normGeneric_nDimension(int idFunction, double *parameters,int nrParameters, double** x, double* y, int nrData, int xDim);
double estimateFunction_nDimension(int idFunction, double *parameters, int nrParameters, double* x, int xDim);
void leastSquares_nDimension(int idFunction, double* parameters, int nrParameters,
                  double** x, double* y, int nrData,int xDim, double* lambda,
                  double* parametersDelta, double* parametersChange,bool isWeighted, double* weights);



#endif // FUNCTIONS_H
