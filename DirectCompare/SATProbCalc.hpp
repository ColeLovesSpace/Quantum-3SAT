#ifndef SolveSAT_HPP
#define SolveSAT_HPP
 
#include <math.h>
 
// extern "C"
// double function(double kx, double ky, double kz, double * eps, double sigma, double G_max);
 
extern "C"
double SolveSAT(int n, double *SAT, int c, int i, int v);

extern "C"
int SolveSATbpp(int n, double * SAT, int c, int i);

extern "C"
double * SolveSATactual(int n, double *SAT, int c, int i, int v);

extern "C"
double * SolveSATgates(int n, double *SAT, int c, int i, int v);

extern "C"
int SolveSATiterations(int n, double *SAT, int c, int i, int v);

// extern "C"
// int SolveSATfast(int n, double *SAT, int c, int i, int v);


 
#endif