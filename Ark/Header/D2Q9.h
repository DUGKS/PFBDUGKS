#ifndef _ARK_D2Q9_H
#define _ARK_D2Q9_H
//
#include <math.h>

#define Q  9

double const

RT = 1.0/3.0,

MaxU = 1.0,

MaSpan = 1.0;//sqrt(3.0*RT)

typedef struct{
  double x;
  double y;
} DmQn;

typedef struct{
	#ifdef _ARK_ALLENCAHN_FLIP
	double *hEq;
	double *hS;
	#endif

	#ifdef _ARK_MOMENTUM_FLIP
	double *fEq;
	double *fS;
	#endif

	#ifdef _ARK_ENERGY_FLIP
	double *gEq;
	double *gS;
	#endif
} DDF_EqS;

extern double const xi_u[Q];

extern double const xi_v[Q];

extern double const omega[Q];

#endif