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

extern double const xi_u[Q];

extern double const xi_v[Q];

extern double const omega[Q];

#endif