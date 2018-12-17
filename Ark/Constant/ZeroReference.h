#ifndef _ZERO_REFERENCE_H
#define _ZERO_REFERENCE_H

#include <math.h>

#define CHLENGTH 256

//----------------------Mesh file---------------
const int 

NL = CHLENGTH;

const double

ChLength = CHLENGTH,

MinL = 1.0,

X_Beg = 0.0,

X_End = CHLENGTH,

Y_Beg = 0.0,

Y_End = CHLENGTH,

Lx = CHLENGTH,

Ly = CHLENGTH;


//-----------------------------reference parameters-------------------------

const double 

T0 = 1.0/3,

R0 = 1.0,

Lambda0 = 1.5,

Rho0 = 1.0,//4.5435E-2,//CS

U0 = 0.0,//W_i*TC_r,

V0 = 0.0,

Ma = 0.0,

p0 = 1,

Re = 0.0,

// Mu0 = U0*ChLength*Rho0/Re,

Mu0 = 1.0/6.0,//3.413333E-3,//0.06827,//,////,//

Nu0 = 1.0/6.0,

Tau0 = 0.5;

// Kn = 16.0*Mu0/(5.0*Rho0*R0*T0)*sqrt(1.0/(4.0*Lambda0*PI));

//------------------------------Atomic species-------------------------------------
// const double 

// Omega0 = 0.5, //hard sphere(HS) = 0.5, variable hard sphere(VHS) = 0.68

// Pr = 2.0/3.0, //Prandtl Number

// nK = 0,//internal degree 0 = single  2 = double

// agEq = (nK + 3.0 - 2.0)/2.0,//used in the Equilibrium of g(x,xi,t)

// Cv = (nK + 3.0)/2.0,

// Gamma = (nK + 5.0)/(nK + 3.0);
double const

TINY = 1.0E-14,

PI = 3.141592653589793;
#endif

