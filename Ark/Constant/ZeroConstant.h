#ifndef _ZERO_CONSTANT_H
#define _ZERO_CONSTANT_H

#include "ZeroFlip.h"
#include "ZeroReference.h"

const double 

	CFL = 0.5,

	Dt = 0.5,// CFL >0 ? CFL*MinL/MaxU : 1E-4

	hDt = 0.25;// 0.5*Dt

#ifdef _ARK_ALLENCAHN_FLIP

const double

	M_Phi = 0.1, tauM = 0.3,

	PhiL = 1.0, PhiV = -1, 

	RhoL = 10.0, RhoV = 1.0, RhoRatio = 10.0,

	aPhi = 4.5,// (RhoL-RhoV)/(PhiL-PhiV);

	bPhi = 5.5,// RhoV-aPhi*PhiV;

	MuL = 10.0/6.0, MuV = 1.0/6.0,

	NuL = 1.0/6.0, NuV = 1.0/6.0,

	centerX = 0.5*CHLENGTH, centerY = 0.5*CHLENGTH,	radius = 0.25*CHLENGTH,

	Gx = 0.0, Gy = 0.0,

	Sigma = 1E-3,

	wI = 8.0,

	Beta = 12.0*1E-3/8.0, // 12.0*Sigma/wI;

	Kappa = 3.0*1E-3*8.0/2.0, // 3.0*Sigma*wI/2;

	Cn = 8.0/CHLENGTH,

	Pe = 0.0;

#endif


const int

MaxLevel = 8,

MinLevel = 6,

End_Step = 0,//PhaseFieldAC::iT + 1000,//log(2.0)/(8.0*PI*PI*Nu0*dt),

ConvergenceInterval = 1000, //SumRho,SumT,independent

ResidualInterval = 1000, //print to screen

writeFileInterval = 1000; //always >= ResidualControl

const double 

RESIDUAL = 1E-10,

AdaptError = 1E-5; 

#endif