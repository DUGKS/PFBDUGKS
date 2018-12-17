#ifndef _ARK_BASILISKDUGKS_H
#define _ARK_BASILISKDUGKS_H

/**
ZeroConstant should be placed at the beginning place since it includes ZeroFlip.h,
a macro file that contains an array of switchs. 
*/
#include "./Constant/ZeroConstant.h"

//varied according to the velocity model used.
#include "./Header/D2Q9.h"
//be in consistent with velocity model
#include "./Header/DDF.h"

/***************************************************************\
*					define macro variables
\***************************************************************/
scalar rho[], phi[], u[], v[], p[], ftau[];

scalar GradX_rho[],GradY_rho[],GradX_phi[],GradY_phi[],ForceX[], ForceY[];

scalar prev_phiu[],prev_phiv[];

scalar phi_1k[],rho_1k[],u_1k[],v_1k[],Error_1k[];
//macroscopic varaibels List
scalar *msvList = 

{rho,phi,u,v,p,GradX_rho,GradY_rho,GradX_phi,GradY_phi,ForceX,ForceY,prev_phiu,prev_phiv};

scalar *adaptList;

scalar *boundaryList;
		
int step = 0; int ThreadNum = 1;

double totalDensity = 0.0;

double residualPer1k = 1.0;

#endif
