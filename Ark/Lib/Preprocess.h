#ifndef _ARK_PREPROCESS_H
#define _ARK_PREPROCESS_H

char caseName[128] = {'\0'};

void SetAdaptList()
{	
	adaptList = list_concat(hTildeList,fTildeList);
	adaptList = list_concat(adaptList,msvList);

	// boundaryList = list_concat(hTildeList,hBPList);
	// boundaryList = list_concat(boundaryList,fTildeList);
	// boundaryList = list_concat(boundaryList,fBPList);
	// boundaryList = list_concat(boundaryList,msvList);

	boundaryList = list_concat(hBPList,fBPList);
	boundaryList = list_concat(boundaryList,msvList);

	scalar *haloScalarList = list_concat((scalar *) hFluxList,(scalar *) fFluxList);
	haloFluxList = vectors_from_scalars(haloScalarList);
	free(haloScalarList);
}
void InterpolationScheme()
{
	for(scalar s in all)
	{
		// s.prolongation = refine_biquadratic;
		// s.refine = s.prolongation;
		//s.restriction = coarsen_quadratic;
	}
}
void update_xEq(Point point, double hEq[], double hS[],double fEq[],double fS[]);

double analyticalPhi(double x,double y)
{
	double phi = 
	0.5*(PhiL+PhiV) + 0.5*(PhiL-PhiV)*		//bubble
	// 0.5*(PhiL+PhiV) - 0.5*(PhiL-PhiV)*		//drop
	tanh
	(
		(sqrt((x-centerX)*(x-centerX) + (y-centerY)*(y-centerY))-radius)*2.0/wI
	);
	// double phi = 
	// 0.5*(PhiL+PhiV) + 0.5*(PhiL-PhiV)*		//bubble
	// // 0.5*(PhiL+PhiV) - 0.5*(PhiL-PhiV)*		//drop
	// tanh
	// (
	// 	(sqrt((y-centerY)*(y-centerY))-radius)*2.0/wI
	// );
	return phi;

}
double analyticalPhiPlate(double x,double y)
{
	double phi = 
	0.5*(PhiL+PhiV) + 0.5*(PhiL-PhiV)*		//bubble
	// 0.5*(PhiL+PhiV) - 0.5*(PhiL-PhiV)*		//drop
	tanh
	(
		(sqrt((y-centerY)*(y-centerY))-radius)*2.0/wI
	);
	return phi;

}

void AC_Bubble()
{
	strcpy(caseName,__func__);
	printSplitLine(sp);
	printf("Case Name : %s\n",caseName);
	printSplitLine(nl);
	//
	foreach()
	{
		int k = 0;
		#ifdef _ARK_ALLENCAHN_FLIP
		phi[] = analyticalPhi(x,y);
		rho[] = aPhi*phi[] + bPhi;
		//phi[] = 1;
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		u[]   = U0;
		v[]   = U0;
		p[]   = p0;
		ftau[] = Tau0;
		#endif

		double hEq[Q] = {0.0}, hS[Q] = {0.0};
		double fEq[Q] = {0.0}, fS[Q] = {0.0};
		update_xEq(point,hEq,hS,fEq,fS);

		k = 0;
		scalar hTilde,fTilde;
		for(hTilde,fTilde in hTildeList,fTildeList)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			hTilde[] = hEq[k];
			#endif

			#ifdef _ARK_MOMENTUM_FLIP
			fTilde[] = fEq[k];
			#endif
			++k;
		}
	}
	#ifdef _ARK_ADAPT_FLIP
	adapt_wavelet({phi},(double[]){AdaptError},MaxLevel,MinLevel,adaptList);
	#endif
}
void AC_Plate()
{
	strcpy(caseName,__func__);
	printSplitLine(sp);
	printf("Case Name : %s\n",caseName);
	printSplitLine(nl);

	foreach()
	{
		int k = 0;
		#ifdef _ARK_ALLENCAHN_FLIP
		phi[] = analyticalPhiPlate(x,y);
		rho[] = aPhi*phi[] + bPhi;
		//phi[] = 1;
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		u[]   = U0;
		v[]   = U0;
		p[]   = p0;
		ftau[] = Tau0;
		#endif

		double hEq[Q] = {0.0}, hS[Q] = {0.0};
		double fEq[Q] = {0.0}, fS[Q] = {0.0};
		update_xEq(point,hEq,hS,fEq,fS);

		k = 0;
		scalar hTilde,fTilde;
		for(hTilde,fTilde in hTildeList,fTildeList)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			hTilde[] = hEq[k];
			#endif

			#ifdef _ARK_MOMENTUM_FLIP
			fTilde[] = fEq[k];
			#endif
			++k;
		}
	}
	#ifdef _ARK_ADAPT_FLIP
	adapt_wavelet({phi},(double[]){AdaptError},MaxLevel,MinLevel,adaptList);
	#endif
}
void ErrorMsgForConstants(const char* nameA, const char* nameB, double A, double B)
{
	printErrorLine(' ');
	printf("%s !=  %s: %g\t%g\n",nameA,nameB,A,B);
	printErrorLine(nl);
	exit(-1);
}
bool isTiny(double a)
{
	return fabs(a) < TINY ? true : false;
}
bool doubleEqual(double a,double b)
{
	return isTiny(a-b);
}
void CheckConstants()
{
	if(!doubleEqual(MinL,ChLength/NL))
	{
		printErrorLine(' ');
		printf("MinL !=  ChLength/NL: %g\t%g\n",MinL,ChLength/NL);
		printErrorLine(nl);
		exit(-1);;
	}
	if(X_End-X_Beg != Lx || Y_End-Y_Beg != Ly)
	{
		printErrorLine(' ');
		printf("Lx or Ly doesn't match : %g\t%g\n",Lx,Ly);
		printErrorLine(nl);
		exit(-1);
	}
	if(Lambda0 != 1.0/(2.0*R0*T0))
	{
		printErrorLine(' ');
		printf("Lambda0 != 1.0/(2*R0*T0) : %g\t%g\n",Lambda0,1.0/(2.0*R0*T0));
		printErrorLine(nl);
		exit(-1);
	}
	if(RT != R0*T0)
	{
		printErrorLine(' ');
		printf("RT != R0*T0 : %g\t%g\n",RT,R0*T0);
		printErrorLine(nl);
		exit(-1);
	}
	if(Nu0 != Mu0/Rho0)
	{
		printErrorLine(' ');
		printf("Nu0 != Mu0/Rho0 : %g\t%g\n",Nu0,Mu0/Rho0);
		printErrorLine(nl);
		exit(-1);
	}
	if(Tau0 != Nu0/RT)
	{
		printErrorLine(' ');
		printf("Tau0 != Nu0/RT : %g\t%g\n",Tau0,Nu0/RT);
		printErrorLine(nl);
		exit(-1);
	}
	if(MaxU != sqrt(3.0*RT))
	{
		printErrorLine(' ');
		printf("MaxU != sqrt(3.0*RT) : %g\t%g\n",MaxU,sqrt(3.0*RT));
		printErrorLine(nl);
		exit(-1);
	}
	if(Dt != CFL*MinL/MaxU)
	{
		printErrorLine(' ');
		printf("Dt != CFL*MinL/MaxU : %g\t%g\n",Dt,CFL*MinL/MaxU);
		printErrorLine(nl);
		exit(-1);
	}
	if(hDt != 0.5*Dt)
	{
		printErrorLine(' ');
		printf("hDt != 0.5*Dt : %g\t%g\n",hDt,0.5*Dt);
		printErrorLine(nl);
		exit(-1);
	}
	#ifdef _ARK_ALLENCAHN_FLIP
	if(!doubleEqual(tauM , M_Phi/RT))
	{
		printErrorLine(' ');
		printf("tauM != M_Phi/RT : %g\t%g\n",tauM,M_Phi/RT);
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(aPhi,(RhoL-RhoV)/(PhiL-PhiV)))
	{
		printErrorLine(' ');
		printf("aPhi != (RhoL-RhoV)/(PhiL-PhiV) : %g\t%g\n",aPhi,(RhoL-RhoV)/(PhiL-PhiV));
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(bPhi,(RhoV-aPhi*PhiV)))
	{
		printErrorLine(' ');
		printf("bPhi != (RhoV-aPhi*PhiV) : %g\t%g\n",bPhi,(RhoV-aPhi*PhiV));
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(MuL , RhoL*NuL))
	{
		printErrorLine(' ');
		printf("MuL != RhoL*NuL : %g\t%g\n",MuL,RhoL*NuL);
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(MuV , RhoV*NuV))
	{
		printErrorLine(' ');
		printf("MuV != RhoV*NuV : %g\t%g\n",MuV,RhoV*NuV);
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(Cn , wI/ChLength))
	{
		printErrorLine(' ');
		printf("Cn != wI/ChLength : %g\t%g\n",Cn,wI/ChLength);
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(Pe , U0*ChLength/M_Phi))
	{
		printErrorLine(' ');
		printf("Pe != U0*ChLength/M_Phi : %g\t%g\n",Pe,U0*ChLength/M_Phi);
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(RhoRatio , RhoL/RhoV))
	{
		printErrorLine(' ');
		printf("RhoRatio != RhoL/RhoV : %g\t%g\n",RhoRatio,RhoL/RhoV);
		printErrorLine(nl);
		exit(-1);
	}
	if(!doubleEqual(Beta , 12.0*Sigma/wI))
	{
		ErrorMsgForConstants("Beta","12.0*Sigma/wI",Beta,12.0*Sigma/wI);
	}
	if(!doubleEqual(Kappa , 3.0*Sigma*wI/2.0))
	{
		ErrorMsgForConstants("Kappa","3.0*Sigma*wI/2.0",Kappa,3.0*Sigma*wI/2);
	}
	#endif
}
#endif