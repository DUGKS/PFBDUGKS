#ifndef _ARK_DUGKSSOLVER_H
#define _ARK_DUGKSSOLVER_H

double const aT = 4.0/3.0, bT = 1.0/3.0;

//-------------------------D2Q9AC.h-------------------------

double calcMu(double phi);

void MicroFlux_2D();

void update_MacroVar();

void update_Source();

void update_xEq(Point point,DDF_EqS xEqxS);

//------------------------Output.h---------------------------

void Output_Flowfield(int step);

void Output_Residual(int step);

void Output_TotalDensity(int step);


//------------------------DUGKSSolver.h---------------------------

void update_xBarPlus();

void update_xTilde();

void update_TotalDensity(int step);

void update_Residual(int step);

void DUGKSSolver()
{
	step = 0;ThreadNum = 24;
	Output_Flowfield(step);

	printSplitLine(sp);
	printf("iteration start    ThreadNum : %d\n",ThreadNum);
	printSplitLine(nl);

	while(residualPer1k > RESIDUAL)
	{
		update_xBarPlus();

		#ifdef _ARK_ADAPT_FLIP
		boundary(boundaryList);
		#endif

		MicroFlux_2D();

		#ifdef _ARK_ADAPT_FLIP
		boundary_flux(haloFluxList);
		#endif

		update_xTilde();
		update_MacroVar();
		update_Source();
		
		++step;
		if(step%ConvergenceInterval == 0)
		{
			update_TotalDensity(step);
		}
		if(step%ResidualInterval == 0)
		{
			update_Residual(step);			
		}

		#ifdef _ARK_ADAPT_FLIP
		if(step % 100 == 0 && step < 1000)
		//if(step < 10)
		{
		boundary(adaptList);
		adapt_wavelet({phi},(double[]){AdaptError},MaxLevel,MinLevel,adaptList);
		}
		#endif
	}
}
void update_xBarPlus()
{
	foreach()
	{
		int k = 0;
		//Equilibrium and Source DF
		double hEq[Q] = {0.0}, hS[Q] = {0.0};
		double fEq[Q] = {0.0}, fS[Q] = {0.0};

		DDF_EqS xEqxS;

		#ifdef _ARK_ALLENCAHN_FLIP
		xEqxS.hEq = hEq;
		xEqxS.hS = hS;
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		xEqxS.fEq = fEq;
		xEqxS.fS = fS;
		#endif
		
		update_xEq(point,xEqxS);
		//update \bar{h}^{+}(x,xi,t)
		#ifdef _ARK_ALLENCAHN_FLIP
		double 
		aBPh = (2.0*tauM - hDt)/(2.0*tauM + Dt),
		bBPh = 1.0 - aBPh,
		cBPh = tauM*bBPh;

		k = 0;
		scalar hBP,hTilde;
		for(hBP, hTilde in hBPList, hTildeList)
		{
			hBP[] = aBPh*hTilde[] + bBPh*hEq[k] + cBPh*hS[k];
			++k;
		}
		#endif
		//update \bar{f}^{+}(x,xi,t)
		#ifdef _ARK_MOMENTUM_FLIP
		double fMu = calcMu(phi[]);	ftau[] = fMu/(rho[]*RT);
		double 
		aBPf = (2.0*ftau[] - hDt)/(2.0*ftau[] + Dt), 
		bBPf = 1.0 - aBPf, 
		cBPf = ftau[]*bBPf;

		k = 0;
		scalar fBP,fTilde;
		for(fBP, fTilde in fBPList, fTildeList)
		{
			fBP[] = aBPf*fTilde[] + bBPf*fEq[k] + cBPf*fS[k];
			++k;
		}
		#endif
	}
}
//
void update_xTilde()
{
	foreach()
	{
		int k = 0;

		#ifdef _ARK_ALLENCAHN_FLIP
		double hFluxSum = 0.0;
		k = 0;
		scalar hBP,hTilde; vector hFlux;
		for(hTilde, hBP, hFlux in hTildeList, hBPList, hFluxList)
		{
			hFluxSum = (hFlux.x[] - hFlux.x[1,0])*xi_u[k]
					 + (hFlux.y[] - hFlux.y[0,1])*xi_v[k];
			hTilde[] = aT*hBP[] - bT*hTilde[] + Dt*hFluxSum/Delta;
			++k;
		}
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		double fFluxSum = 0.0;
		k = 0;
		scalar fBP,fTilde; vector fFlux;
		for(fTilde, fBP, fFlux in fTildeList, fBPList, fFluxList)
		{
			fFluxSum = (fFlux.x[] - fFlux.x[1,0])*xi_u[k]
					 + (fFlux.y[] - fFlux.y[0,1])*xi_v[k];

			fTilde[] = aT*fBP[] - bT*fTilde[] + Dt*fFluxSum/Delta;
			++k;
		}
		#endif	
	}
}
void update_Residual(int step)
{
	#if _OPENMP
	omp_set_num_threads(1);
	#endif

	double dRho = 0.0, SumdRho = 0.0;
	double SumRho = 0.0;

	// #if _OPENMP
	// foreach(reduction(+:SumdRho) reduction(+:SumRho))
	// {
	// 	dRho = phi[] - phi_1k[];
	// 	SumdRho += dRho*dRho;
	// 	SumRho += phi[]*phi[];
	// }
	// foreach()
	// {
	// 	phi_1k[] = phi[];
	// }
	// #else
	foreach()
	{
		Error_1k[] = rho[] - rho_1k[];
		SumdRho += Error_1k[]*Error_1k[];
		SumRho += rho[]*rho[];
		rho_1k[] = rho[];
	}
	//#endif
	residualPer1k = sqrt(SumdRho/(SumRho + 1.0E-30));

	#ifdef _ARK_MOMENTUM_FLIP
	#endif

	Output_Residual(step);
	#ifndef _ARK_SERVER_FLIP
	printf("%d\t\t%.8e\t\t%.8e\t\t%.8e\n",step, residualPer1k, totalDensity, SumRho);
	#endif

	if(step%writeFileInterval == 0)
	{
		Output_Flowfield(step);
	}

	#if _OPENMP
	omp_set_num_threads(ThreadNum);
	#endif
}

void update_TotalDensity(int step)
{
	totalDensity = 0.0;
	#if _OPENMP
	foreach(reduction(+:totalDensity))
	#else
	foreach()
	#endif
	{
		totalDensity = totalDensity + rho[];
	}
	Output_TotalDensity(step);
}
#endif