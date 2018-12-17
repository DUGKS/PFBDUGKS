#ifndef _ARK_D2Q9AC_H
#define _ARK_D2Q9AC_H

const char DmQnName[] = "D2Q9AC"; 

DmQn const xi[Q] = {{0,0},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};


double const xi_u[Q] = {0,1,1,0,-1,-1,-1,0,1};

double const xi_v[Q] = {0,0,1,1,1,0,-1,-1,-1};

int const _BB[Q] = {0,5,6,7,8,1,2,3,4};

const double omega[Q]=
{
	4.0/9.0,
 	1.0/9.0, 1.0/36.0,
 	1.0/9.0, 1.0/36.0,
 	1.0/9.0, 1.0/36.0,
 	1.0/9.0, 1.0/36.0
};

const double Kp = 0.6;

double SourcePhi(double phi)
{
	return -4.0*(phi-PhiL)*(phi-PhiV)/(wI*(PhiL-PhiV));
}
double ChemicalPotential(double phi,double laplacianPhi)
{
	return
	(
	2*Beta*(phi-PhiL)*(phi-PhiV)*(2*phi-PhiL-PhiV)
	-Kappa*laplacianPhi
	);
}
double calcRho(double phi)
{
	return aPhi*phi+bPhi;
}
double resetPhi(double phi)
{
	// if(phi < 0)
	// 	return 0;
	// if(phi > 1)
	// 	return 1;
	if(phi < 0.0001)
		return 0;
	if(phi > 0.9999)
		return 1;
}
inline double calcMu(double phi)
{
	return Tau0*RT*(aPhi*phi+bPhi);
}
void update_MacroVar()
{
	foreach()
	{
		int k = 0;
		//update macroscopic variables
		#ifdef _ARK_ALLENCAHN_FLIP
		phi[] = 0.0;
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		rho[] = 0.0;u[] = 0.0;v[] = 0.0;p[] = 0.0; 
		#endif

		k = 0;
		scalar hTilde,fTilde;
		for(hTilde,fTilde in hTildeList,fTildeList)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			phi[] += hTilde[];
			#endif

			#ifdef _ARK_MOMENTUM_FLIP
			p[] += fTilde[];
			u[] += fTilde[]*xi_u[k];
			v[] += fTilde[]*xi_v[k];
			#endif
			++k;
		}
		//phi[] = resetPhi(phi[]);
		rho[] = calcRho(phi[]);
	}
}
void update_Source()
{
	foreach()
	{
		// gradient of phi and rho
		GradX_phi[] = Grad_x(point,phi)/Delta;
		GradY_phi[] = Grad_y(point,phi)/Delta;
		GradX_rho[] = aPhi*GradX_phi[];
		GradY_rho[] = aPhi*GradY_phi[];

		#ifdef _ARK_MOMENTUM_FLIP
		//CSF
		double laplacianPhi = 0.0,modF = 0.0;
		laplacianPhi = Laplacian(point,phi)/(Delta*Delta);
		modF = ChemicalPotential(phi[],laplacianPhi);
		ForceX[] = modF*GradX_phi[];
		ForceY[] = modF*GradY_phi[];
		//body force, !!pay attention
		ForceX[] += Gx*rho[];
		ForceY[] -= Gy*(rho[] - RhoL);

		u[] += hDt*ForceX[]; v[] += hDt*ForceY[];
		u[] /= rho[]; v[] /= rho[];
		double uu = u[]*u[] + v[]*v[];
		//dynamic pressure
		p[] -= fTilde0[];
		p[] += hDt*(u[]*GradX_rho[] + v[]*GradY_rho[])-rho[]*omega[0]*uu*Lambda0;
		p[] *= Kp;
		#endif
		//
		double L = sqrt(GradX_phi[]*GradX_phi[] + GradY_phi[]*GradY_phi[]);
		double partial_t_phiu = (phi[]*u[] - prev_phiu[])/Dt;
		double partial_t_phiv = (phi[]*v[] - prev_phiv[])/Dt;
		if(L != 0)
		{
			double modPhi = SourcePhi(phi[]);
			GradX_phi[] = GradX_phi[]*modPhi/L + partial_t_phiu/RT;
			GradY_phi[] = GradY_phi[]*modPhi/L + partial_t_phiv/RT;
		}
		else
		{
			GradX_phi[] = partial_t_phiu/RT;
			GradY_phi[] = partial_t_phiv/RT;
		}
		prev_phiu[] = phi[]*u[];
		prev_phiv[] = phi[]*v[];
	}
}
void update_xEq(Point point,double hEq[],double hS[],double fEq[],double fS[])
{
	double u1,GAMMA;
	double uu = u[]*u[] + v[]*v[];
	int k = 0;
	for(k = 0;k < Q;++k)
	{
		u1 = (xi_u[k]*u[] + xi_v[k]*v[])/RT;
		GAMMA = u1 + 0.5*u1*u1 - uu*Lambda0;
		//
		#ifdef _ARK_ALLENCAHN_FLIP
		hEq[k] = omega[k]*(1.0 + u1)*phi[];
		hS[k]  = omega[k]*(xi_u[k]*GradX_phi[] + xi_v[k]*GradY_phi[]);
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		fEq[k] = omega[k]*(p[]/RT + rho[]*GAMMA) - (int)(!k)*p[]/RT;
		fS[k]  = omega[k]*
		(
			(1+GAMMA)*((xi_u[k]-u[])*ForceX[] + (xi_v[k]-v[])*ForceY[])/RT
			+
			GAMMA*((xi_u[k]-u[])*GradX_rho[] + (xi_v[k]-v[])*GradY_rho[])
		);
		#endif
	}
}
void MicroFlux_2D()
{
	foreach_face()
	{		
		//Index
		int k = 0;
		/**
		Discrete Distribution Function; 
		xBh  : xBar at half time; Latex Code : \bar{x}(x,t+h)
		xEqh : Equilibrium DF at half time;
		xSh  : Source DF at half time;
		*/
		double hBh[Q] = {0.0}, hEqh[Q] = {0.0}, hSh[Q] = {0.0};
		double fBh[Q] = {0.0}, fEqh[Q] = {0.0}, fSh[Q] = {0.0};
		//DF bar at half time; Central Scheme
		double GradhBP_x = 1000, GradhBP_y = 1000;
		double GradfBP_x = 1000, GradfBP_y = 1000;

		scalar hBP,fBP;
		k = 0;
		for(hBP,fBP in hBPList,fBPList)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			GradhBP_x = (hBP[]-hBP[-1])/Delta;
			GradhBP_y = 0.25*(hBP[0,1]-hBP[0,-1]+hBP[-1,1]-hBP[-1,-1])/Delta;
			hBh[k] = 0.5*(hBP[]+hBP[-1])-hDt*(xi[k].x*GradhBP_x + xi[k].y*GradhBP_y);
			#endif

			#ifdef _ARK_MOMENTUM_FLIP
			GradfBP_x = (fBP[]-fBP[-1])/Delta;
			GradfBP_y = 0.25*(fBP[0,1]-fBP[0,-1]+fBP[-1,1]-fBP[-1,-1])/Delta;
			fBh[k] = 0.5*(fBP[]+fBP[-1])-hDt*(xi[k].x*GradfBP_x + xi[k].y*GradfBP_y);
			#endif

			++k;
		}
		//macro variables
		double phih = 0.0,GradX_phih = 0.0,GradY_phih = 0.0;
		double rhoh = 0.0,GradX_rhoh = 0.0,GradY_rhoh = 0.0;
		double ph = 0.0;
		double uh = 0.0,vh = 0.0;
		double ForceXh = 0.0,ForceYh = 0.0;
		double uu = 0.0;

		//macroscopic varaibels at faces
		GradX_phih = (GradX_phi[] + GradX_phi[-1])*0.5;
		GradY_phih = (GradY_phi[] + GradY_phi[-1])*0.5;
		//
		ForceXh = (ForceX[]+ForceX[-1])*0.5;
		ForceYh = (ForceY[]+ForceY[-1])*0.5;
		//
		GradX_rhoh = (GradX_rho[] + GradX_rho[-1])*0.5;
		GradY_rhoh = (GradY_rho[] + GradY_rho[-1])*0.5;
		// density
		for(k = 0;k < Q;++k)
		{
			phih += hBh[k];
		}
		rhoh = calcRho(phih);
		//momentum
		for(k = 0;k < Q;++k)
		{
			ph += fBh[k];
			uh += fBh[k]*xi_u[k];
			vh += fBh[k]*xi_v[k];
		}
		//velocity
		uh += 0.5*hDt*ForceXh; vh += 0.5*hDt*ForceYh;
		uh /= rhoh; vh /= rhoh;
		uu = uh*uh + vh*vh;
		//dynamic pressure
		ph -= fBh[0];
		ph += 0.5*hDt*(GradX_rhoh*uh + GradY_rhoh*vh)
			  -omega[0]*rhoh*Lambda0*uu;
		ph *= Kp;

		//Equilibrium DF and Source term
		double u1,GAMMA;
		for(k = 0;k < Q;++k)
		{
			u1 = (xi_u[k]*uh + xi_v[k]*vh)/RT;
			GAMMA = u1 + 0.5*u1*u1 - uu*Lambda0;
		
			#ifdef _ARK_ALLENCAHN_FLIP
			hEqh[k] = omega[k]*phih*(1.0+u1);
			hSh[k]  = omega[k]*(GradX_phih*xi_u[k] + GradY_phih*xi_v[k]);
			#endif
	
			#ifdef _ARK_MOMENTUM_FLIP
			fEqh[k] = omega[k]*(ph/RT + rhoh*GAMMA)-(int)(!k)*ph/RT;
			fSh[k]  = omega[k]*
			(
				(1+GAMMA)*((xi_u[k]-uh)*ForceXh + (xi_v[k]-vh)*ForceYh)/RT
			+	GAMMA*((xi_u[k]-uh)*GradX_rhoh + (xi_v[k]-vh)*GradY_rhoh)
			);
			#endif
		}
		//Original DF
		#ifdef _ARK_ALLENCAHN_FLIP
		double ah = 2.0*tauM/(2.0*tauM + hDt), bh = 1.0 - ah, ch = tauM*bh;
		#endif

		#ifdef _ARK_MOMENTUM_FLIP
		double fMuh = calcMu(phih), ftauh = fMuh/(rhoh*RT);
		double af = 2.0*ftauh/(2.0*ftauh + hDt), bf = 1.0 - af, cf = ftauh*bf;
		#endif
		//
		k = 0;
		vector hFlux,fFlux;
		for(hFlux,fFlux in hFluxList,fFluxList)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			hFlux.x[] = ah*hBh[k] + bh*hEqh[k] + ch*hSh[k];
			#endif

			#ifdef _ARK_MOMENTUM_FLIP
			fFlux.x[] = af*fBh[k] + bf*fEqh[k] + cf*fSh[k];
			#endif

			++k;
		}
	}
}
#endif