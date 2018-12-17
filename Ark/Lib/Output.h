#ifndef _ARK_OUTPUT_H
#define _ARK_OUTPUT_H

extern char caseName[];

extern double totalDensity;

extern double residualPer1k;

void Output_Flowfield(int step)
{
	#if _OPENMP
	omp_set_num_threads(1);
	#endif
	
	int cellNum = 0;
    foreach()
    {
		cellNum = cellNum + 1;
    }
  	//
  	char Name_flowfield[256];memset(Name_flowfield,0,256);
  	sprintf(Name_flowfield,"../FlowField/global/step%d_CHL%g_Ma%g_%s_t%g.dat",
  		step,ChLength,Ma,caseName,step*Dt);

  	FILE *OutFile_flowfield;
  	OutFile_flowfield = fopen(Name_flowfield,"w");
  	fprintf(OutFile_flowfield, "variables = x,y,phi,rho,u,v,p,phi_x,phi_y,Error\n");
  	fprintf(OutFile_flowfield, "Zone T = %s_step%d\n", caseName,step);
  	fprintf(OutFile_flowfield, "N = %d, E = %d\n",4*cellNum, cellNum);
	fprintf(OutFile_flowfield,"DATAPACKING = BLOCK, ZONETYPE = FEQUADRILATERAL\n");
	fprintf(OutFile_flowfield, "VARLOCATION = \
	([1-2]=NODAL,[3-10]=CellCentered)\n");
	foreach()
	{
		fprintf(OutFile_flowfield,"%f\t%f\t%f\t%f\n",
			x-0.5*Delta,x+0.5*Delta,x+0.5*Delta,x-0.5*Delta);	
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%f\t%f\t%f\t%f\n",
			y-0.5*Delta,y-0.5*Delta,y+0.5*Delta,y+0.5*Delta);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",phi[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",rho[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",u[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",v[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",p[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",GradX_phi[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",GradY_phi[]);
	}
	foreach()
	{
		fprintf(OutFile_flowfield,"%.10e\n",Error_1k[]);
	}
	for(int iii = 0;iii < cellNum; ++iii)
	{
		fprintf(OutFile_flowfield, "%d %d %d %d\n",
		0 + iii * 4 + 1, 1 + iii * 4 + 1, 2 + iii * 4 + 1, 3 + iii * 4 + 1);
	}
	fclose(OutFile_flowfield);

	#if _OPENMP
	omp_set_num_threads(ThreadNum);
	#endif
}
void Output_TotalDensity(int step)
{
	char Name_TotalDensity[256];memset(Name_TotalDensity,0,256);
	sprintf(Name_TotalDensity,"../FlowField/Convergence/SumRho_mu%g_Re%g_%s.dat",
		Mu0,Re,caseName);
	FILE *OutFile_TotalDensity = fopen(Name_TotalDensity,"a");
	fprintf(OutFile_TotalDensity, "%.6e\t%.8e\n",step*Dt,totalDensity);
	fclose(OutFile_TotalDensity);
}
void Output_Residual(int step)
{
	char Name_Residual[256];memset(Name_Residual,0,256);
	sprintf(Name_Residual,"../FlowField/Convergence/Residual_mu%g_Re%g_%s.dat",
		Mu0,Re,caseName);

	FILE *OutFile_Residual = fopen(Name_Residual,"a");
	fprintf(OutFile_Residual, "%.6e\t%.8e\n",step*Dt,residualPer1k);
	fclose(OutFile_Residual);
}
#endif