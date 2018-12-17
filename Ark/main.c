#include "grid/quadtree.h"
//3x3 linear refine 
#include "./Header/additional-attributes.h"
//basic declaration
#include "BasiliskDUGKS.h"
//print split line or Error line
#include "./Lib/InforPrint.h"
//gradient scheme
#include "./Lib/GradScheme.h"
//Initialization function
#include "./Lib/Preprocess.h"
//output functions
#include "./Lib/Output.h"
//2D DUGKS solver
#include "./Lib/DUGKSSolver.h"
// specialised functions for macroscopic variables and microflux
#include "./Header/D2Q9AC.h"

int main()
{
	origin(X_Beg,Y_Beg);
	size(ChLength);
	//unrefine ((x < 50) || (x > 200) || (y < 50) || (y > 200)); 
	foreach_dimension()
    	periodic (right);
	init_grid(NL);

	#include "./Constant/ZeroCondition.h"

	DUGKSSolver();
	
	Output_Flowfield(step);
}