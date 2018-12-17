#ifndef _ARK_DISTRIBUTIONFUNC_FLIP
#define _ARK_DISTRIBUTIONFUNC_FLIP

/***************************************************************\
			   DF for phase index in cell fields
\***************************************************************/
scalar 
hTilde0[],hTilde1[],hTilde2[],hTilde3[],hTilde4[],hTilde5[],hTilde6[],hTilde7[],hTilde8[];

scalar
hBP0[],hBP1[],hBP2[],hBP3[],hBP4[],hBP5[],hBP6[],hBP7[],hBP8[];  

scalar *hTildeList = 
{hTilde0,hTilde1,hTilde2,hTilde3,hTilde4,hTilde5,hTilde6,hTilde7,hTilde8};

scalar *hBPList = 
{hBP0,hBP1,hBP2,hBP3,hBP4,hBP5,hBP6,hBP7,hBP8};

/***************************************************************\
			   DF for momentum in cell fields
\***************************************************************/
//auxiliary distribution function(DF),fTilde;
scalar  
fTilde0[],fTilde1[],fTilde2[],fTilde3[],fTilde4[],fTilde5[],fTilde6[],fTilde7[],fTilde8[];
//auxiliary DF, fBarPlus;		
scalar  
fBP0[],fBP1[],fBP2[],fBP3[],fBP4[],fBP5[],fBP6[],fBP7[],fBP8[];

scalar *fTildeList = 
{fTilde0,fTilde1,fTilde2,fTilde3,fTilde4,fTilde5,fTilde6,fTilde7,fTilde8};

scalar *fBPList = 
{fBP0,fBP1,fBP2,fBP3,fBP4,fBP5,fBP6,fBP7,fBP8};

/***************************************************************\
*					DF for flux in  face fields
\***************************************************************/
//micro flux of phase index
face vector 
hFlux0[],hFlux1[],hFlux2[],hFlux3[],hFlux4[],hFlux5[],hFlux6[],hFlux7[],hFlux8[];

vector *hFluxList = 
{hFlux0,hFlux1,hFlux2,hFlux3,hFlux4,hFlux5,hFlux6,hFlux7,hFlux8};

//micro-flux of momentum
face vector 
fFlux0[],fFlux1[],fFlux2[],fFlux3[],fFlux4[],fFlux5[],fFlux6[],fFlux7[],fFlux8[];

vector *fFluxList = 
{fFlux0,fFlux1,fFlux2,fFlux3,fFlux4,fFlux5,fFlux6,fFlux7,fFlux8};

vector *haloFluxList;


#endif