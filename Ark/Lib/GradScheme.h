#ifndef _ARK_GRADSCHEME_H
#define _ARK_GRADSCHEME_H

const double apsi = 4.0,bpsi = 1.0,cpsi = 6.0;


double Grad_x(Point point, scalar s)
{
	return 
	(apsi*(s[1,0]-s[-1,0]) + bpsi*(s[1,1]-s[-1,1] + s[1,-1]-s[-1,-1]))/(cpsi*2.0);
}
double Grad_y(Point point, scalar s)
{
	return 
	(apsi*(s[0,1]-s[0,-1]) + bpsi*(s[1,1]-s[1,-1] + s[-1,1]-s[-1,-1]))/(cpsi*2.0);
}
double Laplacian(Point point, scalar s)
{
	return 
	(
	4.0*(s[1,0]+s[0,1]+s[-1,0]+s[0,-1])
	+   (s[-1,-1] + s[-1,1] + s[1,-1] + s[1,1])
	-   20.0*s[]
	)/6.0;
}
#endif