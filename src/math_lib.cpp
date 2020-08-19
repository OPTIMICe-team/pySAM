#include "math_lib.h"

#define emc 0.57721566490153286060651209008240243104215933593992 // Euler-Mascheroni constant
#define toll 1.e-18

double digammaR(double x)
{
/*
	double old_res=-emc;
	double new_res=old_res+(x-1.)/x;
	double err=abs((old_res-new_res)/old_res);
	double n=0.;
	while(err>toll)
	{
		n++;
		old_res=new_res;
		new_res=old_res+(x-1.)/((n+1.)*(n+x));
		err=abs((old_res-new_res)/old_res);
	}
	return new_res;
*/
	double result=0.;
	double z=x;
	while(z<10.)
	{
		result-=(1./z);
		z++;
	}
	result+=(log(z)-0.5/z-1./(12.*z*z)+1./(120.*z*z*z*z)-1./(252.*z*z*z*z*z*z));
}

double trigammaR(double x)
{
/*
	double tol=1e-6;
	double old_res=1./(x*x);
	double new_res=old_res+1./((x+1.)*(x+1.));
	double err=abs(old_res-new_res);
	double n=1.;
	while(err>toll)
	{
		n++;
		old_res=new_res;
		new_res=old_res+1./((x+n)*(x+n));
		err=abs((old_res-new_res)/old_res);
	}
	return new_res;
*/
	double result=0.;
	double z=x;
	while(z<10.)
	{
		result+=(1./(z*z));
		z++;
	}
	result+=(1./z+1./(2.*z*z)+1./(6.*z*z*z)-1./(30.*z*z*z*z*z)+1./(42.*z*z*z*z*z*z*z));
}
