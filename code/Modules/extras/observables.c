#include "phi4.h"
#include "math.h"

double mag;
int i;

double M()
{
    mag=0.;
    for(i=0;i<V;i++)
	{
	    mag+=phi[i];
	}

    return mag;
}

double Mabs()
{
    return fabs(M());
}

double Mtwo()
{
    return pow(M(),2);
}

double Mfour()
{
    return pow(M(),4);
}
