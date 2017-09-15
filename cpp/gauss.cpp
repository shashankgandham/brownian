#include "parameters.hpp"
#include <cstdlib>

void gauss(double v1, double v2) {
	double x1,x2,z = 2, sd;
	while(z > 1) {
		x1 = 2.0 * ((double)rand()/RAND_MAX) - 1;
		x2 = 2.0 * ((double)rand()/RAND_MAX) - 1;
		z=x1*x1+x2*x2;
	}
	z= sqrt((-2.0*log(z))/z);
	v1=x1*z*sqrt(kbt/mass_fl);
	v2=x2*z*sqrt(kbt/mass_fl);
}
