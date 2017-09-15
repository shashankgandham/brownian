#include "parameters.hpp"
#include <cstdlib>

void initialize_activity() {
	double aaa; //seriously better names
	for(int i = 1; i <= no_of_colloid; i++) {
		ra[3*i - 2] = pos_colloid[3*i - 2] - rand()*llx; 
//And ran1--- variable or function?
		ra[3*i - 1] = pos_colloid[3*i - 1] - rand()*lly;
		ra[3*i] = pos_colloid[3*i] - rand()*llx;
		aaa = sqrt(pow(ra[3*i - 2],2) + pow(ra[3*i -1], 2) + pow(ra[3*i],2));

		for(int j = 0; j <= 2; j++) {
			ra[3*i - j] /= aaa;
			vel_colloid[3*i - j] += ra[3*i - j]*v0;
		}
	}
}
