#include "parameters.hpp"

void tumble(){
	for (int i = 1; i <= no_of_colloid; i++) {
		ra[3*i-2] = mod(pos_colloid[3*i-2] - ran()*lx, lx);
		ra[3*i-1] = mod(pos_colloid[3*i-1]- ran()*ly, ly);
		ra[3*i] = mod(pos_colloid[3*i] - ran()*lz, lz);
		double rax = sqrt(ra[3*i-2]*ra[3*i-2] + ra[3*i-1]*ra[3*i-1] + ra[3*i]*ra[3*i]);
		ra[3*i-2] /= rax, ra[3*i-1] /= rax, ra[3*i] /= rax;
	}
}
