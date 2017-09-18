#include "parameters.hpp"
#include <cstring>
#include <cmath>
void neighbour_list_md() {
	double x1,x2,y1,y2,z1,z2,x12,y12,z12,r2;
	memset(n_neighbour, 0, sizeof(int)*(no_of_colloid+2));
	
	for(int i = 1; i < no_of_colloid; i++) {
		for(int j = i + 1; j <= no_of_colloid; j++) {
			x1=pos_colloid[3*i-2]; y1=pos_colloid[3*i-1]; z1=pos_colloid[3*i];
			x2=pos_colloid[3*j-2]; y2=pos_colloid[3*j-1]; z2=pos_colloid[3*j];
			x12=x1-x2; y12=y1-y2; z12=z1-z2;
			x12 = x12 - llx*round(x12*inv_llx);
			y12 = y12 - lly*round(y12*inv_lly);
			z12 = z12 - llz*round(z12*inv_llz);
			//Seriously someone learn about mod;
			r2=x12*x12+y12*y12+z12*z12;
			if(r2 < neigh_cutoff2){ 
				n_neighbour[i] = n_neighbour[i] + 1;
				neighbour[n_neighbour[i]][i] = j;
			}
		}
	}
}
