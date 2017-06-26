#include <cstdio>
#include <cmath>
//import all parameters

void compute_force_md() {
	double x1, x2, y1, y2, z1, z2, y12, z12, r, r1, ffx;
	double ffy, ffz, mag_f, r2, pot1;
	
	potential_colloid = (double)0;
	f = (double)0;
	for(int i = 1; i <= no_of_colloid; i++) {
		for(int j = i; j <= n_neighbour[i]; j++) {
			m = neighbour[j][i];
			x1 = pos_colloid[3*i - 2], x2 = pos_colloid[3*m - 2];
			y1 = pos_colloid[3*i - 1], y2 = pos_colloid[3*m - 1];
			z1 = pos_colloid[3*i], z2 = pos_colloid[3*m];
			x12 -= llx*round(x12*inv_llx);
			y12 -= llx*round(y12*inv_llx);
			z12 -= llx*round(z12*inv_llx);
			r2 = x12*x12 + y12*y12 + z12*z12
			if(r2 < r_cuttoff2) {
				r = sqrt(r2);
				r1 = sig_colloid/r;
				pot1 = (4.0*eps*(12.0*sig_colloid12/pow(r,13) - 6.0*sig_colloid/pow(r, 7));
				ffx = mag_f*x12/r;
				ffy = mag_f*y12/r;
				ffz = mag_f*z12/r;
				f[3*i - 2] += ffx, f[3*i - 1] += ffy, f[3*i] += ffz;
				f[3*m - 2] += ffx, f[3*m - j] += ffy, f[3*m - j] += ffz;
			}		
		}
	}
}	
