#include "parameters.hpp"

void compute_force_md() {
	double x[4], y[4], z[4], r, r1, ffx,ffy, ffz, mag_f, r2;
	int m;
	potential_colloid = 0;
	memset(f, 0, sizeof f);
	for(int i = 1; i <= no_of_colloid; i++) {
		for(int j = i; j <= n_neighbour[i]; j++) {
			m = neighbour[j][i];
			x[1] = pos_colloid[3*i - 2], x[2] = pos_colloid[3*m - 2];
			y[1] = pos_colloid[3*i - 1], y[2] = pos_colloid[3*m - 1];
			z[1] = pos_colloid[3*i], z[2] = pos_colloid[3*m];
			x[0] = mod(x[1] - x[2], lx), y[0] = mod(y[1] - y[2], lx), z[0] = mod(z[1] - z[2], lz);
			r2 = x[0]*x[0] + y[0]*y[0] + z[0]*z[0];
			if(r2 < pow(r_cutoff, 2)) {
				r = sqrt(r2);
				r1 = sig_colloid/r;

				potential_colloid += 4*eps*(pow(r1, 12) - pow(r1, 6)) - ufc + fc*r;
				mag_f = 4.0*eps*(12.0*pow(sig_colloid,12)/pow(r,13) - 6.0*sig_colloid/pow(r, 7));

				ffx = mag_f*x[0]/r, ffy = mag_f*y[0]/r, ffz = mag_f*z[0]/r;
				f[3*i - 2] += ffx, f[3*i - 1] += ffy, f[3*i] += ffz;
				f[3*m - 2] += ffx, f[3*m - 2] += ffy, f[3*m] += ffz;
			}
		}
	}
}
