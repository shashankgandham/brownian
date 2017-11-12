#include "parameters.hpp"

double compute_force_md() {
	double x, y, z, r, ffx,ffy, ffz, mag_f, potential_colloid = 0;
	double r_cutoff = pow(2, 1.0/6.0)*sig_colloid;
	double fc = 4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
	double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;
	int m;
	memset(f, 0, 3*no_of_colloid + 2);
	for(int i = 1; i <= no_of_colloid; i++) {
		for(int j = i; j <= n_neighbour[i]; j++) {
			m = neighbour[j][i];
			x = mod(pos_colloid[3*i - 2] - pos_colloid[3*m - 2], lx);
			y = mod(pos_colloid[3*i - 1] - pos_colloid[3*m - 1], ly);
			z = mod(pos_colloid[3*i] - pos_colloid[3*m], lz);
			r = sqrt(x*x + y*y + z*z);
			if(r < r_cutoff) {
				potential_colloid += 4*eps*(pow(sig_colloid/r, 12) - pow(sig_colloid/r, 6)) - ufc + fc*r;
				mag_f = 4.0*eps*(12.0*pow(sig_colloid,12)/pow(r,13) - 6.0*sig_colloid/pow(r, 7));
				ffx = mag_f*x/r, ffy = mag_f*y/r, ffz = mag_f*z/r;
				f[3*i - 2] += ffx, f[3*i - 1] += ffy, f[3*i] += ffz;
				f[3*m - 2] += ffx, f[3*m - 1] += ffy, f[3*m] += ffz;
			}
		}
	}
	return potential_colloid;
}
