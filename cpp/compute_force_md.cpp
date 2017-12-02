#include "parameters.hpp"

void compute_force_md() {
	double mag_f, r_cutoff = pow(2, 1.0/6.0)*sig_colloid, r;
	double fc = 4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
	double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;
	point temp, ff;
	potential_colloid = 0, memset(f, 0, no_of_colloid + 2);
	for(int i = 1; i <= no_of_colloid; i++) {
		for(int j = i; j <= n_neighbour[i]; j++) {
			temp = mod(pos_colloid[i] - pos_colloid[neighbour[j][i]], len);
			r = sqrt((temp*temp).sum());
			if(r < r_cutoff) {
				potential_colloid += 4*eps*(pow(sig_colloid/r, 12) - pow(sig_colloid/r, 6)) - ufc + fc*r;
				mag_f = 4.0*eps*(12.0*pow(sig_colloid,12)/pow(r,13) - 6.0*sig_colloid/pow(r, 7));
				ff = (temp*mag_f)/r;
				f[i] += ff, f[neighbour[j][i]] += ff;
			}
		}
	}
}
