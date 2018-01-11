#include "parameters.hpp"

point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *ra, *old_force, len = point(30, 30, 30);
int n = 10, niter = 21000, file = 0, nbin = 300, maxpart = 100, no_of_colloid = 1, nbox, ntab = 32, seed = 77777, **nbr, **up_nbr;
int no_of_fluid = len.prod()*10, *no_neigh, *neigh_fl[10005], *neighbour[256], *n_neighbour, *box_neigh[512], *iv, **box_part, *fluid_no, **cell_part;
double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.04;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid, I_colloid = 0.4*mass_colloid*sigma*sigma*0.25, potential_colloid;

//COMPUTE FORCE
double mag_f, r_cutoff = pow(2, 1.0/6.0)*sig_colloid, r;
double fc = 4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;

double ran() {
	static int im1 = 2147483563, im2 = 2147483399, ia1 = 40014, ia2 = 40692, iq1 = 53668, iq2 = 52774, iy;
	static int imm = im1 - 1, ir1 = 12211, ir2 = 3791, ntab = 32, ndiv = 1 + imm/ntab, idum = 123456789, j, k;
	double eps = 1.2e-7, rnmx = 1 - eps, am = 1.0/im1;
	if(seed <= 0) {
		seed = std::max(-seed,1);
		idum = seed;
		for(j = ntab + 8; j > 1; j--) {
			k = seed/iq1;
			seed = ia1*(seed - k*iq1) - k*ir1;
			if(seed < 0) seed += im1;
			if(j <= ntab) iv[j] = seed;
		}
		iy = iv[1];
	}
	k = seed/iq1;
	seed = ia1*(seed - k*iq1) - k*ir1;
	if(seed < 0) seed += im1;
	k = idum/iq2;
	idum = ia2*(idum - k*iq2) - k*ir2;
	if(idum < 0) idum = idum + im2;
	j = 1 + iy/ndiv;
	iy = iv[j] - idum;
	iv[j] = seed;
	if(iy < 1) iy = iy + imm;
	return std::min(am*iy, rnmx);
}

