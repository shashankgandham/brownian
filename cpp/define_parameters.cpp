#include "parameters.hpp"
int n = 10, niter = 21000, file = 0, lx = 30, ly = 30, lz = 30 , nbin = 300, maxpart = 100, no_of_colloid = 1;
int NTAB = 32, seed = 77777, IY;
int no_of_fluid = lx*ly*lz*10, *no_neigh, *neigh_fl[10000], *neighbour[200], *n_neighbour, *box_neigh[512], *IV;
double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.04;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid;
double *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *ra, **dist, *old_force, *f;

double ran() {
	static int IM1 = 2147483563, IM2 = 2147483399, IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
	static int IMM1 = IM1-1, IR1 = 12211, IR2 = 3791, NTAB = 32, NDIV=1+IMM1/NTAB, IDUM2 = 123456789, J, K;
	double EPS=1.2e-7,RNMX=1.-EPS,AM=1./IM1;
	if(seed <= 0) {
		seed = std::max(-seed,1);
		IDUM2 = seed;
		for(J = NTAB + 8; J > 1; J--) {
			K = seed/IQ1;
			seed =IA1*(seed-K*IQ1)-K*IR1;
			if(seed < 0) seed += IM1;
			if(J <= NTAB) IV[J] = seed;
		}
		IY=IV[1];
	}
	K = seed/IQ1;
	seed = IA1*(seed - K*IQ1) - K*IR1;
	if(seed < 0) seed += IM1;
	K=IDUM2/IQ2;
	IDUM2 = IA2*(IDUM2 - K*IQ2) - K*IR2;
	if(IDUM2 < 0) IDUM2 = IDUM2 + IM2;

	J=1+IY/NDIV;
	IY = IV[J]-IDUM2;
	IV[J] = seed;
	if(IY < 1) IY=IY+IMM1;
	return std::min(AM*IY,RNMX);
}
