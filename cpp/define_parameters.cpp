//Simulation Specifications
#include "parameters.hpp"

//Parameters
int n = 10, niter = 21000, file = 0, seed = 77777, lx = 30, ly = 30, lz = 30 , nbin = 300, maxpart = 100, no_of_colloid = 1;
double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.04;

//Formulae
int no_of_fluid = lx*ly*lz*10;
double dt = ndt/(double)n;
double r_cutoff = pow(2, 1.0/6.0)*sig_colloid, space_limit = 1.3*sig_colloid;
double fc = 4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6))+fc*r_cutoff;
double sigma = 0.80*sig_colloid;
double I_colloid = 0.4*mass_colloid*sigma*sigma*0.25;
double neigh_cutoff = 3.0*sig_colloid, length_cutoff = sigma*0.5 + 2;

//Variables
double *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *ra, **dist, *old_force, *f;
int *no_neigh, *neigh_fl[10000], *neighbour[200], *n_neighbour;

void initialize() {
	n_neighbour = (int *)calloc(sizeof(int), (no_of_colloid + 2));
	no_neigh = (int *)malloc((no_of_colloid + 2)*sizeof(int));
	pos_fl = (double *)malloc((3 * no_of_fluid + 2)*sizeof(double));
	vel_fl = (double *)malloc((3 * no_of_fluid + 2)*sizeof(double));
	f = (double *)malloc((3 * no_of_colloid + 2)*sizeof(double));
	old_force = (double *)malloc((3 * no_of_colloid + 2)*sizeof(double));
	pos_colloid = (double *)malloc((3 * no_of_colloid + 2)*sizeof(double));
	vel_colloid = (double *)malloc((3 * no_of_colloid + 2)*sizeof(double));
	ang_vel_colloid = (double *)malloc((3 * no_of_colloid + 2)*sizeof(double));
	ra = (double *)malloc((3 * no_of_colloid + 2)*sizeof(double));
	dist = (double **)malloc(sizeof(double *)*(no_of_colloid + 2));
	std::srand(seed);

	for(int i = 0; i < 10000; i++)
		neigh_fl[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
	for(int i = 0; i < 200; i++)
		neighbour[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
	for(int i = 0; i < no_of_colloid; i++)
		dist[i] = (double *)malloc(sizeof(double)*(no_of_colloid + 2));
}

double mod(double a, double b) {
	return (fmod(fmod(a, b) + b, b));
}
