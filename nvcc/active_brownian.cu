#include "parameters.h"
point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *ra, *old_force, len = point(30, 30, 30), *dump_vel_fl, *cell_vel;
int n = 10, niter = 21000, file = 0, nbin = 300, maxpart = 100, no_of_colloid = 10, nbox, **nbr, **up_nbr, *cnt, *up_cnt, *fluid_no, *iv, seed = 77777;
int no_of_fluid = len.prod()*10, *no_neigh, **neigh_fl, **neighbour, *n_neighbour, **box_neigh, **box_part, **cell_part, ntab = 32, nn;
double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.04;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid, I_colloid = 0.1*mass_colloid*sigma*sigma, potential_colloid;

int main() {
    double ke_colloid, ke_fluid, ang_ke_colloid, energy_colloid;
    point mom = point(0, 0, 0);
    initialize();
    initialize_colloid();
    initialize_fluid();
} 
