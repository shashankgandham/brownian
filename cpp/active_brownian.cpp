#include "parameters.hpp"

point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *ra, *old_force, *dump_vel_fl, len = point(30, 30, 30);
int n = 10, niter = 21000, file = 0, nbin = 300, maxpart = 100, no_of_colloid = 10, nbox, **nbr, **up_nbr, *cnt, *up_cnt, *fluid_no, *iv, seed = 77777;
int no_of_fluid = len.prod()*10, *no_neigh, **neigh_fl, *neighbour[256], *n_neighbour, *box_neigh[512], **box_part, **cell_part, ntab = 32, nn;
double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.0;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid, I_colloid = 0.1*mass_colloid*sigma*sigma, potential_colloid;

int main(int argv, char *argc[]) {
	clock_t begin = clock();
	no_of_fluid = len.prod()*atoi(argc[1]), no_of_colloid = atoi(argc[2]);
	double ke_colloid, ke_fluid, ang_ke_colloid, energy_colloid;
	point mom = point(0, 0, 0);
	initialize();
	initialize_colloid();
	initialize_fluid();
	create_box();
	neighbour_list_mpcd();
	neighbour_list_md();
	compute_force_md();
	tumble();
	printf(" After Tumble\n");
	for(nn = 1; nn <= niter; nn++) {
		printf("%12d\n", nn);
		rotation_mpcd();
		run();
		for(int l = 1; l <= n; l++) {
			std::copy(f, f + no_of_colloid + 2, old_force);
			update_pos_md();
			neighbour_list_md();
			update_pos_mpcd();
			neighbour_list_mpcd();
			if(!(l%10) && nn > 10000) updown_velocity();
			fluid_colloid_collision();
			update_activity_direction();
			compute_force_md();
			update_velocity_colloid();
		}
		ke_colloid = ke_fluid = ang_ke_colloid = 0;
		for(int i = 1; i <= no_of_colloid; i++) {
			ke_colloid 	   += (vel_colloid[i]*vel_colloid[i]).sum();
			ang_ke_colloid += (ang_vel_colloid[i]*ang_vel_colloid[i]).sum();
			mom 		   += (vel_colloid[i]*mass_colloid);
		//	pos_colloid[i].print();
		//	vel_colloid[i].print();
		//	ang_vel_colloid[i].print();
		}
		for(int i = 1; i <= no_of_fluid; i++) {
			ke_fluid += (vel_fl[i]*vel_fl[i]).sum();
			mom += (vel_fl[i]*mass_fl);
		}
		ke_colloid = 0.5*mass_colloid*ke_colloid;
		ang_ke_colloid = 0.5*I_colloid*ang_ke_colloid;
		energy_colloid = potential_colloid + ke_colloid + ang_ke_colloid;
		ke_fluid = 0.5*ke_fluid*mass_fl;
<<<<<<< HEAD
		printf("%.32lf\n", energy_colloid + ke_fluid);
=======
	//	printf("%.32lf\n", energy_colloid + ke_fluid);
>>>>>>> 11667a0cdeb0c037f9cdd4cb05ecb67d7f95da62
	}
	clock_t end = clock();
	printf("%lf\n", (double)(end - begin)/CLOCKS_PER_SEC);
	return 0;
}
