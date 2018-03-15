#include "parameters.cuh"

int main() {
	point mom;
	double ke_fluid, energy_colloid;
	initialize();
	initialize_rand();
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
		energy_colloid = *potential_colloid;
		energy_colloid += 0.5*mass_colloid*thrust::transform_reduce(thrust::seq, vel_colloid + 1, vel_colloid + no_of_colloid + 1, mod_value(), (double)0, add_double());
		energy_colloid += 0.5*I_colloid*thrust::transform_reduce(thrust::seq, ang_vel_colloid + 1, ang_vel_colloid + no_of_colloid + 1, mod_value(), (double)0, add_double());
		mom = thrust::reduce(thrust::seq, vel_colloid + 1, vel_colloid + no_of_colloid + 1, point(0, 0, 0), add_point())*mass_colloid;
		mom += thrust::reduce(thrust::seq, vel_fl + 1, vel_fl + no_of_fluid + 1, point(0, 0, 0), add_point())*mass_fl;
		ke_fluid = 0.5*mass_fl*thrust::transform_reduce(thrust::seq, vel_fl + 1, vel_fl + no_of_fluid + 1, mod_value(), (double)0, add_double());
		printf("%.32lf\n", ((*mom)*(*mom)).sum());
		printf("%.32lf\n", *energy_colloid);
	}
	return 0;
}