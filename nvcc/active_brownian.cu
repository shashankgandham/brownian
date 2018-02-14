#include "parameters.cuh"
int main() {
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
		cudaDeviceSynchronize();
		for(int i = 1; i <= no_of_colloid; i++) {
			pos_colloid[i].print();
			vel_colloid[i].print();
			ang_vel_colloid[i].print();
		}
	}
	return 0;
}
