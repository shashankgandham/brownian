#include "parameters.cuh"
int main() {
	point mom = point(0, 0, 0);
	double ke_colloid, ke_fluid, ang_ke_colloid, energy_colloid;
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
		mom = point(0, 0, 0);
		ke_colloid = ke_fluid = ang_ke_colloid = 0;
		for(int i = 1; i <= no_of_colloid; i++) {
			ke_colloid 	   += (vel_colloid[i]*vel_colloid[i]).sum();
			ang_ke_colloid += (ang_vel_colloid[i]*ang_vel_colloid[i]).sum();
			mom 		   += (vel_colloid[i]*mass_colloid);
			//pos_colloid[i].print();
		//	vel_colloid[i].print();
		//	ang_vel_colloid[i].print();
		}
		for(int i = 1; i <= no_of_fluid; i++) {
			ke_fluid += (vel_fl[i]*vel_fl[i]).sum();
			mom += (vel_fl[i]*mass_fl);
		}
		ke_colloid = 0.5*mass_colloid*ke_colloid;
		ang_ke_colloid = 0.5*I_colloid*ang_ke_colloid;
		energy_colloid = *potential_colloid + ke_colloid + ang_ke_colloid;
		ke_fluid = 0.5*ke_fluid*mass_fl;
	
		printf("%.32lf\n", energy_colloid);
	}
	return 0;
}
