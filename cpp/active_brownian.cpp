#include "parameters.hpp"

int main() {
	initialize();
	int energy_colloid;
	double mom1_x, mom1_y, mom1_z, mom2_x, mom2_y, mom2_z;

	if(!file) {
		initialize_colloid();
		initialize_fluid();
	}
	create_box();
	neighbour_list_mpcd();
	compute_force_md();
	tumble();
	printf("After Tumble\n");
	for(nn = 1; nn <= niter; nn++) {
		printf("%d\n", nn);
		rotation_mpcd();
		run();
		for(int l = 1; l <= n; l++) {
			qq = (nn  - 1)*10 + l;
			std::copy(f, f + 3*no_of_colloid + 2, old_force);
			update_pos_md();
			neighbour_list_md();
			update_pos_mpcd();
			neighbour_list_mpcd();
			if(!(qq%10) && nn > 10000)
				updown_velocity();
			fluid_colloid_collision();
			update_activity_direction();
			compute_force_md();
			update_velocity_colloid();
		}
		ke_colloid = ke_fluid = ang_ke_colloid = 0;
		mom2_x = mom2_y = mom2_z = 0;
		for(int i = 1; i <= no_of_colloid; i++) {
			for(int j = 0; j <= 2; j++) {
				ke_colloid += vel_colloid[3*i - j]*vel_colloid[3*i - j];
				ang_ke_colloid += ang_vel_colloid[3*i - j]*ang_vel_colloid[3*i - j];
			}
		}
		ke_colloid = 0.5*mass_colloid*ke_colloid;
		ang_ke_colloid = 0.5*I_colloid*ang_ke_colloid;
		energy_colloid = potential_colloid + ke_colloid + ang_ke_colloid;

		for(int i = 1; i <= no_of_fluid; i++) {
			for(int j = 0; j <= 2; j++)
				ke_fluid += vel_fl[3*i - j]*vel_fl[3*i - j];
		}
		ke_fluid = 0.5*ke_fluid*mass_fl;
		mom1_x = mom1_y = mom1_z = 0;
		mom2_x = mom2_y = mom2_z = 0;
		for(int i = 1; i <= no_of_fluid; i++) {
			mom1_x += mass_fl*vel_fl[3*i - 2];
			mom1_y += mass_fl*vel_fl[3*i - 1];
			mom1_y += mass_fl*vel_fl[3*i];
		}
		for(int i = 1; i <= no_of_colloid; i++) {
			mom2_x += mass_colloid*vel_colloid[3*i - 2];
			mom2_y += mass_colloid*vel_colloid[3*i - 1];
			mom2_z += mass_colloid*vel_colloid[3*i];
		}
		mom_x = mom1_x + mom2_x;
		mom_y = mom1_y + mom2_y;
		mom_z = mom1_z + mom2_z;
	}
	return 0;
}
