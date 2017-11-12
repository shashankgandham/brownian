#include "parameters.hpp"

int main() {
	int energy_colloid, nbox;
	double mom_x, mom_y, mom_z, potential_colloid, ke_colloid, ke_fluid, ang_ke_colloid;
	double I_colloid = 0.4*mass_colloid*sigma*sigma*0.25;
	initialize();
	initialize_colloid(I_colloid);
	initialize_fluid();
	nbox = create_box();
	neighbour_list_mpcd(nbox);
	compute_force_md();
	tumble();
	printf("After Tumble\n");
	for(int nn = 1; nn <= niter; nn++) {
		printf("%d\n", nn);
		rotation_mpcd();
		run();
		for(int l = 1; l <= n; l++) {
			std::copy(f, f + 3*no_of_colloid + 2, old_force);
			update_pos_md();
			neighbour_list_md();
			update_pos_mpcd();
			neighbour_list_mpcd(nbox);
			if(!(l%10) && nn > 10000)
				updown_velocity();
			fluid_colloid_collision(I_colloid);
			update_activity_direction();
			potential_colloid = compute_force_md();
			update_velocity_colloid();
		}
		ke_colloid = ke_fluid = ang_ke_colloid = 0;
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
		mom_x = mom_y = mom_z = 0;
		for(int i = 1; i <= no_of_fluid; i++) {
			mom_x += mass_fl*vel_fl[3*i - 2];
			mom_y += mass_fl*vel_fl[3*i - 1];
			mom_y += mass_fl*vel_fl[3*i];
		}
		for(int i = 1; i <= no_of_colloid; i++) {
			mom_x += mass_colloid*vel_colloid[3*i - 2];
			mom_y += mass_colloid*vel_colloid[3*i - 1];
			mom_z += mass_colloid*vel_colloid[3*i];
		}
	}
	return 0;
}
