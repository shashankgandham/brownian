//done
#include "parameters.hpp"
void update_velocity_colloid(){
	int i;
	double dtby2m;

	dtby2m = dt/(mass_colloid*2.0);
	for (i = 1; i <= no_of_colloid; i++){
		vel_colloid[3*i-2]=vel_colloid[3*i-2]+(old_force[3*i-2]+f[3*i-2])*dtby2m; 
        vel_colloid[3*i-1]=vel_colloid[3*i-1]+(old_force[3*i-1]+f[3*i-1])*dtby2m;
        vel_colloid[3*i]=vel_colloid[3*i]+(old_force[3*i]+f[3*i])*dtby2m;
	}
}
