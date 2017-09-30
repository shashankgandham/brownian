#include "parameters.hpp"
#include <cstdio>
#include <cstdlib>

void initialize_colloid() {
	int counter, check, nofp = 0, x, y, z, r, mmm;
	double avr_vel_colloid_x, avr_vel_colloid_y, avr_vel_colloid_z, tx, ty, tz;
	double average_vel_colloid_x,average_vel_colloid_y,average_vel_colloid_z;

	for(int k = 40; k <= (lx-1)*10; k += 50) {
		for(int j = 40; j <= (ly-1)*10; j+= 50) {
			for(int i = 40; i <= (lz-1)*10; i += 50) {
				nofp++;
				if(nofp <= no_of_colloid){
					pos_colloid[3*nofp-2] = i/10.0;
					pos_colloid[3*nofp-1] = j/10.0;
					pos_colloid[3*nofp] = k/10.0;
				}
			}
		}
	}
	counter = avr_vel_colloid_x = avr_vel_colloid_y = avr_vel_colloid_z = 0;
	while(counter < no_of_colloid) {
		tx = ran()*lx, ty = ran()*ly, tz = ran()*lz;
		check = 1;
		for(int j = 1; j <= counter; j++) {
			x = tx - pos_colloid[3*j-2];
			y = ty - pos_colloid[3*j-1];
			z = tz - pos_colloid[3*j];
			if(abs(x) > lx/2.0)
				x = lx - abs(x);
			if(abs(y) > ly/2.0)
				y = ly - abs(y);
			if(abs(z) > lz/2.0)
				z = lz - abs(z);
			//MOD maybe?
			r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
			if(r < space_limit)
				check = 0;

		}
		if(check==1) {
			counter++;
			pos_colloid[3*counter - 2] = tx;
			pos_colloid[3*counter - 1] = ty;
			pos_colloid[3*counter] = tz;
		}
	}
	mmm = sqrt(kbt1/mass_colloid);
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[3*j-2] = (((double)rand()/RAND_MAX)-0.5)*vscale_colloid;
		vel_colloid[3*j-1] = (((double)rand()/RAND_MAX)-0.5)*vscale_colloid;
		vel_colloid[3*j] =   (((double)rand()/RAND_MAX)-0.5)*vscale_colloid;
		avr_vel_colloid_x = avr_vel_colloid_x + vel_colloid[3*j-2];
		avr_vel_colloid_y = avr_vel_colloid_y + vel_colloid[3*j-1];
		avr_vel_colloid_z = avr_vel_colloid_z + vel_colloid[3*j];

	}
	avr_vel_colloid_x = avr_vel_colloid_x/(float)no_of_colloid;
	avr_vel_colloid_y = avr_vel_colloid_y/(float)no_of_colloid;
	avr_vel_colloid_z = avr_vel_colloid_z/(float)no_of_colloid;
	average_vel_colloid_x = 0.0;
	average_vel_colloid_y = 0.0;
	average_vel_colloid_z = 0.0;
	ke_colloid1=0.0;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[3*j-2] = vel_colloid[3*j-2] - avr_vel_colloid_x;
		vel_colloid[3*j-1] = vel_colloid[3*j-1] - avr_vel_colloid_y;
		vel_colloid[3*j] = vel_colloid[3*j] - avr_vel_colloid_z;
		average_vel_colloid_x = average_vel_colloid_x + vel_colloid[3*j-2];
		average_vel_colloid_y = average_vel_colloid_y + vel_colloid[3*j-1];
		average_vel_colloid_z = average_vel_colloid_z + vel_colloid[3*j];
		ke_colloid1 = ke_colloid1 + pow(vel_colloid[3*j-2], 2) + pow(vel_colloid[3*j-1], 2) + pow(vel_colloid[3*j], 2);
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		ang_vel_colloid[3*j-2] = (((double)rand()/RAND_MAX)-0.5)*ang_vscale_colloid;
		ang_vel_colloid[3*j-1] = (((double)rand()/RAND_MAX)-0.5)*ang_vscale_colloid;
		ang_vel_colloid[3*j]   = (((double)rand()/RAND_MAX)-0.5)*ang_vscale_colloid;
	}
	ang_ke_colloid1=0.0;
	for(int j = 1; j <= no_of_colloid; j++)
		ang_ke_colloid1 = ang_ke_colloid1 + pow(ang_vel_colloid[3*j-2],2) + pow(ang_vel_colloid[3*j-1], 2) + pow(ang_vel_colloid[3*j], 2);
}
