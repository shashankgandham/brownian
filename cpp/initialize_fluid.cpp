#include "parameters.hpp"
#include <cstdlib>
#include <cstdio>
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))
void initialize_fluid() {
	int counter_fl,check_fl;
	double x, y, z, ke1_fluid, avr_vel_fl_x, avr_vel_fl_y, avr_vel_fl_z, tx, ty, tz;
	double average_vel_fl_x, average_vel_fl_y, average_vel_fl_z;
	counter_fl = 0, avr_vel_fl_x = avr_vel_fl_y = avr_vel_fl_z = 0;
	while(counter_fl < no_of_fluid) {
		tx = ran()*lx, ty = ran()*ly, tz = ran()*lz;
		check_fl = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			x = mod(tx - pos_colloid[3*j-2], lx);
			y = mod(ty - pos_colloid[3*j-1], ly);
			z = mod(tz - pos_colloid[3*j], lz);
			check_fl = (x*x + y*y + z*z < sigma*0.5)? 0: check_fl;
		}
		if(check_fl) {
			counter_fl++;
			pos_fl[3*counter_fl-2] = tx;
			pos_fl[3*counter_fl-1] = ty;
			pos_fl[3*counter_fl] = tz;
		}
	}
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[3*j-2] = (ran() - 0.5)*vscale_fluid;
		vel_fl[3*j-1] = (ran() - 0.5)*vscale_fluid;
		vel_fl[3*j]   = (ran() - 0.5)*vscale_fluid;
		avr_vel_fl_x += vel_fl[3*j-2];
		avr_vel_fl_y += vel_fl[3*j-1];
		avr_vel_fl_z += vel_fl[3*j];
	}
	avr_vel_fl_x /= no_of_fluid;
	avr_vel_fl_y /= no_of_fluid;
	avr_vel_fl_z /= no_of_fluid;

	average_vel_fl_x = average_vel_fl_y = average_vel_fl_z = ke1_fluid = 0;
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[3*j-1] = vel_fl[3*j-1] - avr_vel_fl_x;
		vel_fl[3*j-2] = vel_fl[3*j-2] - avr_vel_fl_y;
		vel_fl[3*j]   = vel_fl[3*j] - avr_vel_fl_z;
		average_vel_fl_x += vel_fl[3*j-1];
		average_vel_fl_y += vel_fl[3*j-2];
		average_vel_fl_z += vel_fl[3*j];
		ke1_fluid = ke1_fluid + 0.50*mass_fl*(pow(vel_fl[3*j-1], 2) + pow(vel_fl[3*j-2], 2) + pow(vel_fl[3*j], 2));

	}
}
