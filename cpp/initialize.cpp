#include "parameters.hpp"

void initialize_colloid(double I_colloid) {
	int counter, check, nofp = 0, x, y, z, r;
	double avr_vel_colloid_x, avr_vel_colloid_y, avr_vel_colloid_z, tx, ty, tz, space_limit = 1.3*sig_colloid;
	double ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);

	for(int k = 40; k <= (lx-1)*10; k += 50) {
		for(int j = 40; j <= (ly-1)*10; j+= 50) {
			for(int i = 40; i <= (lz-1)*10; i += 50) {
				nofp++;
				if(nofp <= no_of_colloid) {
					pos_colloid[3*nofp-2] = i/10.0;
					pos_colloid[3*nofp-1] = j/10.0;
					pos_colloid[3*nofp]   = k/10.0;
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
			x = (abs(x) > lx/2.0)? lx - abs(x): x;
			y = (abs(y) > ly/2.0)? ly - abs(y): y;
			z = (abs(z) > lz/2.0)? lz - abs(z): z;

			r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
			if(r < space_limit)
				check = 0;
		}
		if(check == 1) {
			counter++;
			pos_colloid[3*counter - 2] = tx;
			pos_colloid[3*counter - 1] = ty;
			pos_colloid[3*counter] = tz;
		}
	}

	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[3*j-2] = (ran() - 0.5)*vscale_colloid;
		vel_colloid[3*j-1] = (ran() - 0.5)*vscale_colloid;
		vel_colloid[3*j]   = (ran() - 0.5)*vscale_colloid;
		avr_vel_colloid_x += vel_colloid[3*j-2];
		avr_vel_colloid_y += vel_colloid[3*j-1];
		avr_vel_colloid_z += vel_colloid[3*j];
	}

	avr_vel_colloid_x /= no_of_colloid;
	avr_vel_colloid_y /= no_of_colloid;
	avr_vel_colloid_z /= no_of_colloid;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[3*j-2] = vel_colloid[3*j-2] - avr_vel_colloid_x;
		vel_colloid[3*j-1] = vel_colloid[3*j-1] - avr_vel_colloid_y;
		vel_colloid[3*j] = vel_colloid[3*j] - avr_vel_colloid_z;
		ang_vel_colloid[3*j-2] = (ran() - 0.5)*ang_vscale_colloid;
		ang_vel_colloid[3*j-1] = (ran() - 0.5)*ang_vscale_colloid;
		ang_vel_colloid[3*j]   = (ran() - 0.5)*ang_vscale_colloid;
	}
}

void initialize_fluid() {
	int counter_fl, check_fl;
	double x, y, z, avr_vel_fl_x, avr_vel_fl_y, avr_vel_fl_z, tx, ty, tz;
	double average_vel_fl_x, average_vel_fl_y, average_vel_fl_z, vscale_fluid = sqrt(12.0*kbt/mass_fl);
	counter_fl = 0, avr_vel_fl_x = avr_vel_fl_y = avr_vel_fl_z = 0;
	while(counter_fl < no_of_fluid) {
		tx = ran()*lx, ty = ran()*ly, tz = ran()*lz;
		check_fl = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			x = mod(tx - pos_colloid[3*j-2], lx);
			y = mod(ty - pos_colloid[3*j-1], ly);
			z = mod(tz - pos_colloid[3*j], lz);
			check_fl = ((x*x + y*y + z*z) < sigma*0.5)? 0: check_fl;
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
		printf("%lf %lf %lf %lf\n", vel_fl[3*j - 2], vel_fl[3*j - 1], vel_fl[3*j], vscale_fluid);
		exit(0);
		avr_vel_fl_x += vel_fl[3*j-2];
		avr_vel_fl_y += vel_fl[3*j-1];
		avr_vel_fl_z += vel_fl[3*j];
	}
	avr_vel_fl_x /= no_of_fluid;
	avr_vel_fl_y /= no_of_fluid;
	avr_vel_fl_z /= no_of_fluid;
	average_vel_fl_x = average_vel_fl_y = average_vel_fl_z = 0;

	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[3*j-2] = vel_fl[3*j-2] - avr_vel_fl_x;
		vel_fl[3*j-1] = vel_fl[3*j-1] - avr_vel_fl_y;
		vel_fl[3*j]   = vel_fl[3*j] - avr_vel_fl_z;
		average_vel_fl_x += vel_fl[3*j-2];
		average_vel_fl_y += vel_fl[3*j-1];
		average_vel_fl_z += vel_fl[3*j];
		printf("%lf %lf %lf\n", vel_fl[3*j - 2], vel_fl[3*j - 1], vel_fl[3*j]);
	}
}
