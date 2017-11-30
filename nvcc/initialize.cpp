#include "parameters.hpp"

void initialize() {
	point **pointers[] = {&pos_fl, &vel_fl, &f, &pos_colloid, &vel_colloid, &ang_vel_colloid, &old_force, &ra};
	n_neighbour = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
	iv 			= (int *)malloc(sizeof(int)*(ntab + 2));
	no_neigh 	= (int *)malloc((no_of_colloid + 2)*sizeof(int));

	for(int i = 0; i < 8; i++) {
		if(i < 2) *pointers[i] = (point *)malloc((no_of_fluid   + 2)*sizeof(point));
		else 	  *pointers[i] = (point *)malloc((no_of_colloid + 2)*sizeof(point));
	}
	for(int i = 0; i <= 10000; i++) {
		if(i <= 200) neighbour[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
		if(i <= 500) box_neigh[i] = (int *)malloc(sizeof(int)*(len.prod()    + 2));
		neigh_fl[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
	}
}

void initialize_colloid() {
	int counter = 0, check, nofp = 0;
	double space_limit = 1.3*sig_colloid, ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
	point avr_vel_colloid = point(0, 0, 0), t, temp;
	for(int k = 40; k <= (len.x-1)*10; k += 50) {
		for(int j = 40; j <= (len.y-1)*10; j+= 50) {
			for(int i = 40; i <= (len.z-1)*10; i += 50) {
				if(nofp < no_of_colloid)
					pos_colloid[++nofp] = point(i, j, k)/10.0;
			}
		}
	}
	while(counter < no_of_colloid) {
		t = point(ran()*len.x, ran()*len.y, ran()*len.z);
		check = 1;
		for(int j = 1; j <= counter; j++) {
			temp = img(t - pos_colloid[j], len);
			temp = point(abs(temp.x), abs(temp.y), abs(temp.z));
			check = ((temp*temp).sum() < space_limit)? 0: check;
		}
		if(check == 1) pos_colloid[++counter] = t;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = point(ran() - 0.5, ran() - 0.5, ran() - 0.5)*vscale_colloid;
		avr_vel_colloid += vel_colloid[j];
	}
	avr_vel_colloid = avr_vel_colloid/no_of_colloid;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j] - avr_vel_colloid;
		ang_vel_colloid[j] = point((ran() - 0.5), ran() - 0.5, ran() - 0.5)*ang_vscale_colloid;
	}
}

void initialize_fluid() {
	int counter = 0, check, count = 0;
	double vscale_fluid = sqrt(12.0*kbt/mass_fl);
	point avr_vel = point(0, 0, 0), t, temp;

	while(counter < no_of_fluid) {
		t = point(ran()*len.x, ran()*len.y, ran()*len.z);
		check = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < sigma*0.5)? 0: check;
		}
		if(check) pos_fl[++counter] = t;
	}
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[j]   = point(ran() - 0.5, ran() - 0.5 , ran() - 0.5)*vscale_fluid;
		avr_vel += vel_fl[j];
	}
	avr_vel = avr_vel/no_of_fluid;
	for(int j = 1; j <= no_of_fluid; j++)
		vel_fl[j] = vel_fl[j] - avr_vel;
}
