#include "parameters.hpp"

void initialize() {
	point **pointers[] = {&pos_fl, &vel_fl, &f, &pos_colloid, &vel_colloid, &ang_vel_colloid, &old_force, &ra};
	n_neighbour = (int * )malloc(sizeof(int)*(no_of_colloid + 2));
	iv 			= (int * )malloc(sizeof(int)*(ntab + 2));
	no_neigh 	= (int * )malloc((no_of_colloid + 2)*sizeof(int));
	fluid_no 	= (int * )calloc(len.prod() + 2, sizeof(int));
	box_part 	= (int **)calloc(maxpart + 2, sizeof(int *));
	cell_part = (int **)malloc(sizeof(int *)*(maxpart  + 2));
	nbr = (int **)malloc(sizeof(int *) * 7005);
	up_nbr = (int **)malloc(sizeof(int *) * 7005);

	for(int i = 0; i < 8; i++)
		*pointers[i] = (point *)malloc((((i < 2)?no_of_fluid:no_of_colloid) + 2)*sizeof(point));

	for(int i = 0; i <= 10000; i++) {
		if(i <= 200)     neighbour[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
		if(i <= 500)     box_neigh[i] = (int *)malloc(sizeof(int)*(len.prod()    + 2));
		if(i < maxpart)  box_part[i]  = (int *)malloc((len.prod() + 1)*sizeof(int));
		if(i <= maxpart) cell_part[i] = (int *)malloc(sizeof(int)*(len.prod() + 2));
		if(i <= 7000)    nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
		if(i <= 7000)   up_nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
		neigh_fl[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
	}
}

void initialize_colloid() {
	int counter = 0, check, nofp = 0, x = len.x, y = len.y, z = len.z;
	double space_limit = 1.3*sig_colloid, ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
	point avr_vel = point(0, 0, 0), t, temp, iter = point(4, 4, 4), lim = len - point(1, 1, 1);

	for(int i = 0; i <= lim.prod(); i += 5, iter.next(lim, 5), nofp++) {
		if(nofp < no_of_colloid) pos_colloid[nofp] = iter;
		else break;
	}
	while(counter < no_of_colloid) {
		t = t.random()*len;
		check = 1;
		for(int j = 1; j <= counter; j++) {
			temp = abs(img(t - pos_colloid[j], len));
			check = ((temp*temp).sum() < space_limit)? 0: check;
		}
		if(check) pos_colloid[++counter] = t;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j].random(0.5)*vscale_colloid;
		avr_vel += vel_colloid[j];
	}
	avr_vel = avr_vel/no_of_colloid;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j] - avr_vel;
		ang_vel_colloid[j] = ang_vel_colloid[j].random(0.5)*ang_vscale_colloid;
	}
}

void initialize_fluid() {
	int counter = 0, check;
	double vscale_fluid = sqrt(12.0*kbt/mass_fl);
	point avr_vel = point(0, 0, 0), t, temp;
	while(counter < no_of_fluid) {
		t = t.random()*len;
		check = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < sigma*0.5)? 0: check;
		}
		if(check) pos_fl[++counter] = t;
	}
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[j] = vel_fl[j].random(0.5)*vscale_fluid;
		avr_vel += vel_fl[j];
	}
	avr_vel = avr_vel/no_of_fluid;
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[j] = vel_fl[j] - avr_vel;
	}
}
