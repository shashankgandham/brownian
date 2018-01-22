#include "parameters.hpp"

void initialize() {
	size_t devPitch;
	cudaError_t s1, s2, s3, s4, s5, s6, s7;
	point **ppointers[]  = {&pos_fl, &vel_fl, &dump_vel_fl, &f, &pos_colloid, &vel_colloid, &ang_vel_colloid, &old_force, &ra};
	int   **ipointers[]  = {&fluid_no, &n_neighbour, &no_neigh, &cnt, &up_cnt};
	int isize[]          = {(int)len.prod(), no_of_colloid};
 	int psize[]          = {no_of_fluid, no_of_colloid};
	s1 = cudaMallocPitch(&box_part,  &devPitch, sizeof(int)*(len.prod() + 2), maxpart + 2);	
	s2 = cudaMallocPitch(&cell_part, &devPitch, sizeof(int)*(len.prod() + 2), maxpart + 2);
	s3 = cudaMallocPitch(&box_neigh, &devPitch, sizeof(int)*(len.prod() + 2), 512);
	s4 = cudaMallocPitch(&nbr,       &devPitch, sizeof(int)*(no_of_colloid + 2), 7002);
	s5 = cudaMallocPitch(&up_nbr,    &devPitch, sizeof(int)*(no_of_colloid + 2), 7002);
	s6 = cudaMallocPitch(&neighbour, &devPitch, sizeof(int)*(no_of_colloid + 2), 256);
	s7 = cudaMallocPitch(&neigh_fl,  &devPitch, sizeof(int)*(no_of_colloid + 2), 10002);
	cell_vel = (point *)malloc((len.prod() + 1)*sizeof(point));
	iv = (int *)calloc(ntab + 2, sizeof(int));
	for(int i = 0; i < 9; i++) {
		if(i < 5) *ipointers[i] = (int   *)calloc(isize[i>0] + 2, sizeof(int)  );
				  *ppointers[i] = (point *)calloc(psize[i>2] + 2, sizeof(point));
	}
}

void initialize_colloid() {
	int counter = 0, check, nofp = 0, x = len.x, y = len.y, z = len.z;
	double space_limit = 1.3*sig_colloid, ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
	point avr_vel = point(0, 0, 0), t, temp, iter = point(4, 4, 4), lim = len - point(1, 1, 1);

	for(int i = 0; i <= lim.prod(); i += 5, iter.next(lim, point(5, 5, 5), point(4, 4, 4)), nofp++) {
		if(nofp < no_of_colloid) pos_colloid[nofp] = iter;
		else break;
	}
    int tempx = 0;
	while(counter < no_of_colloid) {
		t = t.random(point(0, 0, 0), len);
		check = 1;
		for(int j = 1; j <= counter; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < space_limit)? 0: check;
        }
		if(check)
			pos_colloid[++counter] = t;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j].random(point(0.5, 0.5, 0.5))*vscale_colloid;
		avr_vel += vel_colloid[j];
	}
	avr_vel = avr_vel/no_of_colloid;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j] - avr_vel;
		ang_vel_colloid[j] = (t.random(point(0.5, 0.5, 0.5)))*ang_vscale_colloid;
	}
}

void initialize_fluid() {
	int counter = 0, check;
	double vscale_fluid = sqrt(12.0*kbt/mass_fl);
	point avr_vel = point(0, 0, 0), t, temp;
	while(counter < no_of_fluid) {
		t = t.random(point(0, 0, 0), len);
		check = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < sigma*0.5)? 0: check;
		}
		if(check)
			pos_fl[++counter] = t;
	}
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[j] = vel_fl[j].random(point(0.5, 0.5, 0.5))*vscale_fluid;
		avr_vel += vel_fl[j];
	}
	avr_vel = avr_vel/no_of_fluid;
	for(int j = 1; j <= no_of_fluid; j++)
		vel_fl[j] = vel_fl[j] - avr_vel;
}
