#include "parameters.cuh"
#include <cstring>
point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *ra, *old_force, len = point(30, 30, 30), *cell_vel, **rot, *dump_vel_fl, **u;
int n = 10, niter = 21000, file = 0, nbin = 300, maxpart = 100, no_of_colloid = 10, nbox, **nbr, **up_nbr, *cnt, *up_cnt, *fluid_no, *iv, *seed, *iy;
int no_of_fluid = len.prod()*10, *no_neigh, **neigh_fl, **neighbour, *n_neighbour, **box_neigh, **box_part, **cell_part, nn, ran_c = 0, *idum;

double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.04;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid, I_colloid = 0.1*mass_colloid*sigma*sigma, *potential_colloid, *rana, *ranb;

void initialize() {
	point **ppointers[]  = {&pos_fl, &vel_fl, &f, &pos_colloid, &vel_colloid, &ang_vel_colloid, &old_force, &ra};
	int   **ipointers[]  = {&fluid_no, &n_neighbour, &no_neigh, &cnt, &up_cnt};
	int isize[]          = {(int)len.prod(), no_of_colloid };
	int psize[]          = {no_of_fluid, no_of_colloid};
	cudaMallocManaged(&box_part,  (len.prod() + 2)*sizeof(int *));
	cudaMallocManaged(&cell_part, (len.prod() + 2)*sizeof(int *));
	cudaMallocManaged(&rot, (len.prod() + 2)*sizeof(point *));
	cudaMallocManaged(&cell_vel, (len.prod() + 2)*sizeof(point));
	cudaMallocManaged(&dump_vel_fl, (no_of_fluid + 2)*sizeof(point));
	cudaMallocManaged(&nbr, 7005*sizeof(int *));
	cudaMallocManaged(&up_nbr, 7005*sizeof(int *));
	cudaMallocManaged(&box_neigh, sizeof(int *)*512);
	cudaMallocManaged(&neighbour, sizeof(int *)*256);
	cudaMallocManaged(&neigh_fl,  sizeof(int *)*(no_of_colloid + 2));
	cudaMallocManaged(&u,  sizeof(point *)*(no_of_colloid + 2));
	cudaMallocManaged(&iv, sizeof(int)*64);
	cudaMallocManaged(&seed, sizeof(int));
	cudaMallocManaged(&idum, sizeof(int));
	cudaMallocManaged(&iy, sizeof(int));
	cudaMallocManaged(&potential_colloid, sizeof(double));
	cudaMallocManaged(&rana, sizeof(double)*(len.prod() + 2));
	cudaMallocManaged(&ranb, sizeof(double)*(len.prod() + 2));
	*seed = 77777, *idum = 123456789, *iy = 0;
	for(int i = 0; i < 64; i++) iv[i] = 0;
	for(int i = 0; i < 8; i++) {
		if(i < 5)  cudaMallocManaged(ipointers[i], (isize[i>0] + 2)*sizeof(int));
		cudaMallocManaged(ppointers[i], (psize[i>1] + 2)*sizeof(point));
	}
	for(int i = 0; i <= len.prod(); i++) {
		if(i <= 500)       cudaMallocManaged(&box_neigh[i], sizeof(int)*(len.prod()    + 2));
		if(i <= 200)       cudaMallocManaged(&neighbour[i], sizeof(int)*(no_of_colloid + 2));
		if(i <= 7000)      cudaMallocManaged(&nbr[i],       sizeof(int)*(no_of_colloid + 2));
		if(i <= 7000)      cudaMallocManaged(&up_nbr[i],    sizeof(int)*(no_of_colloid + 2));
		if(i <= no_of_colloid)	cudaMallocManaged(&neigh_fl[i],  sizeof(int)*(10000 + 2));
		if(i <= no_of_colloid)	cudaMallocManaged(&u[i],  sizeof(point)*(10000 + 2));
							cudaMallocManaged(&box_part[i],  sizeof(int)*(maxpart    + 2));
						   cudaMallocManaged(&cell_part[i], sizeof(int)*(maxpart    + 2));
						   cudaMallocManaged(&rot[i],		sizeof(point)*4);
	}
}
void initialize_colloid() {
	int counter = 0, check, nofp = 0;
	double space_limit = 1.3*sig_colloid, ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
	point avr_vel = point(0, 0, 0), t, temp, iter = point(4, 4, 4), lim = len - point(1, 1, 1);

	for(int i = 0; i <= lim.prod(); i += 5, iter.next(lim, point(5, 5, 5), point(4, 4, 4)), nofp++) {
		if(nofp < no_of_colloid) pos_colloid[++nofp] = iter;
		else break;
	}
	while(counter < no_of_colloid) {
		t = t.random(iv, seed, idum, iy)*len;
		check = 1;
		for(int j = 1; j <= counter; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < space_limit)? 0: check;
		}
		if(check)
			pos_colloid[++counter] = t;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = (vel_colloid[j].random(iv, seed, idum, iy) - point(0.5, 0.5, 0.5))*vscale_colloid;
		avr_vel += vel_colloid[j];
	}
	avr_vel = avr_vel/no_of_colloid;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j] - avr_vel;
		ang_vel_colloid[j] = (t.random(iv, seed, idum, iy) - point(0.5, 0.5, 0.5))*ang_vscale_colloid;
	}
}

void initialize_fluid() {
	int counter = 0, check;
	double vscale_fluid = sqrt(12.0*kbt/mass_fl);
	point avr_vel = point(0, 0, 0), t, temp;
	while(counter < no_of_fluid) {
		t = t.random(iv, seed, idum, iy)*len;
		check = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < sigma*0.5)? 0: check;
		}
		if(check) 
			pos_fl[++counter] = t;
	}
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[j] = (vel_fl[j].random(iv, seed, idum, iy) - point(0.5, 0.5, 0.5))*vscale_fluid;
		avr_vel += vel_fl[j];
	}
	avr_vel = avr_vel/no_of_fluid;
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[j] = vel_fl[j] - avr_vel;
	}
}