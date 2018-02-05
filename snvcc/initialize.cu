#include "parameters.cuh"

point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *ra, *old_force, len = point(30, 30, 30);
int n = 10, niter = 21000, file = 0, nbin = 300, maxpart = 100, no_of_colloid = 10, nbox, **nbr, **up_nbr, *cnt, *up_cnt, *fluid_no;
int no_of_fluid = len.prod()*10, *no_neigh, *neigh_fl[10005], *neighbour[256], *n_neighbour, *box_neigh[512], **box_part, **cell_part, nn;
double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0.04;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid, I_colloid = 0.1*mass_colloid*sigma*sigma, potential_colloid;

void initialize() {
	point **ppointers[]  = {&pos_fl, &vel_fl, &f, &pos_colloid, &vel_colloid, &ang_vel_colloid, &old_force, &ra};
	int   **ipointers[]  = {&fluid_no, &n_neighbour, &no_neigh, &cnt, &up_cnt};
	int isize[]          = {(int)len.prod(), no_of_colloid };
	int psize[]          = {no_of_fluid, no_of_colloid};

	box_part 	= (int **)calloc((maxpart + 2),sizeof(int *));
	cell_part 	= (int **)calloc((maxpart + 2),sizeof(int *));
	nbr 		= (int **)calloc(7005,sizeof(int *));
	up_nbr 		= (int **)calloc(7005,sizeof(int *));

	for(int i = 0; i < 8; i++) {
		if(i < 5) *ipointers[i] = (int   *)calloc(isize[i>0] + 2, sizeof(int)  );
				  *ppointers[i] = (point *)calloc(psize[i>1] + 2, sizeof(point));
	}
	for(int i = 0; i <= 10000; i++) {
		if(i <= 500)      box_neigh[i] = (int *)calloc(sizeof(int),(len.prod()    + 2));
		if(i <= maxpart)  box_part[i]  = (int *)calloc(sizeof(int),(len.prod()    + 2));
		if(i <= maxpart)  cell_part[i] = (int *)calloc(sizeof(int),(len.prod()    + 2));
		if(i <= 200)      neighbour[i] = (int *)calloc(sizeof(int),(no_of_colloid + 2));
		if(i <= 7000)     nbr[i]       = (int *)calloc(sizeof(int),(no_of_colloid + 2));
		if(i <= 7000)     up_nbr[i]    = (int *)calloc(sizeof(int),(no_of_colloid + 2));
						  neigh_fl[i]  = (int *)calloc(sizeof(int),(no_of_colloid + 2));
	}
}
__global__ void d_initialize_colloid(point *pos_colloid, point *vel_colloid, point *ang_vel_colloid, double sig_colloid, double kbt1, 
									double I_colloid, double mass_colloid, int no_of_colloid, point len) {
	int counter = 0, check, nofp = 0;
	double space_limit = 1.3*sig_colloid, ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
	point avr_vel = point(0, 0, 0), t, temp, iter = point(4, 4, 4), lim = len - point(1, 1, 1);

	for(int i = 0; i <= lim.prod(); i += 5, iter.next(lim, point(5, 5, 5), point(4, 4, 4)), nofp++) {
		if(nofp < no_of_colloid) pos_colloid[++nofp] = iter;
		else break;
	}
	while(counter < no_of_colloid) {
		t = t.random(point(0, 0, 0), len);
		check = 1;
		for(int j = 1; j <= counter; j++) {
			d_img(&temp, t - pos_colloid[j], len);
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
void initialize_colloid() {
	point *d_pos_colloid, *d_vel_colloid, *d_ang_vel_colloid;
	cudaMalloc(&d_pos_colloid, sizeof(point)*(no_of_colloid + 2));
	cudaMalloc(&d_vel_colloid, sizeof(point)*(no_of_colloid + 2));
	cudaMalloc(&d_ang_vel_colloid, sizeof(point)*(no_of_colloid + 2));
	d_initialize_colloid<<<1, 1>>> (d_pos_colloid, d_vel_colloid, d_ang_vel_colloid, sig_colloid, kbt1, 
							I_colloid, mass_colloid, no_of_colloid, len); 
	cudaMemcpy(vel_colloid, d_vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaMemcpy(ang_vel_colloid, d_ang_vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaMemcpy(pos_colloid, d_pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_ang_vel_colloid), cudaFree(d_vel_colloid), cudaFree(d_pos_colloid);
}

__global__ void d_initialize_fluid(point *pos_fl, point *vel_fl, point *pos_colloid, double kbt, double mass_fl, 
								double sigma, int no_of_fluid, int no_of_colloid, point len) {
	int counter = 0, check;
	double vscale_fluid = sqrt(12.0*kbt/mass_fl);
	point avr_vel = point(0, 0, 0), t, temp;
	while(counter < no_of_fluid) {
		t = t.random(point(0, 0, 0), len);
		check = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			d_img(&temp, t - pos_colloid[j], len);
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

void initialize_fluid() {
	point *d_pos_fl, *d_vel_fl, *d_pos_colloid;
	cudaMalloc(&d_pos_fl, sizeof(point)*(no_of_fluid + 2));
	cudaMalloc(&d_vel_fl, sizeof(point)*(no_of_fluid + 2));
	cudaMalloc(&d_pos_colloid, sizeof(point)*(no_of_colloid + 2));
	cudaMemcpy(d_pos_colloid, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	d_initialize_fluid<<<1, 1>>> (d_pos_fl, d_vel_fl, d_pos_colloid, kbt, mass_fl, sigma, no_of_fluid, no_of_colloid, len); 
	cudaMemcpy(vel_fl, d_vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaMemcpy(pos_fl, d_pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_pos_fl), cudaFree(d_vel_fl), cudaFree(d_pos_colloid);
}