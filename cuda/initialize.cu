#include "parameters.cuh"

curandState_t *state;
point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *ra, *old_force, len = point(30, 30, 30), *cell_vel, **rot, *dump_vel_fl, **u, **vc, **om;
int n = 10, niter = 21000, file = 0, nbin = 300, maxpart = 100, no_of_colloid = 10, nbox, **nbr, **up_nbr, *cnt, *up_cnt, *fluid_no, **dp;
int no_of_fluid = len.prod()*10, *no_neigh, **neigh_fl, **neighbour, *n_neighbour, **box_neigh, **box_part, **cell_part, nn, seed = 77777;

double kbt = 1, kbt1 = 1, ndt = 0.1, dv = 0.1, mass_fl = 1.0, mass_colloid = 654.1, sig_colloid = 5.0, eps = 1.0, v0 = 0;
double dt = ndt/(double)n, sigma = 0.80*sig_colloid, I_colloid = 0.1*mass_colloid*sigma*sigma, *potential_colloid, *rana, *ranb;

dim3 thr(512), thrs(32, 32), blk;
void initialize() {
	point **ppointers[]  = {&pos_fl, &vel_fl, &f, &pos_colloid, &vel_colloid, &ang_vel_colloid, &old_force, &ra};
	int   **ipointers[]  = {&fluid_no, &n_neighbour, &no_neigh, &cnt, &up_cnt};
	int isize[]          = {(int)len.prod(), no_of_colloid };
	int psize[]          = {no_of_fluid, no_of_colloid};
	cudaMallocManaged(&box_part,  sizeof(int *)*(len.prod() + 2));
	cudaMallocManaged(&cell_part, sizeof(int *)*(len.prod() + 2));
	cudaMallocManaged(&box_neigh, sizeof(int *)*512);
	cudaMallocManaged(&neighbour, sizeof(int *)*256);
	cudaMallocManaged(&neigh_fl,  sizeof(int *)*(no_of_colloid + 2));
	cudaMallocManaged(&dp,  	  sizeof(int *)*(no_of_colloid + 2));
	cudaMallocManaged(&nbr, 	  sizeof(int *)*7005);
	cudaMallocManaged(&up_nbr,    sizeof(int *)*7005);
	cudaMallocManaged(&u,  		  sizeof(point *)*(no_of_colloid + 2));
	cudaMallocManaged(&rot, 	  sizeof(point *)*(len.prod() + 2));
	cudaMallocManaged(&cell_vel, (len.prod() + 2)*sizeof(point));
	cudaMallocManaged(&dump_vel_fl, (no_of_fluid + 2)*sizeof(point));
	cudaMallocManaged(&potential_colloid, sizeof(double));
	cudaMallocManaged(&vc, (no_of_colloid + 2)*sizeof(point *));
	cudaMallocManaged(&om, (no_of_colloid + 2)*sizeof(point *));
	cudaMallocManaged(&state, sizeof(curandState_t)*(no_of_fluid + 2));
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
		if(i <= no_of_colloid)	cudaMallocManaged(&vc[i],  sizeof(point)*(10000 + 2));
		if(i <= no_of_colloid)	cudaMallocManaged(&u[i],  sizeof(point)*(10000 + 2));
		if(i <= no_of_colloid)	cudaMallocManaged(&om[i],  sizeof(point)*(10000 + 2));
		if(i <= no_of_colloid)	cudaMallocManaged(&dp[i],  sizeof(int)*(512)); //512 = nbox
								cudaMallocManaged(&box_part[i],  sizeof(int)*(maxpart    + 2));
						   		cudaMallocManaged(&cell_part[i], sizeof(int)*(maxpart    + 2));
						   		cudaMallocManaged(&rot[i],		sizeof(point)*4);
	}
}
__global__ void conserv_mom(point *vel, point avr, int no) {
	int i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	if(i <= no) vel[i] = vel[i] - avr;
}

__global__ void initialize_colloid(point *pos_colloid, point *vel_colloid, point *ang_vel_colloid, int no_of_colloid, double sig_colloid, double kbt, double I_colloid, double kbt1, double mass_colloid, point len, curandState_t *state) {
	int counter = 0, check;
	double space_limit = 1.3*sig_colloid, ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid), vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
	point avr_vel = point(0, 0, 0), t, temp;
	while(counter < no_of_colloid) {
		t = t.rand(&state[1])*len;
		check = 1;
		for(int j = 1; j <= counter; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < space_limit)? 0: check;
		}
		if(check) pos_colloid[++counter] = t;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = (vel_colloid[j].rand(&state[1]) - point(0.5, 0.5, 0.5))*vscale_colloid;
		avr_vel += vel_colloid[j];
	}
	avr_vel = avr_vel/no_of_colloid;
	for(int j = 1; j <= no_of_colloid; j++) {
		vel_colloid[j] = vel_colloid[j] - avr_vel;
		ang_vel_colloid[j] = (t.rand(&state[1]) - point(0.5, 0.5, 0.5))*ang_vscale_colloid;
	}
}
void initialize_colloid() {
	initialize_colloid<<<1, 1>>>(pos_colloid, vel_colloid, ang_vel_colloid, no_of_colloid, sig_colloid, kbt, I_colloid, kbt1, mass_colloid, len, state);
}
__global__ void d_initialize_fluid(point *pos_fl, point *vel_fl, point *pos_colloid, int no_of_colloid, int no_of_fluid, double kbt, double mass_fl, double sigma, point len, curandState_t *state) {
	int check = 1, i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	double vscale_fluid = sqrt(12.0*kbt/mass_fl);
	point t, temp;
	while(true) {
		t = t.rand(&state[i])*len;
		check = 1;
		for(int j = 1; j <= no_of_colloid; j++) {
			temp = img(t - pos_colloid[j], len);
			check = (sqrt((temp*temp).sum()) < sigma*0.5)? 0: check;
		}
		if(check) {
			pos_fl[i] = t;
			break;
		}
	}
	vel_fl[i] = (vel_fl[i].rand(&state[i]) - point(0.5, 0.5, 0.5))*vscale_fluid;
}

void initialize_fluid() {
	blk = dim3((no_of_fluid + thr.x - 1)/thr.x);
	point avr_vel;	
	d_initialize_fluid<<<blk, thr>>>(pos_fl, vel_fl, pos_colloid, no_of_colloid, no_of_fluid, kbt, mass_fl, sigma, len, state);
	cudaDeviceSynchronize();
	avr_vel = thrust::reduce(thrust::device, vel_fl + 1, vel_fl + no_of_fluid + 1, point(0, 0, 0), add_point())/no_of_fluid;
	conserv_mom<<<blk, thr>>>(vel_fl, avr_vel, no_of_fluid);
}

__global__ void curand_setup(curandState_t *state, int seed) {
	int i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	curand_init(seed, i, 0, &state[i]);
}
void initialize_rand() {
	blk = dim3((no_of_fluid + thr.x - 1)/thr.x);
	curand_setup<<<blk, thr>>>(state, seed);
}
