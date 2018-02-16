#include "parameters.cuh"
#include <cuda_profiler_api.h>

__global__ void d_tumble(point *ra, point *pos_colloid, point len, int no_of_colloid, int *iv, int *seed, int *idum, int *iy){
	for(int i = 1; i <= no_of_colloid; i++) {
		ra[i] = img(pos_colloid[i] - ra[i].random(iv, seed, idum, iy)*len, len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
	}
}

void tumble() {
	d_tumble<<<1, 1>>>(ra, pos_colloid, len, no_of_colloid, iv, seed, idum, iy);
}

__global__ void d_nbrc(point *ra, point *vel_colloid, point *pos_fl, point *pos_colloid, point len, 
		int *no_neigh, int **nbr, int **neigh_fl, int *cnt, int no_of_colloid, double v0, double sigma) {
	point vector;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
//	int j = blockIdx.y*blockDim.y + threadIdx.y + 1;
//	for(int i = 1; i <= no_of_colloid; i++) {
	if(i <= no_of_colloid) {
		vel_colloid[i] += ra[i]*v0;
		cnt[i] = 0;
		for(int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[i][j]] - pos_colloid[i], len);
			if((vector*vector).sum() <= power(sigma*0.5+0.5, 2) && (vector*vel_colloid[i]).sum() <= 0)
				nbr[atomicAdd(&cnt[i], 1) + 1][i] = neigh_fl[i][j];
		}
	}
}
void d_velc(point *ra, point *vel_fl, int **nbr, int *cnt, int no_of_colloid, double mass_colloid, double mass_fl, double v0) {
	point del; double temp;
	for(int i = 1; i <= no_of_colloid; i++) {
		del = ra[i]*v0;
		for(int j = 1; j <= cnt[i]; j++) {
			temp = mass_colloid/(mass_fl*cnt[i]);
			//atomic Sub Double supported in tesla P100;
			vel_fl[nbr[j][i]] = vel_fl[nbr[j][i]] - del*temp;
		}
	}
}

void run() {
	dim3 thr(32), blk((no_of_colloid + thr.x -1)/thr.x);
	cudaDeviceSynchronize();
	d_nbrc<<<blk, thr>>>(ra, vel_colloid, pos_fl, pos_colloid, len, no_neigh, nbr, neigh_fl, cnt, no_of_colloid, v0, sigma);
	cudaDeviceSynchronize();
	d_velc(ra, vel_fl, nbr, cnt, no_of_colloid, mass_colloid, mass_fl, v0);
}
void updown_velocity() {
	point up_vel = point(0, 0, 0), vector, vel;
	for (int i = 1; i <= no_of_colloid; i++){
		cnt[i] = 0, up_cnt[i] = 0, vel = point(0, 0, 0);
		for (int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[i][j]] - pos_colloid[i], len);

			if((vector*vector).sum() <= pow((sigma*0.5 + 0.5), 2) && (vector*vel_colloid[i]).sum() <= 0.0)
				nbr[++cnt[i]][i] = neigh_fl[i][j];

			if((vector*vector).sum() <= pow((sigma*0.5 + 0.1), 2) && (vector*vel_colloid[i]).sum() <= 0.0)
				up_nbr[++up_cnt[i]][i] = neigh_fl[i][j];
		}
		for (int j = 1; j <= cnt[i]; j++)
			vel += vel_fl[nbr[j][i]];

		for(int j = 1; j <= up_cnt[i]; j++)
			up_vel += vel_fl[up_nbr[j][i]];

		up_vel = (up_cnt[i] > 0)? up_vel/up_cnt[i] - vel_colloid[i]: up_vel;
		vel    = (up_cnt[i] > 0)? vel/cnt[i] - vel_colloid[i]: vel;
	}
}
