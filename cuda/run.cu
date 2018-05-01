#include "parameters.cuh"

__global__ void d_tumble(point *ra, point *pos_colloid, point len, int no_of_colloid, curandState_t *state){
	int i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	if(i <= no_of_colloid)  {
		ra[i] = img(pos_colloid[i] - ra[i].rand(&state[i])*len, len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
	}
}

void tumble() {
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	d_tumble<<<blk, thr>>>(ra, pos_colloid, len, no_of_colloid, state);
}

__global__ void d_nbrc(point *ra, point *vel_colloid, point *pos_fl, point *pos_colloid, point len, 
		int *no_neigh, int **nbr, int **neigh_fl, int *cnt, int no_of_colloid, double v0, double sigma) {
	point vector;
	int j = blockIdx.x*blockDim.x + threadIdx.x + 1;
	int i = blockIdx.y*blockDim.y + threadIdx.y + 1;
	if(i <= no_of_colloid) {
		vel_colloid[i] += ra[i]*v0, cnt[i] = 0;
		if(j <= no_neigh[i]) {
			vector = img(pos_fl[neigh_fl[i][j]] - pos_colloid[i], len);
			if((vector*vector).sum() <= power(sigma*0.5+0.5, 2) && (vector*vel_colloid[i]).sum() <= 0)
				nbr[atomicAdd(&cnt[i], 1) + 1][i] = neigh_fl[i][j];
		}
	}
}
__global__ void d_velc(point *ra, point *vel_fl, int **nbr, int *cnt, int no_of_colloid, double mass_colloid, double mass_fl, double v0) {
	point del; double temp;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_colloid) {
		del = ra[i]*v0, temp = mass_colloid/(mass_fl*cnt[i]);
		for(int j = 1; j <= cnt[i]; j++) 
			vel_fl[nbr[j][i]] -= del*temp;
	}
}

void run() {
	blk = dim3((10000 + thr.x -1)/thr.x, (no_of_colloid + thr.y - 1)/thr.y);
	d_nbrc<<<blk, thr>>>(ra, vel_colloid, pos_fl, pos_colloid, len, no_neigh, nbr, neigh_fl, cnt, no_of_colloid, v0, sigma);
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	d_velc<<<blk, thr>>>(ra, vel_fl, nbr, cnt, no_of_colloid, mass_colloid, mass_fl, v0);
}

__global__ void helper_upd(int no_of_colloid, int *cnt, int *up_cnt, int *no_neigh, point **vel, point **up_vel, 
			int **neigh_fl, point *pos_fl, point *pos_colloid, point *vel_colloid, point *vel_fl, point len, double sigma) {
	
	point vector;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;	
	int j = blockIdx.y*blockDim.y + threadIdx.y + 1;	
	if(i <= no_of_colloid) {
		if(j <= no_neigh[i]) {
			vector = img(pos_fl[neigh_fl[i][j]] - pos_colloid[i], len);

			if((vector*vector).sum() <= pow((sigma*0.5 + 0.5), 2) && (vector*vel_colloid[i]).sum() <= 0.0)
				vel[i][atomicAdd(&cnt[i], 1) + 1] = vel_fl[neigh_fl[i][j]];

			if((vector*vector).sum() <= pow((sigma*0.5 + 0.1), 2) && (vector*vel_colloid[i]).sum() <= 0.0)
				up_vel[i][atomicAdd(&up_cnt[i], 1) + 1] = vel_fl[neigh_fl[i][j]];
		}
	}
}
__global__ void calc_upd(int no_of_colloid, int *cnt, int *up_cnt, point **vel, point **up_vel, point *vel_colloid) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;	
	if(i <= no_of_colloid) {
		vel[i][0]    = thrust::reduce(thrust::device, vel[i], vel[i] + cnt[i] + 1, point(0, 0, 0), add_point());	
		vel[i][0]    = (cnt[i])? vel[i][0]/cnt[i] - vel_colloid[i]: vel[i][0];
		up_vel[i][0] = thrust::reduce(thrust::device, up_vel[i], up_vel[i] + up_cnt[i] + 1, point(0, 0, 0), add_point());	
		up_vel[i][0] = (up_cnt[i])? up_vel[i][0]/up_cnt[i] - vel_colloid[i]: up_vel[i][0];
	}
}
void updown_velocity() {
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	imemset<<<blk, thr>>>(cnt, no_of_colloid);
	imemset<<<blk, thr>>>(up_cnt, no_of_colloid);
	helper_upd<<<blk, thr>>>(no_of_colloid, cnt, up_cnt, no_neigh, vel, up_vel, neigh_fl, pos_fl, 
			pos_colloid, vel_colloid, vel_fl, len, sigma);
	calc_upd<<<blk, thr>>>(no_of_colloid, cnt, up_cnt, vel, up_vel, vel_colloid);
}
