#include "parameters.cuh"

__global__ void d_tumble(point *ra, point *pos_colloid, point len, int no_of_colloid, curandState_t *state){
	int i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	if(i <= no_of_colloid)  {
		ra[i] = img(pos_colloid[i] - ra[i].rand(&state[i])*len, len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
	}
}

void tumble() {
	int thr = 256, blk = (no_of_colloid + thr - 1)/thr;
	d_tumble<<<blk, thr>>>(ra, pos_colloid, len, no_of_colloid, state);
}

__global__ void d_nbrc(point *ra, point *vel_colloid, point *pos_fl, point *pos_colloid, point len, 
		int *no_neigh, int **nbr, int **neigh_fl, int *cnt, int no_of_colloid, double v0, double sigma) {
	point vector;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_colloid) {
		vel_colloid[i] += ra[i]*v0;
		cnt[i] = 0;
		for(int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[i][j]] - pos_colloid[i], len);
			if((vector*vector).sum() <= power(sigma*0.5+0.5, 2) && (vector*vel_colloid[i]).sum() <= 0)
				nbr[++cnt[i]][i] = neigh_fl[i][j];
		}
	}
}
__global__ void d_velc(point *ra, point *vel_fl, int **nbr, int *cnt, int no_of_colloid, double mass_colloid, double mass_fl, double v0) {
	point del; double temp;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_colloid) {
		del = ra[i]*v0, temp = mass_colloid/(mass_fl*cnt[i]);
		for(int j = 1; j <= cnt[i]; j++) {
			vel_fl[nbr[j][i]] -= del*temp;
		}
	}
}
void run() {
	dim3 thr(32), blk((no_of_colloid + thr.x -1)/thr.x);
	d_nbrc<<<blk, thr>>>(ra, vel_colloid, pos_fl, pos_colloid, len, no_neigh, nbr, neigh_fl, cnt, no_of_colloid, v0, sigma);
	d_velc<<<blk, thr>>>(ra, vel_fl, nbr, cnt, no_of_colloid, mass_colloid, mass_fl, v0);
}
__global__ void d_updown_velocity(int no_of_colloid, int *cnt, int *up_cnt, int **neigh_fl, point *pos_fl, point *pos_colloid, point *vel_colloid) {
	point up_vel = point(0, 0, 0), vector, vel;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;	
	if(i <= no_of_colloid) {
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

void updown_velocity() {
	int thr = 256, blk = (no_of_colloid + thr - 1)/thr;
	d_updown_velocity<<<blk, thr>>>(no_of_colloid, cnt, up_cnt, neigh_fl, pos_fl, pos_colloid, vel_colloid)
}