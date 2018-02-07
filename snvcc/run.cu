#include "parameters.cuh"

__global__ void d_tumble(point *ra, point *pos_colloid, point len, int no_of_colloid){
	for (int i = 1; i <= no_of_colloid; i++) {
		ra[i] = img(pos_colloid[i] - ra[i].random(point(0, 0, 0), len), len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
    }
}

void tumble() {
	point *d_ra, *d_pos_colloid;
	cudaMalloc(&d_ra, sizeof(point)*(no_of_colloid + 2));
	cudaMalloc(&d_pos_colloid, sizeof(point)*(no_of_colloid + 2));
	cudaMemcpy(d_pos_colloid, pos_colloid, sizeof(point)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ra, ra, sizeof(point)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	d_tumble<<<1, 1>>>(d_ra, d_pos_colloid, len, no_of_colloid);
	cudaMemcpy(pos_colloid, d_pos_colloid, sizeof(point)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(ra, d_ra, sizeof(point)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaFree(d_ra), cudaFree(d_pos_colloid);
}

__global__ void d_run(point *ra, point *vel_colloid, point *vel_fl, point *pos_fl, point *pos_colloid, point len, 
				int *no_neigh, int **nbr, int **neigh_fl, int *cnt, int *up_cnt, int no_of_colloid, double mass_fl, 
				double v0, double mass_colloid, double sigma) {
	point vector, del;
	double temp;
	for(int i = 1; i <= no_of_colloid; i++) {
		vel_colloid[i] += ra[i]*v0, del = ra[i]*v0;
		cnt[i] = up_cnt[i] = 0;
		for(int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[j][i]] - pos_colloid[i], len);
			if((vector*vector).sum() <= power(sigma*0.5+0.5, 2) && (vector*vel_colloid[i]).sum() <= 0)
				nbr[++cnt[i]][i] = neigh_fl[j][i];
 		}
		for(int j = 1; j <= cnt[i]; j++) {
			temp = mass_colloid/(mass_fl*cnt[i]);
			vel_fl[nbr[j][i]] = vel_fl[nbr[j][i]] - del*temp;
		}
	}
}

void run() {
	point *d_ra, *d_vel_colloid, *d_vel_fl, *d_pos_fl, *d_pos_colloid;
	int *d_no_neigh, **d_nbr, **d_neigh_fl, *d_cnt, *d_up_cnt, *h_nbr[7005], *h_neigh_fl[10005];
	cudaMalloc(&d_ra, sizeof(point)*sizeof(no_of_colloid + 2));
	cudaMalloc(&d_vel_colloid, sizeof(point)*(no_of_colloid + 2));
	cudaMalloc(&d_pos_colloid, sizeof(point)*(no_of_colloid + 2));
	cudaMalloc(&d_vel_fl, sizeof(point)*(no_of_fluid + 2));
	cudaMalloc(&d_pos_fl, sizeof(point)*(no_of_fluid + 2));
	cudaMalloc(&d_no_neigh, sizeof(int)*(no_of_colloid + 2));
	cudaMalloc(&d_cnt, sizeof(int)*(no_of_colloid + 2));
	cudaMalloc(&d_up_cnt, sizeof(int)*(no_of_colloid + 2));
	cudaMalloc(&d_nbr, sizeof(int *)*7005);
	cudaMalloc(&d_neigh_fl, sizeof(int *)*10005);

	cudaMemcpy(d_ra, ra, sizeof(point)*sizeof(no_of_colloid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vel_colloid, vel_colloid, sizeof(point)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pos_colloid, pos_colloid, sizeof(point)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vel_fl, vel_fl, sizeof(point)*(no_of_fluid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pos_fl, pos_fl, sizeof(point)*(no_of_fluid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_no_neigh, no_neigh, sizeof(int)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_cnt, cnt, sizeof(int)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_up_cnt, up_cnt, sizeof(int)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	for(int i = 0; i <= 10000; i++) {
		if(i <= 7000) {
			cudaMalloc(&h_nbr[i], sizeof(int)*(no_of_colloid + 2));
			cudaMemcpy(h_nbr[i], nbr[i], sizeof(int)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
		}
		cudaMalloc(&h_neigh_fl[i], sizeof(int)*(no_of_colloid + 2));
		cudaMemcpy(h_neigh_fl[i], neigh_fl[i], sizeof(int)*(no_of_colloid + 2), cudaMemcpyHostToDevice);
	}
	cudaMemcpy(d_nbr, h_nbr, sizeof(int *)*7005, cudaMemcpyHostToDevice);
	cudaMemcpy(d_neigh_fl, h_neigh_fl, sizeof(int *)*10005, cudaMemcpyHostToDevice);
	d_run<<<1, 1>>>(d_ra, d_vel_colloid, d_vel_fl, d_pos_fl, d_pos_colloid, len, d_no_neigh, d_nbr, d_neigh_fl, 
		d_cnt, d_up_cnt, no_of_colloid, mass_fl, v0, mass_colloid, sigma);

	cudaMemcpy(ra, d_ra, sizeof(point)*sizeof(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_colloid, d_vel_colloid, sizeof(point)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(pos_colloid, d_pos_colloid, sizeof(point)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_fl, d_vel_fl, sizeof(point)*(no_of_fluid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(pos_fl, d_pos_fl, sizeof(point)*(no_of_fluid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(no_neigh, d_no_neigh, sizeof(int)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(cnt, d_cnt, sizeof(int)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	cudaMemcpy(up_cnt, d_up_cnt, sizeof(int)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
	for(int i = 0; i <= 10000; i++) {
		if(i <= 7000) {
			cudaMemcpy(nbr[i], h_nbr[i], sizeof(int)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
			cudaFree(nbr[i]);
		}
		cudaMemcpy(neigh_fl[i], h_neigh_fl[i], sizeof(int)*(no_of_colloid + 2), cudaMemcpyDeviceToHost);
		cudaFree(neigh_fl[i]);
	}
	cudaFree(d_ra), cudaFree(d_pos_fl), cudaFree(d_vel_fl), cudaFree(d_cnt), cudaFree(d_up_cnt), cudaFree(d_nbr);
	cudaFree(d_pos_colloid), cudaFree(d_vel_colloid), cudaFree(d_no_neigh), cudaFree(d_neigh_fl);
}
void updown_velocity() {
	point up_vel = point(0, 0, 0), vector, vel;
	for (int i = 1; i <= no_of_colloid; i++){
		cnt[i] = 0, up_cnt[i] = 0, vel = point(0, 0, 0);
		for (int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[j][i]] - pos_colloid[i], len);
			
			if((vector*vector).sum() <= pow((sigma*0.5 + 0.5), 2) && (vector*vel_colloid[i]).sum() <= 0.0)
				nbr[++cnt[i]][i] = neigh_fl[j][i];

			if((vector*vector).sum() <= pow((sigma*0.5 + 0.1), 2) && (vector*vel_colloid[i]).sum() <= 0.0)
				up_nbr[++up_cnt[i]][i] = neigh_fl[j][i];
		}
		for (int j = 1; j <= cnt[i]; j++)
			vel += vel_fl[nbr[j][i]];

		for(int j = 1; j <= up_cnt[i]; j++)
			up_vel += vel_fl[up_nbr[j][i]];

		up_vel = (up_cnt[i] > 0)? up_vel/up_cnt[i] - vel_colloid[i]: up_vel;
		vel    = (up_cnt[i] > 0)? vel/cnt[i] - vel_colloid[i]: vel;
	}
}