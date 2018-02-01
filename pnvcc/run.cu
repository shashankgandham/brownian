#include "parameters.cuh"

__global__ void d_tumble(int no_of_colloid, point *d_ra, point *d_pos, point len){
	for (int i = 1; i <= no_of_colloid; i++) {
		//d_img(&d_pos[i],d_pos[i] - d_ra[i].random(point(0, 0, 0), len), len);
		//d_ra[i] = d_ra[i]/sqrt((d_ra[i]*d_ra[i]).sum());
		point c;
		d_random(&c, point(0, 0, 0), len, seed ,iv);
    }
}

void tumble(){
	/*for (int i = 1; i <= no_of_colloid; i++) {
		ra[i] = img(pos_colloid[i] - ra[i].random(point(0, 0, 0), len), len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
		ra[i].print();
    }*/
    
    point *d_ra, *d_pos;
    cudaMalloc(&d_ra, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&d_pos, (no_of_colloid + 2)*sizeof(point));
    cudaMemcpy(d_ra, ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_pos, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	
  	d_tumble<<<1,1>>>(no_of_colloid, d_ra, d_pos, len);

  	cudaMemcpy(ra, d_ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(d_ra); cudaFree(d_pos); 
  	for (int i = 1; i <= no_of_colloid; i++)
  		ra[i].print();
}

void run() {
	point vector, del;
	double temp;
	for(int i = 1; i <= no_of_colloid; i++) {
		vel_colloid[i] += ra[i]*v0, del = ra[i]*v0;
		cnt[i] = up_cnt[i] = 0;
		for(int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[j][i]] - pos_colloid[i], len);
			if((vector*vector).sum() <= pow(sigma*0.5+0.5, 2) && (vector*vel_colloid[i]).sum() <= 0)
				nbr[++cnt[i]][i] = neigh_fl[j][i];
		}
		for(int j = 1; j <= cnt[i]; j++) {
			temp = mass_colloid/(mass_fl*cnt[i]);
			vel_fl[nbr[j][i]] = vel_fl[nbr[j][i]] - del*(mass_colloid/(mass_fl*cnt[i]));
		}
	}
}

void updown_velocity(){
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
