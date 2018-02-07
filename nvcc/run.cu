#include "parameters.cuh"

__global__ void d_tumble(point *ra, point *pos_colloid, point len, int no_of_colloid){
	for (int i = 1; i <= no_of_colloid; i++) {
		ra[i] = img(pos_colloid[i] - ra[i].random(point(0, 0, 0), len), len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
    }
}

void tumble() {
	d_tumble<<<1, 1>>>(ra, pos_colloid, len, no_of_colloid);
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
	d_run<<<1, 1>>>(ra, vel_colloid, vel_fl, pos_fl, pos_colloid, len, no_neigh, nbr, neigh_fl, 
		cnt, up_cnt, no_of_colloid, mass_fl, v0, mass_colloid, sigma);
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