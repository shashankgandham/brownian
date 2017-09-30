#include "parameters.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
void run() {
	int jj, mm, cnt[no_of_colloid], **nbr, up_cnt[no_of_colloid], **up_nbr;
	double vector_x, vector_y, vector_z, vector, dot, xx, delx, dely, delz;

	nbr = (int **)malloc(sizeof(int *) * 7005);
	up_nbr = (int **)malloc(sizeof(int *) * 7005);
	for(int i = 0; i <= 7000; i++) {
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
		up_nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 2));
	}

	for(int i = 1; i <= no_of_colloid; i++){
		vel_colloid[3*i - 2] += v0*ra[3*i - 2];
		vel_colloid[3*i - 1] += v0*ra[3*i - 1];
		vel_colloid[3*i] += v0*ra[3*i];
		delx = v0*ra[3*i - 2], dely = v0*ra[3*i - 1], delz = v0*ra[3*i];

		cnt[i] = up_cnt[i] = 0;
		for(int j = 1; j <= no_neigh[i]; j++) {
			jj = neigh_fl[j][i];

			vector_x = mod(pos_fl[3*jj-2] - pos_colloid[3*i-2], lx);
			vector_y = mod(pos_fl[3*jj-1] - pos_colloid[3*i-1], ly);
			vector_z = mod(pos_fl[3*jj] - pos_colloid[3*i], lz);
			vector = vector_x*vector_x + vector_y*vector_y + vector_z*vector_z;

			dot = vector_x*vel_colloid[3*i-2] + vector_y*vel_colloid[3*i-1] + vector_z*vel_colloid[3*i];
			if(vector <= pow(sigma*0.5+0.5, 2) && dot <= 0)
				nbr[++cnt[i]][i] = jj;
		}
		for(int j = 1; j <= cnt[i]; j++) {
			jj = nbr[j][i];
			xx = vel_fl[3*jj-2];
			vel_fl[3*jj - 2] =vel_fl[3*jj-2] - delx*mass_colloid/(mass_fl*float(cnt[i]));
			vel_fl[3*jj - 1] =vel_fl[3*jj-1] - dely*mass_colloid/(mass_fl*float(cnt[i]));
			vel_fl[3*jj] = vel_fl[3*jj] - delz*mass_colloid/(mass_fl*float(cnt[i]));
		}
	}

	for(int i = 0; i <= 7000; i++)
		free(nbr[i]), free(up_nbr[i]);
	free(nbr), free(up_nbr);
}
