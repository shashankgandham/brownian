#include "parameters.hpp"

void updown_velocity(){
	int jj, mm, cnt[no_of_colloid], **nbr, up_cnt[no_of_colloid], **up_nbr;
	coord up_vel = coord(0, 0, 0), vector, vel;

	nbr = (int **)malloc(sizeof(int *)*7005);
	up_nbr = (int **)malloc(sizeof(int *)*7005);
	for(int i = 0; i <= 7000; i++) {
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));
		up_nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));
	}

	for (int i = 1; i <= no_of_colloid; i++){
		cnt[i] = 0, up_cnt[i] = 0, vel = coord(0, 0, 0);
		for (int j = 1; j <= no_neigh[i]; j++) {
			jj = neigh_fl[j][i];
			vector = img(pos_fl[jj] - pos_colloid[i], len);

			if((vector*vector).sum() <= (sigma*0.5 + 0.5)*(sigma*0.5 + 0.5) && (vector*vel_colloid[i]).sum() <= 0.0)
				nbr[++cnt[i]][i] = jj;

			if((vector*vector).sum() <= (sigma*0.5 + 0.1)*(sigma*0.5 + 0.1) && (vector*vel_colloid[i]).sum() <= 0.0)
				up_nbr[++up_cnt[i]][i] = jj;
		}
		for (int j = 1; j <= cnt[i]; j++)
			vel += vel_fl[nbr[j][i]];

		for(int j = 1; j <= up_cnt[i]; j++)
			up_vel += vel_fl[nbr[j][i]];

		up_vel = (up_cnt[i] > 0)? up_vel/up_cnt[i] - vel_colloid[i]: up_vel;
		vel    = (up_cnt[i] > 0)? vel/cnt[i] - vel_colloid[i]: vel;
	}
	for(int i = 1; i <= 7000; i++)
		free(nbr[i]), free(up_nbr[i]);
	free(nbr), free(up_nbr);
}
