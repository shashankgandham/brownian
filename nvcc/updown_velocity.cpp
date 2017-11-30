#include "parameters.hpp"

void updown_velocity(){
	int cnt[no_of_colloid], **nbr, up_cnt[no_of_colloid], **up_nbr;
	point up_vel = point(0, 0, 0), vector, vel;

	nbr = (int **)malloc(sizeof(int *)*7005);
	up_nbr = (int **)malloc(sizeof(int *)*7005);
	for(int i = 0; i <= 7000; i++) {
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));
		up_nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));
	}

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
	for(int i = 1; i <= 7000; i++)
		free(nbr[i]), free(up_nbr[i]);
	free(nbr), free(up_nbr);
}
