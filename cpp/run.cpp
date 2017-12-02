#include "parameters.hpp"

void tumble(){
	double r;
	for (int i = 1; i <= no_of_colloid; i++) {
		ra[i] = img(pos_colloid[i] - ra[i].random()*len, len);
		ra[i] = ra[i]/sqrt((ra[i]*ra[i]).sum());
	}
}

void run() {
	int cnt[no_of_colloid], **nbr, up_cnt[no_of_colloid];
	point vector, del;

	nbr    = (int **)malloc(sizeof(int *) * 7005);
	for(int i = 0; i <= 7000; i++)
		nbr[i]    = (int *)malloc(sizeof(int)*(no_of_colloid + 2));

	for(int i = 1; i <= no_of_colloid; i++) {
		vel_colloid[i] += ra[i]*v0, del = ra[i]*v0;
		cnt[i] = up_cnt[i] = 0;
		for(int j = 1; j <= no_neigh[i]; j++) {
			vector = img(pos_fl[neigh_fl[j][i]] - pos_colloid[i], len);
			if((vector*vector).sum() <= pow(sigma*0.5+0.5, 2) && (vector*vel_colloid[i]).sum() <= 0)
				nbr[++cnt[i]][i] = neigh_fl[j][i];
		}
		for(int j = 1; j <= cnt[i]; j++) {
			vel_fl[nbr[j][i]] = vel_fl[nbr[j][i]] - del*mass_colloid/(mass_fl*cnt[i]);
		}
	}

	for(int i = 0; i <= 7000; i++)
		free(nbr[i]);
	free(nbr);
}
