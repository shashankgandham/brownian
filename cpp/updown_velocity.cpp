#include "parameters.hpp"

void updown_velocity(){
	double up_velx, up_vely, up_velz, vector_x, vector_y, vector_z, vector, dot, velx , vely, velz;
	int ii, jj, mm, cnt[no_of_colloid], **nbr, up_cnt[no_of_colloid], **up_nbr;

	nbr = (int **)malloc(sizeof(int *)*7001);
	for(int i = 0; i <= 7000; i++)
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));

	up_nbr = (int **)malloc(sizeof(int *)*7001);
	for(int i = 0; i <= 7000; i++)
		up_nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));

	for (int i = 1; i <= no_of_colloid; i++){
		cnt[i] = 0, up_cnt[i]=0;
		for (int j = 1; j <= no_neigh[i]; j++) {
			jj = neigh_fl[j][i];

			vector_x = mod(pos_fl[3*jj-2] - pos_colloid[3*i-2], lx);
			vector_y = mod(pos_fl[3*jj-1]-pos_colloid[3*i-1], ly);
			vector_z = mod(pos_fl[3*jj]-pos_colloid[3*i], lz);
			vector = vector_x*vector_x + vector_y*vector_y + vector_z*vector_z;

			dot = vector_x*vel_colloid[3*i-2] + vector_y*vel_colloid[3*i-1] + vector_z*vel_colloid[3*i];

			if(vector <= (sigma*0.5+0.5)*(sigma*0.5+0.5) && dot <= 0.0){
				cnt[i] = cnt[i] + 1;
				nbr[cnt[i]][i] = jj;
			}

			if(vector <= (sigma*0.5+0.1)*(sigma*0.5+0.1) && dot >= 0.0){
				up_cnt[i]=up_cnt[i]+1;
				up_nbr[up_cnt[i]][i]=jj;
			}
		}
		velx = vely = velz = 0;

		if (cnt[i] > 0 && up_cnt[i] > 0){
			for (int j = 1; j <= cnt[i]; j++){
				mm = nbr[j][i];
				velx += vel_fl[3*mm-2];
				vely += vel_fl[3*mm-1];
				velz += vel_fl[3*mm];
			}
			up_velx = up_vely = up_velz = 0;

			for(int j = 1; j <= up_cnt[i]; j++){
				mm = up_nbr[j][i];
				up_velx += vel_fl[3*mm - 2];
				up_vely += vel_fl[3*mm - 1];
				up_velz += vel_fl[3*mm];
			}

			up_velx /= up_cnt[i], up_vely /= up_cnt[i], up_velz /= up_cnt[i];
			velx /= cnt[i], vely /= cnt[i], velz /= cnt[i];

			up_velx -= vel_colloid[3*i - 2];
			up_vely -= vel_colloid[3*i - 1];
			up_velz -= vel_colloid[3*i];

			velx -= vel_colloid[3*i - 2];
			vely -= vel_colloid[3*i - 1];
			velz -= vel_colloid[3*i];

		}
	}
	for(int i = 1; i <= 7000; i++)
		free(nbr[i]);
	for(int i = 1; i <= 7000; i++)
		free(up_nbr[i]);
	free(nbr);
	free(up_nbr);
}
