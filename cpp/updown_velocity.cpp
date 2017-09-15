//done
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "parameters.hpp"
void updown_velocity(){
	double up_velx2,up_vely2,up_velz2;
	double vector_x,vector_y,vector_z,vector,dot;
	double velx2,vely2,velz2,velx12,vely12,velz12;
	int i,ii,j,jj,mm,cnt[no_of_colloid],**nbr;
	int up_cnt[no_of_colloid],**up_nbr;

	nbr = (int **)malloc(sizeof(int *)*7001);
	for(int i = 0; i <= 7000; i++) {
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));
	}
	up_nbr = (int **)malloc(sizeof(int *)*7001);
	for(int i = 0; i <= 7000; i++) {
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid + 1));
	}
	for (i = 1; i <= no_of_colloid; i++){
		cnt[i]=0;
		up_cnt[i]=0;
		for (j = 1; j <= no_neigh[i]; j++){
			jj=neigh_fl[j][i];

			vector_x=pos_fl[3*jj-2]-pos_colloid[3*i-2];
			vector_y=pos_fl[3*jj-1]-pos_colloid[3*i-1];
			vector_z=pos_fl[3*jj]-pos_colloid[3*i];

			vector_x =vector_x - llx*round(vector_x*inv_llx);
			vector_y =vector_y - lly*round(vector_y*inv_lly);
			vector_z =vector_z - llz*round(vector_z*inv_llz);

			vector=vector_x*vector_x+vector_y*vector_y+vector_z*vector_z;

			dot=vector_x*vel_colloid[3*i-2]+vector_y*vel_colloid[3*i-1]+vector_z*vel_colloid[3*i];
			if(vector <= (sigma*0.5+0.5)*(sigma*0.5+0.5) && dot <= 0.0){
				cnt[i]=cnt[i]+1;   
				nbr[cnt[i]][i]=jj;  
			}

			if(vector <= (sigma*0.5+0.1)*(sigma*0.5+0.1) && dot >= 0.0){
				up_cnt[i]=up_cnt[i]+1;
				up_nbr[up_cnt[i]][i]=jj; 
			}
		}
		velx2=0.0;
		vely2=0.0;
		velz2=0.0;

		if (cnt[i]>0 && up_cnt[i]>0){
			for (ii = 1; ii <= cnt[i]; i++){
				mm=nbr[ii][i];
				velx2=velx2+vel_fl[3*mm-2];
				vely2=vely2+vel_fl[3*mm-1];
				velz2=velz2+vel_fl[3*mm];
			}
			up_velx2=0.0;
			up_vely2=0.0;
			up_velz2=0.0;

			for(ii=1;ii <= up_cnt[i]; ii++){
				mm=up_nbr[ii][i];
				up_velx2=up_velx2+vel_fl[3*mm-2];
				up_vely2=up_vely2+vel_fl[3*mm-1];
				up_velz2=up_velz2+vel_fl[3*mm];
			}

			up_velx2=up_velx2/double(up_cnt[i]); 
			up_vely2=up_vely2/double(up_cnt[i]);
			up_velz2=up_velz2/double(up_cnt[i]);

			velx2=velx2/double(cnt[i]);
			vely2=vely2/double(cnt[i]);
			velz2=velz2/double(cnt[i]);

			//write(172,fmt='(9g25.15)') qq,velx2,up_velx2,vely2,up_vely2,velz2,up_velz2; 

			up_velx2=up_velx2-vel_colloid[3*i-2]; 
			up_vely2=up_vely2-vel_colloid[3*i-1]; 
			up_velz2=up_velz2-vel_colloid[3*i];

			velx2=velx2-vel_colloid[3*i-2];     
			vely2=vely2-vel_colloid[3*i-1];
			velz2=velz2-vel_colloid[3*i];

			//write(173,fmt='(9g25.15)') qq,velx2,up_velx2,vely2,up_vely2,velz2,up_velz2
		} 
	}
	for(int i = 1; i <= 7000; i++) {
		free(nbr[i]);
	}
	for(int i = 1; i <= 7000; i++) {
		free(up_nbr[i]); 
	}
	free(nbr);
	printf("up 1\n");
	free(up_nbr);
	printf("up 2\n");
}

