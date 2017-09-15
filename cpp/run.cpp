#include "parameters.hpp"
#include <cmath>
#include <cstdlib>
void run() {    
	int i,ii,j,jj,mm,cnt[no_of_colloid],**nbr, up_cnt[no_of_colloid],**up_nbr;
	double cos_b1,cos_b2,cos_b3,sin_b1,sin_b2,sin_b3, b1,b2,b3,bb1,bb2,bb3,ran1,v01,v02,v03;
	double m11,m12,m13,m21,m22,m23,m31,m32,m33, vector_x,vector_y,vector_z,vector,dot;
	double velx1,vely1,velz1,velx2,vely2,velz2,velx12,vely12,velz12, vv, theta,change,ke_before,ke_after,xx,vv1,vv2;
	double up_velx2,up_vely2,up_velz2,aaa,delx,dely,delz, dumx,dumy,dumz;
	delx=0.0;dely=0.0;delz=0.0;
	
	nbr = (int **)malloc(sizeof(int) * 7000);
	for(int i = 0; i < 7000; i++) 
		nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid));
	
	up_nbr = (int **)malloc(sizeof(int) * 7000);
	for(int i = 0; i < 7000; i++) 
		up_nbr[i] = (int *)malloc(sizeof(int)*(no_of_colloid));
	
	for(int i = 1; i <= no_of_colloid; i++){
		dumx=vel_colloid[3*i-2];
		dumy=vel_colloid[3*i-1];
		dumz=vel_colloid[3*i];

		vel_colloid[3*i-2]=vel_colloid[3*i-2]+v0*ra[3*i-2];
		vel_colloid[3*i-1]=vel_colloid[3*i-1]+v0*ra[3*i-1];
		vel_colloid[3*i]=vel_colloid[3*i]+v0*ra[3*i];

		delx=vel_colloid[3*i-2]-dumx;
		dely=vel_colloid[3*i-1]-dumy;
		delz=vel_colloid[3*i]-dumz;
		cnt[i]=0;up_cnt[i]=0;
		
		for(int j=1; j <= no_neigh[i]; j++) {
			jj=neigh_fl[j][i];
			
			vector_x=pos_fl[3*jj-2]-pos_colloid[3*i-2];
			vector_y=pos_fl[3*jj-1]-pos_colloid[3*i-1];
			vector_z=pos_fl[3*jj]-pos_colloid[3*i];

			vector_x =vector_x - llx*int(vector_x*inv_llx);
			vector_y =vector_y - lly*int(vector_y*inv_lly);
			vector_z =vector_z - llz*int(vector_z*inv_llz);

			vector=vector_x*vector_x+vector_y*vector_y+vector_z*vector_z;

			dot=vector_x*vel_colloid[3*i-2]+vector_y*vel_colloid[3*i-1]+vector_z*vel_colloid[3*i];

			if(vector <= pow(sigma*0.5+0.5, 2) && dot <= 0) {;
			; //write something & condition check
				cnt[i]=cnt[i]+1;
				nbr[cnt[i]][i]=jj; 
			}
		}
		for(int ii=1; ii <= cnt[i]; ii++) {
			mm=nbr[ii][i];
			xx=vel_fl[3*mm-2];
			vel_fl[3*mm-2]=vel_fl[3*mm-2]-delx*mass_colloid/(mass_fl*float(cnt[i]));
			vel_fl[3*mm-1]=vel_fl[3*mm-1]-dely*mass_colloid/(mass_fl*float(cnt[i]));
			vel_fl[3*mm]=vel_fl[3*mm]-delz*mass_colloid/(mass_fl*float(cnt[i]));
		}
	}

	for(int i = 0; i < 7000; i++) 
		free(nbr[i]);
	free(nbr);
	
	for(int i = 0; i < 7000; i++) 
		free(up_nbr[i]);
	free(up_nbr);
		
}
