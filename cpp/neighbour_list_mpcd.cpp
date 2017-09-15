#include "parameters.hpp"
#include <cstring>
#include <cstdio>
#include <cstdlib>
void neighbour_list_mpcd() {
	int i,j,k,ii,box_no,fluid_no[lx*ly*lz], **box_part;
	int p,l,lxly,mm, bx1,bx2,by1,by2,bz1,bz2;
	double x1,x2,y1,y2,z1,z2,sigmaby2;
	int di,idx,idy,idz,int_sigma, x1c,y1c,z1c,cbox;

	box_part = (int **)malloc(sizeof(int *)*maxpart);
	for(int i = 0; i < maxpart; i++) 
		box_part[i] = (int *)malloc(sizeof(int)*(lx*ly*lz));

	sigmaby2 = sigma*0.50;
	memset(fluid_no, 0, sizeof fluid_no);
	memset(box_part, 0, sizeof box_part); 
	lxly= lx*ly;
	i=1,no_of_fluid; //what does this mean??
	box_no = 1 + int(pos_fl[3*i-2])+lx*int(pos_fl[3*i-1])+lxly*int(pos_fl[3*i]);
	fluid_no[box_no] = fluid_no[box_no] + 1;
	j = fluid_no[box_no];
	box_part[j][box_no] = i;
	for(int j = 1; j <= no_of_colloid; j++) {
		no_neigh[j] = 0;
		x1=int(pos_colloid[3*j-2]);
		y1=int(pos_colloid[3*j-1]);
		z1=int(pos_colloid[3*j]);
		cbox=1+x1+y1*lx+z1*lx*ly;
		for(int k = 1; k <= nbox; k++) {
			if(box_neigh[k][cbox] > 125000) 
				//write(*,*)'error',cbox,j,nn   write to file again
				mm = box_neigh[k][cbox];
			for(int i=1; i <= fluid_no[mm]; i++) {
				ii = box_part[i][mm];
				no_neigh[j] = no_neigh[j]+1;
				neigh_fl[no_neigh[j]][j]=ii;
			} 
		} 
	}
	
	for(int i = 1; i < maxpart; i++)
		free(box_part[i]);
	free(box_part); 
}
