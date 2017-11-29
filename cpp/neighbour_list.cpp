#include "parameters.hpp"

void neighbour_list_md() {
	double x, y, z, r, neigh_cutoff = 3.0*sig_colloid;
	memset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));

	for(int i = 1; i < no_of_colloid; i++) {
		for(int j = i + 1; j <= no_of_colloid; j++) {
			x = img(pos_colloid[3*i - 2] - pos_colloid[3*j - 2], lx);
			y = img(pos_colloid[3*i - 1] - pos_colloid[3*j - 1], ly);
			z = img(pos_colloid[3*i - 0] - pos_colloid[3*j - 0], lz);
			r = x*x + y*y + z*z;
			if(r < pow(neigh_cutoff,2)) {
				n_neighbour[i] = n_neighbour[i] + 1;
				neighbour[n_neighbour[i]][i] = j;
			}
		}
	}
}

void neighbour_list_mpcd() {
	int box_no, *fluid_no, **box_part, mm, cbox;
	double x, y, z;
	fluid_no = (int *)calloc(lx*ly*lz + 1, sizeof(int));
	box_part = (int **)calloc(maxpart + 1, sizeof(int *));

	for(int i = 0; i < maxpart; i++)
		box_part[i] = (int *)calloc(lx*ly*lz + 1, sizeof(int));

	for(int	i = 1; i <= no_of_fluid; i++) {
		box_no = 1 + int(pos_fl[3*i-2]) + lx*int(pos_fl[3*i-1]) + lx*ly*int(pos_fl[3*i]);
		box_part[++fluid_no[box_no]][box_no] = i;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		no_neigh[j] = 0;
		x = int(pos_colloid[3*j-2]);
		y = int(pos_colloid[3*j-1]);
		z = int(pos_colloid[3*j]);
		cbox = 1 + x + y*lx + z*lx*ly;
		for(int k = 1; k <= nbox; k++) {
			mm = box_neigh[k][cbox];
			for(int i = 1; i <= fluid_no[mm]; i++) {
				no_neigh[j] = no_neigh[j] + 1;
				neigh_fl[no_neigh[j]][j] = box_part[i][mm];
			}
		}
	}

	for(int i = 0; i < maxpart; i++)
		free(box_part[i]);
	free(box_part);
	free(fluid_no);
}
