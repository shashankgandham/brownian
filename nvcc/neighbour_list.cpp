#include "parameters.hpp"

void neighbour_list_md() {
	double neigh_cutoff = 3.0*sig_colloid;
	point temp;
	memset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));
	for(int i = 1; i < no_of_colloid; i++) {
		for(int j = i + 1; j <= no_of_colloid; j++) {
			temp = img(pos_colloid[i] - pos_colloid[j], len);
			if((temp*temp).sum() < pow(neigh_cutoff,2)) {
				n_neighbour[i]++;
				neighbour[n_neighbour[i]][i] = j;
			}
		}
	}
}

void neighbour_list_mpcd() {
	int box_no, *fluid_no, **box_part, mm, cbox;
	point temp;
	fluid_no = (int *)calloc(len.prod() + 1, sizeof(int));
	box_part = (int **)calloc(maxpart + 1, sizeof(int *));

	for(int i = 0; i < maxpart; i++)
		box_part[i] = (int *)calloc(len.prod() + 1, sizeof(int));

	for(int	i = 1; i <= no_of_fluid; i++) {
		box_no = 1 + int(pos_fl[i].x) + len.x*int(pos_fl[i].y) + len.x*len.y*int(pos_fl[i].z);
		box_part[++fluid_no[box_no]][box_no] = i;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		no_neigh[j] = 0;
		temp = point(int(pos_colloid[j].x), int(pos_colloid[j].y), int(pos_colloid[j].z));
		cbox = 1 + temp.x + temp.y*len.x + temp.z*len.x*len.y;
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
	free(box_part), free(fluid_no);
}
