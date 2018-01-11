#include "parameters.hpp"

void create_box() {
	int tbox, box;
	point temp, iter = point(1, 0, 0);
	for(int i = 1; i <= len.prod(); i++, iter.next(len, 1), nbox = 0) {
		for(int z = 0; z <= 6; z++) {
			for(int y = 0; y <= 6; y++) {
				for(int x = 0; x <= 6; x++) {
					box = iter.cell(len);
					tbox = mod(iter - point(3 - x, 3 - y, 3 - z), len).cell(len);
					if(tbox != box) box_neigh[++nbox][box] = tbox;
				}
			}
		}
	}
}

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
	int box_no, mm, cbox;
	point temp;
	memset(fluid_no, 0, sizeof(int)*(len.prod() + 2));
	for(int	i = 1; i <= no_of_fluid; i++) {
		box_no = 1 + pos_fl[i].cell(len);
		box_part[++fluid_no[box_no]][box_no] = i;
	}
	for(int j = 1; j <= no_of_colloid; j++) {
		no_neigh[j] = 0;
		cbox = 1 + pos_colloid[j].cell(len);
		for(int k = 1; k <= nbox; k++) {
			mm = box_neigh[k][cbox];
			for(int i = 1; i <= fluid_no[mm]; i++) {
				no_neigh[j] = no_neigh[j] + 1;
				neigh_fl[no_neigh[j]][j] = box_part[i][mm];
			}
		}
	}
}
