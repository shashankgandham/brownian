#include "parameters.hpp"

inline point cmod(point a, point b) {  
    if(a.x <= 0) a.x += b.x;
    if(a.y <  0) a.y += b.y;
    if(a.z <  0) a.z += b.z;
    if(a.x >  b.x) a.x -= b.x;
    if(a.y >= b.y) a.y -= b.y;
    if(a.z >= b.z) a.z -= b.z;
    return a;
}
void create_box() {
	int tbox, box;
	point temp, iter = point(1, 0, 0);
	for(int i = 1; i <= len.prod(); i++, iter.next(len, point(1, 1, 1), point(1, 0, 0))) {
		nbox = 0, box = iter.cell(len);
		for(int z = -3; z <= 3; z++) {
			for(int y = -3; y <= 3; y++) {
				for(int x = -3; x <= 3; x++) {
					temp = cmod(iter + point(x, y, z), len);
					tbox = temp.cell(len);
					if(tbox != box) {
                        box_neigh[++nbox][box] = tbox;
                    }
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
				neighbour[++n_neighbour[i]][i] = j;
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
				neigh_fl[++no_neigh[j]][j] = box_part[i][mm];
			}
		}
	}
}
