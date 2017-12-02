#include "parameters.hpp"

void create_box() {
	int tbox, box;
	point temp;
	for(int l = 1; l <= len.z; l++) {
		for(int m = 1; m <= len.y; m++) {
			for(int mn = 1; mn <= len.x; mn++) {
				nbox = 0;
				for(int kk = l; kk <= l + 6; kk++) {
					for(int pp = m; pp <= m + 6; pp++) {
						for(int ii = mn; ii <= mn + 6; ii++) {
							temp = mod(point(ii - 3, pp - 3, kk - 3), len);
							box = (l - 1)*len.x*len.y + (m-1)*len.x + mn;
							tbox = (temp.z - 1)*len.x*len.y + (temp.y - 1)*len.x + temp.x;
							if(tbox != box) box_neigh[++nbox][box] = tbox;
						}
					}
				}
			}
		}
	}
}
