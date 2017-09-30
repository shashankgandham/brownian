#include "parameters.hpp"

int create_box(int *box_neigh[512]) {
	int temp, box, k, j, i, nbox;

	for(int l = 1; l <= lz; l++) {
		for(int m = 1; m <= ly; m++) {
			for(int mn = 1; mn <= lx; mn++) {
				box = (l - 1)*lx*ly + (m-1)*lx + mn;
				nbox = 0;
				for(int kk = l; kk <= l + 6; kk++) {
					k = mod(kk - 3, lz);
					for(int pp = m; pp <= m + 6; pp++) {
						j = mod(pp - 3, ly);
						for(int ii = mn; ii <= mn + 6; ii++) {
							i = mod(ii - 3, lx);
							temp = (k - 1)*lx*ly + (j-1)*lx + i;
							if(temp /= box) {
								nbox++;
								box_neigh[nbox][box] = (k - 1)*lx*ly + (j - 1)*lx + 1;
							}
						}
					}
				}
			}
		}
	}
	return nbox;
}
