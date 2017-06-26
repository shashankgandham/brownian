// This function looks like it can be seriously optimized by reducing the number of loops;

void create_box() {
	int alx, aly, alz, temp, box, k, j, i;
	for(int l = 1; l <= lz; l++) {
		for(int m = 1; m <= ly; m++) {
			for(int mn = 1; mn <= lx; mn++) {
				box = (l - 1)*lx*ly + (m-1)*lx + mn;
				nbox = 0;
				for(kk = l; kk <= l + 6; kk++) {
					k = kk - 3;
					if(kk <= 0) 
						k += lz; 
					if(k > lz)
						k = k - lz;		
					//Use Mod maybe?
					for(int pp = m; pp <= m + 6; pp++) {
						j = pp-3;
						if(j <= 0)
							j += ly;
						else if( j > ly)
							j -= ly;
						//Mod?
						for(int ii = mn; ii <= mn + 6; ii++) {
							i = ii - 3;
							if( i <= 0) 
								i += lx;
							else
								i -= lx;
							//Serious considerations for mod;
							temp = (k - 1)*lx*ly + (j-1)*lx + i
							if(temp /= box) {
								nbox++;
								box[nbox][box] = (k - 1)*lx*ly + (j - 1)*lx + 1;
							}
						}
					}
					
				}
			}
		}
	}
}
