#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <climits>
#include <cuda_profiler_api.h>
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define _USE_MATH_DEFINES

inline CUDA_CALLABLE_MEMBER double ran(int *iv, int *seed, int *idum, int *iy, int debug = 0) {
	static int im1 = 2147483563, im2 = 2147483399, ia1 = 40014, ia2 = 40692, iq1 = 53668, iq2 = 52774, j, k;
	static int imm, ir1 = 12211, ir2 = 3791, ntab = 32, ndiv;
	double eps, rnmx, am; eps = 1.2e-7, rnmx = 1 - eps, am = 1.0/im1; imm = im1 - 1, ndiv = 1 + imm/ntab;
	if(*seed <= 0) {
		*seed = max(-*seed,1);
		*idum = *seed;
		for(j = ntab + 8; j >= 1; j--) {
			k = *seed/iq1;
			*seed = ia1*(*seed - k*iq1) - k*ir1;
			if(*seed < 0) *seed += im1;
			if(j <= ntab) iv[j] = *seed;
		}
		*iy = iv[1];
	}
	k = *seed/iq1;
	*seed = ia1*(*seed - k*iq1) - k*ir1;
	if(*seed < 0) *seed += im1;
	k = *idum/iq2;
	*idum = ia2*(*idum - k*iq2) - k*ir2;
	if(*idum < 0) *idum = *idum + im2;
	j = 1 + *iy/ndiv;
	*iy = iv[j] - *idum;
	iv[j] = *seed;
	if(*iy < 1) *iy = *iy + imm;
	if(debug) {
		printf("%d %d %d\n", *seed, *idum, *iy);
	}
	return min(am*(*iy), rnmx);

}

struct point {
	double x, y, z;
	CUDA_CALLABLE_MEMBER point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z){}
	CUDA_CALLABLE_MEMBER point operator+ (const point &op) { return point(x + op.x, y + op.y, z + op.z); }
	CUDA_CALLABLE_MEMBER point operator- (const point &op) { return point(x - op.x, y - op.y, z - op.z); }
	CUDA_CALLABLE_MEMBER point operator/ (const double &op){ return point(x / op,   y / op,   z / op);   }
	CUDA_CALLABLE_MEMBER point operator/ (const point &op) { return point(x / op.x, y / op.y, z / op.z); }
	CUDA_CALLABLE_MEMBER point operator* (const point &op) { return point(x * op.x, y * op.y, z * op.z); }
	CUDA_CALLABLE_MEMBER point operator* (const double &op){ return point(x * op, y * op, z * op); }
	CUDA_CALLABLE_MEMBER point operator+=(const point &op) { x += op.x, y += op.y, z += op.z; return *this; }
	CUDA_CALLABLE_MEMBER point operator-=(const point &op) { x -= op.x, y -= op.y, z -= op.z; return *this; }
	CUDA_CALLABLE_MEMBER bool operator==(const point &op) { return (x == op.x && y == op.y && z == op.z); }

	CUDA_CALLABLE_MEMBER double sum()  { return (x + y + z); }
	CUDA_CALLABLE_MEMBER double prod() { return (x * y * z); }
	CUDA_CALLABLE_MEMBER int cell(point len) { return int(x) + len.x*int(y) + len.x*len.y*int(z); }
	CUDA_CALLABLE_MEMBER void print()  { printf("%36.32lf %35.32lf %35.32lf\n", x, y, z); }
	CUDA_CALLABLE_MEMBER point random(int *iv, int *seed, int *idum, int *iy, int debug = 0){
		x = ran(iv, seed, idum, iy, debug); y = ran(iv, seed, idum, iy, debug); z = ran(iv, seed, idum, iy, debug);
		return *this;
	}
	CUDA_CALLABLE_MEMBER void next(point len, point inc = point(1, 1, 1), point start = point(1, 1, 1)) {
		x += inc.x;
		if(x > len.x) y += inc.y, x = start.x;
		if(y > len.y) z += inc.z, y = start.y;
	}
};
inline CUDA_CALLABLE_MEMBER point round(point a) { return point(round(a.x), round(a.y), round(a.z));}
inline CUDA_CALLABLE_MEMBER point mod(point a, point b) { return a - b*round((a - b/2)/b);} 
inline CUDA_CALLABLE_MEMBER point img(point a, point b) { return a - b*round(a/b);}
inline CUDA_CALLABLE_MEMBER double power(double x, int r) { double ans = 1; for(int i = 1; i <=r; i++) ans *= x; return ans; }

extern point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *old_force, *ra, len, *cell_vel, *del_v, **rot, *dump_vel_fl, **u, **vc, **om;
extern int n, niter, file, nbin, no_of_fluid, maxpart, no_of_colloid, nbox, **nbr, **up_nbr, *cnt, *up_cnt, nn, *iv, *seed, *iy, **dp;
extern int **neighbour, *n_neighbour, *no_neigh, **neigh_fl, **box_neigh, **box_part, *fluid_no, **cell_part, ran_c, *idum;
extern double kbt, kbt1, ndt, dt, mass_colloid, sig_colloid, eps, v0, sigma, dv, mass_fl, I_colloid, *potential_colloid, *rana, *ranb;

void create_box(), compute_force_md(), fluid_colloid_collision(), initialize(), initialize_fluid(), initialize_colloid();
void neighbour_list_md(), neighbour_list_mpcd(), rotation_mpcd(), run(), tumble(), updown_velocity();
void update_velocity_colloid(), update_pos_md(), update_pos_mpcd() ,update_activity_direction();	
