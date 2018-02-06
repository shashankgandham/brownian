#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <climits>
#define max(a, b) (a>b?a:b)
#define min(a, b) (a<b?a:b)
#define CUDA_CALLABLE_MEMBER __host__ __device__

inline CUDA_CALLABLE_MEMBER double ran() {
	static int im1 = 2147483563, im2 = 2147483399, ia1 = 40014, ia2 = 40692, iq1 = 53668, iq2 = 52774;
	static int imm, ir1 = 12211, ir2 = 3791, ndiv, idum = 123456789, seed = 77777, iy, j, k, ntab = 32, iv[32];
	double eps, rnmx, am; eps = 1.2e-7, rnmx = 1 - eps, am = 1.0/im1; ndiv = 1 + imm/ntab; imm = im1 - 1;
	if(seed <= 0) {
		seed = max(-seed,1);
		idum = seed;
		for(j = ntab + 8; j >= 1; j--) {
			k = seed/iq1;
			seed = ia1*(seed - k*iq1) - k*ir1;
			if(seed < 0) seed += im1;
			if(j <= ntab) iv[j] = seed;
		}
		iy = iv[1];
	}
	k = seed/iq1;
	seed = ia1*(seed - k*iq1) - k*ir1;
	if(seed < 0) seed += im1;
	k = idum/iq2;
	idum = ia2*(idum - k*iq2) - k*ir2;
	if(idum < 0) idum = idum + im2;
	j = 1 + iy/ndiv;
	iy = iv[j] - idum;
	iv[j] = seed;
	if(iy < 1) iy = iy + imm;
	return min(am*iy, rnmx);
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
	CUDA_CALLABLE_MEMBER point random(point dec = point(0, 0 ,0), point mul = point(1, 1, 1)) {
		x = ran(); y = ran(); z = ran();
		*this = (*this)*mul - dec;
		return *this;
	}
	CUDA_CALLABLE_MEMBER void next(point len, point inc = point(1, 1, 1), point start = point(1, 1, 1)) {
		x += inc.x;
		if(x > len.x) y += inc.y, x = start.x;
		if(y > len.y) z += inc.z, y = start.y;
	}
};
inline CUDA_CALLABLE_MEMBER void d_round(point *c, point a) { *c = point(round(a.x), round(a.y), round(a.z));}
inline CUDA_CALLABLE_MEMBER void d_mod(point *c, point a, point b) { d_round(c, (a - b/2)/b); *c = a - b*(*c);} 
inline CUDA_CALLABLE_MEMBER void d_img(point *c, point a, point b) {d_round(c, a/b); *c = a - b*(*c);}
inline CUDA_CALLABLE_MEMBER double power(double x, int r) { double ans = 1; for(int i = 1; i <=r; i++) ans *= x; return ans; }
extern point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *old_force, *ra, len;
extern int n, niter, file, nbin, no_of_fluid, maxpart, no_of_colloid, nbox, **nbr, **up_nbr, *cnt, *up_cnt, ntab, nn;
extern int *neighbour[256], *n_neighbour, *no_neigh, *neigh_fl[10005], *box_neigh[512], **box_part, *fluid_no, **cell_part;
extern double kbt, kbt1, ndt, dt, mass_colloid, sig_colloid, eps, v0, sigma, dv, mass_fl, I_colloid, potential_colloid;

void create_box(), compute_force_md(), fluid_colloid_collision(), initialize(), initialize_fluid(), initialize_colloid();
void neighbour_list_md(), neighbour_list_mpcd(), rotation_mpcd(), run(), tumble(), updown_velocity();
void update_velocity_colloid(), update_pos_md(), update_pos_mpcd() ,update_activity_direction();	