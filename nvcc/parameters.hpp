#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

double ran();
struct point{
	double x, y, z;
	inline point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z){}
	inline point operator+ (const point &op) { return point(x + op.x, y + op.y, z + op.z); }
	inline point operator- (const point &op) { return point(x - op.x, y - op.y, z - op.z); }
	inline point operator/ (const double &op){ return point(x / op,   y / op,   z / op);   }
	inline point operator/ (const point &op) { return point(x / op.x, y / op.y, z / op.z); }
	inline point operator* (const point &op) { return point(x * op.x, y * op.y, z * op.z); }
	inline point operator* (const double &op){ return point(x * op, y * op, z * op); }
	inline point operator+=(const point &op) { x += op.x, y += op.y, z += op.z; return *this; }
	inline point operator-=(const point &op) { x -= op.x, y -= op.y, z -= op.z; return *this; }

	inline double sum()  { return (x + y + z); }
	inline double prod() { return (x * y * z); }
	inline void print()  { printf("%.16lf %0.16lf %0.16lf\n", x, y, z); }
	inline point random(double dec = 0) {
		x = ran() - dec; y = ran() - dec; z = ran() - dec;
		return *this;
	}
	inline void next(point len, int inc) {
		x += inc;
		if(x == len.x + 1) y += inc, x = 1;
		if(y == len.y + 1) z += inc, y = 1;
	}
	inline int cell(point len) {
		return int(x) + len.x*int(y) + len.x*len.y*int(z);
	}
};

inline point abs(point a)   { return ((point(abs(a.x), abs(a.y), abs(a.z)) )); }
inline point round(point a) { return ((point(round(a.x), round(a.y), round(a.z)))); }
inline point mod(point a, point b) { return ((a - b*(round((a - b/2)/b)))); }
inline point img(point a, point b) { return ((a - b*round(a/b)));    }


extern int n, niter, file, nbin, no_of_fluid, maxpart, no_of_colloid, ntab, nbox, **nbr, **up_nbr;
extern int *neighbour[256], *n_neighbour, *no_neigh, *neigh_fl[10005], *box_neigh[512], *iv, **box_part, *fluid_no, **cell_part;
extern double kbt, kbt1, ndt, dt, mass_colloid, sig_colloid, eps, v0, sigma, dv, mass_fl, I_colloid, potential_colloid;
extern point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *old_force, *ra, len;
extern double mag_f, r_cutoff, r, fc, ufc;

void create_box(), compute_force_md(), fluid_colloid_collision(), initialize(), initialize_fluid(), initialize_colloid();
void neighbour_list_md(), neighbour_list_mpcd(), rotation_mpcd(), run(), tumble(), updown_velocity();
void update_velocity_colloid(), update_pos_md(), update_pos_mpcd() ,update_activity_direction();
