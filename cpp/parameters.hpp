#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

double ran();
struct point{
	double x, y, z;
	point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z){}
	point operator+ (const point &op) { return point(x + op.x, y + op.y, z + op.z); }
	point operator- (const point &op) { return point(x - op.x, y - op.y, z - op.z); }
	point operator/ (const double &op){ return point(x / op,   y / op,   z / op);   }
	point operator* (const point &op) { return point(x * op.x, y * op.y, z * op.z); }
	point operator* (const double &op){ return point(x * op, y * op, z * op); }
	point operator+=(const point &op) { x += op.x, y += op.y, z += op.z; return *this; }
	point operator-=(const point &op) { x -= op.x, y -= op.y, z -= op.z; return *this; }

	double sum()  { return (x + y + z); }
	double prod() { return (x * y * z); }
	void print() { printf("%.16lf %0.16lf %0.16lf\n", x, y, z); }

	point random(double dec = 0) {
		x = ran() - dec; y = ran() - dec; z = ran() - dec;
		return point (x, y, z);
	}
};

extern int n, niter, file, nbin, no_of_fluid, maxpart, no_of_colloid, ntab, nbox;
extern int *neighbour[200], *n_neighbour, *no_neigh, *neigh_fl[10000], *box_neigh[512], *iv;
extern double kbt, kbt1, ndt, dt, mass_colloid, sig_colloid, eps, v0, sigma, dv, mass_fl, I_colloid, potential_colloid;
extern point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *old_force, *ra, len;

void create_box(), compute_force_md(), fluid_colloid_collision(), initialize(), initialize_fluid(), initialize_colloid();
void neighbour_list_md(), neighbour_list_mpcd(), rotation_mpcd(), run(), tumble(), updown_velocity();
void update_velocity_colloid(), update_pos_md(), update_pos_mpcd() ,update_activity_direction();

point img(point, point), mod(point, point);
