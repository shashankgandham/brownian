//Simulation Specifications
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#define ran() (((double)rand()/RAND_MAX))
#define mod(a, b) ((fmod(fmod(a, b) + b, b)))

extern int n, niter, file, lx, ly, lz, nbin, no_of_fluid, maxpart, no_of_colloid;
extern int *neighbour[200], *n_neighbour, *no_neigh, *neigh_fl[10000], *box_neigh[512];
extern double kbt, kbt1, ndt, dt, mass_colloid, sig_colloid, eps, v0, sigma, dv, mass_fl;
extern double *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *ra, **dist, *old_force, *f;


int create_box();
void neighbour_list_md();
void neighbour_list_mpcd(int);
void rotation_mpcd();
void initialize_fluid();
void initialize_colloid(double);
void fluid_colloid_collision(double);
double compute_force_md();
void tumble();
void run();
void updown_velocity();
void update_velocity_colloid();
void update_pos_md();
void update_pos_mpcd();
void update_activity_direction();
void initialize();
