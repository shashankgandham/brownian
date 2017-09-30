//Simulation Specifications
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#define ran() (((double)rand()/RAND_MAX))

extern int n, niter, *size_cluster, n_cluster, *identify, file;
extern double kbt, kbt1, ndt, dt, pi, **dist;
extern int nn, qq, lx, ly, lz;
extern double ang_ke_colloid, ang_ke_colloid1;
extern int nbin, no_of_fluid, maxpart, no_of_colloid;
extern double dv, mass_fl, vscale_fluid, *pos_fl, *vel_fl, ke_fluid, scale_fac_mpcd;
extern double mass_colloid, vscale_colloid, sig_colloid, eps, r_cutoff, space_limit, fc, ufc;
extern double scale_fac, ang_scale_fac, potential_colloid, ke_colloid;
extern int *neighbour[200], *n_neighbour, *no_neigh, *neigh_fl[10000], *box_neigh[500], nbox;
extern double *f, *old_force, *pos_colloid, *vel_colloid;
extern double *ang_vel_colloid, mom_x, mom_y, mom_z, sigma;
extern double x_pos, y_pos, z_pos, dx, dy, dz, ke_colloid1, *ra, *old_ra, *rra, v0;
extern double I_colloid, ang_vscale_colloid, neigh_cutoff, length_cutoff;


void create_box();
void neighbour_list_md();
void neighbour_list_mpcd();
void rotation_mpcd();
void initialize_fluid();
void initialize_colloid();
void fluid_colloid_collision();
void compute_force_md();
void tumble();
void run();
void initialize_activity();
void updown_velocity();
void update_velocity_colloid();
void update_pos_md();
void update_pos_mpcd();
void update_activity_direction();
void initialize();
double mod(double, double);
