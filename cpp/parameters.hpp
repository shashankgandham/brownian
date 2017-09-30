//Simulation Specifications
#include <cmath>

extern int n, niter, *size_cluster, n_cluster, *identify, file_colloid;
extern double kbt, kbt1, ndt, dt, pi, **dist;
extern int thermostat, step_thermo, zzz, zzzz, nn, qq, lx, ly, lz, dt2;
extern double llx, lly, llz, llxby2, llyby2, llzby2, inv_llx, inv_lly, inv_llz, ang_ke_colloid, ang_ke_colloid1;
extern int nbin, *mb_vel, no_of_fluid, maxpart, no_of_colloid;
extern double dv, mass_fl, vscale_fluid, *pos_fl, *vel_fl, ke_fluid, scale_fac_mpcd;
extern double mass_colloid, vscale_colloid, sig_colloid, eps, r_cutoff, space_limit, fc, ufc;
extern double scale_fac, ang_scale_fac, *f_x, *f_y, *f_z, potential_colloid, ke_colloid;
extern int *neighbour[200], *n_neighbour, *no_neigh, *neigh_fl[10000], *box_neigh[500], nbox;
extern double *f, r_cutoff2, sig_colloid12, sig_colloid6, sig_colloidby2_square, *old_force, *pos_colloid, *vel_colloid;
extern double *ang_vel_colloid, rfx, rfy, rfz, rcx, rcy, rcz, mom_x, mom_y, mom_z, sigma, sigma_square; //actual diameter of colloid
extern double rsx, rsy, rsz, x_pos, y_pos, z_pos, dx, dy, dz, ke_colloid1, *ra, *old_ra, *rra, v0;
extern double I_colloid, ang_vscale_colloid, neigh_cutoff, neigh_cutoff2, length_cutoff, length_cutoff2, u1, u2, u3;


void create_box();
void neighbour_list_md();
void neighbour_list_mpcd();
void rotation_mpcd();
void initialize_fluid();
void initialize_colloid();
void fluid_colloid_collision();
void gauss(double, double);
void compute_force_md();
void tumble();
void run();
void initialize_activity();
void updown_velocity();
void update_velocity_colloid();
void update_pos_md();
void update_pos_mpcd();
void update_activity_direction();
void stochastic_reflection();
void initialize();
double mod(double, double);
