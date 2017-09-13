//Simulation Specifications
#include <cmath>
#include <cstdlib>
#include "parameters.hpp"

int n = 10, niter = 21000;
double kbt = 1, kbt1 = 1, ndt = 0.1, dt = ndt/(double)n; 
// Check Initialization of dt works??????
double pi = 3.14159265;
/* 
   kbt is system temperature after equilibrium and kbt1 is to initiate
   process at different temperature
 */
int thermostat = 1, step_thermo = 2;
int lx = 30, ly = 30, lz = 30, dt2 = dt*dt;
double llx = lx, lly = ly, llz = lz;
double llxby2 = llx/2, llyby2 = lly/2, llzby2 = llz/2;
double inv_llx = 1/llx, inv_lly = 1/lly, inv_llz = 1/llz;
double ang_ke_colloid, ang_ke_colloid1;
int zzz = -4280145, zzzz = 77777, nn, qq;

int nbin = 300; //For Maxwell velocity
double dv = 0.1; //Distribution
int *mb_vel = (int *)malloc(nbin + 2);
//Fluid Specifications

double mass_fl = 1.0, vscale_fluid = sqrt(12.0*kbt/mass_fl);
int no_of_fluid = lx*ly*lz*10, maxpart = 100; 

double *pos_fl = (double *)malloc(3 * no_of_fluid + 2), ke_fluid, scale_fac_mpcd;
double *vel_fl = (double *)malloc(3 * no_of_fluid + 2);
double *f_x = (double *)malloc(no_of_colloid);
double *f_y = (double *)malloc(no_of_colloid);
double *f_z = (double *)malloc(no_of_colloid);
int *n_neighbour = (int *)malloc(no_of_colloid);
double *f = (double *)malloc(3 * no_of_colloid);
double *old_force = (double *)malloc(3 * no_of_colloid);
double *pos_colloid = (double *)malloc(3 * no_of_colloid);
double *vel_colloid = (double *)malloc(3 * no_of_colloid);
double *ang_vel_colloid = (double *)malloc(3 * no_of_colloid);
int *no_neigh = (int *)malloc(no_of_colloid);
double *ra = (double *)malloc(3 * no_of_colloid);
double *old_ra = (double *)malloc(3 * no_of_colloid);
double *rra = (double *)malloc(3 * no_of_colloid);
int *size_cluster = (int *)malloc(no_of_colloid);
int *identify = (int *)malloc(no_of_colloid);

//Colloid Specifications
int no_of_colloid=1;
double mass_colloid = 654.10, vscale_colloid = sqrt(12.0*kbt1/mass_colloid);
double sig_colloid = 5.0, eps = 1.0;
double r_cutoff=pow(2, 1.0/6.0)*sig_colloid, space_limit = 1.3*sig_colloid;
double fc=4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
double ufc= 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6))+fc*r_cutoff;
double scale_fac, ang_scale_fac, potential_colloid, ke_colloid;
int *neighbour[200];
double r_cutoff2 = pow(r_cutoff, 2);
double sig_colloid12 = pow(sig_colloid, 12), sig_colloid6 = pow(sig_colloid, 6);
double sig_colloidby2_square = pow(sig_colloid, 2)/4.0;
double rfx, rfy, rfz, rcx, rcy, rcz, mom_x, mom_y, mom_z;
double sigma=0.80*sig_colloid, sigma_square = pow(sigma,2); //actual diameter of colloid
double rsx, rsy, rsz, x_pos, y_pos, z_pos;
double dx, dy, dz, ke_colloid1;
/*
   sig_colloid is minimum distance between two colloids, r_cutoff is 
   range of L-J Force,space llimit is initial sepration between any two 
   colloids, space_limit must be less than r_cutoff eps is epsilon in L-J
 */
double I_colloid = 0.4*mass_colloid*sigma*sigma*0.25; //MOMENT OF INERTIA
//Check for precedency and compare with fortran value!!!
double ang_vscale_colloid = sqrt(12.0*kbt1/I_colloid);
double neigh_cutoff=3.0*sig_colloid, neigh_cutoff2=neigh_cutoff*neigh_cutoff; //cutoff for costructing neighbour list
double length_cutoff = sigma*0.5 + 2; //2 more than 100*v*dt
double length_cutoff2 = pow(length_cutoff, 2);
double u1, u2, u3;
int *neigh_fl[10000]; //seg_fault maybe?
int *box_neigh[500], nbox;
double v0 = 0.04;
//Cluster Specifications
int n_cluster;
double **dist;
int file_colloid=0; //set 1 if initial data to be read from file

void initialize() {
	for(int i = 0; i < 10000; i++) 
		neigh_fl[i] = (int *)malloc(sizeof(int)*no_of_colloid);
	
	for(int i = 0; i < 500; i++)
		box_neigh[i] = (int *)malloc(sizeof(int)*lx*ly*lz);
	for(int i = 0; i < 200; i++) 
		neighbour[i] = (int *)malloc(sizeof(int)*no_of_colloid);
	dist = (double **)malloc(sizeof(double *)*no_of_colloid);
	for(int i = 0; i < no_of_colloid; i++) {
		dist[i] = (double *)malloc(sizeof(double)*no_of_colloid);
	}
}
