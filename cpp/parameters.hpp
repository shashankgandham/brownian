//Simulation Specifications
extern const int n = 10, niter = 21000;
extern const double kbt = 1, kbt1 = 1, ndt 0.1, dt = ndt/(double)n; 
// Check Initialization of dt works??????
extern const double pi = 3.14159265;
/* 
	kbt is system temperature after equilibrium and kbt1 is to initiate
	process at different temperature
*/
extern const int thermostat = 1, step_thermo = 2;
extern const int lx = 30, ly = 30, lz = 30, dt2 = dt*dt;
extern const double llx = lx, lly = ly, llz = lz;
extern const double llxby2 = llx/2, llyby2 = lly/2, llzby2 = llz/2;
extern const double inv_llx = 1/llx, inv_lly = 1/lly, inv_llz = 1/llz;
extern double ang_ke_colloid, ang_ke_colloid1;
extern int zzz = -4280145, zzzz = 77777, nn, qq;

extern const int nbin = 300; //For Maxwell velocity
extern double dv = 0.1; //Distribution
extern int mb_vel[nbin + 2];

//Fluid Specifications

extern const double mass_fl = 1.0, vscale_fluid = sqrt(12.0*kbt/mass_fl);
extern const int no_of_fluid = lx*ly*lz*10, maxpart = 100; 
extern double pos_fl[3*no_of_fluid + 2] vel_fl[3*no_of_fluid + 2], ke_fluid, scale_fac_mpcd;

//Colloid Specifications
extern const int no_of_colloid=1;
extern const double mass_colloid = 654.10, vscale_colloid = sqrt(12.0*kbT1/mass_colloid);
extern const double sig_colloid = 5.0, eps = 1.0;
extern const double r_cutoff=pow(2, 1.0/6.0)*sig_colloid, space_limit = 1.3*sig_colloud;
extern const double fc=4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cut_off,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cut_off, 7));
extern const double ufc= 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_solloid/r_cutoff, 6) + fc*r_cutoff;
extern double scale fac, ang_scale_fac, f_x[no_of_colloid], f_y[no_of_colloid],f_z[no_of_colloid], potential_colloid, ke_colloid;
extern int neighbour[200][no_of_colloid], n_neighbour[no_of_colloid];
extern double f[3*no_of_colloid], r_cutoff2 = pow(r_cutoff, 2);
extern double sig_colloid12 = pow(sig_colloid, 12), sig_colloid6 = pow(sig_colloid, 6);
extern const double sig_colloidby2_square = pow(sig_colloid, 2)/4.0;
extern double old_force[3*no_of_colloid], pos_colloid[3*no_of_colloid], vel_colloid[3*no_of_colloid];
extern double ang_vel_colloid[3*no_of_colloid];
extern double rfx, rfy, rfz, rcx, rcy, rcz, mom_x, mom_y, mom_z;
extern const double sigma=0.80*sig_colloid, sigma_square = pow(sigma,2); //actual diameter of colloid
extern double rsx, rsy, ssz, x_pos, y_pos, z_pos;
extern double dx, dy, dz, ke_colloid1;
/*
sig_colloid is minimum distance between two colloids, r_cutoff is 
range of L-J Force,space llimit is initial sepration between any two 
colloids, space_limit must be less than r_cutoff eps is epsilon in L-J
*/
extern const double I_colloid = 0.4*mass_colloid*pow(sigma, 2)*0.25 //MOMENT OF INERTIA
//Check for precedency and compare with fortran value!!!
extern const double ang_vscale_colloid = sqrt(12.0*kbT1/I_colloid);
extern const double neigh_cutoff=3.0*sig_colloid, neigh_cutoff2=neigh_cutoff*neigh_cutoff; //cutoff for costructing neighbour list

extern const double length_cutoff = signma*0.5 + 2; //2 more than 100*v*dt
extern const double length_cutoff2 = pow(length_cutoff, 2);
extern double u1, u2, u3
extern int no_neigh[no_of_colloid], neigh_fl[10000][no_colloid]; //seg_fault maybe?
extern int box_neigh[500][lx*ly*lz], nbox;
extern double ra[3*no_colloid], old_ra[3*no_of_colloid], rra[3*no_of_colloid];
extern const double v0 = 0.04;

//Cluster Specifications
extern int size_cluster[no_of_colloud], n_cluster, identify[no_of_colloid];
extern double dist[no_of_colloid][no_of_colloid];
extern const int file_colloid=0; //set 1 if initial data to be read from file
