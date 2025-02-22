#include "parameters.cuh"
	
double sig_colloid12 = pow(sig_colloid, 12), sig_colloid6 = pow(sig_colloid, 6);
double r_cutoff = pow(2, 1.0/6.0)*sig_colloid, r, fc = 0; //4.0*eps*(12.0*(sig_colloid12/pow(r_cutoff,13)) - 6.0*(sig_colloid6/pow(r_cutoff, 7)));
double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;

__global__ void d_compute_force_md(point *f, int *n_neighbour, int **neighbour, point *pos_colloid, double sig_colloid, double sig_colloid12, double sig_colloid6, double r_cutoff, double fc, double ufc, double eps, double *potential_colloid, point len, int no_of_colloid) {
	point temp, ff;
	double t1, t2, mag_f = 0, r;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_colloid) {
		f[i] = point(0, 0, 0);
		for(int j = 1; j <= n_neighbour[i]; j++) {
			temp = img(pos_colloid[i] - pos_colloid[neighbour[j][i]], len);
			r = sqrt((temp*temp).sum());
			if(r < r_cutoff) {
				*potential_colloid += 4*eps*(power(sig_colloid/r, 12) - power(sig_colloid/r, 6)) - ufc + fc*r;
				t1 = sig_colloid12/power(r,13), t2 = sig_colloid6/power(r, 7);
				mag_f = 4.0*eps*(12.0*t1 - 6.0*t2) - fc;
				ff = (temp*mag_f)/r;
				f[i] += ff, f[neighbour[j][i]] -= ff;
			}   
		}
	}
}
__global__ void d_update_activity_direction(point *ang_vel_colloid, point *ra, double dt, int no_of_colloid) {
	point m[4], b, sb, cb;
	int i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	if(i <= no_of_colloid) {
		b  = ang_vel_colloid[i]*dt;
		sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));
		m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
		m[2] =  point((sb.x*sb.y)*cb.z + cb.x*sb.z, (-sb.x*sb.y)*sb.z + cb.x*cb.z, -sb.x*cb.y);
		m[3] =  point((-cb.x*sb.y)*cb.z + sb.x*sb.z, (cb.x*sb.y)*sb.z + sb.x*cb.z, cb.x*cb.y);
		ra[i] = point((m[1]*ra[i]).sum(), (m[2]*ra[i]).sum(), (m[3]*ra[i]).sum());
	}
}
__global__ void d_update_pos_md(point *pos_colloid, point *vel_colloid, point *f, point *old_force, double dt, double mass_colloid, point len, int no_of_colloid) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_colloid) {
		old_force[i] = f[i];
		pos_colloid[i] = mod(pos_colloid[i] + vel_colloid[i]*dt + f[i]*0.5*dt*dt/mass_colloid, len);
	}
}
__global__ void d_update_pos_mpcd(point *pos_fl, point *vel_fl, double dt, point len, int no_of_fluid) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;	
	if(i <= no_of_fluid) pos_fl[i] = mod(pos_fl[i] + vel_fl[i]*dt, len);
}
__global__ void d_update_vel_colloid(point *vel_colloid, point *old_force, point *f, double dt, double mass_colloid, int no_of_colloid) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;	
	if(i <= no_of_colloid) vel_colloid[i] += (old_force[i] + f[i])*dt/(mass_colloid*2); 
} 

void compute_force_md() {
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	d_compute_force_md<<<blk, thr>>>(f, n_neighbour, neighbour, pos_colloid, sig_colloid, sig_colloid12, sig_colloid6, r_cutoff, fc, ufc, eps, potential_colloid, len, no_of_colloid);
}
void update_activity_direction() {
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	d_update_activity_direction<<<blk, thr>>>(ang_vel_colloid, ra, dt, no_of_colloid);
}
void update_pos_md() {
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	d_update_pos_md<<<blk, thr>>>(pos_colloid, vel_colloid, f, old_force, dt, mass_colloid, len, no_of_colloid);
}
void update_pos_mpcd() {
	blk = dim3((no_of_fluid + thr.x - 1)/thr.x);
	d_update_pos_mpcd<<<blk, thr>>>(pos_fl, vel_fl, dt, len, no_of_fluid);
}
void update_velocity_colloid() {
	blk = dim3((no_of_colloid + thr.x - 1)/thr.x);
	d_update_vel_colloid<<<blk, thr>>>(vel_colloid, old_force, f, dt, mass_colloid, no_of_colloid);
}
