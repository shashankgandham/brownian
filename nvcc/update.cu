#include "parameters.cuh"

double dtb2 = dt/(mass_colloid*2);

__global__ void d_compute_force_md(point *f, int *n_neighbour, int **neighbour, point *pos_colloid, double sig_colloid, double eps, double *potential_colloid, point len, int no_of_colloid) {
	double mag_f, r_cutoff, fc, ufc, sig_colloid12, sig_colloid6, r;
	r_cutoff = power(2, 1.0/6.0)*sig_colloid, r, fc = 4.0*eps*(12.0*(power(sig_colloid,12)/power(r_cutoff,13)) - 6.0*(power(sig_colloid, 6)/power(r_cutoff, 7)));
	ufc = 4.0*eps*(power(sig_colloid/r_cutoff, 12) - power(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;
	sig_colloid12 = power(sig_colloid, 12), sig_colloid6 = power(sig_colloid, 6);
	point temp, ff;
	double t1, t2;
	for(int i = 1; i <= no_of_colloid; i++) {
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
	for(int i = 1; i <= no_of_colloid; i++) {
		b  = ang_vel_colloid[i]*dt;
		sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));
		m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
		m[2] =  point((sb.x*sb.y)*cb.z + cb.x*sb.z, (-sb.x*sb.y)*sb.z + cb.x*cb.z, -sb.x*cb.y);
		m[3] =  point((-cb.x*sb.y)*cb.z + sb.x*sb.z, (cb.x*sb.y)*sb.z + sb.x*cb.z, cb.x*cb.y);
		ra[i] = point((m[1]*ra[i]).sum(), (m[2]*ra[i]).sum(), (m[3]*ra[i]).sum());
	}
}
__global__ void d_update_pos_md(point *pos_colloid, point *vel_colloid, point *f, double dt, double mass_colloid, point len, int no_of_colloid) {
	double dt2= dt*dt, ddt = 0.5*dt2/mass_colloid;
	for(int i = 1; i <= no_of_colloid; i++)
		pos_colloid[i] = mod(pos_colloid[i] + vel_colloid[i]*dt + f[i]*ddt, len);
}
__global__ void d_update_pos_mpcd(point *pos_fl, point *vel_fl, double dt, point len, int no_of_fluid) {
	for(int i = 1; i <= no_of_fluid; i++)
	pos_fl[i] = mod(pos_fl[i] + vel_fl[i]*dt, len);
}
__global__ void d_update_vel_colloid(point *vel_colloid, point *old_force, point *f, double dtb2, int no_of_colloid){
	for(int i = 1; i <= no_of_colloid; i++) 
		vel_colloid[i] += old_force[i] + f[i]*dtb2; 
}

void compute_force_md() { 
	double *d_potential_colloid;
	cudaMalloc(&d_potential_colloid, sizeof(double));
	cudaMemcpy(d_potential_colloid, &potential_colloid, sizeof(double), cudaMemcpyHostToDevice);	
	d_compute_force_md<<<1, 1>>>(f, n_neighbour, neighbour, pos_colloid, sig_colloid, eps, d_potential_colloid, len, no_of_colloid);
	cudaMemcpy(&potential_colloid, d_potential_colloid, sizeof(double), cudaMemcpyHostToDevice);	

}

void update_activity_direction() {
	d_update_activity_direction<<<1, 1>>>(ang_vel_colloid, ra, dt, no_of_colloid);
}
void update_pos_md() {
	d_update_pos_md<<<1, 1>>>(pos_colloid, vel_colloid, f, dt, mass_colloid, len, no_of_colloid);
}

void update_pos_mpcd() {
	d_update_pos_mpcd<<<1, 1>>>(pos_fl, vel_fl, dt, len, no_of_fluid);
}

void update_velocity_colloid() {
	d_update_vel_colloid<<<1, 1>>>(vel_colloid, old_force, f, dtb2, no_of_colloid);
}
