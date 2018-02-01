#include "parameters.cuh"

double dtb2 = dt/(mass_colloid*2);

__device__ double power(double x, int r) {
	double ans = 1;
	for(int i = 1; i <=r; i++)
		ans *= x;
	return ans;
}

//PARALLELIZED
__global__ void d_compute_force_md(point *f, int *n_neighbour, int *neighbour[], int *pos_colloid, double sig_colloid, double eps, double *potential_colloid, point len) {
	double mag_f, r_cutoff, fc, ufc, sig_colloid12, sig_colloid6, r;
	r_cutoff = power(2, 1.0/6.0)*sig_colloid, r, fc = 4.0*eps*(12.0*(power(sig_colloid,12)/power(r_cutoff,13)) - 6.0*(power(sig_colloid, 6)/power(r_cutoff, 7)));
	ufc = 4.0*eps*(power(sig_colloid/r_cutoff, 12) - power(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;
	sig_colloid12 = power(sig_colloid, 12), sig_colloid6 = power(sig_colloid, 6);
	point temp, ff;
	double t1, t2;
	int i = blockIdx.x + 1, j = blockIdx.y + 1;
	if(j <= n_neighbour[i]) {
		d_img(&temp, pos_colloid[i] - pos_colloid[neighbour[j][i]], len);
		r = sqrt((temp*temp).sum());
		if(r < r_cutoff) {
			potential_colloid += 4*eps*(power(sig_colloid/r, 12) - power(sig_colloid/r, 6)) - ufc + fc*r;
			t1 = sig_colloid12/power(r,13), t2 = sig_colloid6/power(r, 7);
			mag_f = 4.0*eps*(12.0*t1 - 6.0*t2) - fc;
			ff = (temp*mag_f)/r;
			f[i] += ff, f[neighbour[j][i]] -= ff;
		}   
	}

}
__global__ void d_update_activity_direction(point *ang_vel_colloid, point *ra, double dt) {
	point m[4], b, sb, cb;
	int i = blockIdx.x + 1;
	b  = ang_vel_colloid[i]*dt;
	sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));
	m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
	m[2] =  point((sb.x*sb.y)*cb.z + cb.x*sb.z, (-sb.x*sb.y)*sb.z + cb.x*cb.z, -sb.x*cb.y);
	m[3] =  point((-cb.x*sb.y)*cb.z + sb.x*sb.z, (cb.x*sb.y)*sb.z + sb.x*cb.z, cb.x*cb.y);

	ra[i] = point((m[1]*ra[i]).sum(), (m[2]*ra[i]).sum(), (m[3]*ra[i]).sum());
}
__global__ void d_update_pos_md(point *pos_colloid, point *vel_colloid, point *f, double dt, double mass_colloid, point len) {
	double dt2= dt*dt, ddt = 0.5*dt2/mass_colloid;
	int i = blockIdx.x + 1;
	d_mod(&pos_colloid[i], pos_colloid[i] + vel_colloid[i]*dt + f[i]*ddt, len);
}
__global__ void d_update_pos_mpcd(point *pos_fl, point *vel_fl, double dt, point len) {
	int i = (blockIdx.x * blockDim.x) + threadIdx.x + 1;
	d_mod(&pos_fl[i], pos_fl[i] + vel_fl[i]*dt, len);
}
__global__ void d_update_vel_colloid(point *d_vel, point *d_old_force, point *d_f, double dtb2){
	int i = blockIdx.x + 1;
	d_vel[i] += d_old_force[i] + d_f[i]*dtb2; 
}

void compute_force_md() {
	point *d_pos_colloid; size_t Pitch;
	int *d_n_neighbour, **d_neighbour;
	potential_colloid = 0;
	memset(f, 0, sizeof(int)*(no_of_colloid + 2));
	cudaMalloc(&d_pos_colloid, (no_of_colloid + 2)*sizeof(point));
	cudaMalloc(&d_n_neighbour, (no_of_colloid + 2)*sizeof(int)); 
	cudaMallocPitch(&d_neighbour, &Pitch, (no_of_colloid + 2)*sizeof(int)); 
	cudaMemcpy(d_ang_vel, ang_vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ra, ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	
	d_compute_force_md<<<no_of_colloid, 1>>>(d_f, d_n_neighbour, d_neighbour, d_pos_colloid, sig_colloid, eps, d_potential_colloid, len);
	cudaMemcpy(ra, d_ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_ang_vel), cudaFree(d_ra);

}

void update_activity_direction() {
	point *d_ang_vel, *d_ra;
	cudaMalloc(&d_ang_vel, (no_of_colloid + 2)*sizeof(point));
	cudaMalloc(&d_ra, (no_of_colloid + 2)*sizeof(point)); 
	cudaMemcpy(d_ang_vel, ang_vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ra, ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	

	d_update_activity_direction<<<no_of_colloid, 1>>>(d_ang_vel, d_ra, dt);
	cudaMemcpy(ra, d_ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_ang_vel), cudaFree(d_ra);
}
void update_pos_md() {
	point *d_vel, *d_pos, *d_f;

	cudaMalloc(&d_vel, (no_of_colloid + 2)*sizeof(point));
	cudaMalloc(&d_pos, (no_of_colloid + 2)*sizeof(point));
	cudaMalloc(&d_f,   (no_of_colloid + 2)*sizeof(point)); 
	cudaMemcpy(d_pos, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vel, vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_f, f, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	

	d_update_pos_md<<<no_of_colloid, 1>>>(d_pos, d_vel, d_f, dt, mass_colloid, len);
	cudaMemcpy(pos_colloid, d_pos, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_pos), cudaFree(d_vel), cudaFree(d_f);
}

void update_pos_mpcd() {
	point *d_pos, *d_vel;
	cudaMalloc(&d_vel, (no_of_fluid + 2)*sizeof(point));
	cudaMalloc(&d_pos, (no_of_fluid + 2)*sizeof(point)); 
	cudaMemcpy(d_vel, vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pos, pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	

	d_update_pos_mpcd<<<no_of_fluid/1000, 1000>>>(d_pos, d_vel, dt, len);
	cudaMemcpy(pos_fl, d_pos, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_vel), cudaFree(d_pos);
}

void update_velocity_colloid() {
	point *d_vel, *d_old_force, *d_f;
	cudaMalloc(&d_vel, (no_of_colloid + 2)*sizeof(point));
	cudaMalloc(&d_old_force, (no_of_colloid + 2)*sizeof(point));
	cudaMalloc(&d_f, (no_of_colloid + 2)*sizeof(point)); 
	cudaMemcpy(d_old_force, old_force, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vel, vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_f, f, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	

	d_update_vel_colloid<<<no_of_colloid, 1>>>(d_vel, d_old_force, d_f, dtb2);
	cudaMemcpy(vel_colloid, d_vel, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaFree(d_old_force), cudaFree(d_vel), cudaFree(d_f);
}
