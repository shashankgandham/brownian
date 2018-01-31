#include "parameters.cuh"

double mag_f, r_cutoff = pow(2, 1.0/6.0)*sig_colloid, r;
double fc = 4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;
double sig_colloid12 = pow(sig_colloid, 12), sig_colloid6 = pow(sig_colloid, 6);
double dtb2 = dt/(mass_colloid*2);

double power(double x, int r) {
    double ans = 1;
    for(int i = 1; i <=r; i++)
        ans *= x;
    return ans;
}
void compute_force_md() {
	point temp, ff;
	potential_colloid = 0;
    memset(f, 0, sizeof(int)*(no_of_colloid + 2));
    double t1, t2;
	for(int i = 1; i <= no_of_colloid; i++) {
		for(int j = 1; j <= n_neighbour[i]; j++) {
			temp = img(pos_colloid[i] - pos_colloid[neighbour[j][i]], len);
			r = sqrt((temp*temp).sum());
			if(r < r_cutoff) {
				potential_colloid += 4*eps*(pow(sig_colloid/r, 12) - pow(sig_colloid/r, 6)) - ufc + fc*r;
				t1 = sig_colloid12/power(r,13), t2 = sig_colloid6/power(r, 7);
                mag_f = 4.0*eps*(12.0*t1 - 6.0*t2) - fc;
                ff = (temp*mag_f)/r;
				f[i] += ff, f[neighbour[j][i]] -= ff;
            }
		}   
	}
}
//PARALLELIZED
__global__ void d_update_activity_direction(point *ang_vel_colloid, point *ra, double dt) {
	point m[4], b, sb, cb;
	b  = ang_vel_colloid[blockIdx.x]*dt;
	sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));
	m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
	m[2] =  point(sb.x*sb.y*cb.z + cb.x*sb.z, -sb.x*sb.y*sb.z + cb.x*cb.z, -sb.x*cb.y);
	m[3] =  point(-cb.x*sb.y*cb.z + sb.x*sb.z, cb.x*sb.y*sb.z + sb.x*cb.z, cb.x*cb.y);
	ra[blockIdx.x] = point((m[1]*ra[blockIdx.x]).sum(), (m[2]*ra[blockIdx.x]).sum(), (m[3]*ra[blockIdx.x]).sum());
}
__global__ void d_update_pos_md(point *pos_colloid, point *vel_colloid, point *f, double dt, double mass_colloid, point len) {
	double dt2= dt*dt, ddt = 0.5*dt2/mass_colloid;
	d_mod(&pos_colloid[blockIdx.x], pos_colloid[blockIdx.x] + vel_colloid[blockIdx.x]*dt + f[blockIdx.x]*ddt, len);
}
__global__ void d_update_pos_mpcd(point *pos_fl, point *vel_fl, double dt, point len) {
  	d_mod(&pos_fl[blockIdx.x], pos_fl[blockIdx.x] + vel_fl[blockIdx.x]*dt, len);
}
__global__ void d_update_vel_colloid(point *d_vel, point *d_old_force, point *d_f, double dtb2){
	d_vel[blockIdx.x] += d_old_force[blockIdx.x] + d_f[blockIdx.x]*dtb2; 
}

void update_activity_direction() {
	point *d_ang_vel, *d_ra;
    cudaMalloc(&d_ang_vel, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&d_ra,   	   (no_of_colloid + 2)*sizeof(point)); 
  	cudaMemcpy(d_ang_vel, ang_vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_ra, 	  ra			 , (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	
  	
  	d_update_activity_direction<<<no_of_colloid + 1, 1>>>(d_ang_vel, d_ra, dt);
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
  	cudaMemcpy(d_f, 			f, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	
  	
  	d_update_pos_md<<<no_of_colloid + 1, 1>>>(d_pos, d_vel, d_f, dt, mass_colloid, len);
  	cudaMemcpy(pos_colloid, d_pos, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(d_pos), cudaFree(d_vel), cudaFree(d_f);
	for(int i = 1; i <=10; i++)
		pos_colloid[i].print();
	printf("\n");
	exit(0);

}

void update_pos_mpcd() {
	point *d_pos, *d_vel;
    cudaMalloc(&d_vel, (no_of_fluid + 2)*sizeof(point));
    cudaMalloc(&d_pos, (no_of_fluid + 2)*sizeof(point)); 
  	cudaMemcpy(d_vel, vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_pos, pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	
  	
  	d_update_pos_mpcd<<<no_of_fluid + 1,1>>>(d_pos, d_vel, dt, len);
  	cudaMemcpy(vel_fl, d_vel, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
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
  	
  	d_update_vel_colloid<<<no_of_colloid + 1, 1>>>(d_vel, d_old_force, d_f, dtb2);
  	cudaMemcpy(vel_colloid, d_vel, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(d_old_force), cudaFree(d_vel), cudaFree(d_f);
}
