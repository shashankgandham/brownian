#include "parameters.cuh"
//update pos md is parallelized
//update_vel_colloid is parallelized
double mag_f, r_cutoff = pow(2, 1.0/6.0)*sig_colloid, r;
double fc = 4.0*eps*(12.0*(pow(sig_colloid,12)/pow(r_cutoff,13)) - 6.0*(pow(sig_colloid, 6)/pow(r_cutoff, 7)));
double ufc = 4.0*eps*(pow(sig_colloid/r_cutoff, 12) - pow(sig_colloid/r_cutoff, 6)) + fc*r_cutoff;
double sig_colloid12 = pow(sig_colloid, 12), sig_colloid6 = pow(sig_colloid, 6);

double power(double x, int r) {
    double ans = 1;
    for(int i = 1; i <=r; i++)
        ans *= x;
    return ans;
}
void compute_force_md() {
	point temp, ff;
	potential_colloid = 0;
    for(int i = 0; i <= no_of_colloid; i++) f[i] = 0;
    double t1, t2, t3, rp;
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

__global__ void d_update_activity(int no_of_colloid, point *d_ang_vel, point *d_ra, double dt){
	point m[4], b, sb, cb, temp;
	for (int i = 1; i <= no_of_colloid; i++){
		b  = d_ang_vel[i]*dt;
		sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));
		m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
		temp.x = sb.x*sb.y, temp.y = -sb.x*sb.y;
		m[2] =  point(temp.x*cb.z + cb.x*sb.z, temp.y*sb.z + cb.x*cb.z, -sb.x*cb.y);
		temp.x = -cb.x*sb.y, temp.y = cb.x*sb.y;
		m[3] =  point(temp.x*cb.z + sb.x*sb.z, temp.y*sb.z + sb.x*cb.z, cb.x*cb.y);
		temp = d_ra[i];
		d_ra[i] = point((m[1]*temp).sum(), (m[2]*temp).sum(), (m[3]*temp).sum());
		
	}
}

//FLOATING POINT ERROR
void update_activity_direction() {
	point m[4], b, sb, cb, temp;
	for (int i = 1; i <= no_of_colloid; i++) {
		b  = ang_vel_colloid[i]*dt;
		sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));
		m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
		temp.x = sb.x*sb.y, temp.y = -sb.x*sb.y;
		m[2] =  point(temp.x*cb.z + cb.x*sb.z, temp.y*sb.z + cb.x*cb.z, -sb.x*cb.y);
		temp.x = cb.x*sb.y, temp.y = cb.x*sb.y;
		m[3] =  point(sb.x*sb.z - temp.x*cb.z, temp.y*sb.z + sb.x*cb.z, cb.x*cb.y);
		temp = ra[i];
		ra[i] = point((m[1]*temp).sum(), (m[2]*temp).sum(), (m[3]*temp).sum());
		
	}

	/*	
    point *d_ang_vel, *d_ra;

    cudaMalloc(&d_ang_vel, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&d_ra, (no_of_colloid + 2)*sizeof(point));
    
    cudaMemcpy(d_ang_vel, ang_vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_ra, ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	
  	d_update_activity<<<1,1>>>(no_of_colloid, d_ang_vel, d_ra, dt);

  	cudaMemcpy(ra, d_ra, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(d_ra); cudaFree(d_ang_vel);

  	for (int i = 1; i <= no_of_colloid; i++)
  		ra[i].print();
  	*/
}

__global__ void d_update_md(point *cuda_pos, point *cuda_vel, point *cuda_f, double dt, double mass_colloid, point len){
	double dt2 = dt * dt, ddt = 0.5*dt2/mass_colloid;
	d_mod(&cuda_pos[blockIdx.x], cuda_pos[blockIdx.x] += cuda_vel[blockIdx.x]*dt + cuda_f[blockIdx.x]*ddt, len);

}

//PARALLELIZED
void update_pos_md() {
	point *cuda_pos, *cuda_vel, *cuda_f, *mod;
    cudaMalloc(&mod, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&cuda_pos, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&cuda_vel, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&cuda_f, (no_of_colloid + 2)*sizeof(point)); 
  	cudaMemcpy(cuda_pos, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(cuda_vel, vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(cuda_f, f, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	
  	
  	d_update_md<<<no_of_colloid + 1,1>>>(cuda_pos,cuda_vel,cuda_f,dt,mass_colloid, len);
  	cudaMemcpy(pos_colloid, cuda_pos, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(cuda_pos); cudaFree(cuda_vel); cudaFree(cuda_f);	
}


__global__ void d_update_mpcd(point *d_pos, point *d_vel, double dt, point len){
	d_mod(&d_pos[blockIdx.x], d_pos[blockIdx.x] += d_vel[blockIdx.x]*dt, len);	
}

//ERROR
void update_pos_mpcd() {
	/*for (int i = 1; i <= no_of_fluid; i++) {
		pos_fl[i] = mod(pos_fl[i] += vel_fl[i] * dt, len);
		pos_fl[i].print();
	}
	*/
	point *d_pos, *d_vel;
	
	
    cudaMalloc(&d_pos, (no_of_fluid + 2)*sizeof(point));
    cudaMalloc(&d_vel, (no_of_fluid + 2)*sizeof(point));
    cudaMemcpy(d_pos, pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_vel, vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	
  	d_update_mpcd<<<(no_of_fluid +1)/4096,4096>>>(d_pos,d_vel,dt,len);
  	cudaDeviceSynchronize();

  	cudaMemcpy(pos_fl, d_pos, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(d_pos); cudaFree(d_vel); 
  	
}

__global__ void d_update_vel_colloid(point *d_vel, point *d_old_force, point *d_f, double dtb2){
	d_vel[blockIdx.x] += d_old_force[blockIdx.x] + d_f[blockIdx.x]*dtb2; 
}

//PARALLELIZED
void update_velocity_colloid() {
    double dtb2 = dt/(mass_colloid*2);
	
    point *d_vel, *d_old_force, *d_f;
    cudaMalloc(&d_vel, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&d_old_force, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&d_f, (no_of_colloid + 2)*sizeof(point)); 
  	cudaMemcpy(d_old_force, old_force, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_vel, vel_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_f, f, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);  	
  	
  	d_update_vel_colloid<<<no_of_colloid + 1,1>>>(d_vel, d_old_force, d_f, dtb2);

  	cudaMemcpy(vel_colloid, d_vel, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
  	cudaFree(d_old_force); cudaFree(d_vel); cudaFree(d_f);
 }