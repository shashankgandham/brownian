#include "parameters.cuh"

inline __device__ point crossmul(point a, point b) {
	return point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline __device__ point stochastic_reflection(point rf, point rs, double mass_fl, double kbt, point len, int *iv, int *seed, int *idum, int *iy) {
	double m_beta = mass_fl/kbt, random_e = power(1 - ran(iv, seed, idum, iy), 2), val, v[4], x[4], z = 2;
	point un, ut, t, n;
	n = rs/sqrt((rs*rs).sum());
	val = sqrt(-log(random_e)/m_beta);

	un = n*val;
	t = img((t.random(iv, seed, idum, iy)*len - rf), len);
	ut = crossmul(un, t);
	ut = ut/sqrt((ut*ut).sum());
	while(z > 1) {
		x[1] = 2.0 * ran(iv, seed, idum, iy) - 1, x[2] = 2.0 * ran(iv, seed, idum, iy) - 1;
		z = x[1]*x[1] + x[2]*x[2];
	}
	z = sqrt((-2.0*log(z))/z);
	v[1] = x[1]*z*sqrt(kbt/mass_fl); v[2] = x[2]*z*sqrt(kbt/mass_fl);
	return ut*v[1] + un;
}
__global__ void d_fluid_colloid_collision(int *no_neigh, point *pos_colloid, point *pos_fl, point *vel_colloid, 
		point *ang_vel_colloid, point *dump_vel_fl, double mass_colloid, double I_colloid,
		double mass_fl, double dt, point *vel_fl, point len, double sigma, int no_of_colloid,
		double kbt, int **neigh_fl, int *iv, int *seed, int *idum, int *iy) {
	point rr, rs, u, omega, vc;
	for(int j = 1; j <= no_of_colloid; j++) {
		vc = omega = point(0, 0, 0);
		for(int i = 1; i <= no_neigh[j]; i++) {
			int l = neigh_fl[i][j];
			rr = img(pos_colloid[j] - pos_fl[l], len);
			if((rr*rr).sum() <= pow(sigma, 2)*0.25) {
				pos_fl[l] = mod(pos_fl[l] - vel_fl[l]*dt* 0.5, len);
				rs = img(pos_fl[l] - pos_colloid[j], len);
				u  = stochastic_reflection(pos_fl[l], rs, mass_fl, kbt, len, iv, seed, idum, iy);
				vel_fl[l] = u + vel_colloid[j] + crossmul(ang_vel_colloid[j], rs);
				vc += (dump_vel_fl[l] - vel_fl[l]);
				u = (dump_vel_fl[l] - vel_fl[l]);
				omega += crossmul(rs, (dump_vel_fl[l] - vel_fl[l]));
				pos_fl[l] = mod(pos_fl[l] + vel_fl[l]*dt*0.5, len);
			}
		}
		vel_colloid[j] 	   += vc*mass_fl/mass_colloid;
		ang_vel_colloid[j] += omega*mass_fl/I_colloid;
	}
}

void fluid_colloid_collision() {
	point *dump_vel_fl;
	cudaMallocManaged(&dump_vel_fl, sizeof(point)*(no_of_fluid + 2));
	cudaMemcpy(dump_vel_fl, vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	d_fluid_colloid_collision<<<1, 1>>> (no_neigh, pos_colloid, pos_fl, vel_colloid, ang_vel_colloid, dump_vel_fl, 
			mass_colloid, I_colloid, mass_fl, dt, vel_fl, len, sigma, no_of_colloid, kbt, neigh_fl, iv, seed, idum, iy);
}
