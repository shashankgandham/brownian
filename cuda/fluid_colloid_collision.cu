#include "parameters.cuh"

inline CUDA_CALLABLE_MEMBER point crossmul(point a, point b) {
	return point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline __device__ point stochastic_reflection(point rf, point rs, double mass_fl, double kbt, point len, curandState_t *state) {
	double m_beta = mass_fl/kbt, random_e = power(1 - curand_uniform_double(state), 2), val, v[4], x[4], z = 2;
	point un, ut, t, n;
	n = rs/sqrt((rs*rs).sum());
	val = sqrt(-log(random_e)/m_beta);

	un = n*val;
	t = img((t.rand(state)*len - rf), len);
	ut = crossmul(un, t);
	ut = ut/sqrt((ut*ut).sum());
	while(z > 1) {
		x[1] = 2.0 * curand_uniform_double(state) - 1, x[2] = 2.0 * curand_uniform_double(state) - 1;
		z = x[1]*x[1] + x[2]*x[2];
	}
	z = sqrt((-2.0*log(z))/z);
	v[1] = x[1]*z*sqrt(kbt/mass_fl); v[2] = x[2]*z*sqrt(kbt/mass_fl);
	return ut*v[1] + un;
}

__global__ void d_fluid_colloid_collision(int *no_neigh, point *pos_colloid, point *pos_fl, 
		point *vel_colloid, point *ang_vel_colloid, point *dump_vel_fl, double mass_colloid, 
		double I_colloid, double mass_fl, double dt, point *vel_fl, point len, double sigma, 
		int no_of_colloid, double kbt, int **neigh_fl, point **vc, point **om, curandState_t *state) {

	point rr, rs, uu;
	int j = blockIdx.y*blockDim.y + threadIdx.y + 1, i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(j <= no_of_colloid) {
		vc[j][0] = om[j][0] = point(0, 0, 0);
		if(i <= no_neigh[j]) {
			vc[j][i] = om[j][i] = point(0, 0, 0);
			int l = neigh_fl[j][i];
			rr = img(pos_colloid[j] - pos_fl[l], len);
			if((rr*rr).sum() <= pow(sigma, 2)*0.25) {
				pos_fl[l] = mod(pos_fl[l] - vel_fl[l]*dt* 0.5, len);
				rs = img(pos_fl[l] - pos_colloid[j], len);
				uu  = stochastic_reflection(pos_fl[l], rs, mass_fl, kbt, len, &state[i]);
				vel_fl[l] = uu + vel_colloid[j] + crossmul(ang_vel_colloid[j], rs);
				pos_fl[l] = mod(pos_fl[l] + vel_fl[l]*dt*0.5, len);
				point t1 = dump_vel_fl[l] - vel_fl[l], t2;
				t2  = crossmul(rs, t1);
				vc[j][i] = t1; om[j][i] = t2;
			}
		}
	}
}
__global__ void d_dump(point *dump_vel_fl, point *vel_fl, int no_of_fluid) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_fluid) dump_vel_fl[i] = vel_fl[i];
}
__global__ void update_fcc(point **vc, point **om, point *vel_colloid, point *ang_vel_colloid, 
		int *no_neigh, int no_of_colloid, double mass_colloid, double mass_fl, double I_colloid) {

	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if(j <= no_of_colloid) {
		vc[j][0] = thrust::reduce(thrust::device, vc[j], vc[j] + no_neigh[j] + 1, point(0, 0, 0), add_point());
		om[j][0] = thrust::reduce(thrust::device, om[j], om[j] + no_neigh[j] + 1, point(0, 0, 0), add_point());
		vel_colloid[j] += vc[j][0]*mass_fl/mass_colloid;
		ang_vel_colloid[j] += om[j][0]*mass_fl/I_colloid;
	}
}
void fluid_colloid_collision() {
	blk = dim3((no_of_fluid + thr.x - 1)/thr.x);
	d_dump<<<blk, thr>>> (dump_vel_fl, vel_fl, no_of_fluid);
	blk = dim3((10000 + thrs.x - 1)/thrs.x, (no_of_colloid + thrs.y - 1)/thrs.y);
	d_fluid_colloid_collision<<<blk, thrs>>>(no_neigh, pos_colloid, pos_fl, vel_colloid, 
			ang_vel_colloid, dump_vel_fl, mass_colloid, I_colloid, mass_fl, dt, vel_fl, 
			len, sigma, no_of_colloid, kbt, neigh_fl, vc, om, state);
	blk = dim3((no_of_colloid + thr.x -1)/thr.x);
	update_fcc<<<blk, thr>>>(vc, om, vel_colloid, ang_vel_colloid, no_neigh, no_of_colloid, 
			mass_colloid, mass_fl, I_colloid);
}
