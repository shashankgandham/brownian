#pragma once
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <curand_kernel.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

extern curandState_t *state;
extern dim3 thr, thrs, blk;
struct point {
	double x, y, z;
	CUDA_CALLABLE_MEMBER point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z){}
	CUDA_CALLABLE_MEMBER point operator+ (const point &op) { return point(x + op.x, y + op.y, z + op.z); }
	CUDA_CALLABLE_MEMBER point operator- (const point &op) { return point(x - op.x, y - op.y, z - op.z); }
	CUDA_CALLABLE_MEMBER point operator/ (const double &op){ return point(x / op,   y / op,   z / op);   }
	CUDA_CALLABLE_MEMBER point operator/ (const point &op) { return point(x / op.x, y / op.y, z / op.z); }
	CUDA_CALLABLE_MEMBER point operator* (const point &op) { return point(x * op.x, y * op.y, z * op.z); }
	CUDA_CALLABLE_MEMBER point operator* (const double &op){ return point(x * op, y * op, z * op); }
	CUDA_CALLABLE_MEMBER point operator+=(const point &op) { x += op.x, y += op.y, z += op.z; return *this; }
	CUDA_CALLABLE_MEMBER point operator-=(const point &op) { x -= op.x, y -= op.y, z -= op.z; return *this; }
	CUDA_CALLABLE_MEMBER bool operator==(const point &op) { return (x == op.x && y == op.y && z == op.z); }

	CUDA_CALLABLE_MEMBER double sum()  { return (x + y + z); }
	CUDA_CALLABLE_MEMBER double prod() { return (x * y * z); }
	CUDA_CALLABLE_MEMBER int cell(point len) { return int(x) + len.x*int(y) + len.x*len.y*int(z); }
	CUDA_CALLABLE_MEMBER void print()  { printf("%36.32lf %35.32lf %35.32lf\n", x, y, z); }
	 __device__ point rand(curandState_t *state){
		x = curand_uniform_double(state); y = curand_uniform_double(state); z = curand_uniform_double(state);
		return *this;
	}
	CUDA_CALLABLE_MEMBER void next(point len, point inc = point(1, 1, 1), point start = point(1, 1, 1)) {
		x += inc.x;
		if(x > len.x) y += inc.y, x = start.x;
		if(y > len.y) z += inc.z, y = start.y;
	}
};
struct add_point: public thrust::binary_function<point &, point &, point &> {
	CUDA_CALLABLE_MEMBER point operator()(const point &a, const point &b) {
		return point(a.x+b.x, a.y+b.y, a.z+b.z);
	}
};
struct add_double: public thrust::binary_function<double &, double &, double &> {
	CUDA_CALLABLE_MEMBER double operator()(const double &a, const double &b) {
		return a + b;
	}
};
struct mod_value: public thrust::unary_function<point &, double &> {
	CUDA_CALLABLE_MEMBER double operator()(const point &b) {
		return b.x*b.x + b.y*b.y + b.z*b.z;
	}
};

inline CUDA_CALLABLE_MEMBER point round(point a) { return point(round(a.x), round(a.y), round(a.z));}
inline CUDA_CALLABLE_MEMBER point mod(point a, point b) { return a - b*round((a - b/2)/b);} 
inline CUDA_CALLABLE_MEMBER point img(point a, point b) { return a - b*round(a/b);}
inline CUDA_CALLABLE_MEMBER double power(double x, int r) { double ans = 1; for(int i = 1; i <=r; i++) ans *= x; return ans; }

extern point *pos_colloid, *pos_fl, *vel_colloid, *vel_fl, *ang_vel_colloid, *f, *old_force, *ra, len, *cell_vel, *del_v, **rot, *dump_vel_fl, **u, **vc, **om;
extern int n, niter, file, nbin, no_of_fluid, maxpart, no_of_colloid, nbox, **nbr, **up_nbr, *cnt, *up_cnt, nn, **dp;
extern int **neighbour, *n_neighbour, *no_neigh, **neigh_fl, **box_neigh, **box_part, *fluid_no, **cell_part;
extern double kbt, kbt1, ndt, dt, mass_colloid, sig_colloid, eps, v0, sigma, dv, mass_fl, I_colloid, *potential_colloid;

void create_box(), compute_force_md(), fluid_colloid_collision(), initialize(), initialize_fluid(), initialize_colloid();
void neighbour_list_md(), neighbour_list_mpcd(), rotation_mpcd(), run(), tumble(), updown_velocity();
void update_velocity_colloid(), update_pos_md(), update_pos_mpcd() ,update_activity_direction(), initialize_rand();	
