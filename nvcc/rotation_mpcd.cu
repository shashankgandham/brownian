#include "parameters.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

__global__ void cellpart_sync(int **cell_part, int *fluid_no, point len) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= len.prod())
		thrust::sort(thrust::seq, cell_part[i] + 1, cell_part[i] + fluid_no[i] + 1);
}
__global__ void d_cellpart(int **cell_part, int *fluid_no, int no_of_fluid, point *pos_fl, point rr, point len) {
	int cell_no, i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	point temp;
	if(i <= no_of_fluid) {
		cell_no = mod(pos_fl[i] + rr, len).cell(len) + 1;
		cell_part[cell_no][atomicAdd(&(fluid_no[cell_no]), 1) + 1] = i;
	}
}
__global__ void  d_cellvel(point *cell_vel, point *vel_fl, int **cell_part, int *fluid_no, point len) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= len.prod() && fluid_no[i] > 1) {
		for(int	j = 1; j <= fluid_no[i]; j++)
				cell_vel[i] += vel_fl[cell_part[i][j]]/fluid_no[i];
	}

}
__global__ void  d_velfl(point *cell_vel, point *vel_fl, int **cell_part, int *fluid_no, point *rot[4], point len) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	point del_v;
	int k;
	if(i <= len.prod() && fluid_no[i] > 1) {
		for(int	j = 1; j <= fluid_no[i]; j++) {
			k = cell_part[i][j];
				del_v = vel_fl[k] - cell_vel[i];
				vel_fl[k] = cell_vel[i] + point((rot[i][1]*del_v).sum(), (rot[i][2]*del_v).sum(), (rot[i][3]*del_v).sum());
		}
	}

}


void d_rotation_mpcd(point *vel_fl, point *pos_fl, point *cell_vel, point **rot, int *fluid_no, int **cell_part, int no_of_fluid, 
		point len, double kbt, double mass_fl, double *rana, double *ranb) {
	double r[4], ir[4], theta, phi, rho, ct, st, ict;
	point del_v, temp;
//	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	for(int i = 1; i <= len.prod(); i++) {
//	if(i <= len.prod()) {
		if (fluid_no[i] > 1) {
			rho = 2*rana[i] - 1, phi = 2.0*M_PI*ranb[i], theta = M_PI*(130/180.0);
			r[1] = cos(phi)*sqrt(1 - rho*rho), r[2] = sin(phi)*sqrt(1 - rho*rho), r[3] = rho;
			ct = cos(theta), st = sin(theta), ict = 1 - cos(theta);
			ir[1] = ict*r[1], ir[2] = ict*r[2], ir[3] = ict*r[3];
			rot[i][1] = point(ir[1]*r[1] + ct, 	 ir[1]*r[2] - st*r[3], ir[1]*r[3] + st*r[2]);
			rot[i][2] = point(ir[2]*r[1] + st*r[3], ir[2]*r[2] + ct, 	   ir[2]*r[3] - st*r[1]);
			rot[i][3] = point(ir[3]*r[1] - st*r[2], ir[3]*r[2] + st*r[1], ir[3]*r[3] + ct);
		}
	}
}
__global__ void d_rotate(int *fluid_no, int**cell_part, point *vel_fl, point *cell_vel, point len, double mass_fl, double kbt) {
	point del_v;
	double var, scale_fac_mpcd;
	int k, i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= len.prod()) {
		var = 0.0;
		if(fluid_no[i] > 1) {
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[i][j];
				del_v = vel_fl[k] - cell_vel[i];
				var += (del_v*del_v).sum();
			}
			scale_fac_mpcd = sqrt(3.0 * (fluid_no[i] - 1) * kbt/(mass_fl * var));
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[i][j];
				del_v = vel_fl[k] - cell_vel[i];
				vel_fl[k] = cell_vel[i] + del_v*scale_fac_mpcd;
			}
		}
	}
}

void rotation_mpcd() {
	point rr;
	int thr = 256, blk = (no_of_fluid + thr -1)/thr;
	
	cudaDeviceSynchronize();
	rr = rr.random(iv, seed, idum, iy) - point(0.5, 0.5, 0.5);
	cudaMemset(cell_vel, 0, (len.prod() + 2)*sizeof(point));
	cudaMemset(fluid_no, 0, (len.prod() + 2)*sizeof(int));
	d_cellpart<<<blk, thr>>>(cell_part, fluid_no, no_of_fluid, pos_fl, rr, len);
	
	blk = (len.prod() + thr - 1)/thr;
	cellpart_sync<<<blk, thr>>> (cell_part, fluid_no, len);
	d_cellvel<<<blk, thr>>>(cell_vel, vel_fl, cell_part, fluid_no, len);

	cudaDeviceSynchronize();
	for(int  i = 1; i <= len.prod(); i++) {
		if(fluid_no[i] > 1) {
			rana[i] = ran(iv, seed, idum, iy);
			ranb[i] = ran(iv, seed, idum, iy);
		}
	}
 
	d_rotation_mpcd(vel_fl, pos_fl, cell_vel, rot, fluid_no, cell_part, no_of_fluid, len, kbt, mass_fl, rana, ranb);
	d_velfl<<<blk, thr>>>(cell_vel, vel_fl, cell_part, fluid_no, rot, len);
	d_rotate<<<blk, thr>>> (fluid_no, cell_part, vel_fl, cell_vel, len, mass_fl, kbt);
}
