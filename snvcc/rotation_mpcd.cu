#include "parameters.cuh"

__device__ double d_asin(double x) { return asin(x); }

__global__ void d_rotation_mpcd(point *vel_fl, point *pos_fl, int *fluid_no, int **cell_part, int no_of_fluid, 
						point len, double kbt, double mass_fl) {
	int k, cell_no;
	double r[4], ir[4], theta, phi, rho, var, scale_fac_mpcd, ct, st, ict;
	point *cell_vel, del_v, rr, rot[4], temp;
	cell_vel = (point *)malloc((len.prod() + 2)*sizeof(point));
	memset(fluid_no, 0, (len.prod() + 2)*sizeof(int));
	rr.random(point(0.5, 0.5, 0.5));
	for(int i = 1; i <= no_of_fluid; i++) {
		temp = mod(pos_fl[i] + rr, len);
		cell_no = 1 + temp.cell(len);
		cell_part[++fluid_no[cell_no]][cell_no] = i;
	}
	memset(cell_vel, 0, sizeof cell_vel);
	for(int i = 1; i <= len.prod(); i++) {
		if (fluid_no[i] > 1) {
			for(int j = 1; j <= fluid_no[i]; j++)
				cell_vel[i] += vel_fl[cell_part[j][i]]/fluid_no[i];
			rho = 2*ran() - 1, phi = 4.0*d_asin(1.0)*ran(), theta = 2*d_asin(1)*(130/180.0);
			r[1] = cos(phi)*sqrt(1 - rho*rho), r[2] = sin(phi)*sqrt(1 - rho*rho), r[3] = rho;
			ct = cos(theta), st = sin(theta), ict = 1 - cos(theta);
			ir[1] = ict*r[1], ir[2] = ict*r[2], ir[3] = ict*r[3];
			rot[1] = point(ir[1]*r[1] + ct, 	 ir[1]*r[2] - st*r[3], ir[1]*r[3] + st*r[2]);
			rot[2] = point(ir[2]*r[1] + st*r[3], ir[2]*r[2] + ct, 	   ir[2]*r[3] - st*r[1]);
			rot[3] = point(ir[3]*r[1] - st*r[2], ir[3]*r[2] + st*r[1], ir[3]*r[3] + ct);
			for(int j = 1;j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				del_v = vel_fl[k] - cell_vel[i];
				vel_fl[k] = cell_vel[i] + point((rot[1]*del_v).sum(), (rot[2]*del_v).sum(), (rot[3]*del_v).sum());
			}
		}
	}
	for(int i = 1; i <= len.prod(); i++) {
		var = 0.0;
		if(fluid_no[i] > 1) {
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				del_v = vel_fl[k] - cell_vel[i];
				var += (del_v*del_v).sum();
			}
			scale_fac_mpcd = sqrt(3.0 * (fluid_no[i] - 1) * kbt/(mass_fl * var));
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				del_v = vel_fl[k] - cell_vel[i];
				vel_fl[k] = cell_vel[i] + del_v*scale_fac_mpcd;
			}
		}
	}
}

void rotation_mpcd() {
	point *d_vel_fl, *d_pos_fl;
	int *d_fluid_no, **d_cell_part, *h_cell_part[128];
	cudaMalloc(&d_vel_fl, (no_of_fluid + 2)*sizeof(point));
	cudaMalloc(&d_pos_fl, (no_of_fluid + 2)*sizeof(point));
	cudaMalloc(&d_fluid_no, (len.prod() + 2)*sizeof(int));
	cudaMalloc(&d_cell_part, (maxpart + 2)*sizeof(int *));
	for(int i = 0; i <= maxpart; i++) {
		cudaMalloc(&h_cell_part[i], sizeof(int)*(len.prod() + 2));
		cudaMemcpy(h_cell_part[i], cell_part[i], sizeof(int)*(len.prod() + 2), cudaMemcpyHostToDevice);
	}
	cudaMemcpy(d_cell_part, h_cell_part, (maxpart + 2)*sizeof(int *), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pos_fl, pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vel_fl, vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fluid_no, fluid_no, (len.prod() + 2)*sizeof(point), cudaMemcpyHostToDevice);
	d_rotation_mpcd<<<1, 1>>> (d_vel_fl, d_pos_fl, d_fluid_no, d_cell_part, no_of_fluid, len, kbt, mass_fl);
	cudaMemcpy(pos_fl, d_pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaMemcpy(vel_fl, d_vel_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	cudaMemcpy(fluid_no, d_fluid_no, (len.prod() + 2)*sizeof(point), cudaMemcpyDeviceToHost);
	for(int i = 0; i <= maxpart; i++) {
		cudaMemcpy(cell_part[i], h_cell_part[i], sizeof(int)*(len.prod() + 2), cudaMemcpyDeviceToHost);
		cudaFree(h_cell_part[i]);
	}
	cudaFree(d_pos_fl), cudaFree(d_vel_fl), cudaFree(fluid_no), cudaFree(d_cell_part);
}