#include "parameters.cuh"

__device__ double d_asin(double x) { return asin(x); }

__global__ void d_rotation_mpcd(point *vel_fl, point *pos_fl, int *fluid_no, int **cell_part, int no_of_fluid, 
						point len, double kbt, double mass_fl, int *iv, int *seed, int *idum, int *iy) {
	int k, cell_no;
	double r[4], ir[4], theta, phi, rho, var, scale_fac_mpcd, ct, st, ict;
	point *cell_vel, del_v, rr, rot[4], temp;
	cell_vel = (point *)malloc((len.prod() + 2)*sizeof(point));
	memset(fluid_no, 0, (len.prod() + 2)*sizeof(int));
	rr = rr.random(iv, seed, idum, iy) - point(0.5, 0.5, 0.5);
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
			rho = 2*ran(iv, seed, idum, iy) - 1, phi = 4.0*d_asin(1.0)*ran(iv, seed, idum, iy), theta = 2*d_asin(1)*(130/180.0);
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
	d_rotation_mpcd<<<1, 1>>> (vel_fl, pos_fl, fluid_no, cell_part, no_of_fluid, len, kbt, mass_fl, iv, seed, idum, iy);
}
