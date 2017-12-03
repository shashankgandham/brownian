#include "parameters.hpp"

void rotation_mpcd() {
	int j, k, **cell_part, cell_no, *fluid_no;
	double r[4], theta, phi, rho, var, scale_fac_mpcd, ct, st, ict;
	point *tmp_pos, cell_vel[(int)len.prod() + 1], del_v, del_v1[no_of_fluid + 1], rr, rot[4];

	tmp_pos   = (point *)malloc(sizeof(point)*(no_of_fluid + 2));
	fluid_no  = (int *)malloc(sizeof(int) *(len.prod() + 2));
	cell_part = (int **)malloc(sizeof(int *)*(maxpart  + 2));
	memset(fluid_no, 0, (len.prod() + 2)*sizeof(int));
	std::copy(pos_fl, pos_fl + no_of_fluid + 2, tmp_pos);

	for(int i = 0; i <= maxpart; i++)
		cell_part[i] = (int *)malloc(sizeof(int)*(len.prod() + 2));

	rr.random(0.5);
	for(int i = 1; i <= no_of_fluid; i++) {
		tmp_pos[i] = mod(tmp_pos[i] + rr, len);
		cell_no = 1 + int(tmp_pos[i].x) + len.x*int(tmp_pos[i].y) + len.x*len.y*int(tmp_pos[i].z);
		cell_part[++fluid_no[cell_no]][cell_no] = i;
	}
	memset(cell_vel, 0, sizeof cell_vel);
	for(int i = 1; i <= len.prod(); i++) {
		if (fluid_no[i] > 1) {
			for(int j = 1; j <= fluid_no[i]; j++)
				cell_vel[i] += vel_fl[cell_part[j][i]]/fluid_no[i];

			rho = 2*ran() - 1, phi = 4.0*asin(1.0)*ran(), theta = 2*asin(1)*(130/180.0);
			r[1] = cos(phi)*sqrt(1 - rho*rho), r[2] = sin(phi)*sqrt(1 - rho*rho), r[3] = rho;
			ct = cos(theta), st = sin(theta), ict = 1 - cos(theta);

			rot[1] = point(ict*r[1]*r[1] + ct, ict*r[1]*r[2] - st*r[3], ict*r[1]*r[3] + st*r[2]);
			rot[2] = point(ict*r[2]*r[1] + st*r[3], ict*r[2]*r[2] + ct, ict*r[2]*r[3] - st*r[1]);
			rot[3] = point(ict*r[3]*r[1] - st*r[2], ict*r[3]*r[2] + st*r[1], ict*r[3]*r[3] + ct);

			for(int j = 1;j <= fluid_no[i]; j++) {
				del_v = vel_fl[cell_part[j][i]] - cell_vel[i];
				vel_fl[cell_part[j][i]] = cell_vel[i] + point((rot[1]*del_v).sum(), (rot[2]*del_v).sum(), (rot[3]*del_v).sum());
			}
		}
	}
	for(int i = 1; i <= len.prod(); i++) {
		var = 0.0;
		if(fluid_no[i] > 1) {
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				del_v1[k] = vel_fl[k] - cell_vel[i];
				var += (del_v1[k]*del_v1[k]).sum();
			}
			scale_fac_mpcd = sqrt(3.0 * (fluid_no[i] - 1) * kbt/(mass_fl * var));
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				vel_fl[k] = cell_vel[i] + del_v1[k]*scale_fac_mpcd;
			}
		}
	}
	for(int i = 0; i <= maxpart; i++)
		free(cell_part[i]);
	free(tmp_pos), free(cell_part), free(fluid_no);
}
