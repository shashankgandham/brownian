#include "parameters.hpp"
#include <assert.h>

void rotation_mpcd() {
	int j, k, **cell_part, cell_no, *fluid_no;
	double cell_velx[lx*ly*lz], cell_vely[lx*ly*lz], cell_velz[lx*ly*lz], r[4], theta, phi, rho, rot[4][4], rr[4];
	double del_vx, del_vy, del_vz, *tmp_pos, var, del_vx1[no_of_fluid], del_vy1[no_of_fluid], del_vz1[no_of_fluid];
	double scale_fac_mpcd;

	tmp_pos = (double *)malloc(sizeof(double)*(3*no_of_fluid + 2));
	fluid_no = (int *)malloc(sizeof(int)*(lx*ly*lz + 1));
	cell_part = (int **)malloc(sizeof(int *)*(maxpart + 2));
	memset(fluid_no, 0, (lx*ly*lz + 2)*sizeof(int));
	std::copy(pos_fl, pos_fl + 3*no_of_fluid + 2, tmp_pos);

	for(int i = 0; i <= maxpart; i++)
		cell_part[i] = (int *)malloc(sizeof(int)*(lx*ly*lz + 2));

	for(int i = 1; i <= 3; i++)
		rr[i] = ran() - 0.5;

	for(int i=1; i <= no_of_fluid; i++) {
		tmp_pos[3*i-2] = mod(tmp_pos[3*i - 2] + rr[1], lx);
		tmp_pos[3*i-1] = mod(tmp_pos[3*i-1] + rr[2], ly);
		tmp_pos[3*i] = mod(tmp_pos[3*i] + rr[3], lz);
	}
	for(int i=1; i <= no_of_fluid; i++) {
		cell_no = 1 + int(tmp_pos[3*i-2]) + lx*int(tmp_pos[3*i-1]) + lx*ly*int(tmp_pos[3*i]);
		fluid_no[cell_no] = fluid_no[cell_no] + 1;
		j = fluid_no[cell_no];
		cell_part[j][cell_no] = i;
	}
	memset(cell_velx, 0, sizeof cell_velx);
	memset(cell_vely, 0, sizeof cell_vely);
	memset(cell_velz, 0, sizeof cell_velz);
	for(int i = 1; i <= lx*ly*lz; i++) {
		if (fluid_no[i] > 1) {
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				cell_velx[i] += vel_fl[3*k-2]/fluid_no[i];
				cell_vely[i] += vel_fl[3*k-1]/fluid_no[i];
				cell_velz[i] += vel_fl[3*k]/fluid_no[i];
			}
			rho = 2*ran() - 1, phi = 4.0*asin(1.0)*ran();
			r[1] = cos(phi)*sqrt(1-rho*rho);
			r[2] = sin(phi)*sqrt(1-rho*rho);
			r[3] = rho;
			theta = 2*asin(1)*(130/180.0);
			rot[1][1] = (1 - cos(theta))*r[1]*r[1] + cos(theta);
			rot[1][2] = (1 - cos(theta))*r[1]*r[2] - sin(theta)*r[3];
			rot[1][3] = (1 - cos(theta))*r[1]*r[3] + sin(theta)*r[2];
			rot[2][1] = (1 - cos(theta))*r[1]*r[2] + sin(theta)*r[3];
			rot[2][2] = (1 - cos(theta))*r[2]*r[2] + cos(theta);
			rot[2][3] = (1 - cos(theta))*r[2]*r[3] - sin(theta)*r[1];
			rot[3][1] = (1 - cos(theta))*r[1]*r[3] - sin(theta)*r[2];
			rot[3][2] = (1 - cos(theta))*r[2]*r[3] + sin(theta)*r[1];
			rot[3][3] = (1 - cos(theta))*r[3]*r[3] + cos(theta);
			for(int j = 1;j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				del_vx = vel_fl[3*k-2] - cell_velx[i];
				del_vy = vel_fl[3*k-1] - cell_vely[i];
				del_vz = vel_fl[3*k] - cell_velz[i];

				vel_fl[3*k-2] = cell_velx[i] + rot[1][1]*del_vx + rot[1][2]*del_vy + rot[1][3]*del_vz;
				vel_fl[3*k-1] = cell_vely[i] + rot[2][1]*del_vx + rot[2][2]*del_vy + rot[2][3]*del_vz;
				vel_fl[3*k]   = cell_velz[i] + rot[3][1]*del_vx + rot[3][2]*del_vy + rot[3][3]*del_vz;
			}
		}
	}
	for(int i = 1; i <= lx*ly*lz; i++) {
		var = 0.0;
		if(fluid_no[i] > 1) {
			for(int j =1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				del_vx1[k] = vel_fl[3*k - 2] - cell_velx[i];
				del_vy1[k] = vel_fl[3*k - 1] - cell_vely[i];
				del_vz1[k] = vel_fl[3*k]   - cell_velz[i];
				var += del_vx1[k]*del_vx1[k] + del_vy1[k]*del_vy1[k] + del_vz1[k]*del_vz1[k];
			}
			scale_fac_mpcd = sqrt(3.0 * (fluid_no[i] - 1) * kbt/(mass_fl * var));
			for(int j = 1; j <= fluid_no[i]; j++) {
				k = cell_part[j][i];
				vel_fl[3*k-2] = cell_velx[i] + del_vx1[k]*scale_fac_mpcd;
				vel_fl[3*k-1] = cell_vely[i] + del_vy1[k]*scale_fac_mpcd;
				vel_fl[3*k]   = cell_velz[i] + del_vz1[k]*scale_fac_mpcd;
			}
		}
	}
	for(int i = 0; i <= maxpart; i++)
		free(cell_part[i]);
	free(tmp_pos);
	free(cell_part);
	free(fluid_no);
}
