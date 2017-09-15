//done
#include "parameters.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>

void rotation_mpcd() { 
	int i,j,k,l1,l2,m, fluid_no[lx*ly*lz], **cell_part, cell_no;
	double cell_velx[lx*ly*lz], cell_vely[lx*ly*lz], cell_velz[lx*ly*lz],r1, r2, r3, theta, phi, rho, prob;
	double rot11,rot12,rot13,rot21,rot22,rot23,rot31,rot32,rot33, del_vx, del_vy, del_vz, mx1,my1,mz1;
	double tmp_pos[3*no_of_fluid + 2], rr1,rr2,rr3,var,del_vx1[no_of_fluid],del_vy1[no_of_fluid],del_vz1[no_of_fluid];
	cell_part = (int **)malloc(sizeof(int *)*maxpart);
	for(int i = 0; i < maxpart; i++) {
		cell_part[i] = (int *)malloc(sizeof(int)*(lx*ly*lz));
	}

	memset(fluid_no, 0, sizeof fluid_no);
	rr1 = rand() - 0.5;
	rr2 = rand() - 0.5;
	rr3 = rand() - 0.5;

	std::copy(pos_fl, pos_fl + sizeof pos_fl, tmp_pos);

	for(int i=1; i <= no_of_fluid; i++) {
		tmp_pos[3*i-2] = tmp_pos[3*i-2] + rr1;
		tmp_pos[3*i-1] =tmp_pos[3*i-1] + rr2;
		tmp_pos[3*i] = tmp_pos[3*i] + rr3;
		tmp_pos[3*i-2] = tmp_pos[3*i-2] - llx*round((tmp_pos[3*i-2]-llxby2)*inv_llx);
		tmp_pos[3*i-1] = tmp_pos[3*i-1] - lly*round((tmp_pos[3*i-1]-llyby2)*inv_lly);
		tmp_pos[3*i] = tmp_pos[3*i] - llz*round((tmp_pos[3*i]-llzby2)*inv_llz);
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
				cell_velx[i] = cell_velx[i] + vel_fl[3*k-2]/float(fluid_no[i]);
				cell_vely[i] = cell_vely[i] + vel_fl[3*k-1]/float(fluid_no[i]);
				cell_velz[i] = cell_velz[i] + vel_fl[3*k]/float(fluid_no[i]);
			}
			rho = 2.0*rand()-1.0;
			phi = 4.0*asin(1.0)*rand();
			r1 = cos(phi)*sqrt(1-rho*rho);
			r2 = sin(phi)*sqrt(1-rho*rho);
			r3 = rho;
			theta = 2.0*asin(1.0)*(130.0/180.0);    
			rot11 = (1.0-cos(theta))*r1*r1 + cos(theta);
			rot12 = (1.0-cos(theta))*r1*r2 - sin(theta)*r3;
			rot13 = (1.0-cos(theta))*r1*r3 + sin(theta)*r2;
			rot21 = (1.0-cos(theta))*r1*r2 + sin(theta)*r3;
			rot22 = (1.0-cos(theta))*r2*r2 + cos(theta);
			rot23 = (1.0-cos(theta))*r2*r3 - sin(theta)*r1;
			rot31 = (1.0-cos(theta))*r1*r3 - sin(theta)*r2;
			rot32 = (1.0-cos(theta))*r2*r3 + sin(theta)*r1;
			rot33 = (1.0-cos(theta))*r3*r3 + cos(theta);
			for(int j = 1;j <= fluid_no[i]; j++) {
				k= cell_part[j][i];
				del_vx = vel_fl[3*k-2] - cell_velx[i];
				del_vy = vel_fl[3*k-1] - cell_vely[i];
				del_vz = vel_fl[3*k] - cell_velz[i];

				vel_fl[3*k-2] = cell_velx[i] + rot11*del_vx + rot12*del_vy + rot13*del_vz;
				vel_fl[3*k-1] = cell_vely[i] + rot21*del_vx + rot22*del_vy + rot23*del_vz;
				vel_fl[3*k]   = cell_velz[i] + rot31*del_vx + rot32*del_vy + rot33*del_vz;
			}
		}
	}
	if (nn % 1 == 0){
		for(int i = 1; i <= lx*ly*lz; i++) {
			var = 0.0;
			if(fluid_no[i] > 1) { 
				for(int j =1; j <= fluid_no[i]; j++) {
					k = cell_part[j][i];
					del_vx1[k] = vel_fl[3*k - 2] - cell_velx[i];
					del_vy1[k] = vel_fl[3*k - 1] - cell_vely[i];
					del_vz1[k] = vel_fl[3*k]   - cell_velz[i];
					var = var + del_vx1[k]*del_vx1[k] + del_vy1[k]*del_vy1[k] + del_vz1[k]*del_vz1[k];
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
	}
}
