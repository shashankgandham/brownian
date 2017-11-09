#include "parameters.hpp"

void update_activity_direction(){
	double sb[4], cb[4], b[4], raa[4], m[4][4];
	for (int i = 1; i <= no_of_colloid; i++) {
		b[1] = ang_vel_colloid[3*i-2]*dt;
		b[2] = ang_vel_colloid[3*i-1]*dt;
		b[3] = ang_vel_colloid[3*i]*dt;
		for(int j = 1; j <= 3; j++)
			sb[j] = sin(b[j]), cb[j] = cos(b[j]);

		m[1][1] =  cb[2]*cb[3], m[1][2] = -cb[2]*sb[3], m[1][3] =  sb[2];
		m[2][1] =  sb[1]*sb[2]*cb[3] + cb[1]*sb[3];
		m[2][2] = -sb[1]*sb[2]*sb[3] + cb[1]*cb[3];
		m[3][1] = -cb[1]*sb[2]*cb[3] + sb[1]*sb[3];
		m[3][2] =  cb[1]*sb[2]*sb[3] + sb[1]*cb[3];
		m[2][3] = -sb[1]*cb[2], m[3][3] =  cb[1]*cb[2];

		raa[0] = ra[3*i-2], raa[1] = ra[3*i-1], raa[2] = ra[3*i];
		ra[3*i-2] = m[1][1]*raa[0] + m[1][2]*raa[1] + m[1][3]*raa[2];
		ra[3*i-1] = m[2][1]*raa[0] + m[2][2]*raa[1] + m[2][3]*raa[2];
		ra[3*i] = m[3][1]*raa[0] + m[3][2]*raa[1] + m[3][3]*raa[2];
	}
}

void update_pos_md() {
	for(int i = 1; i <= no_of_colloid; i++) {
		pos_colloid[3*i - 2] += vel_colloid[3*i - 2]*dt + f[3*i -2]*0.05*dt*dt/mass_colloid;
		pos_colloid[3*i - 1] += vel_colloid[3*i - 1]*dt + f[3*i - 1]*0.05*dt*dt/mass_colloid;
		pos_colloid[3*i] += vel_colloid[3*i]*dt + f[3*i]*0.05*dt*dt/mass_colloid;
		pos_colloid[3*i - 2] = mod(pos_colloid[3*i - 2], lx);
		pos_colloid[3*i - 1] = 	mod(pos_colloid[3*i - 1], ly);
		pos_colloid[3*i] = 	mod(pos_colloid[3*i], lz);
	}
}

void update_pos_mpcd(){
	check_velfl();
	for (int i = 1; i <= no_of_fluid; i++) {
		pos_fl[3*i-2] = mod(pos_fl[3*i-2] + vel_fl[3*i-2]*dt, lx);
		pos_fl[3*i-1] = mod(pos_fl[3*i-1] + vel_fl[3*i-1]*dt, ly);
		pos_fl[3*i]   = mod(pos_fl[3*i] + vel_fl[3*i]*dt, lz);
		if(!(pos_fl[3*i - 2] > 0 && pos_fl[3*i - 2] <= 30))
			printf("%d --> %lf %lf\n", i, pos_fl[3*i - 2], vel_fl[3*i - 2]);
		if(!(pos_fl[3*i - 1] > 0 && pos_fl[3*i - 1] <= 30))
			printf("%d --> %lf %lf\n", i, pos_fl[3*i - 1], vel_fl[3*i - 1]);
		if(!(pos_fl[3*i - 0] > 0 && pos_fl[3*i - 0] <= 30))
			printf("%d --> %lf %lf\n", i, pos_fl[3*i - 0], vel_fl[3*i - 0]);
	}
}

void update_velocity_colloid() {
	for (int i = 1; i <= no_of_colloid; i++){
		vel_colloid[3*i-2] += (old_force[3*i-2]+f[3*i-2])*dt/(mass_colloid*2.0);
		vel_colloid[3*i-1] += (old_force[3*i-1]+f[3*i-1])*dt/(mass_colloid*2.0);
		vel_colloid[3*i]   += (old_force[3*i]+f[3*i])*dt/(mass_colloid*2.0);
	}
}
