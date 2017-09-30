#include "parameters.hpp"

void update_activity_direction(){
	double b[4], raa[4], m[4][4];
	for (int i = 1; i <= no_of_colloid; i++) {
		b[1] = ang_vel_colloid[3*i-2]*dt;
		b[2] = ang_vel_colloid[3*i-1]*dt;
		b[3] = ang_vel_colloid[3*i]*dt;

		m[1][1] =  cos(b[2])*cos(b[3]);
		m[1][2] = -cos(b[2])*sin(b[3]);
		m[1][3] =  sin(b[2]);
		m[2][1] =  sin(b[1])*sin(b[2])*cos(b[3]) + cos(b[1])*sin(b[3]);
		m[2][2] = -sin(b[1])*sin(b[2])*sin(b[3]) + cos(b[1])*cos(b[3]);
		m[2][3] = -sin(b[1])*cos(b[2]);
		m[3][1] = -cos(b[1])*sin(b[2])*cos(b[3]) + sin(b[1])*sin(b[3]);
		m[3][2] =  cos(b[1])*sin(b[2])*sin(b[3]) + sin(b[1])*cos(b[3]);
		m[3][3] =  cos(b[1])*cos(b[2]);

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
	for (int i = 1; i <= no_of_fluid; i++) {
		pos_fl[3*i-2] = mod(pos_fl[3*i-2] + vel_fl[3*i-2]*dt, lx);
		pos_fl[3*i-1] = mod(pos_fl[3*i-1] + vel_fl[3*i-1]*dt, ly);
		pos_fl[3*i]   = mod(pos_fl[3*i] + vel_fl[3*i]*dt, lz);
	}
}

void update_velocity_colloid() {
	for (int i = 1; i <= no_of_colloid; i++){
		vel_colloid[3*i-2] += (old_force[3*i-2]+f[3*i-2])*dt/(mass_colloid*2.0);
		vel_colloid[3*i-1] += (old_force[3*i-1]+f[3*i-1])*dt/(mass_colloid*2.0);
		vel_colloid[3*i]   += (old_force[3*i]+f[3*i])*dt/(mass_colloid*2.0);
	}
}
