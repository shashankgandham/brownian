#include "parameters.hpp"

void update_activity_direction(){
	coord raa, m[4], b, sb, cb;
	for (int i = 1; i <= no_of_colloid; i++) {
		b = ang_vel_colloid[i]*dt;
		sb.x = sin(b.x), cb.x = cos(b.x);
		sb.y = sin(b.y), cb.x = cos(b.y);
		sb.z = sin(b.z), cb.x = cos(b.z);

		m[1].x =  cb.y*cb.z;
		m[1].y = -cb.y*sb.z;
		m[1].z =  sb.y;
		m[2].x =  sb.x*sb.y*cb.z + cb.x*sb.z;
		m[2].y = -sb.x*sb.y*sb.z + cb.x*cb.z;
		m[2].z = -sb.x*cb.y;
		m[3].x = -cb.x*sb.y*cb.z + sb.x*sb.z;
		m[3].y =  cb.x*sb.y*sb.z + sb.x*cb.z;
		m[3].z =  cb.x*cb.y;
		ra[i] = coord((m[1]*ra[i]).sum(), (m[2]*ra[i]).sum(), (m[3]*ra[i]).sum());
	}
}

void update_pos_md() {
	for(int i = 1; i <= no_of_colloid; i++) {
		pos_colloid[i] +=  vel_colloid[i]*dt + f[i]*0.05*dt*dt/mass_colloid;
		pos_colloid[i]  =  mod(pos_colloid[i], len);
	}
}

void update_pos_mpcd() {
	for (int i = 1; i <= no_of_fluid; i++)
		pos_fl[i] = mod(pos_fl[3*i] + vel_fl[3*i]*dt, len);
}

void update_velocity_colloid() {
	for (int i = 1; i <= no_of_colloid; i++)
		vel_colloid[i] += (old_force[i] + f[i])*dt/(mass_colloid*2.0);
}
