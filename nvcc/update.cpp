#include "parameters.hpp"

void update_activity_direction(){
	point m[4], b, sb, cb;
	for (int i = 1; i <= no_of_colloid; i++) {
		b  = ang_vel_colloid[i]*dt;
		sb = point(sin(b.x), sin(b.y), sin(b.z)), cb = point(cos(b.x), cos(b.y), cos(b.z));

		m[1] =  point(cb.y*cb.z, -cb.y*sb.z, sb.y);
		m[2] =  point(sb.x*sb.y*cb.z + cb.x*sb.z, -sb.x*sb.y*sb.z + cb.x*cb.z, -sb.x*cb.y);
		m[3] =  point(-cb.x*sb.y*cb.z + sb.x*sb.z, cb.x*sb.y*sb.z + sb.x*cb.z, cb.x*cb.y);
		ra[i] = point((m[1]*ra[i]).sum(), (m[2]*ra[i]).sum(), (m[3]*ra[i]).sum());
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
		pos_fl[i] = mod(pos_fl[i] + vel_fl[i]*dt, len);
}

void update_velocity_colloid() {
	for (int i = 1; i <= no_of_colloid; i++)
		vel_colloid[i] += (old_force[i] + f[i])*dt/(mass_colloid*2.0);
}
