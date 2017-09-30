#include "parameters.hpp"

void stochastic_reflection(double *u, double rfx, double rfy, double rfz, double rsx, double rsy, double rsz) {
	double nx, ny, nz, unx, uny, unz, m_beta, random_e, ut, utx, uty, utz, txx, tyy, tzz, val, v[2], den, x[4], z = 2;

	m_beta = mass_fl/kbt;
	den = sqrt(rsx*rsx + rsy*rsy + rsz*rsz);
	random_e = pow(1 - ran(), 2);
	val = sqrt(-log(random_e)/m_beta);

	unx = val*rsx/den;
	uny = val*rsy/den;
	unz = val*rsz/den;

	txx = mod(ran()*lx - rfx, lx);
	tyy = mod(ran()*ly - rfy, ly);
	tzz = mod(ran()*lz - rfz, lz);

	utx = uny*tzz - tyy*unz;
	uty = unz*txx - unx*tzz;
	utz = unx*tyy - txx*uny;

	ut = sqrt(utx*utx + uty*uty + utz*utz);
	utx /= ut, uty /= ut, utz /= ut;

	while(z > 1) {
		x[1] = 2.0 * ran() - 1, x[2] = 2.0 * ran() - 1;
		z = x[1]*x[1] + x[2]*x[2];
	}
	z = sqrt((-2.0*log(z))/z);
	v[1] = x[1]*z*sqrt(kbt/mass_fl);
	v[2] = x[2]*z*sqrt(kbt/mass_fl);

	utx = v[1]*utx;
	uty = v[1]*uty;
	utz = v[1]*utz;

	u[1] = unx + utx;
	u[2] = uny + uty;
	u[3] = unz + utz;
}

void fluid_colloid_collision(double I_colloid) {
	double dump_vel_fl[3*no_of_fluid], vc[4], omega[3], v[4], rr, rrx, rry, rrz, u[4];
	double rfx, rfy, rfz, rsx, rsy, rsz, rcx, rcy, rcz;

	for (int i = 3; i <= 3*no_of_fluid; i++)
		dump_vel_fl[i] = vel_fl[i];

	for (int j = 1; j <= no_of_colloid; j++){
		vc[1] = vc[2] = vc[3] = 0;
		omega[1] = omega[2] = omega[3] = 0;

		for (int i = 1; i <= no_neigh[j]; i++){
			int lll = neigh_fl[i][j];
			rrx = mod(pos_colloid[3*j-2] - pos_fl[3*lll-2], lz);
			rry = mod(pos_colloid[3*j-1] - pos_fl[3*lll-1], ly);
			rrz = mod(pos_colloid[3*j] - pos_fl[3*lll], lz);
			rr = rrx*rrx + rry*rry + rrz*rrz;

			if(rr <= pow(sigma, 2)*0.25) {
				pos_fl[3*lll - 2] = mod(pos_fl[3*lll - 2] - vel_fl[3*lll - 2]*dt*0.5, lx);
				pos_fl[3*lll - 1] = mod(pos_fl[3*lll - 1] - vel_fl[3*lll - 1]*dt*0.5, ly);
				pos_fl[3*lll] = mod(pos_fl[3*lll] - vel_fl[3*lll]*dt*0.5, lz);

				rcx = pos_colloid[3*j-2], rcy = pos_colloid[3*j-1], rcz = pos_colloid[3*j];
				rfx = pos_fl[3*lll-2], rfy = pos_fl[3*lll-1], rfz = pos_fl[3*lll];
				rsx = mod(rfx - rcx, lx), rsy = mod(rfy - rcy, ly), rsz = mod(rfz - rcz, lz);

				stochastic_reflection(u, rfx, rfy, rfz, rsx, rsy, rsz);

				/*velocity of fluid in streaming step */
				vel_fl[3*lll-2] = u[1] + vel_colloid[3*j-2] + ang_vel_colloid[3*j-1]*rsz - ang_vel_colloid[3*j]*rsy;
				vel_fl[3*lll-1] = u[2] + vel_colloid[3*j-1] - ang_vel_colloid[3*j-2]*rsz + ang_vel_colloid[3*j]*rsx;
				vel_fl[3*lll] = u[3] + vel_colloid[3*j] + ang_vel_colloid[3*j-2]*rsy - ang_vel_colloid[3*j-1]*rsx;

				v[1] = dump_vel_fl[3*lll-2] - vel_fl[3*lll-2];
				v[2] = dump_vel_fl[3*lll-1] -vel_fl[3*lll-1];
				v[3] = dump_vel_fl[3*lll] - vel_fl[3*lll];

				/*calculation for linear and angular velocity of the colloidal particle*/
				vc[1] += v[1], vc[2] += v[2], vc[3] += v[3];
				omega[1] +=  rsy*v[3] - v[2]*rsz;
				omega[2] += -rsx*v[3] + v[1]*rsz;
				omega[3] +=  rsx*v[2] - v[1]*rsy;

				/*updating position for next dt/2 time */
				pos_fl[3*lll - 2] = mod(pos_fl[3*lll - 2] + vel_fl[3*lll - 2]*dt*0.5, lx);
				pos_fl[3*lll - 1] = mod(pos_fl[3*lll - 1] + vel_fl[3*lll - 1]*dt*0.5, ly);
				pos_fl[3*lll] = mod(pos_fl[3*lll] + vel_fl[3*lll]*dt*0.5, lz);
			}
		}
		vel_colloid[3*j-2] += vel_colloid[3*j-2] + vc[1]*mass_fl/mass_colloid;
		vel_colloid[3*j-1] += vel_colloid[3*j-1] + vc[2]*mass_fl/mass_colloid;
		vel_colloid[3*j] += vel_colloid[3*j] + vc[3]*mass_fl/mass_colloid;

		/*angular velocity of colloidal particle */

		ang_vel_colloid[3*j-2] += omega[1]*mass_fl/I_colloid;
		ang_vel_colloid[3*j-1] += omega[2]*mass_fl/I_colloid;
		ang_vel_colloid[3*j] += omega[3]*mass_fl/I_colloid;
	}
}
