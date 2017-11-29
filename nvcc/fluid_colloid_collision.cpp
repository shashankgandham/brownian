#include "parameters.hpp"

coord crossmul(coord a, coord b) {
	return coord(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

void stochastic_reflection(coord *u, coord rf, coord rs) {
	double m_beta = mass_fl/kbt, random_e = pow(1 - ran(), 2), val, v[4], den, x[4], z = 2;
	coord un, ut, t;
	den = sqrt((rs*rs).sum());
	val = sqrt(-log(random_e)/m_beta);

	un = (rs*val)/den;

	t = img(len*ran() - rf, len);
	ut = crossmul(un, t);
	ut = ut/sqrt((ut*ut).sum());

	while(z > 1) {
		x[1] = 2.0 * ran() - 1, x[2] = 2.0 * ran() - 1;
		z = x[1]*x[1] + x[2]*x[2];
	}

	z = sqrt((-2.0*log(z))/z);
	v[1] = x[1]*z*sqrt(kbt/mass_fl);
	v[2] = x[2]*z*sqrt(kbt/mass_fl);

	ut = ut*v[1];
	*u = ut + un;
}

void fluid_colloid_collision() {
	coord rr, rf, rs, rc, dump_vel_fl[no_of_fluid], u, v, omega, vc;

	for (int i = 1; i < no_of_fluid; i++)
		dump_vel_fl[i] = vel_fl[i];

	for (int j = 1; j <= no_of_colloid; j++) {
		vc = omega = coord(0, 0, 0);
		for (int i = 1; i <= no_neigh[j]; i++) {
			int l = neigh_fl[i][j];
			rr = img(pos_colloid[j] - pos_fl[l], len);
			if((rr*rr).sum() <= pow(sigma, 2)*0.25) {
				pos_fl[l] = mod(pos_fl[l] - vel_fl[l]*dt*0.5, len);
				rc = pos_colloid[j], rf = pos_fl[l];
				rs = img(rf - rc, len);
				stochastic_reflection(&u, rf, rs);
				vel_fl[l] = u + vel_colloid[j] + crossmul(ang_vel_colloid[j], rs);

				v   = dump_vel_fl[l] - vel_fl[l];
				vc += v;

				omega += crossmul(rs, v);

				pos_fl[l] = mod(pos_fl[l] + vel_fl[l]*dt*0.5, len);
			}
		}
		vel_colloid[j] 	   += vc*mass_fl/mass_colloid;
		ang_vel_colloid[j] += omega*mass_fl/I_colloid;
	}
}
