#include "parameters.hpp"

inline point crossmul(point a, point b) {
	return point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline point stochastic_reflection(point rf, point rs) {
	double m_beta = mass_fl/kbt, random_e = pow(1 - ran(), 2), val, v[4], x[4], z = 2;
	point un, ut, t, n;
	n = rs/sqrt((rs*rs).sum());
	val = sqrt(-log(random_e)/m_beta);

	un = n*val;
	t = img(t.random(rf, len), len);
	ut = crossmul(un, t);
	ut = ut/sqrt((ut*ut).sum());
	while(z > 1) {
		x[1] = 2.0 * ran() - 1, x[2] = 2.0 * ran() - 1;
		z = x[1]*x[1] + x[2]*x[2];
	}
	z = sqrt((-2.0*log(z))/z);
	v[1] = x[1]*z*sqrt(kbt/mass_fl); v[2] = x[2]*z*sqrt(kbt/mass_fl);
	return ut*v[1] + un;
}

void fluid_colloid_collision() {
	point rr, rs, dump_vel_fl[no_of_fluid + 1], u, omega, vc;
	std::copy(vel_fl, vel_fl + no_of_fluid + 1, dump_vel_fl);

	for (int j = 1; j <= no_of_colloid; j++) {
		vc = omega = point(0, 0, 0);
		for (int i = 1; i <= no_neigh[j]; i++) {
			int l = neigh_fl[i][j];
			rr = img(pos_colloid[j] - pos_fl[l], len);
			if((rr*rr).sum() <= pow(sigma, 2)*0.25) {
				pos_fl[l] = mod(pos_fl[l] - vel_fl[l]*dt*0.5, len);
				rs = img(pos_fl[l] - pos_colloid[j], len);
				u  = stochastic_reflection(pos_fl[l], rs);
				vel_fl[l] = u + vel_colloid[j] + crossmul(ang_vel_colloid[j], rs);
                vc += (dump_vel_fl[l] - vel_fl[l]);
				u = (dump_vel_fl[l] - vel_fl[l]);
				omega += crossmul(rs, (dump_vel_fl[l] - vel_fl[l]));
				pos_fl[l] = mod(pos_fl[l] + vel_fl[l]*dt*0.5, len);
			}
		}
		vel_colloid[j] 	   += vc*mass_fl/mass_colloid;
		ang_vel_colloid[j] += omega*mass_fl/I_colloid;
	}
}
