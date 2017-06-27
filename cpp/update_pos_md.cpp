void update_pos_md() {
	int i;
	double ddt;
	ddt = 0.05*dt2/mass_colloid;
	for(int i = 1; i <= no_of_colloid; i++) {
		dx = vel_colloid[3*i - 2]*dt + f[3*i -2]*ddt;
		dy = vel_colloid[3*i - 1]*dt + f[3*i - 1]*ddt;
		dz = vel_colloid[3*i]*dt + f[3*i]*ddt;
		pos_colloid[3*i - 2] += dx;
		pos_colloid[3*i - 1] += dy;
		pos_colloid[3*i] += dz;
		if(pos_colloid[3*i - 2] > llx) 
			pos_colloid[3*i - 2] -= llx;
		if(pos_colloid[3*i - 1] > lly) 
			pos_colloid[3*i - 1] -= lly;
		if(pos_colloid[3*i] > lly) 
			pos_colloid[3*i] -= llz;
	}
}
