#include "parameters.hpp"
#include <cstdlib>

void initialize_fluid() {
	int i,j,k,counter_fl,check_fl;
	double ran1,x12,y12,z12,r,ke1_fluid,zz1;
	double avr_vel_fl_x,avr_vel_fl_y,avr_vel_fl_z,tx,ty,tz;
	double average_vel_fl_x,average_vel_fl_y,average_vel_fl_z;
	counter_fl = 0, avr_vel_fl_x = 0.0, avr_vel_fl_y = 0.0, avr_vel_fl_z = 0;
	std::srand(zzzz);	
	while(counter_fl < no_of_fluid) {
		tx = rand()*llx;
		ty = rand()*lly;
		tz = rand()*llz;
		check_fl =1;
		for(int j=1; j <= no_of_colloid; j++) {
			x12 = tx - pos_colloid[3*j-2];
			y12 = ty - pos_colloid[3*j-1];
			z12 = tz - pos_colloid[3*j];
			x12 = x12 - llx*round(x12*inv_llx);
			y12 = y12 - lly*round(y12*inv_lly);
			z12 = z12 - llz*round(z12*inv_llz);
			//Does anybody know how mod works???
			r = sqrt(pow(x12, 2) + pow(y12, 2)+ pow(z12, 2));
			if(r < sigma*0.5) {
				check_fl = 0;
			}
		}
		if(check_fl==1) {
			counter_fl = counter_fl + 1; 
			pos_fl[3*counter_fl-2] = tx;        		
			pos_fl[3*counter_fl-1] = ty;
			pos_fl[3*counter_fl] = tz;
		}
	}
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[3*j-2] = (rand()-0.5)*vscale_fluid;
		vel_fl[3*j-1] = (rand()-0.5)*vscale_fluid;
		vel_fl[3*j]   = (rand()-0.5)*vscale_fluid;
		avr_vel_fl_x = avr_vel_fl_x + vel_fl[3*j-2];
		avr_vel_fl_y = avr_vel_fl_y + vel_fl[3*j-1];
		avr_vel_fl_z = avr_vel_fl_z + vel_fl[3*j];
	}
	avr_vel_fl_x = avr_vel_fl_x/(float)no_of_fluid;
	avr_vel_fl_y = avr_vel_fl_y/(float)no_of_fluid;
	avr_vel_fl_z = avr_vel_fl_z/(float)no_of_fluid;

	average_vel_fl_x = 0.0;
	average_vel_fl_y = 0.0;
	average_vel_fl_z = 0.0;
	ke1_fluid=0.0;
	for(int j = 1; j <= no_of_fluid; j++) {
		vel_fl[3*j-1] = vel_fl[3*j-1] - avr_vel_fl_x;
		vel_fl[3*j-2] = vel_fl[3*j-2] - avr_vel_fl_y;
		vel_fl[3*j]   = vel_fl[3*j] - avr_vel_fl_z;
		average_vel_fl_x = average_vel_fl_x + vel_fl[3*j-1];
		average_vel_fl_y = average_vel_fl_y + vel_fl[3*j-2];
		average_vel_fl_z = average_vel_fl_z + vel_fl[3*j];
		ke1_fluid = ke1_fluid + 0.50*mass_fl*(pow(vel_fl[3*j-1], 2) + pow(vel_fl[3*j-2], 2) + pow(vel_fl[3*j], 2));

	}
	for(int i = 1; i <= no_of_fluid; i++) {
		//Write to file bitch
		/*
			write(222,*) pos_fl[3*i-1],pos_fl[3*i-2],pos_fl[3*i]; 
			write(223,*) vel_fl[3*i-1],vel_fl[3*i-2],vel_fl[3*i]; 
		*/
	}
}
