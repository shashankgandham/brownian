//done
#include <cstdio>
#include <cmath>
#include "parameters.hpp"
void update_activity_direction(){
	int i;
	double b1,b2,b3, raa[3*no_of_colloid];
	double m11,m12,m13,m21,m22,m23,m31,m32,m33;

	for (i = 1; i <= no_of_colloid; i++){
		b1=ang_vel_colloid[3*i-2]*dt;
		b2=ang_vel_colloid[3*i-1]*dt;
		b3=ang_vel_colloid[3*i]*dt;

		/*!FIRST MAKE A 3D ROTATION MATRIX      */
		m11=cos(b2)*cos(b3);
		m12=-cos(b2)*sin(b3);
		m13=sin(b2);
		m21=sin(b1)*sin(b2)*cos(b3)+cos(b1)*sin(b3);
		m22=-sin(b1)*sin(b2)*sin(b3)+cos(b1)*cos(b3);
		m23=-sin(b1)*cos(b2);
		m31=-cos(b1)*sin(b2)*cos(b3)+sin(b1)*sin(b3);
		m32=cos(b1)*sin(b2)*sin(b3)+sin(b1)*cos(b3);
		m33=cos(b1)*cos(b2);

		/*!UPDATE THE DIRECTION OF ACTIVITY*/
		raa[3*i-2]=ra[3*i-2];
		raa[3*i-1]=ra[3*i-1];
		raa[3*i]=ra[3*i];

		ra[3*i-2]=m11*raa[3*i-2]+m12*raa[3*i-1]+m13*raa[3*i];
		ra[3*i-1]=m21*raa[3*i-2]+m22*raa[3*i-1]+m23*raa[3*i];
		ra[3*i]=m31*raa[3*i-2]+m32*raa[3*i-1]+m33*raa[3*i];
	}
}
