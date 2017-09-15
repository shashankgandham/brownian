//done
#include "parameters.hpp"
#include <cmath>
#include <cstdlib>

void stochastic_reflection() {
	double nx,ny,nz,unx,uny,unz,xx,m_beta,random_e,ux,uy,uz, uuu,utx,uty,txx,tyy,tzz,ttt,utz,ran1,val,v1,v2, den;
	m_beta=mass_fl/kbt;
	den = sqrt(rsx*rsx + rsy*rsy + rsz*rsz);
	nx = (rsx)/den;
	ny = (rsy)/den;
	nz = (rsz)/den;
	xx=((double)rand()/RAND_MAX);
	random_e=(1-xx)*(1-xx);
	val = sqrt(-log(random_e)/m_beta);
	unx=val*nx;
	uny=val*ny;
	unz=val*nz;

	txx=((double)rand()/RAND_MAX)*llx;
	tyy=((double)rand()/RAND_MAX)*lly;
	tzz=((double)rand()/RAND_MAX)*llz;

	txx=txx-rfx;
	tyy=tyy-rfy;
	tzz=tzz-rfz;
	txx = txx - llx*round(txx*inv_llx);
	tyy = tyy - lly*round(tyy*inv_lly);
	tzz = tzz - llz*round(tzz*inv_llz);
	
	utx=uny*tzz-tyy*unz;
	uty=unz*txx-unx*tzz;
	utz=unx*tyy-txx*uny;

	uuu = sqrt(utx*utx+uty*uty+utz*utz);

	utx=utx/uuu;
	uty=uty/uuu;
	utz=utz/uuu;
	gauss(v1, v2);
	utx=v1*utx; 
	uty=v1*uty;
	utz=v1*utz;

	u1=unx+utx;
	u2=uny+uty;
	u3=unz+utz;
}
