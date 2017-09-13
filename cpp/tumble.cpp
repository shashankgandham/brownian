//done
#include <cstdlib>
#include <cmath>
#include "parameters.hpp"

void tumble(){
  int i;
  double b1,b2,b3,aaa;
  std::srand(zzzz);
  
  for (i = 1; i <= no_of_colloid; i++){
    b1=rand()*llx;
    b2=rand()*lly;
    b3=rand()*llz;
  
    ra[3*i-2]= pos_colloid[3*i-2]-b1;
    ra[3*i-1]= pos_colloid[3*i-1]-b2;
    ra[3*i]= pos_colloid[3*i]-b3;

    /* minimum image convention*/
    ra[3*i-2] = ra[3*i-2] - llx*round(ra[3*i-2]*inv_llx);
    ra[3*i-1] = ra[3*i-1] - llx*round(ra[3*i-1]*inv_llx);
    ra[3*i] = ra[3*i] - llx*round(ra[3*i]*inv_llx);

    aaa=sqrt(ra[3*i-2]*ra[3*i-2]+ra[3*i-1]*ra[3*i-1]+ra[3*i]*ra[3*i]);

    /*RA IS THE UNIT VECTOR TOWARDS THE ACTIVE VELOCITY*/
    ra[3*i-2]=ra[3*i-2]/aaa; 
    ra[3*i-1]=ra[3*i-1]/aaa;
    ra[3*i]=ra[3*i]/aaa;
  }
}
