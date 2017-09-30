//done
#include <cstdlib>
#include <cmath>
#include "parameters.hpp"

void tumble(){
  double b[4],aaa;
  for (int i = 1; i <= no_of_colloid; i++){
    b[1] = ran()*lx, b[2] = ran()*ly, b[3] = ran()*lz;

    ra[3*i-2] = mod(pos_colloid[3*i-2] - b[1], lx);
    ra[3*i-1] = mod(pos_colloid[3*i-1]- b[2], ly);
    ra[3*i] = mod(pos_colloid[3*i] - b[3], lz);

    aaa = sqrt(ra[3*i-2]*ra[3*i-2] + ra[3*i-1]*ra[3*i-1] + ra[3*i]*ra[3*i]);

    ra[3*i-2] /= aaa, ra[3*i-1] /= aaa, ra[3*i] /= aaa;
  }
}
