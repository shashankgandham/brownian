//done
#include <cmath>
#include "parameters.hpp"
void update_pos_mpcd(){
  int i;
  for (i = 1; i <= no_of_fluid; i++){
    pos_fl[3*i-2]=pos_fl[3*i-2]+vel_fl[3*i-2]*dt;
    pos_fl[3*i-1]=pos_fl[3*i-1]+vel_fl[3*i-1]*dt;
    pos_fl[3*i]  =pos_fl[3*i]+vel_fl[3*i]*dt;

    pos_fl[3*i-2] = pos_fl[3*i-2] - llx*round((pos_fl[3*i-2]-llxby2)*inv_llx);
    pos_fl[3*i-1] = pos_fl[3*i-1] - lly*round((pos_fl[3*i-1]-llyby2)*inv_lly);
    pos_fl[3*i] = pos_fl[3*i] - llz*round((pos_fl[3*i]-llzby2)*inv_llz);
  }
}
