#include <cstdio>
#include <cmath>
#include "parameters.hpp"

void fluid_colloid_collision() {
  double dump_vel_fl[3*no_of_fluid],rr;
  double vc1,vc2,vc3,omega1,omega2,omega3,v1,v2,v3,rrx,rry,rrz;
  int i,j,lll;

  /*find  how many fluid particles collide with colloid particles  */

  for (i = 3; i <= 3*no_of_fluid; i++)
    dump_vel_fl[i] = vel_fl[i];

  for (j = 1; j <= no_of_colloid; j++){
    vc1=0.0;
    vc2=0.0;
    vc3=0.0;
    omega1=0.0;
    omega2=0.0;
    omega3=0.0;      
    for (i = 1; i <= no_neigh[j]; i++){
      lll = neigh_fl[i][j];
      rrx=pos_colloid[3*j-2]-pos_fl[3*lll-2];
      rry=pos_colloid[3*j-1]-pos_fl[3*lll-1];
      rrz=pos_colloid[3*j]-pos_fl[3*lll];

      rrx = rrx - llx*round(rrx*inv_llx);
      rry = rry - lly*round(rry*inv_lly);
      rrz = rrz - llz*round(rrz*inv_llz);

      rr=rrx*rrx+rry*rry+rrz*rrz;
      
      if(rr <= sigma_square*0.25){
        pos_fl[3*lll-2]=pos_fl[3*lll-2]-vel_fl[3*lll-2]*dt*0.5;
        pos_fl[3*lll-1]=pos_fl[3*lll-1]-vel_fl[3*lll-1]*dt*0.5;
        pos_fl[3*lll]=pos_fl[3*lll]-vel_fl[3*lll]*dt*0.5;
      

        pos_fl[3*lll-2] = pos_fl[3*lll-2] - llx*round((pos_fl[3*lll-2]-llxby2)*inv_llx);
        pos_fl[3*lll-1] = pos_fl[3*lll-1] - lly*round((pos_fl[3*lll-1]-llyby2)*inv_lly);
        pos_fl[3*lll] = pos_fl[3*lll] - llz*round((pos_fl[3*lll]-llzby2)*inv_llz);

        rcx=pos_colloid[3*j-2];
        rcy=pos_colloid[3*j-1];
        rcz=pos_colloid[3*j];

        rfx=pos_fl[3*lll-2];
        rfy=pos_fl[3*lll-1];
        rfz=pos_fl[3*lll];

        /*shortest distance vector  of fluid particle from colloid centre*/
        rsx=rfx-rcx; 
        rsy=rfy-rcy;
        rsz=rfz-rcz;

        rsx =rsx - llx*round(rsx*inv_llx);
        rsy =rsy - lly*round(rsy*inv_lly);
        rsz =rsz - llz*round(rsz*inv_llz);

        /*!for stick boundary condition*/
        stochastic_reflection();
           
        /*velocity of fluid in streaming step */
        vel_fl[3*lll-2]=u1+vel_colloid[3*j-2]+ang_vel_colloid[3*j-1]*rsz-ang_vel_colloid[3*j]*rsy; 
        vel_fl[3*lll-1]=u2+vel_colloid[3*j-1]-ang_vel_colloid[3*j-2]*rsz+ang_vel_colloid[3*j]*rsx; 
        vel_fl[3*lll]=u3+vel_colloid[3*j]+ang_vel_colloid[3*j-2]*rsy-ang_vel_colloid[3*j-1]*rsx;

        v1=dump_vel_fl[3*lll-2]-vel_fl[3*lll-2];
        v2=dump_vel_fl[3*lll-1]-vel_fl[3*lll-1];
        v3=dump_vel_fl[3*lll]-vel_fl[3*lll];

        /*calculation for linear and angular velocity of the colloidal particle*/
        vc1=vc1+v1;                               
        vc2=vc2+v2;                              
        vc3=vc3+v3;                               
        omega1=omega1+rsy*v3-v2*rsz;   
        omega2=omega2-rsx*v3+v1*rsz;
        omega3=omega3+rsx*v2-v1*rsy;

        /*updating position for next dt/2 time */
        pos_fl[3*lll-2]=pos_fl[3*lll-2]+vel_fl[3*lll-2]*dt*0.5;
        pos_fl[3*lll-1]=pos_fl[3*lll-1]+vel_fl[3*lll-1]*dt*0.5;
        pos_fl[3*lll]=pos_fl[3*lll]+vel_fl[3*lll]*dt*0.5;

        pos_fl[3*lll-2] = pos_fl[3*lll-2] - llx*round((pos_fl[3*lll-2]-llxby2)*inv_llx);
        pos_fl[3*lll-1] = pos_fl[3*lll-1] - lly*round((pos_fl[3*lll-1]-llyby2)*inv_lly);
        pos_fl[3*lll] = pos_fl[3*lll] - llz*round((pos_fl[3*lll]-llzby2)*inv_llz);
        }
    }
    vel_colloid[3*j-2]=vel_colloid[3*j-2]+vc1*mass_fl/mass_colloid; 
    vel_colloid[3*j-1]=vel_colloid[3*j-1]+vc2*mass_fl/mass_colloid;
    vel_colloid[3*j]=vel_colloid[3*j]+vc3*mass_fl/mass_colloid;

    /*angular velocity of colloidal particle */

    ang_vel_colloid[3*j-2]=ang_vel_colloid[3*j-2]+omega1*mass_fl/I_colloid;
    ang_vel_colloid[3*j-1]=ang_vel_colloid[3*j-1]+omega2*mass_fl/I_colloid;
    ang_vel_colloid[3*j]=ang_vel_colloid[3*j]+omega3*mass_fl/I_colloid;
  }
}       
