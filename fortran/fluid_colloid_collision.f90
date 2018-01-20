subroutine fluid_colloid_collision
	use all_parameters
	implicit none
	double precision  ::dump_vel_fl(3*no_of_fluid),rr, cx, cy, cz, temp1, temp2, temp3
	double precision  ::vc1,vc2,vc3,omega1,omega2,omega3,v1,v2,v3,rrx,rry,rrz
	integer ::i,j,lll

	do i=1,3*no_of_fluid            !ei gota calculation ta () Gomper er paper ar padding er paper dekhe kora.
		dump_vel_fl(i)=vel_fl(i)    !ar anno kono paper dekhar darkar nei.baki gulo bhul bhal.prochur somoi nosto korbe
   	enddo                          !kintu kajer kaaj kichu hobe na

	do j=1,no_of_colloid
		vc1=0.0d0;vc2=0.0d0;vc3=0.0d0
		omega1=0.0d0;omega2=0.0d0;omega3=0.0d0
		do i=1,no_neigh(j)    ! n_neighbour_fl(j)
			lll=neigh_fl(i,j)
            rrx=pos_colloid(3*j-2)-pos_fl(3*lll-2)
            rry=pos_colloid(3*j-1)-pos_fl(3*lll-1)
            rrz=pos_colloid(3*j)-pos_fl(3*lll)
           	rrx = rrx - llx*anint(rrx*inv_llx)
           	rry = rry - lly*anint(rry*inv_lly)
           	rrz = rrz - llz*anint(rrz*inv_llz)

			rr=rrx**2+rry**2+rrz**2
			if(rr.le.sigma_square*0.250d0) then
            	pos_fl(3*lll-2)=pos_fl(3*lll-2)-vel_fl(3*lll-2)*dt*0.5  !!!???? How?
             	pos_fl(3*lll-1)=pos_fl(3*lll-1)-vel_fl(3*lll-1)*dt*0.5
             	pos_fl(3*lll)=pos_fl(3*lll)-vel_fl(3*lll)*dt*0.5

            	pos_fl(3*lll-2) = pos_fl(3*lll-2) - llx*anint((pos_fl(3*lll-2)-llxby2)*inv_llx)!pbc using minimum image
            	pos_fl(3*lll-1) = pos_fl(3*lll-1) - lly*anint((pos_fl(3*lll-1)-llyby2)*inv_lly)
            	pos_fl(3*lll) = pos_fl(3*lll) - llz*anint((pos_fl(3*lll)-llzby2)*inv_llz)

             	rcx=pos_colloid(3*j-2)
             	rcy=pos_colloid(3*j-1)
             	rcz=pos_colloid(3*j)

             	rfx=pos_fl(3*lll-2)
             	rfy=pos_fl(3*lll-1)
             	rfz=pos_fl(3*lll)

             	rsx=rfx-rcx !shortest distance vector  of fluid particle from colloid centre
             	rsy=rfy-rcy
             	rsz=rfz-rcz

                rsx =rsx - llx*anint(rsx*inv_llx)
                rsy =rsy - lly*anint(rsy*inv_lly)
                rsz =rsz - llz*anint(rsz*inv_llz)


            	call stochastic_reflection    !for stick boundary condition
                
				cx = ang_vel_colloid(3*j-1)*rsz-ang_vel_colloid(3*j)*rsy
				cy = ang_vel_colloid(3*j)*rsx -ang_vel_colloid(3*j-2)*rsz
				cz = ang_vel_colloid(3*j-2)*rsy-ang_vel_colloid(3*j-1)*rsx
           		vel_fl(3*lll-2)=u1+vel_colloid(3*j-2)+cx
           		vel_fl(3*lll-1)=u2+vel_colloid(3*j-1)+cy
           		vel_fl(3*lll)=u3+vel_colloid(3*j)+cz
       !         if(nn==94) then
        !       write(*, fmt='(3F36.32)') vel_fl(3*lll-2), vel_fl(3*lll-1), vel_fl(3*lll)
        !      endif
                v1=dump_vel_fl(3*lll-2)-vel_fl(3*lll-2)
          		v2=dump_vel_fl(3*lll-1)-vel_fl(3*lll-1)
          		v3=dump_vel_fl(3*lll)-vel_fl(3*lll)

          		vc1=vc1+v1                               !calculation for linear and
          		vc2=vc2+v2                               !angular velocity of the
          		vc3=vc3+v3                               !colloidal particle
          		temp1 = rsy*v3-rsz*v2; temp2 = rsz*v1-rsx*v3; temp3=rsx*v2-rsy*v1

				omega1=omega1+temp1   !! ??? NOT UNDERSTOOD
          		omega2=omega2+temp2
          		omega3=omega3+temp3
			 	pos_fl(3*lll-2)=pos_fl(3*lll-2)+vel_fl(3*lll-2)*dt*0.5
             	pos_fl(3*lll-1)=pos_fl(3*lll-1)+vel_fl(3*lll-1)*dt*0.5
             	pos_fl(3*lll)=pos_fl(3*lll)+vel_fl(3*lll)*dt*0.5

            	pos_fl(3*lll-2) = pos_fl(3*lll-2) - llx*anint((pos_fl(3*lll-2)-llxby2)*inv_llx)!pbc using minimum image
            	pos_fl(3*lll-1) = pos_fl(3*lll-1) - lly*anint((pos_fl(3*lll-1)-llyby2)*inv_lly)
            	pos_fl(3*lll) = pos_fl(3*lll) - llz*anint((pos_fl(3*lll)-llzby2)*inv_llz)
			endif
		enddo
		vel_colloid(3*j-2)=vel_colloid(3*j-2)+vc1*mass_fl/mass_colloid
		vel_colloid(3*j-1)=vel_colloid(3*j-1)+vc2*mass_fl/mass_colloid
		vel_colloid(3*j)=vel_colloid(3*j)+vc3*mass_fl/mass_colloid

		cx = omega1*mass_fl/I_colloid
		cy = omega2*mass_fl/I_colloid
		cz = omega3*mass_fl/I_colloid
        ang_vel_colloid(3*j-2)=ang_vel_colloid(3*j-2)+cx
        ang_vel_colloid(3*j-1)=ang_vel_colloid(3*j-1)+cy
        ang_vel_colloid(3*j)=ang_vel_colloid(3*j)+cz
    enddo
end subroutine fluid_colloid_collision
