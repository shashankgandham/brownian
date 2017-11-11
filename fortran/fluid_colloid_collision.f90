subroutine fluid_colloid_collision
	use all_parameters
	implicit none
	real*8  ::dump_vel_fl(3*no_of_fluid),rr
	real*8  ::vc1,vc2,vc3,omega1,omega2,omega3,v1,v2,v3,rrx,rry,rrz
	integer ::i,j,lll

!	find  how many fluid particles collide with colloid particles
	do i = 1, 3*no_of_fluid
		dump_vel_fl(i)=vel_fl(i)
	enddo

	do j=1,no_of_colloid
		vc1=0.0d0;vc2=0.0d0;vc3=0.0d0
		omega1=0.0d0;omega2=0.0d0;omega3=0.0d0

		do i=1,no_neigh(j)        ! n_neighbour_fl(j)
			lll=neigh_fl(i,j)
!			do lll=1,no_of_fluid
			rrx=pos_colloid(3*j-2)-pos_fl(3*lll-2)
			rry=pos_colloid(3*j-1)-pos_fl(3*lll-1)
			rrz=pos_colloid(3*j)-pos_fl(3*lll)

			rrx = rrx - llx*anint(rrx*inv_llx)
			rry = rry - lly*anint(rry*inv_lly)
			rrz = rrz - llz*anint(rrz*inv_llz)

!			if(abs(rrx).gt.llxby2) rrx=(llx-abs(rrx))
!			if(abs(rry).gt.llyby2) rry=(lly-abs(rry))
!			if(abs(rrz).gt.llzby2) rrz=(llz-abs(rrz))
			rr = rrx**2 + rry**2 + rrz**2
			if(rr.le.sigma_square*0.250d0) then
				pos_fl(3*lll-2)=pos_fl(3*lll-2)-vel_fl(3*lll-2)*dt*0.5
				pos_fl(3*lll-1)=pos_fl(3*lll-1)-vel_fl(3*lll-1)*dt*0.5
				pos_fl(3*lll)=pos_fl(3*lll)-vel_fl(3*lll)*dt*0.5

				pos_fl(3*lll-2) = pos_fl(3*lll-2) - llx*anint((pos_fl(3*lll-2)-llxby2)*inv_llx)!pbc using minimum image
				pos_fl(3*lll-1) = pos_fl(3*lll-1) - lly*anint((pos_fl(3*lll-1)-llyby2)*inv_lly)
				pos_fl(3*lll) = pos_fl(3*lll) - llz*anint((pos_fl(3*lll)-llzby2)*inv_llz)

!				if(pos_fl(3*i-2).gt.llx)   pos_fl(3*i-2)=pos_fl(3*i-2)-llx    ! !PBC
!				if(pos_fl(3*i-2).le.0.0d0) pos_fl(3*i-2)=pos_fl(3*i-2)+llx
!				if(pos_fl(3*i-1).gt.lly)   pos_fl(3*i-1)=pos_fl(3*i-1)-lly
!				if(pos_fl(3*i-1).le.0.0d0) pos_fl(3*i-1)=pos_fl(3*i-1)+lly
!				if(pos_fl(3*i).gt.llz)     pos_fl(3*i)=pos_fl(3*i)-llz
!				if(pos_fl(3*i).le.0.0d0)   pos_fl(3*i)=pos_fl(3*i)+llz
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

!               if(abs(rsx).gt.llxby2) rsx=(llx-abs(rsx))*(rcx-rfx)/abs(rfx-rcx)
!               if(abs(rsy).gt.llyby2) rsy=(lly-abs(rsy))*(rcy-rfy)/abs(rfy-rcy)
!               if(abs(rsz).gt.llzby2) rsz=(llz-abs(rsz))*(rcz-rfz)/abs(rfz-rcz)
				call stochastic_reflection    !for stick boundary condition
!				write(730,*) rsx**2+rsy**2+rsz**2
				!!velocity of fluid in streaming step!!
				vel_fl(3*lll-2)=u1+vel_colloid(3*j-2)+ang_vel_colloid(3*j-1)*rsz-ang_vel_colloid(3*j)*rsy
				vel_fl(3*lll-1)=u2+vel_colloid(3*j-1)-ang_vel_colloid(3*j-2)*rsz+ang_vel_colloid(3*j)*rsx
				vel_fl(3*lll)=u3+vel_colloid(3*j)+ang_vel_colloid(3*j-2)*rsy-ang_vel_colloid(3*j-1)*rsx

				v1=dump_vel_fl(3*lll-2)-vel_fl(3*lll-2)
				v2=dump_vel_fl(3*lll-1)-vel_fl(3*lll-1)
				v3=dump_vel_fl(3*lll)-vel_fl(3*lll)

				vc1=vc1+v1                               !calculation for linear and
				vc2=vc2+v2                               !angular velocity of the
				vc3=vc3+v3                               !colloidal particle
				omega1=omega1+rsy*v3-v2*rsz   !! ??? NOT UNDERSTOOD
				omega2=omega2-rsx*v3+v1*rsz
				omega3=omega3+rsx*v2-v1*rsy
!				write(734,*) omega1,omega2,omega3
!				!!updating position for next dt/2 time!!

				pos_fl(3*lll-2)=pos_fl(3*lll-2)+vel_fl(3*lll-2)*dt*0.5
				pos_fl(3*lll-1)=pos_fl(3*lll-1)+vel_fl(3*lll-1)*dt*0.5
				pos_fl(3*lll)=pos_fl(3*lll)+vel_fl(3*lll)*dt*0.5
				pos_fl(3*lll-2) = pos_fl(3*lll-2) - llx*anint((pos_fl(3*lll-2)-llxby2)*inv_llx)!pbc using minimum image
				pos_fl(3*lll-1) = pos_fl(3*lll-1) - lly*anint((pos_fl(3*lll-1)-llyby2)*inv_lly)
				pos_fl(3*lll) = pos_fl(3*lll) - llz*anint((pos_fl(3*lll)-llzby2)*inv_llz)

!				if(pos_fl(3*i-2).gt.llx)   pos_fl(3*i-2)=pos_fl(3*i-2)-llx    ! PBC
!				if(pos_fl(3*i-2).le.0.0d0) pos_fl(3*i-2)=pos_fl(3*i-2)+llx
!				if(pos_fl(3*i-1).gt.lly)   pos_fl(3*i-1)=pos_fl(3*i-1)-lly
!				if(pos_fl(3*i-1).le.0.0d0) pos_fl(3*i-1)=pos_fl(3*i-1)+lly
!				if(pos_fl(3*i).gt.llz)     pos_fl(3*i)=pos_fl(3*i)-llz
!				if(pos_fl(3*i).le.0.0d0)   pos_fl(3*i)=pos_fl(3*i)+llz

!				counter(j)=counter(j)+1
			endif
		enddo

!		write(333,*) vc1,vc2,vc3
!		!!velocity of colloidal particle!!
		vel_colloid(3*j-2)=vel_colloid(3*j-2)+vc1*mass_fl/mass_colloid
		vel_colloid(3*j-1)=vel_colloid(3*j-1)+vc2*mass_fl/mass_colloid
		vel_colloid(3*j)=vel_colloid(3*j)+vc3*mass_fl/mass_colloid

!		!!angular velocity of colloidal particle!!
		ang_vel_colloid(3*j-2)=ang_vel_colloid(3*j-2)+omega1*mass_fl/I_colloid
		ang_vel_colloid(3*j-1)=ang_vel_colloid(3*j-1)+omega2*mass_fl/I_colloid
		ang_vel_colloid(3*j)=ang_vel_colloid(3*j)+omega3*mass_fl/I_colloid
	enddo
!	write(052,*) dsqrt(rsx**2+rsy**2+rsz**2)
!	write(033,*) omega1,omega2,omega3
end subroutine fluid_colloid_collision
