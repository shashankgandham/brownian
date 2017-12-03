subroutine initialize_colloid
	use all_parameters
	implicit none
	integer :: i,j,k,counter,check,mm
	integer :: nofp
	real*8  :: ran1,x12,y12,z12,r,mmm
	real*8  :: avr_vel_colloid_x,avr_vel_colloid_y,avr_vel_colloid_z,tx,ty,tz
	real*8  :: average_vel_colloid_x,average_vel_colloid_y,average_vel_colloid_z

	nofp=0
	do k=40,(lz-1)*10,50
		do j=40,(ly-1)*10,50
			do i=40,(lx-1)*10,50
				nofp=nofp+1
				if(nofp.le.no_of_colloid)then
					pos_colloid(3*nofp-2)=dfloat(i)/10d0  !initialize in cubic crystal
					pos_colloid(3*nofp-1)=dfloat(j)/10d0  !positions
					pos_colloid(3*nofp)=dfloat(k)/10d0
				endif
			enddo
		enddo
	enddo

	counter = 0!no_of_colloid  !set counter =0 to use this process of initializing position
	avr_vel_colloid_x = 0.0d0
	avr_vel_colloid_y = 0.0d0
	avr_vel_colloid_z = 0.0d0
	do while(counter < no_of_colloid)
		tx = ran1(zzzz)*llx
		ty = ran1(zzzz)*lly
		tz = ran1(zzzz)*llz

		check =1
		do j = 1, counter
			x12 = tx - pos_colloid(3*j-2)
			y12 = ty - pos_colloid(3*j-1)
			z12 = tz - pos_colloid(3*j)

			if(abs(x12)>llxby2) x12 = llx - abs(x12)
			if(abs(y12)>llyby2) y12 = lly - abs(y12)
			if(abs(z12)>llzby2) z12 = llz - abs(z12)
			r = dsqrt(x12**2+y12**2+z12**2)

			if(r.lt.space_limit) then
				check =0
			endif
		enddo

		if(check==1) then
			counter = counter +1
			pos_colloid(3*counter-2) = tx
			pos_colloid(3*counter-1) = ty
			pos_colloid(3*counter) = tz
		endif
!		write(*,*) 'counter',counter
	enddo

!	pos_colloid(1)=0.50d0
!	pos_colloid(2)=0.50d0
!   pos_colloid(3)=0.50d0
!   pos_colloid(4)=5.50d0
!   pos_colloid(5)=5.50d0
!   pos_colloid(6)=5.50d0
!   pos_colloid(7)=10.50d0
!   pos_colloid(8)=10.50d0
!   pos_colloid(9)=10.50d0
!   pos_colloid(4)=15.0d0
!   pos_colloid(5)=15.0d0
!   pos_colloid(6)=15.0d0
	mmm = dsqrt(kbt1/mass_colloid)

!	do i=1,no_of_colloid
!		vel_colloid(3*i-2)=mmm
!		vel_colloid(3*i-1)=mmm
!		vel_colloid(3*i)=mmm
!	enddo
!	vel_colloid(4)=0.0d0
!	vel_colloid(5)=0.0d0
!	vel_colloid(6)=0.0d0
!	vel_colloid(7)=0.0d0
!	vel_colloid(8)=0.0d0
!	vel_colloid(9)=0.0d0
!	vel_colloid(4)=-mmm
!	vel_colloid(5)=-mmm
!	vel_colloid(6)=-mmm

	do j = 1, no_of_colloid
		vel_colloid(3*j-2) = (ran1(zzzz)-0.5d0)*vscale_colloid
		vel_colloid(3*j-1) = (ran1(zzzz)-0.5d0)*vscale_colloid
		vel_colloid(3*j) = (ran1(zzzz)-0.5d0)*vscale_colloid
		avr_vel_colloid_x = avr_vel_colloid_x + vel_colloid(3*j-2)
		avr_vel_colloid_y = avr_vel_colloid_y + vel_colloid(3*j-1)
		avr_vel_colloid_z = avr_vel_colloid_z + vel_colloid(3*j)

!		if((vel_colloid(3*j-2)**2 + vel_colloid(3*j-1)**2  +vel_colloid(3*j)**2) > 10.0d0) write(*,*) vel_colloid(3*i-2),i,'largev'
	enddo

!	write(*,*) 'no of particles initialized', counter,N
	avr_vel_colloid_x = avr_vel_colloid_x/dfloat(no_of_colloid)
	avr_vel_colloid_y = avr_vel_colloid_y/dfloat(no_of_colloid)
	avr_vel_colloid_z = avr_vel_colloid_z/dfloat(no_of_colloid) !NOTE

	average_vel_colloid_x = 0.0d0
	average_vel_colloid_y = 0.0d0
	average_vel_colloid_z = 0.0d0

	ke_colloid1 = 0.0d0
	do j = 1, no_of_colloid
		vel_colloid(3*j-2) = vel_colloid(3*j-2) - avr_vel_colloid_x
		vel_colloid(3*j-1) = vel_colloid(3*j-1) - avr_vel_colloid_y
		vel_colloid(3*j) = vel_colloid(3*j) - avr_vel_colloid_z
!		write(*,*) vel_colloid(3*j-2)
		average_vel_colloid_x = average_vel_colloid_x + vel_colloid(3*j-2)
		average_vel_colloid_y = average_vel_colloid_y + vel_colloid(3*j-1)
		average_vel_colloid_z = average_vel_colloid_z + vel_colloid(3*j)
		ke_colloid1 = ke_colloid1 + vel_colloid(3*j-2)**2 + vel_colloid(3*j-1)**2 + vel_colloid(3*j)**2
	enddo

!	!!INITIALIZING THE ANGULAR VELOCITY OF COLLOID!!

	do j=1,no_of_colloid
		ang_vel_colloid(3*j-2) = (ran1(zzzz)-0.5d0)*ang_vscale_colloid
		ang_vel_colloid(3*j-1) = (ran1(zzzz)-0.5d0)*ang_vscale_colloid
		ang_vel_colloid(3*j)   = (ran1(zzzz)-0.5d0)*ang_vscale_colloid
	enddo

	ang_ke_colloid1=0.0d0
	do j = 1, no_of_colloid
		ang_ke_colloid1 = ang_ke_colloid1 + ang_vel_colloid(3*j-2)**2 + ang_vel_colloid(3*j-1)**2 + ang_vel_colloid(3*j)**2
	enddo
!	do j=1,N
!		write(1,fmt='(8g25.15)') pos_colloid(3*j-2),pos_colloid(3*j-1),pos_colloid(3*j),vel_colloid(3*j-2),vel_colloid(3*j-1),vel_colloid(3*j)
!   enddo
!   close(1)
end subroutine
