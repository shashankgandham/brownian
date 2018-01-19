subroutine rotation_mpcd      !rotation step for all fluid particles
	use all_parameters
	implicit none
    integer :: i,j,k,l1,l2,m
    integer :: fluid_no(lx*ly*lz),cell_part(maxpart,lx*ly*lz)
    double precision :: cell_velx(lx*ly*lz), cell_vely(lx*ly*lz),ct
    double precision :: r1, r2, r3, theta, ran1, phi, rho, prob, st
    double precision :: rot11,rot12,rot13,rot21,rot22,rot23,rot31
	double precision :: rot32,rot33,ir1, ir2, ir3
	double precision :: del_vx,del_vy,del_vz,mx1,my1,mz1, ict
	double precision :: tmp_pos(3*no_of_fluid), cell_velz(lx*ly*lz)
	double precision :: rr1,rr2,rr3,var,del_vx1(no_of_fluid),del_vy1(no_of_fluid),del_vz1(no_of_fluid)
	integer::cell_no

	fluid_no = 0

	rr1 = ran1(zzzz) - 0.5d0
	rr2 = ran1(zzzz) - 0.5d0
	rr3 = ran1(zzzz) - 0.5d0
	tmp_pos = pos_fl
	do i=1,no_of_fluid
    	tmp_pos(3*i-2)= tmp_pos(3*i-2) + rr1
        tmp_pos(3*i-1)= tmp_pos(3*i-1) + rr2
        tmp_pos(3*i)  = tmp_pos(3*i)   + rr3
        tmp_pos(3*i-2)= tmp_pos(3*i-2) - llx*anint((tmp_pos(3*i-2)-llxby2)*inv_llx)
        tmp_pos(3*i-1)= tmp_pos(3*i-1) - lly*anint((tmp_pos(3*i-1)-llyby2)*inv_lly)
        tmp_pos(3*i)= tmp_pos(3*i) - llz*anint((tmp_pos(3*i)-llzby2)*inv_llz)
    enddo

	do i=1,no_of_fluid
		cell_no = 1 + int(tmp_pos(3*i-2)) + lx*int(tmp_pos(3*i-1)) + lx*ly*int(tmp_pos(3*i))
		fluid_no(cell_no) = fluid_no(cell_no) + 1
        j = fluid_no(cell_no)
        cell_part(j,cell_no) = i
    end do

	cell_velx = 0.0d0; cell_vely = 0.0d0; cell_velz = 0.0d0

	do i=1,lx*ly*lz
    	if (fluid_no(i).gt.1) then
			do j=1,fluid_no(i)
           		k= cell_part(j,i)
           		cell_velx(i) = cell_velx(i) + vel_fl(3*k-2)/dfloat(fluid_no(i))
           		cell_vely(i) = cell_vely(i) + vel_fl(3*k-1)/dfloat(fluid_no(i))
           		cell_velz(i) = cell_velz(i) + vel_fl(3*k)/dfloat(fluid_no(i))
        	end do
        	rho = 2.0d0*ran1(zzzz)-1.0d0
        	phi = 4.0d0*asin(1.0d0)*ran1(zzzz)
        	r1 = cos(phi)*dsqrt(1-rho*rho)
        	r2 = sin(phi)*dsqrt(1-rho*rho)
        	r3 = rho
        	theta = 2.0d0*asin(1.0d0)*(130.0/180.0)    !!Why muliply by 2
			st = sin(theta)
			ct = cos(theta)
			ict = 1.0d0 - ct
			ir1 = ict*r1; ir2 = ict*r2; ir3 = ict*r3
        	rot11 = ir1*r1 + ct; rot12 = ir1*r2- st*r3; rot13 = ir1*r3 + st*r2
        	rot21 = ir2*r1 + st*r3; rot22 = ir2*r2 + ct; rot23 = ir2*r3 - st*r1
        	rot31 = ir3*r1 - st*r2; rot32 = ir3*r2 + st*r1; rot33 = ir3*r3 + ct
	    	do j=1,fluid_no(i)
           		k= cell_part(j,i)
           		del_vx = vel_fl(3*k-2) - cell_velx(i)
           		del_vy = vel_fl(3*k-1) - cell_vely(i)
           		del_vz = vel_fl(3*k)   - cell_velz(i)
		   		vel_fl(3*k-2) = rot11*del_vx + rot12*del_vy + rot13*del_vz
           		vel_fl(3*k-1) = rot21*del_vx + rot22*del_vy + rot23*del_vz
         	 	vel_fl(3*k)   = rot31*del_vx + rot32*del_vy + rot33*del_vz
		   		vel_fl(3*k-2) = cell_velx(i) + vel_fl(3*k-2)
		   		vel_fl(3*k-1) = cell_vely(i) + vel_fl(3*k-1)
           		vel_fl(3*k)   = cell_velz(i) + vel_fl(3*k)
			end do
		end if
	end do
    if (mod(nn,1)==0) then
		do i=1,lx*ly*lz
        	var = 0.0d0
         	if(fluid_no(i).gt.1) then
           		do j=1,fluid_no(i)
              		k=cell_part(j,i)
              		del_vx1(k) = vel_fl(3*k-2) - cell_velx(i)
              		del_vy1(k) = vel_fl(3*k-1) - cell_vely(i)
              		del_vz1(k) = vel_fl(3*k)   - cell_velz(i)
              		ict = del_vx1(k)**2 + del_vy1(k)**2 + del_vz1(k)**2
					var = var + ict
				end do
          		scale_fac_mpcd = dsqrt(3.0d0*(fluid_no(i)-1)*kbt/(mass_fl*var))
          		do j=1,fluid_no(i)
          			k=cell_part(j,i)
            		vel_fl(3*k-2) = cell_velx(i) + del_vx1(k)*scale_fac_mpcd
            		vel_fl(3*k-1) = cell_vely(i) + del_vy1(k)*scale_fac_mpcd
            		vel_fl(3*k)   = cell_velz(i) + del_vz1(k)*scale_fac_mpcd
         		end do
        	endif
		end do
	end if
end subroutine rotation_mpcd
